"""One Claude account's subscription token, minted in GitHub Actions for the control panel.

The Tokens page's Add a Claude account (https://openspineconsortium.com/
symphony/control/tokens.html) used to open a Codespace, which took minutes
to start. This is the same sign-in on an Actions runner, which starts in
seconds: .github/workflows/token.yml installs Claude Code and runs this
script, which runs `claude setup-token` behind a pseudo-terminal as
scripts/tokens.py does. With no browser to open, claude prints its sign-in
link and waits for the code claude.ai shows after the sign-in. The script
posts the link on the issue the page opened, the page shows it as a
button, the person signs in and the page posts the code back on the issue
as the person; the script hands it to claude, takes the token out of what
claude prints (tokens.py's Masker, so no line of claude's reaches the log
with it), seals it for the secrets public key the page gave
(tools/sealbox.py) and posts only the sealed value. The page stores that
value as the slot's secret with GitHub's secrets API, which opens it; the
token itself is never on the page, in the log or in the issue.

The code alone mints nothing: claude's sign-in is PKCE, so only the claude
that printed the link, holding its verifier, can exchange it, and only
once. Only a comment by the person who started the run counts as the code,
and the script deletes it once read.

For someone who may write to the repository but may not write its
organization's secrets, which GitHub leaves to the organization's owners,
the run stores the token itself (the owner's choice of 2026-10-18): the
workflow's job in the environment token-store mints the SymphonyBoard App's
installation token, and the script puts the token into the first free slot
of the organization's pool, never one that holds a token, adds the slot to
SYMPHONY_CLAUDE_ACCOUNTS, marks it in SYMPHONY_CLAUDE_SETTINGS as the
person's own and posts only the slot's number.

  python3 tools/token_run.py         with ISSUE, KEY, KEY_ID, ACTOR, GITHUB_REPOSITORY and GITHUB_TOKEN set,
                                     or STORE_ORG and TOKEN_FILE in place of KEY and KEY_ID
  python3 tools/token_run.py store   then, with APP_TOKEN and REPO_ID too, once claude has gone
  python3 tools/token_run.py look    before the sign-in: lends the person's own token here if they have one
  python3 tools/token_run.py manage  lists, revokes or deletes the person's own token (MANAGE, TARGET_REPO_ID)
  python3 tools/token_run.py fail    where a step before the script failed, FAIL_WHY saying why

A run GitHub started on an issue (GITHUB_EVENT_NAME issues) takes the issue's
author as its actor and the issue's first lines as its request: anyone who
may open an issue can lend a token, with no write access and no page of
ours, and the run closes the issue when it is done (the owner, 2026-10-19:
"it can start on an issue anyone can open").
"""

import base64
import codecs
import json
import os
import re
import select
import sys
import time
import urllib.error
import urllib.request

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "scripts"))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import sealbox  # noqa: E402
import tokens  # noqa: E402

LINK_MARK = "<!-- symphony-token: sign-in -->"
CODE_MARK = "<!-- symphony-token: code -->"
SEALED_MARK = "<!-- symphony-token: sealed -->"
STORED_MARK = "<!-- symphony-token: stored -->"
# what the run says once it has lent, revoked or deleted a person's token, or
# listed where it is lent: the slot, then the ids of the repositories it is lent
# to, never their names, which the issue's readers may not be allowed to see
LENT_MARK = "<!-- symphony-token: lent -->"
FAILED_MARK = "<!-- symphony-token: failed -->"
# what an issue opened for a run says on its first line, since GitHub takes a
# workflow dispatch only from someone who may write and an issue from anyone
# who may read: lend, with "organization" on the second line for a lend with
# no say, which goes where the organization's standing list says and names
# no repository to the lender; manage, with the action and a repository id
LEND_MARK = "<!-- symphony-token: lend -->"
MANAGE_MARK = "<!-- symphony-token: manage -->"
# the organization variable that lists, by id, the repositories a token lent
# with no say runs on, which an owner sets once from the Tokens page
LEND_TO_VARIABLE = "SYMPHONY_LEND_TO"
ORGANIZATION = "organization"
SIGN_IN_LINK = re.compile(r"https://claude\.(?:ai|com)/[A-Za-z0-9/._-]*oauth/authorize\?[A-Za-z0-9%&=._~+-]+")
CODE_SHAPE = re.compile(r"^[A-Za-z0-9_.~#-]{8,512}$")
ISSUE_SHAPE = re.compile(r"^[1-9][0-9]{0,8}$")
KEY_ID_SHAPE = re.compile(r"^[A-Za-z0-9_-]{1,64}$")
WAIT_FOR_LINK = 90
WAIT_FOR_CODE = 600
WAIT_FOR_TOKEN = 60
# claude asks for the code under the link, and a terminal switched to raw
# mode after a key arrived throws that key away, so the code waits for the
# question; past this long it is typed anyway, as the question may be worded
# otherwise in a later claude.
WAIT_FOR_PROMPT = 30
PROMPT = "Paste code here"
# claude's screen takes a chunk of input holding more than one key as a
# paste, a line break in it included, so the code and the Enter that sends
# it go in two writes with this pause between.
PASTE_SETTLE = 1.0
POLL_EVERY = 3
TERMINAL_TYPE = "xterm-256color"
# The organization's pool, which the run writes itself for someone who may
# not write it as themselves (the owner's choice of 2026-10-18): a free slot
# only, never one that holds a token, kept as the person's own.
POOL_SLOTS = tuple(str(n) for n in range(1, 33))
SLOT_PREFIX = "SYMPHONY_CLAUDE_"
ACCOUNTS_VARIABLE = "SYMPHONY_CLAUDE_ACCOUNTS"
SETTINGS_VARIABLE = "SYMPHONY_CLAUDE_SETTINGS"
ORG_SHAPE = re.compile(r"^[A-Za-z0-9](?:[A-Za-z0-9]|-(?=[A-Za-z0-9])){0,38}$")
# What claude must never see: the run's tokens, which the child would inherit.
WITHHELD = ("GITHUB_TOKEN", "APP_TOKEN")
# The App's token is minted only after claude has gone, since a process can
# read its parent's environment however the child's is trimmed.
API = "https://api.github.com"


class RunError(Exception):
    """Why no token was sealed, in words the page shows; never a token."""


def missing_permission(error, org):
    """What to do when GitHub refused the App's token for want of a
    permission, or an empty string. GitHub answers 404, not 403, on an
    organization's variables when the installation lacks Variables, and says
    in X-Accepted-GitHub-Permissions which permission the call needed
    (spine-program's Add Token, 2026-10-20), so the person is told which one
    and where an owner accepts it."""
    if error.code not in (403, 404) or not ORG_SHAPE.match(org or ""):
        return ""
    needed = (error.headers.get("X-Accepted-GitHub-Permissions", "") if error.headers else "").strip()
    if not needed:
        return ""
    words = ", ".join(part.strip().replace("organization_actions_", "organization ").replace("_", " ").replace("=", " ")
                      for part in needed.split(",") if part.strip())
    return (f"; the call needs {words}, which SymphonyBoard's installation on {org} does not hold yet: an owner of {org} "
            f"accepts it at https://github.com/organizations/{org}/settings/installations, then press Add Token again")


class GitHub:
    """The few issue calls the run makes, on the run's own GITHUB_TOKEN."""

    def __init__(self, repo, token, issue):
        self.repo, self.token, self.issue = repo, token, issue

    def call(self, method, path, body=None):
        request = urllib.request.Request(
            API + path, method=method, data=None if body is None else json.dumps(body).encode(),
            headers={"Authorization": f"Bearer {self.token}", "Accept": "application/vnd.github+json",
                     "X-GitHub-Api-Version": "2022-11-28", "Content-Type": "application/json"})
        with urllib.request.urlopen(request, timeout=30) as response:
            text = response.read().decode()
        return json.loads(text) if text else None

    def comment(self, body):
        return self.call("POST", f"/repos/{self.repo}/issues/{self.issue}/comments", {"body": body})

    def comments(self):
        # every page, since on a public issue anyone can post a hundred comments
        # before the author pastes the code
        out, page = [], 1
        while True:
            got = self.call("GET", f"/repos/{self.repo}/issues/{self.issue}/comments?per_page=100&page={page}") or []
            out += got
            if len(got) < 100 or page >= 20:
                return out
            page += 1

    def delete_comment(self, comment_id):
        self.call("DELETE", f"/repos/{self.repo}/issues/comments/{comment_id}")

    def close(self, done):
        self.call("PATCH", f"/repos/{self.repo}/issues/{self.issue}",
                  {"state": "closed", "state_reason": "completed" if done else "not_planned"})


def code_from(comments, actor):
    """The code the person who started the run posted, with its comment's id."""
    for item in comments:
        login = (item.get("user") or {}).get("login", "")
        body = (item.get("body") or "").strip()
        if login.lower() != actor.lower():
            continue
        # the page posts the code under its mark; a person on the issue page pastes it bare
        code = body[len(CODE_MARK):].strip() if body.startswith(CODE_MARK) else body
        if CODE_SHAPE.match(code):
            return code, item["id"]
    return None, None


class Claude:
    """claude setup-token behind a pseudo-terminal, its output read through
    tokens.py's Masker, so a token is kept and never passed on."""

    def __init__(self, command=tokens.SETUP_TOKEN):
        import fcntl
        import pty
        import struct
        import termios
        self.pid, self.master = pty.fork()
        if self.pid == 0:
            for name in ("COLUMNS", "LINES") + WITHHELD:
                os.environ.pop(name, None)
            # GitHub's runners set no TERM, and claude then draws only the cells
            # that change, with cursor moves between, which no token pattern finds
            os.environ["TERM"] = os.environ.get("TERM") or TERMINAL_TYPE
            fcntl.ioctl(0, termios.TIOCSWINSZ, struct.pack("HHHH", *tokens.TERMINAL_SIZE, 0, 0))
            try:
                os.execvp(command[0], list(command))
            except OSError:
                os._exit(127)
        self.masker = tokens.Masker()
        self.decoder = codecs.getincrementaldecoder("utf-8")("replace")
        self.seen = ""
        self.ended = False

    def read(self, seconds):
        """What claude printed within the time, masked; nothing once it ended."""
        if self.ended:
            return ""
        ready, _, _ = select.select([self.master], [], [], seconds)
        if not ready:
            return ""
        try:
            data = os.read(self.master, 4096)
        except OSError:
            data = b""
        if not data:
            self.ended = True
            return self.masker.feed(self.decoder.decode(b"", final=True), final=True)
        shown = self.masker.feed(self.decoder.decode(data))
        self.seen = (self.seen + shown)[-20000:]
        return shown

    def type(self, text):
        os.write(self.master, text.encode())

    def submit(self, code, clock=time.monotonic):
        """The code, then Enter on its own, as a person pastes and presses it."""
        self.type(code)
        settled = clock() + PASTE_SETTLE
        while clock() < settled and not self.ended:
            self.read(0.1)
        self.type("\r")

    def shows(self, words):
        """Whether the screen holds the words, terminal codes and line
        breaks between them folded to single spaces."""
        plain = " ".join(re.sub(tokens.CODE, " ", self.seen).split())
        return " ".join(words.split()) in plain

    def link(self):
        found = SIGN_IN_LINK.findall(self.seen)
        return found[0] if found else None

    def finish(self):
        """claude's exit status, waiting a moment and then stopping it."""
        for _ in range(20):
            done, status = os.waitpid(self.pid, os.WNOHANG)
            if done:
                return os.waitstatus_to_exitcode(status)
            time.sleep(0.1)
        os.kill(self.pid, 15)
        _, status = os.waitpid(self.pid, 0)
        return os.waitstatus_to_exitcode(status)


def wait_for_link(claude, clock=time.monotonic):
    deadline = clock() + WAIT_FOR_LINK
    while clock() < deadline and not claude.ended:
        claude.read(1)
        link = claude.link()
        if link:
            return link
    raise RunError("claude setup-token printed no sign-in link")


def wait_for_prompt(claude, clock=time.monotonic):
    deadline = clock() + WAIT_FOR_PROMPT
    while clock() < deadline and not claude.ended and not claude.shows(PROMPT):
        claude.read(0.5)


def wait_for_code(github, claude, actor, clock=time.monotonic, sleep=time.sleep):
    deadline = clock() + WAIT_FOR_CODE
    while clock() < deadline:
        if claude.ended:
            raise RunError("claude setup-token ended before the code came")
        code, comment_id = code_from(github.comments(), actor)
        if code:
            return code, comment_id
        claude.read(POLL_EVERY)
    raise RunError("no code came within ten minutes; press Add a Claude account again")


def claude_said(claude, code):
    """The last words claude printed, as a line for the page: terminal codes
    and blanks folded, the code taken out, and only what the mask let
    through, so never a token."""
    plain = re.sub(tokens.CODE, " ", claude.seen)
    if code:
        plain = plain.replace(code, "<code>")
    words = " ".join(plain.split())
    return words[-300:]


def wait_for_token(claude, clock=time.monotonic):
    deadline = clock() + WAIT_FOR_TOKEN
    while clock() < deadline and not claude.ended:
        claude.read(1)
        if claude.masker.tokens:
            break
    while not claude.ended and claude.read(0.5):
        pass
    found = claude.masker.tokens
    if len(found) != 1:
        raise RunError("claude printed no token after the code: the code may be wrong or used already; press Add a Claude account again")
    return tokens.whole_token(found[0])


class Pool:
    """The organization's slots, written with the App's installation token,
    which may write the organization's secrets and variables."""

    def __init__(self, org, token, repo_id, call=None, targets=None, keeping=""):
        self.org, self.repo_id = org, repo_id
        # where a stored or lent token goes: this repository, or for a lend with
        # no say the organization's standing list, which may name none yet
        self.targets = [repo_id] if targets is None else list(targets)
        self.keeping = keeping
        self.base = f"/orgs/{org}/actions"
        self.api = call or GitHub("", token, "").call

    def read(self, path):
        try:
            return self.api("GET", path)
        except urllib.error.HTTPError as error:
            if error.code == 404:
                return None
            raise

    def taken(self):
        names, page = set(), 1
        while True:
            got = (self.api("GET", f"{self.base}/secrets?per_page=100&page={page}") or {}).get("secrets") or []
            names.update(item.get("name", "") for item in got)
            if len(got) < 100 or page >= 50:
                return names
            page += 1

    def free_slot(self):
        taken = self.taken()
        for slot in POOL_SLOTS:
            if SLOT_PREFIX + slot not in taken:
                return slot
        raise RunError("all 32 slots of the organization's pool hold a token; an owner frees one with its X")

    def put_secret(self, slot, token):
        key = self.api("GET", f"{self.base}/secrets/public-key")
        sealed = base64.b64encode(sealbox.seal(token.encode(), sealbox.key_from_base64(key["key"]))).decode()
        self.api("PUT", f"{self.base}/secrets/{SLOT_PREFIX}{slot}", {
            "encrypted_value": sealed, "key_id": key["key_id"],
            "visibility": "selected", "selected_repository_ids": self.targets})

    def write_variable(self, name, change, reach=True):
        # reach False changes the value alone, as a delete or a mark does, so
        # no repository is added to the variable's list by a run that lends nothing
        have = self.read(f"{self.base}/variables/{name}")
        value = change(have["value"] if have else "")
        if not have:
            self.api("POST", f"{self.base}/variables", {"name": name, "value": value,
                     "visibility": "selected", "selected_repository_ids": self.targets})
            return
        if value != have["value"]:
            self.api("PATCH", f"{self.base}/variables/{name}", {"name": name, "value": value})
        if reach and have.get("visibility") == "selected":
            for target in self.targets:
                self.api("PUT", f"{self.base}/variables/{name}/repositories/{target}")

    def add_to_accounts(self, slot):
        # an owner's page may write the list at the same moment, so it is read
        # back and written once more where that write took the slot away
        for _ in range(3):
            self.write_variable(ACCOUNTS_VARIABLE, lambda text: with_slot(text, slot))
            have = self.read(f"{self.base}/variables/{ACCOUNTS_VARIABLE}") or {}
            if slot in re.split(r"[,\s]+", have.get("value", "")):
                return
        raise RunError(f"the token is in slot {slot}, but {ACCOUNTS_VARIABLE} would not keep it; an owner adds it")

    def store(self, token, actor):
        """The token in the first free slot, the slot on the accounts list and
        the person named as its owner; the slot's number back."""
        step = "read which slots are free"
        try:
            if isinstance(self, Keeper) and self.own_slot(actor):
                raise RunError("you made a token in another run a moment ago; press Add Token again to lend that one here")
            slot = self.free_slot()
            step = f"store the token in slot {slot}"
            # asked once more just before the write, so a slot another run took
            # since the list is never written over
            if self.read(f"{self.base}/secrets/{SLOT_PREFIX}{slot}") is not None:
                raise RunError(f"slot {slot} was taken while the run looked; press Add Token again")
            # marked as the person's before the secret goes in, so no earlier
            # owner's mark is ever left on this token
            step = f"mark slot {slot} as yours in {SETTINGS_VARIABLE}"
            uid = self.read(f"/users/{actor}") or {}
            self.write_variable(SETTINGS_VARIABLE, lambda text: with_owner(text, slot, actor, uid.get("id"), self.keeping))
            step = f"store the token in slot {slot}"
            self.put_secret(slot, token)
            step = f"add slot {slot} to {ACCOUNTS_VARIABLE}"
            self.add_to_accounts(slot)
        except urllib.error.HTTPError as error:
            # the page waits for a word from the run, so GitHub's refusal becomes one
            raise RunError(f"GitHub answered {error.code} when the run tried to {step} in {self.org}'s pool"
                           + missing_permission(error, self.org)) from None
        return slot


class Keeper(Pool):
    """A person's own token in the pool: one slot whose settings name them as its
    owner, lent to the repositories they choose and taken back whenever they
    like (the owner's choice of 2026-10-19). Only the slot's owner gets here,
    since every call first finds the slot by the person's login."""

    def settings(self):
        have = self.read(f"{self.base}/variables/{SETTINGS_VARIABLE}") or {}
        try:
            data = json.loads(have.get("value") or "{}")
        except json.JSONDecodeError:
            return {}
        return data if isinstance(data, dict) else {}

    def own_slot(self, actor, trusted=True):
        """The slot whose settings name the actor. A slot stored before the
        account's id was kept matches by login alone, so a run anyone may start
        (trusted False) refuses it, since a login passes to whoever takes it
        after a rename; a run started by someone with write access, as the
        page's dispatch is, keeps the id then, and the slot is safe from there."""
        taken = self.taken()
        for slot, setting in sorted(self.settings().items(), key=lambda item: POOL_SLOTS.index(item[0]) if item[0] in POOL_SLOTS else 99):
            owner = setting.get("owner", "") if isinstance(setting, dict) else ""
            if slot in POOL_SLOTS and owner.lower() == actor.lower() and SLOT_PREFIX + slot in taken:
                # a login can pass to someone else once its account is renamed or
                # deleted, so where the slot kept the account's id that must match too
                kept = setting.get("uid")
                if kept is None:
                    if not trusted:
                        raise RunError(f"slot {slot} is yours by login alone, from before the pool kept account ids; press Add Token on the control panel once, which keeps yours, and then lend from an issue")
                    self.keep_uid(slot, actor)
                    return slot
                if str(kept) == self.user_id(actor):
                    return slot
        return ""

    def keep_uid(self, slot, actor):
        uid = (self.read(f"/users/{actor}") or {}).get("id")
        if not isinstance(uid, int):
            return
        def marked(text):
            data = json.loads(text) if (text or "").strip() else {}
            if isinstance(data, dict) and isinstance(data.get(slot), dict):
                data[slot]["uid"] = uid
                return json.dumps(data, separators=(",", ":"))
            return text
        self.write_variable(SETTINGS_VARIABLE, marked, reach=False)

    def user_id(self, login):
        return str((self.read(f"/users/{login}") or {}).get("id", ""))

    def secret_path(self, slot):
        return f"{self.base}/secrets/{SLOT_PREFIX}{slot}"

    def visibility(self, slot):
        return (self.read(self.secret_path(slot)) or {}).get("visibility", "")

    def selected(self, slot):
        seen = self.visibility(slot)
        if seen == "private":
            # private repositories only: a public one would never get the token
            raise RunError(f"slot {slot} is on in {self.org}'s private repositories, which an owner set; an owner switches it to selected repositories")
        return seen == "selected"

    def lent(self, slot):
        if not self.selected(slot):
            return ["all"]
        return self.lent_ids(slot)

    def lent_ids(self, slot):
        ids, page = [], 1
        while True:
            got = (self.api("GET", f"{self.secret_path(slot)}/repositories?per_page=100&page={page}") or {}).get("repositories") or []
            ids += [str(r.get("id")) for r in got if r.get("id")]
            if len(got) < 100 or page >= 50:
                return ids
            page += 1

    def reach_of_others(self, slot):
        """What the slots other than this one reach, read once a run: whether
        any is on in every repository or in the private ones, and the ids the
        others' lists name."""
        wide, ids = False, set()
        for name in self.taken():
            other = name[len(SLOT_PREFIX):] if name.startswith(SLOT_PREFIX) else ""
            if other not in POOL_SLOTS or other == slot:
                continue
            seen = self.visibility(other)
            if seen in ("all", "private"):
                wide = True
            elif seen == "selected":
                ids.update(self.lent_ids(other))
        return wide, ids

    def unlend_variables(self, repo_ids, slot):
        """The two variables taken off each repository no slot reaches any more,
        so its workflows read no lender's name once no token of theirs runs
        there; a slot on in every repository or the private ones may reach any."""
        if not repo_ids:
            return
        wide, others = self.reach_of_others(slot)
        gone = [i for i in repo_ids if str(i) not in others]
        if wide or not gone:
            return
        for name in (ACCOUNTS_VARIABLE, SETTINGS_VARIABLE):
            have = self.read(f"{self.base}/variables/{name}")
            if not have or have.get("visibility") != "selected":
                continue
            for repo_id in gone:
                try:
                    self.api("DELETE", f"{self.base}/variables/{name}/repositories/{repo_id}")
                except urllib.error.HTTPError as error:
                    if error.code != 404:
                        raise

    def lend_all(self, slot):
        for target in self.targets:
            self.lend(slot, target)

    def lend(self, slot, repo_id):
        if self.selected(slot):
            self.api("PUT", f"{self.secret_path(slot)}/repositories/{repo_id}")
        # the workers there read the slot's name and settings from these two
        for name in (ACCOUNTS_VARIABLE, SETTINGS_VARIABLE):
            have = self.read(f"{self.base}/variables/{name}")
            if have and have.get("visibility") == "private":
                raise RunError(f"{name} is on in {self.org}'s private repositories only; an owner switches it to selected repositories or all")
            if have and have.get("visibility") == "selected":
                self.api("PUT", f"{self.base}/variables/{name}/repositories/{repo_id}")

    def revoke(self, slot, repo_id):
        if not self.selected(slot):
            raise RunError(f"slot {slot} is on in every repository of {self.org}, which an owner set; an owner narrows it")
        try:
            self.api("DELETE", f"{self.secret_path(slot)}/repositories/{repo_id}")
        except urllib.error.HTTPError as error:
            if error.code != 404:
                raise
        self.unlend_variables([repo_id], slot)

    def kept_by_organization(self, slot):
        setting = self.settings().get(slot)
        return isinstance(setting, dict) and setting.get("keeping") == ORGANIZATION

    def mark_keeping(self, slot):
        """A lend with no say hands the slot to the organization for good: its
        lender is told no repository from then on, and only delete stays theirs."""
        if not self.keeping:
            return
        def marked(text):
            data = json.loads(text) if (text or "").strip() else {}
            if isinstance(data, dict) and isinstance(data.get(slot), dict):
                data[slot]["keeping"] = self.keeping
                return json.dumps(data, separators=(",", ":"))
            return text
        self.write_variable(SETTINGS_VARIABLE, marked)

    def delete(self, slot):
        # the secret goes first and the owner's mark after it: a run cut short
        # between the two leaves a mark on an empty slot, which the next store
        # into that slot writes over and which own_slot never returns, since it
        # needs the secret too; the other order would leave a token nobody owns,
        # which its lender could no longer delete. A store may take the freed
        # slot between the two, so the mark goes only while it is still this
        # token's (its id), and the slot leaves the accounts list only while no
        # secret is in it again. The repositories it reached are read first, to
        # take the variables off the ones no slot reaches now
        mine = (self.settings().get(slot) or {}).get("id")
        reached = self.lent_ids(slot) if self.visibility(slot) == "selected" else []
        self.api("DELETE", self.secret_path(slot))
        self.write_variable(SETTINGS_VARIABLE, lambda text: without(text, slot, mine), reach=False)
        if self.read(self.secret_path(slot)) is None:
            self.write_variable(ACCOUNTS_VARIABLE, lambda text: ",".join(w for w in re.split(r"[,\s]+", text or "") if w and w != slot), reach=False)
        # the deleted secret is gone from the list, and a store that refilled
        # the slot since reaches what it reaches, so no slot is skipped here
        self.unlend_variables(reached, None)


def without(text, slot, mine=None):
    try:
        data = json.loads(text) if (text or "").strip() else {}
    except json.JSONDecodeError:
        return text
    # a mark another token's store wrote since is that token's, and stays
    # (a slot made before ids has none, so a mark with an id is a new store's)
    theirs = isinstance(data, dict) and isinstance(data.get(slot), dict) and data[slot].get("id") != mine
    if isinstance(data, dict) and not theirs:
        data.pop(slot, None)
        return json.dumps(data, separators=(",", ":"))
    return text


def with_slot(text, slot):
    words = [w for w in re.split(r"[,\s]+", text or "") if w]
    return ",".join(words if slot in words else words + [slot])


def with_owner(text, slot, actor, uid=None, keeping=""):
    # a variable that is not an object by slot is left for an owner to mend,
    # since writing over it would take every other slot's settings away
    try:
        data = json.loads(text) if (text or "").strip() else {}
    except json.JSONDecodeError:
        data = None
    if not isinstance(data, dict):
        raise RunError(f"the token is in slot {slot}, but {SETTINGS_VARIABLE} is not JSON by slot, so it is not marked as yours yet; an owner mends it")
    data[slot] = {"owner": actor, "id": os.urandom(4).hex()}
    if isinstance(uid, int):
        data[slot]["uid"] = uid
    if keeping:
        data[slot]["keeping"] = keeping
    return json.dumps(data, separators=(",", ":"))


def settings(environ, sealing=None):
    issue, key, key_id = environ.get("ISSUE", ""), environ.get("KEY", ""), environ.get("KEY_ID", "")
    actor, repo, token = environ.get("ACTOR", ""), environ.get("GITHUB_REPOSITORY", ""), environ.get("GITHUB_TOKEN", "")
    if not ISSUE_SHAPE.match(issue):
        raise RunError("the run was not given the page's issue")
    if not actor or not repo or not token:
        raise RunError("the run lacks its actor, repository or token")
    if sealing is None:
        sealing = not environ.get("STORE_ORG")
    if not sealing:
        return GitHub(repo, token, issue), None, "", actor
    if not KEY_ID_SHAPE.match(key_id):
        raise RunError("the run was not given the secrets key's id")
    return GitHub(repo, token, issue), sealbox.key_from_base64(key), key_id, actor


def pool_of(environ):
    org, app_token, repo_id = environ.get("STORE_ORG", ""), environ.get("APP_TOKEN", ""), environ.get("REPO_ID", "")
    if not ORG_SHAPE.match(org) or not app_token or not re.fullmatch(r"[1-9][0-9]{0,15}", repo_id):
        raise RunError("the run was not given the organization, the App's token or the repository's id to store the token with")
    return Pool(org, app_token, int(repo_id))


def request(environ):
    """What the issue that started this run asks: its kind (lend, manage or
    nothing) and the words of its second line. A run a dispatch started asks
    nothing here; its inputs say."""
    path = environ.get("GITHUB_EVENT_PATH", "")
    if environ.get("GITHUB_EVENT_NAME") != "issues" or not path:
        return "", []
    try:
        with open(path) as event:
            body = str((json.load(event).get("issue") or {}).get("body") or "")
    except (OSError, ValueError, AttributeError):
        return "", []
    lines = body.replace("\r\n", "\n").split("\n")
    kinds = {LEND_MARK: "lend", MANAGE_MARK: "manage"}
    kind = kinds.get(lines[0].strip(), "")
    return kind, (lines[1].split()[:2] if kind and len(lines) > 1 else [])


def started_by_issue(environ):
    return environ.get("GITHUB_EVENT_NAME") == "issues"


def no_say(environ):
    """A lend with no say: the issue's second line says organization."""
    kind, words = request(environ)
    return kind == "lend" and words[:1] == [ORGANIZATION]


def standing_list(pool):
    have = pool.read(f"{pool.base}/variables/{LEND_TO_VARIABLE}") or {}
    ids = [w for w in re.split(r"[,\s]+", have.get("value", "") or "") if w]
    bad = [w for w in ids if not re.fullmatch(r"[1-9][0-9]{0,15}", w)]
    if bad:
        raise RunError(f"{LEND_TO_VARIABLE} holds something that is not a repository id; an owner mends it")
    return [int(w) for w in ids]


def asked(environ):
    """The action and the repository id a manage run was asked for: the
    dispatch's inputs, or the issue's second line."""
    if environ.get("MANAGE"):
        return environ["MANAGE"], environ.get("TARGET_REPO_ID", "")
    kind, words = request(environ)
    if kind != "manage":
        return "", ""
    return (words + ["", ""])[0], (words + ["", ""])[1]


def finish_issue(github, environ, done):
    """A run an issue started closes it, since no page will; a dispatch's issue
    is the page's to close."""
    if not started_by_issue(environ):
        return
    try:
        github.close(done)
    except Exception:
        pass


def held_file(environ):
    path = environ.get("TOKEN_FILE", "")
    if not path:
        raise RunError("the run was not given the file that holds the token between its steps")
    return path


def hold(token, path):
    """The token on the runner's disk between the sign-in and the store, readable
    by the runner's user alone, since the App's token is minted only once claude
    has gone and no process of claude's can read it."""
    descriptor = os.open(path, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o600)
    with os.fdopen(descriptor, "w") as held:
        held.write(token)


def take_held(path):
    try:
        with open(path) as held:
            return held.read().strip()
    except OSError:
        raise RunError("the sign-in left no token for the store step") from None
    finally:
        if os.path.exists(path):
            os.remove(path)


def run(environ, start=Claude):
    github, key, key_id, actor = settings(environ)
    holding = held_file(environ) if environ.get("STORE_ORG") else ""
    claude = start()
    try:
        github.comment(f"{LINK_MARK}\n{wait_for_link(claude)}")
        print("Posted the sign-in link; waiting for the code.", flush=True)
        code, comment_id = wait_for_code(github, claude, actor)
        print(f"::add-mask::{code}", flush=True)
        github.delete_comment(comment_id)
        wait_for_prompt(claude)
        claude.submit(code)
        try:
            token = wait_for_token(claude)
        except RunError as error:
            said = claude_said(claude, code)
            print(f"claude's last words: {said}", flush=True)
            raise RunError(f"{error}. claude said: {said}") from None
        if holding:
            hold(token, holding)
            print("Held the token for the store step.", flush=True)
            return
        sealed = base64.b64encode(sealbox.seal(token.encode(), key)).decode()
        github.comment(f"{SEALED_MARK}\n{key_id}\n{sealed}")
        print("Posted the sealed token.", flush=True)
    finally:
        claude.finish()


def store(environ, pool=None):
    """The token the sign-in held, into the organization's pool with the App's token."""
    github, _, _, actor = settings(environ)
    pool = pool or keeper_of(environ)
    token = take_held(held_file(environ))
    slot = pool.store(token, actor)
    github.comment(f"{STORED_MARK}\n{slot}")
    finish_issue(github, environ, True)
    print(f"Stored the token in the organization's slot {slot}.", flush=True)


def fail(environ):
    """A word on the issue for a run that stopped before the script could say why,
    as when the App's token is refused or the job is skipped, so the page does not
    wait its fifteen minutes; none when the run said so already."""
    # the key is no matter here, and the step that failed may have been given none
    github, *_ = settings(environ, sealing=False)
    said = [c.get("body", "") for c in github.comments() if (c.get("user") or {}).get("login") == "github-actions[bot]"]
    if not any(body.startswith(FAILED_MARK) for body in said):
        github.comment(f"{FAILED_MARK}\n{environ.get('FAIL_WHY') or 'the run stopped before the token was saved; see its log under Actions, workflow Token'}")
    finish_issue(github, environ, False)


def lend_target(environ, pool, actor):
    """The repository a dispatch asks the token lent to in place of this one,
    the input lend_to: another repository of the same organization, one with no
    workflow Token of its own, as the Tokens page asks when the repository it
    shows has none and runs the workflow in one of the organization's that has
    it. A run an issue started takes none, since its author may hold no role,
    and the id is checked against the organization and the lender's right to
    write there, so a mistyped one lends nowhere."""
    wanted = (environ.get("LEND_TO") or "").strip()
    if not wanted or started_by_issue(environ):
        return None
    if not re.fullmatch(r"[1-9][0-9]{0,15}", wanted):
        raise RunError("lend_to is not a repository id; nothing was lent")
    found = pool.read(f"/repositories/{wanted}")
    owner = str(((found or {}).get("owner") or {}).get("login") or "")
    if not found or owner.lower() != pool.org.lower():
        raise RunError(f"lend_to names no repository of {pool.org} that the App reaches; nothing was lent")
    # the lender may write there, as a dispatch on that repository would have
    # needed, so nobody plants a token in a repository that is not theirs
    had = pool.read(f"/repositories/{wanted}/collaborators/{actor}/permission") or {}
    if had.get("permission") not in ("admin", "write") and had.get("role_name") not in ("admin", "maintain", "write"):
        raise RunError(f"lend_to names a repository {actor} may not write to; nothing was lent")
    return int(wanted)


def keeper_of(environ):
    pool = pool_of(environ)
    keeper = Keeper(pool.org, environ.get("APP_TOKEN", ""), pool.repo_id)
    if no_say(environ):
        # lent where the organization's standing list says, nowhere yet when
        # it names none, and kept as the organization's to place
        keeper.targets = standing_list(keeper)
        keeper.keeping = ORGANIZATION
    try:
        target = lend_target(environ, keeper, environ.get("ACTOR", ""))
    except urllib.error.HTTPError as error:
        raise RunError(f"GitHub answered {error.code} when the run read the repository lend_to names; nothing was lent") from None
    if target is not None:
        keeper.targets = [target]
    return keeper


def revoke_app_token(environ, call=None):
    """The App's token given up at once, so a step that comes after, as claude's
    sign-in does, finds only a dead one where the action kept it."""
    try:
        (call or GitHub("", environ.get("APP_TOKEN", ""), "").call)("DELETE", "/installation/token")
    except Exception:
        # claude must never start beside a live App token, so the step fails
        # and the sign-in after it is skipped
        raise RunError("the run could not give up the App's token before the sign-in, so it stopped; press Add Token again") from None


def say_lent(github, keeper, slot):
    # a token in the organization's keeping names no repository to its lender,
    # who may not be allowed to see them
    where = ORGANIZATION if keeper.kept_by_organization(slot) else ",".join(keeper.lent(slot))
    github.comment(f"{LENT_MARK}\n{slot}\n{where}")


def look(environ, keeper=None, call=None):
    """Before any sign-in: a person who already has a token in the pool has it
    lent to this repository, and no new one is made. Says found=true to the
    workflow so the sign-in steps are skipped."""
    try:
        github, _, _, actor = settings(environ, sealing=False)
        keeper = keeper or keeper_of(environ)
        slot = keeper.own_slot(actor, trusted=not started_by_issue(environ))
        if slot:
            # a token the organization keeps is placed by its owners alone: a
            # lend of it to a repository its lender picks lends nothing more,
            # while a lend with no say places it where the standing list says
            if keeper.keeping or not keeper.kept_by_organization(slot):
                keeper.lend_all(slot)
            keeper.mark_keeping(slot)
            say_lent(github, keeper, slot)
            finish_issue(github, environ, True)
            print(f"Lent slot {slot}, already the person's own, where this run was asked to.", flush=True)
    except urllib.error.HTTPError as error:
        org = environ.get("STORE_ORG", "")
        raise RunError(f"GitHub answered {error.code} when the run looked for your token in {org}'s pool"
                       + missing_permission(error, org)) from None
    finally:
        revoke_app_token(environ, call)
    output = environ.get("GITHUB_OUTPUT")
    if output:
        with open(output, "a") as found:
            found.write(f"found={'true' if slot else 'false'}\n")
    return slot


MANAGE_ACTIONS = ("list", "revoke", "delete")


def manage(environ, keeper=None):
    """A person's own token: where it is lent, taken back from one repository,
    or deleted. A person without a token in the pool is told so."""
    github, _, _, actor = settings(environ, sealing=False)
    action, target = asked(environ)
    if action not in MANAGE_ACTIONS:
        raise RunError(f"the run was asked to {action[:20]!r} a token; it lists, revokes or deletes one")
    if action == "revoke" and not re.fullmatch(r"[1-9][0-9]{0,15}", target):
        raise RunError("the run was not given the repository to take the token back from")
    keeper = keeper or keeper_of(environ)
    try:
        slot = keeper.own_slot(actor, trusted=not started_by_issue(environ))
        if not slot:
            raise RunError(f"you have no token in {keeper.org}'s pool; Add Token makes one")
        if action == "revoke":
            if keeper.kept_by_organization(slot):
                raise RunError(f"your token is in {keeper.org}'s keeping, which you chose with a lend with no say, so its owners place it; you may still delete it")
            keeper.revoke(slot, int(target))
        if action == "delete":
            keeper.delete(slot)
            github.comment(f"{LENT_MARK}\n{slot}\n")
        else:
            say_lent(github, keeper, slot)
        finish_issue(github, environ, True)
    except urllib.error.HTTPError as error:
        raise RunError(f"GitHub answered {error.code} when the run tried to {action} your token in {keeper.org}'s pool"
                       + missing_permission(error, keeper.org)) from None


STAGES = {"sign-in": run, "store": store, "fail": fail, "look": look, "manage": manage}


def main(argv=None):
    stage = (argv if argv is not None else sys.argv[1:] or ["sign-in"])[0]
    if stage not in STAGES:
        print(f"token_run: no stage {stage!r}; the stages are {', '.join(STAGES)}", file=sys.stderr)
        return 2
    try:
        STAGES[stage](os.environ)
    except Exception as error:
        # every failure reaches the page as a word, but only a RunError's own
        # words, which never hold a token; anything else by its kind alone
        known = isinstance(error, (RunError, sealbox.SealError, tokens.TokensError))
        why = str(error) if known else f"the run failed ({type(error).__name__}); see its log under Actions, workflow Token"
        print(f"token_run: {why}", file=sys.stderr)
        if stage != "fail":
            try:
                github, *_ = settings(os.environ, sealing=False)
                github.comment(f"{FAILED_MARK}\n{why}")
                finish_issue(github, os.environ, False)
            except Exception:
                pass
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
