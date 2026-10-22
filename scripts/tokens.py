#!/usr/bin/env python3
"""Store each of your Claude accounts' subscription tokens as an Actions secret, without ever seeing one.

The Claude runs this repository starts in GitHub Actions (@claude, the
workers in Actions and the coordinator workflow, docs/setup/claude-action.md)
run on the secret CLAUDE_CODE_OAUTH_TOKEN, and a benchmark run in Actions
(OpenSpineConsortium/coordbench's bench.yml) spreads its agents over a pool
of accounts, one secret SYMPHONY_CLAUDE_<slot> each, slots 1 to 8, with the
slots in use listed in the variable SYMPHONY_CLAUDE_ACCOUNTS. Each token
comes from `claude setup-token`, a browser sign-in that prints a token good
for about a year; by hand, each then goes through a clipboard into a
settings page, once per account. This wizard does that part. It runs
`claude setup-token` for you behind a pseudo-terminal, lets the sign-in
through to your browser and your keyboard, takes the token out of what
claude prints before it reaches the screen, and hands it to `gh secret set`
over a pipe, which encrypts it on your machine and sends only the sealed
value. The token is never on the screen, on a command line, in a file, in
a log or in this script's output; what is printed is names, dates and
counts. An API key is refused on sight: runs bill nothing (docs/project.md),
and only a subscription token has the shape claude setup-token prints.

  python3 scripts/tokens.py                   the wizard: the state, then one account at a time
  python3 scripts/tokens.py status            the state alone: which slots are set, since when, renew by when
  python3 scripts/tokens.py add [--slot N] [--label WORD] [--also-main] [--paste]
  python3 scripts/tokens.py main [--paste]    CLAUDE_CODE_OAUTH_TOKEN alone, from one more sign-in
  python3 scripts/tokens.py remove N          a slot's secret and its label, off the list of slots
  python3 scripts/tokens.py order 3,1,2       the order the bench deals its jobs to the slots
  python3 scripts/tokens.py profiles          what this machine remembers for each slot
  python3 scripts/tokens.py forget N          a slot's browser profile and cswap account, forgotten here
  python3 scripts/tokens.py show N|main [--copy]  a token back from this machine's copy, on Windows
  python3 scripts/tokens.py share coordbench    the slots seen by more repositories, with no new sign-in
  python3 scripts/tokens.py usage             cswap's usage of each slot's account, for the control panel

Where the secrets go: --org ORG puts them on the organization, seen by the
repositories --repos names (gh's "selected" visibility; --visibility all or
private makes every repository see them); --repo OWNER/NAME puts them on one
repository. With neither, the clone's origin decides: a repository of an
organization takes the organization level with that repository on the list,
a personal one its own level. A secret that exists keeps its visibility and
its list, and the list grows by what --repos names; nothing is taken off a
list here. Needs the GitHub CLI signed in (gh auth login) as an admin of the
organization, or of the repository, and Claude Code (claude) on this
machine; `add --paste` and `main --paste` take the token from a hidden
prompt instead, for a machine without claude. A label is
a word of your own so that status can say which account a slot is; it is
stored as the variable SYMPHONY_CLAUDE_<slot>_LABEL, which anyone who reads
the repository's workflows can see, so never an e-mail address. Refuses to
run for a repository in crowd mode, which holds no Claude token
(docs/setup/crowdsourcing.md). Everything that can fail is tried before
the sign-in, since a token minted and not stored is lost with the process.

Where this machine has Chrome, Edge or Brave, the wizard opens each slot's
sign-in in a browser profile of the slot's own (add --profile NAME does
the same), so a profile signed in to claude.ai as that account once makes
each later token one Approve; where claude-swap (cswap) keeps a list of
accounts, the wizard shows it and takes an account by its number, its
alias offered as the label. It remembers the profile and the cswap number per
slot in tokens.json in this machine's own settings folder, never in the
repository; cswap's credentials are never opened, since a token minted for
Actions is a separate, inference-only one. A Windows shell has no
pseudo-terminal, so there claude keeps the console's keyboard and its
output comes to the wizard on a pipe, through the same mask; the profile
is opened first, on a blank page.

On Windows each token stored is also kept in Windows Credential Manager,
under symphony/<place>/<secret>, encrypted under your Windows sign-in, so
that `show N` gives it back on demand, printed or with --copy on the
clipboard; GitHub never gives a secret back. remove deletes the copy with
the secret, --no-copy keeps none, and elsewhere no copy is kept.
Exit 0 when done, 1 on a refusal, 2 on a usage error made before any call.
"""

import argparse
import codecs
import getpass
import json
import os
import re
import select
import shutil
import struct
import subprocess
import sys
import tempfile
import time
from contextlib import contextmanager
from datetime import datetime, timedelta, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
PROJECT = ROOT / "docs" / "project.md"
CROWD_ON = re.compile(r"^- Crowd mode: on\.", re.MULTILINE)

GH_TIMEOUT = 600
MAIN_SECRET = "CLAUDE_CODE_OAUTH_TOKEN"
SLOT_PREFIX = "SYMPHONY_CLAUDE_"
LABEL_SUFFIX = "_LABEL"
ACCOUNTS_VARIABLE = SLOT_PREFIX + "ACCOUNTS"
SLOTS = tuple(range(1, 9))
# A list of slots is read as coordbench's bench/cooperbench/actions/plan.py
# reads SYMPHONY_CLAUDE_ACCOUNTS: split on commas and blanks, with a real run
# refused on any word that is not one of these.
SLOT_WORDS = {str(slot): slot for slot in SLOTS}
LIST_SEPARATORS = re.compile(r"[,\s]+")
CROWD_VARIABLE = "CLAUDE_CROWD"
TOKEN_LIFE = timedelta(days=365)
RENEW_SOON = timedelta(days=30)
SETUP_TOKEN = ("claude", "setup-token")
VISIBILITIES = ("selected", "all", "private")
# Rows and columns of the terminal claude gets: wide, so it never wraps the
# token over two lines, where no pattern would find it.
TERMINAL_SIZE = (50, 400)
# A held-back tail longer than this is a run of letters, not a token's start;
# a captured run longer than it went on past where a token ends.
HOLD_LIMIT = 400
# A token claude prints is far longer than this; a captured one shorter is
# a piece of one, cut by a terminal narrower than the one claude was given.
TOKEN_SHORTEST = 60
# The shapes are assembled from pieces so that this file never looks like a
# secret to tools/check_docs.py, whose own shapes are built the same way.
# Anything of Anthropic's key shape is kept off the screen; only the
# subscription token's kind is kept.
SECRET_HEAD = "sk-" + "ant-"
TOKEN_HEAD = SECRET_HEAD + "oat"
API_KEY_HEAD = SECRET_HEAD + "api"
SECRET_SHAPE = re.compile(re.escape(SECRET_HEAD) + r"[A-Za-z0-9_-]{20,}")
TOKEN_SHAPE = re.compile(re.escape(TOKEN_HEAD) + r"[0-9]{2}-[A-Za-z0-9_-]{20,}")
# A tail of what claude printed that may be the start of a token whose rest
# is still to come: two characters or more of the head, or the head itself
# with whatever followed it so far.
TOKEN_START = re.compile(
    "(?:" + "|".join(re.escape(SECRET_HEAD[:n]) for n in range(2, len(SECRET_HEAD) + 1))
    + "|" + re.escape(SECRET_HEAD) + r"[A-Za-z0-9_-]*)$"
)
# What the screen shows in a token's place; neither matches a secret's shape.
MASK = TOKEN_HEAD + "**-" + "*" * 16 + " (kept by the wizard, not shown)"
OTHER_MASK = SECRET_HEAD + "****" + " (not a token from claude setup-token; not shown)"
# The rest of a token cut over two lines stands alone on the line after
# the masked part: a run of a token's characters, past line ends, blanks
# and terminal codes, that ends with its line, where a word of claude's is
# followed by a blank or by the code that places the next word. A code cut
# by the chunk's end is waited for.
CODE = r"(?:\x1b\[[0-9;?<>=]*[@-~]|\x1b\][^\x07\x1b]*(?:\x07|\x1b\\)|\x1b[^\[\]])"
CONTINUATION = re.compile(r"^(?:[ \t\r\n]|" + CODE + r")*([A-Za-z0-9_-]*)(" + CODE + r"*)")
PARTIAL_CODE = re.compile(r"\x1b(?:\[[0-9;?<>=]*|\][^\x07\x1b]*\x1b?)?")
# What claude turns on in the terminal and may leave on when it is killed.
TERMINAL_RESET = "\x1b[?2004l\x1b[?1004l\x1b[?1000l\x1b[?1002l\x1b[?1003l\x1b[?1006l\x1b[?2026l\x1b[<u\x1b[?25h\x1b[0m"
LABEL_SHAPE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9 ._-]{0,39}$")
REPO_NAME = re.compile(r"^[A-Za-z0-9._-]+$")
FULL_NAME = re.compile(r"^[A-Za-z0-9][A-Za-z0-9-]*/[A-Za-z0-9._-]+$")
ORIGIN = re.compile(r"github\.com[:/]([A-Za-z0-9][A-Za-z0-9-]*)/([A-Za-z0-9._-]+?)(?:\.git)?/?$")


class TokensError(Exception):
    """Something the wizard could not do, or refuses to do; the message says
    what to do instead, and never holds a token."""


class UsageError(TokensError):
    """Options that cannot go together, found before gh is asked anything."""


# The GitHub CLI.

class Gh:
    """gh's own commands, each run as a process. A value for a secret goes in
    on the pipe to gh's standard input and nowhere else; gh never prints it,
    so neither does anything here."""

    def run(self, args, stdin=None):
        try:
            return subprocess.run(
                ["gh"] + list(args), input=stdin, capture_output=True,
                text=True, encoding="utf-8", timeout=GH_TIMEOUT,
            )
        except FileNotFoundError:
            raise TokensError("gh is not installed: https://cli.github.com, then run gh auth login") from None
        except subprocess.TimeoutExpired:
            raise TokensError(f"gh {args[0]} gave no answer in {GH_TIMEOUT} seconds") from None

    def must(self, args, stdin=None):
        """gh's output, or the error gh named."""
        done = self.run(args, stdin)
        if done.returncode != 0:
            raise TokensError(f"gh {' '.join(args[:2])}: {first_line(done.stderr) or 'failed'}")
        return done.stdout

    def json(self, args):
        text = self.must(args)
        return json.loads(text) if text.strip() else []

    def read(self, args):
        """A read that may find nothing: (found, output). Any failure but a
        404 is an error, since a guess would store a secret in the wrong
        place or in a crowd repository."""
        done = self.run(args)
        if done.returncode == 0:
            return True, done.stdout
        if "(HTTP 404)" in done.stderr:
            return False, ""
        raise TokensError(f"gh {' '.join(args[:2])}: {first_line(done.stderr) or 'failed'}")


def first_line(text):
    return (text or "").strip().splitlines()[0].strip() if (text or "").strip() else ""


# Where the secrets live.

class Target:
    """An organization, whose secrets a list of repositories sees (or every
    repository, with visibility all or private), or one repository."""

    def __init__(self, level, owner, repos=(), visibility="selected", name=None):
        self.level = level
        self.owner = owner
        self.repos = sorted(set(repos))
        self.visibility = visibility
        self.name = name

    @property
    def flags(self):
        if self.level == "org":
            return ["--org", self.owner]
        return ["--repo", self.full_name]

    @property
    def full_name(self):
        return f"{self.owner}/{self.name}"

    @property
    def api(self):
        if self.level == "org":
            return f"orgs/{self.owner}/actions"
        return f"repos/{self.full_name}/actions"

    def repositories(self):
        """The repositories whose crowd mode is checked: the ones named, or
        the one the secrets go on."""
        if self.level == "org":
            return [f"{self.owner}/{repo}" for repo in self.repos]
        return [self.full_name]

    def describe(self):
        if self.level == "repo":
            return f"the repository {self.full_name}"
        if self.visibility == "selected":
            return f"the organization {self.owner}, seen by " + (", ".join(self.repos) or "no repository yet")
        return f"the organization {self.owner}, seen by {'every' if self.visibility == 'all' else 'every private'} repository"


def origin_repository(root=ROOT):
    """OWNER/NAME from the clone's origin, or None when there is none."""
    try:
        done = subprocess.run(["git", "-C", str(root), "remote", "get-url", "origin"],
                              capture_output=True, text=True, encoding="utf-8")
    except OSError:
        return None
    match = ORIGIN.search(done.stdout.strip()) if done.returncode == 0 else None
    return (match.group(1), match.group(2)) if match else None


def check_options(args, root=ROOT):
    """The place the options name, as far as it is known without gh:
    (repos, visibility, origin), or why the options cannot go together."""
    repos = parse_repos(args.repos)
    visibility = args.visibility or "selected"
    if args.repo and (repos or args.visibility):
        raise UsageError("--repo puts the secrets on one repository: --repos and --visibility are for --org")
    if repos and visibility != "selected":
        raise UsageError(f"--repos names the repositories that see a selected secret, and --visibility {visibility} "
                         "makes every repository see it: drop one")
    origin = None if args.repo else origin_repository(root)
    if args.org and not repos and visibility == "selected":
        if origin and origin[0].lower() == args.org.lower():
            repos = [origin[1]]
        else:
            raise UsageError("--org needs the repositories that see the secrets: --repos NAME,NAME (or --visibility all)")
    if not args.org and not args.repo and not origin:
        raise UsageError("this clone has no GitHub origin: say where the secrets go with --org ORG --repos NAME or --repo OWNER/NAME")
    return repos, visibility, origin


def target_from(args, gh, root=ROOT):
    """The level and the place the options name, or the origin's."""
    repos, visibility, origin = check_options(args, root)
    if args.repo:
        owner, name = args.repo.split("/")
        return Target("repo", owner, name=name)
    if args.org:
        return Target("org", args.org, repos, visibility)
    owner, name = origin
    kind = gh.must(["api", f"users/{owner}", "--jq", ".type"]).strip()
    if kind == "Organization":
        return Target("org", owner, repos or [name], visibility)
    if repos or args.visibility:
        raise TokensError(f"{owner} is a personal account, whose repository {name} takes its own secrets: --repos and --visibility are for an organization")
    return Target("repo", owner, name=name)


def parse_repos(text):
    repos = [part.strip() for part in (text or "").replace(" ", ",").split(",") if part.strip()]
    for repo in repos:
        if not REPO_NAME.match(repo):
            raise UsageError(f"--repos takes repository names without their owner, not {repo}")
    return repos


# The token, from claude or from a hidden prompt.

class Masker:
    """What claude prints, with every token taken out before it reaches the
    screen and kept here. A token can arrive in pieces, so a tail that may
    be a token's start is held back until the rest comes, and the last
    characters on the screen are kept for matching, since a lone "s" is let
    through and the token it began is then masked from its second
    character on."""

    def __init__(self):
        self.shown = ""
        self.held = ""
        self.tokens = []
        # What followed the last mask, while it may be a cut token's rest.
        self.after = None
        self.cut_over_lines = False

    def feed(self, chunk, final=False):
        """The text to put on the screen for this chunk. A token that runs
        to the chunk's end may go on in the next one, so it is held back
        whole, as a start is, until the text after it shows where it ends;
        one that has run on past any token's length is masked as it is,
        since it is never shown, and refused later by its length."""
        text = self.shown + self.held + chunk
        shown = len(self.shown)
        pieces, position = [], shown
        for match in SECRET_SHAPE.finditer(text):
            if match.end() == len(text) and not final and len(match.group(0)) <= HOLD_LIMIT:
                break
            found = match.group(0)
            if TOKEN_SHAPE.fullmatch(found) and found not in self.tokens:
                self.tokens.append(found)
            pieces.append((text[position:max(match.start(), shown)], "shown"))
            # A run past the hold limit is refused by its length later, so
            # what follows it is not watched for a cut token's rest.
            kind = "mask" if len(found) <= HOLD_LIMIT else "mask of a run"
            pieces.append((MASK if TOKEN_SHAPE.fullmatch(found) else OTHER_MASK, kind))
            position = max(match.end(), shown)
        cut = len(text)
        if not final:
            start = TOKEN_START.search(text)
            if start and len(start.group(0)) <= HOLD_LIMIT:
                cut = max(start.start(), position)
        self.held = text[cut:]
        pieces.append((text[position:cut], "shown"))
        self.watch(pieces, final)
        emit = "".join(piece for piece, _ in pieces)
        if emit:
            self.shown = (self.shown + emit)[-12:]
        return emit

    def watch(self, pieces, final):
        """A run of a token's characters that stands alone on the line
        after a masked token is the rest of one cut over two lines by a
        terminal narrower than the one claude was given: the piece kept
        is then no token, and the rest is already on the screen, so it is
        noted here and the mint refused. A word of claude's is followed by
        a blank or by the code that places the next word, never by its
        line's end."""
        for piece, kind in pieces:
            if kind == "mask":
                self.after = ""
            elif kind == "mask of a run":
                self.after = None
            elif self.after is not None and piece:
                self.after += piece
                self.decide(final)
        if final and self.after:
            self.decide(final)

    def decide(self, final):
        run = CONTINUATION.match(self.after)
        word, tail = run.group(1), self.after[run.end():]
        if not word:
            if tail and not PARTIAL_CODE.fullmatch(tail):
                self.after = None
        elif tail[:1] in ("\r", "\n") or (final and not tail):
            self.cut_over_lines = True
            self.after = None
        elif tail and not PARTIAL_CODE.fullmatch(tail):
            self.after = None
        if self.after is not None and len(self.after) > HOLD_LIMIT:
            self.after = None


def pty_available():
    """A pseudo-terminal exists on Linux, macOS and WSL, and not in a
    Windows shell, which takes the pasted route."""
    try:
        import pty  # noqa: F401
        import termios  # noqa: F401
    except ImportError:
        return False
    return True


def pump(master, stdin_fd, on_output):
    """Keystrokes to claude, claude's output to on_output, until claude's
    side of the terminal closes."""
    stdin_open = True
    while True:
        wanted = [master] + ([stdin_fd] if stdin_open else [])
        ready, _, _ = select.select(wanted, [], [])
        if master in ready:
            try:
                data = os.read(master, 4096)
            except OSError:
                data = b""
            if not data:
                return
            on_output(data)
        if stdin_fd in ready:
            data = os.read(stdin_fd, 4096)
            if data:
                os.write(master, data)
            else:
                stdin_open = False


def mint(stdin_fd=None, write=None, command=SETUP_TOKEN, env=None):
    """Run claude setup-token behind a pseudo-terminal: the sign-in goes
    through to the browser and the keyboard, everything claude prints
    reaches the screen but the token, and the tokens seen come back with
    claude's exit status and whether one was cut over two lines. env is
    added to claude's environment: BROWSER, which claude opens the sign-in
    link with."""
    import fcntl
    import pty
    import termios
    import tty
    if stdin_fd is None:
        stdin_fd = sys.stdin.fileno()
    if write is None:
        def write(text):
            sys.stdout.write(text)
            sys.stdout.flush()
    pid, master = pty.fork()
    if pid == 0:
        # The width claude wraps at is the terminal's, set here before it
        # starts; a COLUMNS or LINES inherited from the shell would tell a
        # library another width, and a token wrapped over two lines is
        # one no pattern finds whole.
        for name in ("COLUMNS", "LINES"):
            os.environ.pop(name, None)
        os.environ.update(env or {})
        fcntl.ioctl(0, termios.TIOCSWINSZ, struct.pack("HHHH", *TERMINAL_SIZE, 0, 0))
        try:
            os.execvp(command[0], list(command))
        except OSError:
            os._exit(127)
    masker = Masker()
    decoder = codecs.getincrementaldecoder("utf-8")("replace")
    saved = termios.tcgetattr(stdin_fd) if os.isatty(stdin_fd) else None
    if saved is not None:
        tty.setraw(stdin_fd)
    try:
        pump(master, stdin_fd, lambda data: write(masker.feed(decoder.decode(data))))
    except KeyboardInterrupt:
        os.kill(pid, 15)
        raise
    finally:
        if saved is not None:
            # Keys typed after claude ended are dropped, not fed to the
            # wizard's next question.
            termios.tcsetattr(stdin_fd, termios.TCSAFLUSH, saved)
        os.close(master)
        write(masker.feed(decoder.decode(b"", final=True), final=True))
        if saved is not None:
            # The modes claude turns on (bracketed paste, focus and mouse
            # events, its keyboard protocol) are turned off again, since a
            # claude killed mid-screen leaves them on.
            write(TERMINAL_RESET)
        write("\n")
        _, status = os.waitpid(pid, 0)
    return masker.tokens, os.waitstatus_to_exitcode(status), masker.cut_over_lines


def minted_token(mint_with=None, **options):
    """One token from claude setup-token, or why there is none."""
    try:
        tokens, code, cut = (mint_with or mint)(**options)
    except KeyboardInterrupt:
        raise TokensError("stopped; nothing was stored") from None
    if cut:
        raise TokensError("claude's output went on with a run of a token's characters on the line after the token, as a token "
                          "cut over two lines by a narrow terminal does: nothing was stored, and the screen may show part of "
                          "it; run claude setup-token yourself in a wide window and use --paste")
    if len(tokens) == 1:
        return whole_token(tokens[0])
    if tokens:
        raise TokensError("claude printed more than one token, so none was stored; run it again")
    if code is None or code == 127:
        raise TokensError("claude could not be started: install Claude Code (https://code.claude.com/docs/en/setup), or use --paste")
    raise TokensError(f"claude setup-token ended (exit {code}) without printing a token: the sign-in did not finish, or claude prints it "
                      "in a form this wizard does not know; run it again, or run claude setup-token yourself and use --paste")


def whole_token(token):
    """A token of a token's length: a shorter one is a piece cut by a
    terminal, a longer one a run that went on past where a token ends."""
    if len(token) < TOKEN_SHORTEST:
        raise TokensError(f"the token seen is {len(token)} characters, shorter than any claude setup-token prints, so it was "
                          "cut: nothing was stored; run claude setup-token yourself in a wide window and use --paste")
    if len(token) > HOLD_LIMIT:
        raise TokensError("what claude printed ran on past a token's length with no break, so no token could be told apart "
                          "from it: nothing was stored; run claude setup-token yourself and use --paste")
    return token


def hidden_input(prompt):
    """getpass, or a refusal where it could not hide the input: it warns
    so and would echo, in a window that is no console."""
    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter("error", getpass.GetPassWarning)
        try:
            return getpass.getpass(prompt)
        except getpass.GetPassWarning:
            raise TokensError("this window cannot hide what you type, so nothing was pasted: run the wizard from a terminal "
                              "(cmd, PowerShell or Windows Terminal; in Git Bash, through winpty), or let it run claude "
                              "setup-token itself") from None


def pasted_token(ask=None):
    """The token from a hidden prompt, as gh's own prompt would take it."""
    ask = ask or hidden_input
    text = (ask("Paste the token claude setup-token printed, then Enter (nothing is shown): ") or "").strip()
    drop_typed_ahead()
    return check_token(text)


def drop_typed_ahead():
    """A paste with a line break in it leaves the rest in a Windows
    console's input, where the next question would read and echo it."""
    try:
        import msvcrt
    except ImportError:
        return
    while msvcrt.kbhit():
        msvcrt.getwch()


def check_token(text):
    """The token in what was pasted: claude's own line around it, such as
    its export line, is let through; an API key, or anything else, is not."""
    if not text:
        raise TokensError("nothing was pasted; nothing was stored")
    if API_KEY_HEAD in text:
        raise TokensError("that is an API key, which this project never uses (docs/project.md): make a subscription token with claude setup-token")
    found = TOKEN_SHAPE.findall(text)
    if len(set(found)) != 1:
        raise TokensError(f"that is not a token from claude setup-token (one begins with {TOKEN_HEAD}); nothing was stored")
    return whole_token(found[0])


def claude_command():
    """Where claude is: on PATH, or where Claude Code's own installer puts
    it behind a shell alias, which PATH does not show."""
    found = shutil.which(SETUP_TOKEN[0])
    if found:
        return found
    local = Path.home() / ".claude" / "local" / "claude"
    return str(local) if local.is_file() and os.access(local, os.X_OK) else None


def get_token(args, write=None, stdin_fd=None, ask=None, opener=None):
    """The token, minted behind the terminal where one exists, pasted where
    asked or where the machine has none. An opener sends claude's sign-in
    link to one browser profile."""
    if args.paste:
        return pasted_token(ask)
    claude = claude_command()
    if not pty_available():
        if claude is not None and sys.platform == "win32" and not console_keyboard():
            print("This window is no console claude can read its keyboard from (Git Bash's own window, or an IDE), so the "
                  "token is pasted: run the wizard from cmd, PowerShell or Windows Terminal to have it caught by itself.")
            return pasted_token(ask)
        if claude is not None and sys.platform == "win32":
            if opener is not None:
                open_profile(opener)
            show_terminal_codes()
            print("Starting claude setup-token. Approve the sign-in in your browser; the token stays out of sight.")
            if opener is not None:
                print(f"  If the sign-in opens in another profile, copy the link claude shows into {opener.browser.name}'s profile {opener.profile}.")
            return minted_token(mint_with=piped_mint, write=write, command=(claude,) + SETUP_TOKEN[1:])
        print("This shell has no pseudo-terminal, so the token is pasted: run claude setup-token in another window first.")
        return pasted_token(ask)
    if claude is None:
        print("claude is not on this machine's PATH, so the token is pasted: run claude setup-token where it is, then paste it here.")
        return pasted_token(ask)
    print("Starting claude setup-token. Approve the sign-in in your browser; the token stays out of sight.")
    with launcher(opener) as env:
        return minted_token(write=write, stdin_fd=stdin_fd, command=(claude,) + SETUP_TOKEN[1:], env=env)


def piped_mint(command, write=None, env=None, stdin=None, grace=0.5):
    """claude setup-token where there is no pseudo-terminal, as in a Windows
    shell: its keyboard is this console, which it shares, and its output
    comes here on a pipe and goes to the screen through the same mask, so
    the token is caught as behind a terminal. Measured on Linux with
    Claude Code 2.1.302: claude writes to a pipe as to a terminal, with its
    link opened and no line wrapped. The pipe is read in a thread and the
    run ends a moment after claude does, since a program claude started,
    such as a browser, may hold the pipe open long after."""
    if write is None:
        def write(text):
            sys.stdout.write(text)
            sys.stdout.flush()
    try:
        process = subprocess.Popen(list(command), stdin=stdin, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                   env={**os.environ, **(env or {})})
    except OSError:
        return [], None, False
    import queue
    import threading
    chunks = queue.Queue()

    def read_pipe():
        while True:
            try:
                data = os.read(process.stdout.fileno(), 4096)
            except OSError:
                data = b""
            chunks.put(data)
            if not data:
                return

    threading.Thread(target=read_pipe, daemon=True).start()
    masker = Masker()
    decoder = codecs.getincrementaldecoder("utf-8")("replace")
    ended = None
    try:
        while True:
            try:
                data = chunks.get(timeout=0.1)
            except queue.Empty:
                data = None
            if data == b"":
                break
            if data:
                write(masker.feed(decoder.decode(data)))
            if process.poll() is not None:
                ended = ended or time.monotonic()
                if time.monotonic() - ended > grace:
                    break
    except KeyboardInterrupt:
        process.kill()
        raise
    finally:
        write(masker.feed(decoder.decode(b"", final=True), final=True))
        if sys.stdout.isatty():
            # Modes claude turns on stay on in the console if it dies
            # mid-screen, and keys typed meanwhile would answer the
            # wizard's next question.
            write(TERMINAL_RESET)
        write("\n")
        drop_typed_ahead()
        code = process.wait()
    return masker.tokens, code, masker.cut_over_lines


def console_keyboard():
    """Whether this process reads a Windows console, which claude can take
    for its keyboard; a pipe, as in Git Bash's own window or an IDE, it
    cannot, and better found before the sign-in than in it."""
    if not (sys.stdin and sys.stdin.isatty()):
        return False
    try:
        import ctypes
        kernel = ctypes.windll.kernel32
        kernel.GetStdHandle.restype = ctypes.c_void_p
        mode = ctypes.c_uint32()
        return bool(kernel.GetConsoleMode(ctypes.c_void_p(kernel.GetStdHandle(-10)), ctypes.byref(mode)))
    except (AttributeError, OSError):
        return False


def show_terminal_codes():
    """A Windows console shows claude's colours and redraws only once asked
    to read terminal codes, as Windows Terminal does by itself."""
    if sys.platform != "win32":
        return
    try:
        import ctypes
        kernel = ctypes.windll.kernel32
        kernel.GetStdHandle.restype = ctypes.c_void_p
        handle = ctypes.c_void_p(kernel.GetStdHandle(-11))
        mode = ctypes.c_uint32()
        if kernel.GetConsoleMode(handle, ctypes.byref(mode)):
            kernel.SetConsoleMode(handle, mode.value | 0x0004)
    except (AttributeError, OSError):
        pass


def open_profile(opener, settle=2.0):
    """In a Windows shell no launcher goes through BROWSER: Claude Code
    hands a batch file its link through cmd.exe, which cuts it at its =
    and & signs. The profile is opened first instead, on a blank page so
    no second sign-in shows, and the browser sends the link claude opens
    next to the profile used last."""
    try:
        subprocess.Popen([opener.browser.path, f"--profile-directory={opener.profile}", "about:blank"])
    except OSError as error:
        raise TokensError(f"{opener.browser.name} could not be started ({error.strerror or error}); answer n at the profile question") from None
    time.sleep(settle)


# The browser profile and cswap's accounts, remembered on this machine alone.
#
# claude setup-token opens its sign-in link with the program BROWSER names,
# handing it the link as its one argument, on every platform (Claude Code
# 2.1.302's own opener). So a slot's sign-in can go to one profile of a
# Chromium browser through a small launcher made for the one run, and a
# profile signed in to claude.ai as the slot's account once makes every
# later token one Approve. What is remembered is a profile's folder name
# and a cswap account's number per slot, in a file in this machine's own
# settings folder, never in the repository and never an address.

PROFILE_SHAPE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9 ._-]{0,39}$")
SETTINGS_VERSION = 1
READ_LIMIT = 4_000_000
CSWAP_NUMBER = range(1, 1000)


class Browser:
    """A Chromium browser: its program, and the folder whose Local State
    lists its profiles where that is known."""

    def __init__(self, name, path, data=None):
        self.name = name
        self.path = str(path)
        self.data = Path(data) if data else None


class Opener:
    """Where one sign-in link goes: a browser's program and a profile's folder."""

    def __init__(self, browser, profile):
        check_profile(profile)
        check_program(browser.path)
        self.browser = browser
        self.profile = profile


class CswapAccount:
    def __init__(self, number, email, alias, organization):
        self.number = number
        self.email = email
        self.alias = alias
        self.organization = organization


def system_name(env=None, version_file="/proc/version"):
    """win32, darwin, wsl or linux: WSL is Linux whose browsers are
    Windows programs under /mnt/c."""
    env = os.environ if env is None else env
    if sys.platform in ("win32", "darwin"):
        return sys.platform
    if env.get("WSL_DISTRO_NAME"):
        return "wsl"
    try:
        with open(version_file, encoding="utf-8", errors="replace") as handle:
            return "wsl" if "microsoft" in handle.read(1000).lower() else "linux"
    except OSError:
        return "linux"


def browser_candidates(system, env, home):
    """Each Chromium browser where its installer puts it: (name, program,
    profiles folder)."""
    if system == "win32":
        local = env.get("LOCALAPPDATA") or str(home / "AppData" / "Local")
        programs = [env.get("ProgramFiles") or r"C:\Program Files", env.get("ProgramFiles(x86)") or r"C:\Program Files (x86)", local]
        found = []
        for name, folder, program in (("Chrome", r"Google\Chrome", "chrome.exe"), ("Edge", r"Microsoft\Edge", "msedge.exe"),
                                      ("Brave", r"BraveSoftware\Brave-Browser", "brave.exe")):
            data = local + "\\" + folder + "\\User Data"
            found += [(name, base + "\\" + folder + "\\Application\\" + program, data) for base in programs]
        return found
    if system == "darwin":
        support = home / "Library" / "Application Support"
        return [("Chrome", "/Applications/Google Chrome.app/Contents/MacOS/Google Chrome", support / "Google" / "Chrome"),
                ("Edge", "/Applications/Microsoft Edge.app/Contents/MacOS/Microsoft Edge", support / "Microsoft Edge"),
                ("Brave", "/Applications/Brave Browser.app/Contents/MacOS/Brave Browser", support / "BraveSoftware" / "Brave-Browser")]
    if system == "wsl":
        # The Windows user's profiles folder is not read from WSL; a
        # profile is named by its folder instead of picked from a list.
        return [("Chrome", base + "/Google/Chrome/Application/chrome.exe", None)
                for base in ("/mnt/c/Program Files", "/mnt/c/Program Files (x86)")] + [
                ("Edge", base + "/Microsoft/Edge/Application/msedge.exe", None)
                for base in ("/mnt/c/Program Files (x86)", "/mnt/c/Program Files")]
    config = home / ".config"
    return [("Chrome", "google-chrome", config / "google-chrome"), ("Chrome", "google-chrome-stable", config / "google-chrome"),
            ("Chromium", "chromium", config / "chromium"), ("Chromium", "chromium-browser", config / "chromium"),
            ("Edge", "microsoft-edge", config / "microsoft-edge"), ("Brave", "brave-browser", config / "BraveSoftware" / "Brave-Browser")]


def find_browser(system=None, env=None, home=None, which=shutil.which, is_file=os.path.isfile):
    """The browser SYMPHONY_BROWSER names, else the first Chromium browser
    installed, else None."""
    env = os.environ if env is None else env
    system = system or system_name(env)
    home = Path(home) if home else Path.home()
    named = env.get("SYMPHONY_BROWSER", "").strip()
    if named:
        return Browser("your browser", named)
    for name, program, data in browser_candidates(system, env, home):
        path = program if os.path.isabs(program) or "\\" in program else which(program)
        if path and is_file(path):
            return Browser(name, path, data)
    return None


def read_json_file(path):
    """A small JSON file, or None where it is missing, too big or not JSON."""
    try:
        with open(path, "rb") as handle:
            data = handle.read(READ_LIMIT + 1)
        if len(data) > READ_LIMIT:
            return None
        return json.loads(data.decode("utf-8"))
    except (OSError, ValueError):
        return None


def browser_profiles(browser):
    """The browser's own profiles as (folder, name), from its Local State,
    or none where that is not readable."""
    if browser is None or browser.data is None:
        return []
    state = read_json_file(browser.data / "Local State")
    cache = state.get("profile", {}).get("info_cache") if isinstance(state, dict) and isinstance(state.get("profile"), dict) else None
    if not isinstance(cache, dict):
        return []
    found = []
    for folder, info in cache.items():
        if isinstance(folder, str) and PROFILE_SHAPE.fullmatch(folder):
            name = info.get("name") if isinstance(info, dict) else None
            found.append((folder, name if isinstance(name, str) and name.isprintable() else folder))
    return sorted(found, key=lambda each: (each[0] != "Default", each[0]))


def check_profile(profile):
    if not PROFILE_SHAPE.fullmatch(profile or ""):
        raise TokensError("a browser profile is a folder name such as Default, \"Profile 2\" or \"Symphony 1\": up to 40 letters, "
                          "digits, spaces, dots or dashes")


def check_program(path):
    """A program path goes into a launcher between double quotes, so one
    that could end the quotes or expand in a batch file is refused."""
    if not path or not path.isascii() or any(mark in path for mark in ('"', "%", "\n", "\r", "!", "`", "$")):
        raise TokensError("the browser's path holds a letter outside ASCII, a quote, %, !, $ or a backquote, which the launcher "
                          "cannot carry; set SYMPHONY_BROWSER to a path without them, or answer n at the profile question")


def launcher_text(opener):
    """A launcher that hands the link it is given to the profile and returns
    at once, since claude waits on the program it starts. The hangup is
    ignored before the browser is forked, not by nohup after: the browser
    shares claude's terminal, and a hangup sent as that terminal closes
    reached it in the moment before nohup could ignore it."""
    import shlex
    flag = f"--profile-directory={opener.profile}"
    return (f"#!/bin/sh\ntrap '' HUP\n"
            f'{shlex.quote(opener.browser.path)} {shlex.quote(flag)} "$1" </dev/null >/dev/null 2>&1 &\n')


@contextmanager
def launcher(opener):
    """BROWSER for claude's environment, pointing at a launcher in a folder
    of the run's own that is gone when the sign-in is over; none without
    an opener, so claude opens the default browser as it does by itself.
    Linux, macOS and WSL only: Windows opens the profile first instead
    (open_profile)."""
    if opener is None:
        yield {}
        return
    folder = tempfile.mkdtemp(prefix="symphony-browser-")
    try:
        path = Path(folder) / "open"
        path.write_bytes(launcher_text(opener).encode("ascii", "strict"))
        path.chmod(0o700)
        yield {"BROWSER": str(path)}
    finally:
        shutil.rmtree(folder, ignore_errors=True)


def settings_path(system=None, env=None, home=None):
    """Where this machine keeps what the wizard remembers: the user's own
    settings folder, never the repository."""
    env = os.environ if env is None else env
    system = system or system_name(env)
    home = Path(home) if home else Path.home()
    if system == "win32":
        base = Path(env.get("APPDATA") or home / "AppData" / "Roaming")
    else:
        xdg = env.get("XDG_CONFIG_HOME", "")
        base = Path(xdg) if xdg and os.path.isabs(xdg) else home / ".config"
    return base / "symphony" / "tokens.json"


def read_settings(path):
    data = read_json_file(path) if path else None
    if not isinstance(data, dict) or data.get("version") != SETTINGS_VERSION or not isinstance(data.get("slots"), dict):
        return {"version": SETTINGS_VERSION, "slots": {}}
    return data


def place_key(target):
    return f"{target.level}:{target.full_name if target.level == 'repo' else target.owner}"


def remembered(settings, target, slot):
    """What this machine remembers for a slot: a profile's folder and a
    cswap account's number, each only where it still has its shape."""
    entry = settings["slots"].get(place_key(target), {})
    entry = entry.get(str(slot)) if isinstance(entry, dict) else None
    entry = entry if isinstance(entry, dict) else {}
    profile = entry.get("profile")
    number = entry.get("cswap")
    return {"profile": profile if isinstance(profile, str) and PROFILE_SHAPE.fullmatch(profile) else None,
            "cswap": number if type(number) is int and number in CSWAP_NUMBER else None}


def remember(path, settings, target, slot, profile, number):
    """Keep a slot's profile and cswap account for the next run; a file
    that cannot be written costs only the remembering."""
    if path is None:
        return
    entry = {key: value for key, value in (("profile", profile), ("cswap", number)) if value is not None}
    slots = settings["slots"].setdefault(place_key(target), {})
    if not isinstance(slots, dict):
        slots = settings["slots"][place_key(target)] = {}
    if entry == {key: value for key, value in remembered(settings, target, slot).items() if value is not None}:
        return
    if entry:
        slots[str(slot)] = entry
    else:
        slots.pop(str(slot), None)
    try:
        Path(path).parent.mkdir(parents=True, exist_ok=True)
        part = Path(str(path) + ".part")
        part.write_text(json.dumps(settings, indent=2, sort_keys=True) + "\n", encoding="ascii")
        part.chmod(0o600)
        os.replace(part, path)
    except OSError as error:
        print(f"  (not remembered for next time: {path} could not be written, {error.strerror or error})")


def cswap_roots(system, env, home):
    """Where claude-swap keeps its list of accounts: XDG's data folder on
    Linux and WSL, with the older ~/.claude-swap-backup it moves from, and
    that older folder alone on Windows and macOS."""
    roots = []
    if system in ("linux", "wsl"):
        xdg = env.get("XDG_DATA_HOME", "")
        xdg = os.path.expanduser(xdg) if xdg else ""
        roots.append(Path(xdg) / "claude-swap" if xdg and os.path.isabs(xdg) else home / ".local" / "share" / "claude-swap")
    roots.append(home / ".claude-swap-backup")
    return roots


def cswap_accounts(system=None, env=None, home=None):
    """claude-swap's accounts in its own order, read from its sequence.json
    alone: no credential of cswap's is opened, and nothing is fetched."""
    env = os.environ if env is None else env
    system = system or system_name(env)
    home = Path(home) if home else Path.home()
    for root in cswap_roots(system, env, home):
        data = read_json_file(root / "sequence.json")
        if not isinstance(data, dict):
            continue
        records = data.get("accounts") if isinstance(data.get("accounts"), dict) else {}
        found, seen = [], set()
        listed = data.get("sequence")
        for number in listed if isinstance(listed, list) else []:
            record = records.get(str(number))
            if type(number) is not int or number not in CSWAP_NUMBER or not isinstance(record, dict) or number in seen:
                continue
            seen.add(number)
            fields = [record.get(key) for key in ("email", "alias", "organizationName")]
            found.append(CswapAccount(number, *[value if isinstance(value, str) and value.isprintable() else "" for value in fields]))
        return found
    return []


def account_line(account):
    """How an account is shown on this machine's screen; the address is
    shown here to tell accounts apart and is never stored or sent."""
    name = f"{account.alias} ({account.email})" if account.alias else account.email or "an account with no address"
    return f"  {account.number}: {name}" + (f" [{account.organization}]" if account.organization else "")


def chosen_account(accounts, known, ask):
    """The cswap account a slot is for, by its number in cswap's list, or
    None where cswap has none or the person skips it."""
    if not accounts:
        return None
    numbers = {account.number: account for account in accounts}
    print("cswap's accounts on this machine:")
    for account in accounts:
        print(account_line(account))
    default = known if known in numbers else None
    shown = f"{default}, {numbers[default].alias or numbers[default].email}, as last time" if default else "Enter to skip"
    answer = (ask(f"Which of them is this slot for? [{shown}; its number] ") or "").strip()
    if not answer:
        return numbers.get(default)
    if answer.isascii() and answer.isdigit() and len(answer) <= 3 and int(answer) in numbers:
        return numbers[int(answer)]
    raise TokensError("an account is one of the numbers cswap's list shows")


def label_from(account):
    """A cswap alias as the slot's label where it has a label's shape; an
    address never is one."""
    alias = account.alias if account else ""
    return alias if alias and LABEL_SHAPE.fullmatch(alias) else ""


def chosen_profile(browser, slot, known, ask):
    """The browser profile a slot's sign-in opens in: one of the browser's
    own by its number, a folder by its name, a new one named for the slot
    by default, or none, for the default browser."""
    if browser is None:
        return None
    profiles = browser_profiles(browser)
    if profiles:
        print(f"{browser.name}'s profiles:")
        for number, (folder, name) in enumerate(profiles, 1):
            print(f"  {number}: {name}" + (f" ({folder})" if name != folder else ""))
    default = known or f"Symphony {slot}"
    answer = (ask(f"Which {browser.name} profile opens this slot's sign-in? [{default}; "
                  + ("a number above, " if profiles else "") + "a profile folder's name, or n for your default browser] ") or "").strip()
    if not answer:
        return default
    if answer.lower() in ("n", "no", "none"):
        return None
    if answer.isascii() and answer.isdigit() and len(answer) <= 3 and 1 <= int(answer) <= len(profiles):
        return profiles[int(answer) - 1][0]
    check_profile(answer)
    return answer


class Local:
    """What this machine offers the wizard: a browser, cswap's accounts and
    a settings file."""

    def __init__(self, browser=None, accounts=(), settings=None):
        self.browser = browser
        self.accounts = list(accounts)
        self.settings = settings

    @classmethod
    def discover(cls):
        return cls(find_browser(), cswap_accounts(), settings_path())


# The copy on this machine.

# Credential Manager's kind for a credential a program finds by its name,
# the persistence that keeps it on this machine alone, never roaming with
# the Windows profile to another, and the error a missing name gives.
CRED_TYPE_GENERIC = 1
CRED_PERSIST_LOCAL_MACHINE = 2
ERROR_NOT_FOUND = 1168


def copy_name(target, name):
    """The name a copy goes under: the place and the secret, so one slot
    of two places never meets the other."""
    return f"symphony/{place_key(target)}/{name}"


def credential_type():
    """Windows' CREDENTIALW, field for field; DWORD is 32 bits and every
    string a pointer to UTF-16, on Windows as anywhere ctypes builds it."""
    import ctypes

    class FileTime(ctypes.Structure):
        _fields_ = [("low", ctypes.c_uint32), ("high", ctypes.c_uint32)]

    class Credential(ctypes.Structure):
        _fields_ = [("Flags", ctypes.c_uint32), ("Type", ctypes.c_uint32),
                    ("TargetName", ctypes.c_wchar_p), ("Comment", ctypes.c_wchar_p),
                    ("LastWritten", FileTime), ("CredentialBlobSize", ctypes.c_uint32),
                    ("CredentialBlob", ctypes.POINTER(ctypes.c_ubyte)), ("Persist", ctypes.c_uint32),
                    ("AttributeCount", ctypes.c_uint32), ("Attributes", ctypes.c_void_p),
                    ("TargetAlias", ctypes.c_wchar_p), ("UserName", ctypes.c_wchar_p)]

    return Credential


class CredentialManager:
    """Windows Credential Manager, the user's own vault: Windows encrypts
    each entry under the user's sign-in, only a program running as that
    user reads one back, and the Control Panel's Credential Manager (under
    Windows Credentials) lists and deletes them by name. The value is
    UTF-16, as cmdkey and PowerShell write and read a generic credential."""

    def __init__(self, advapi=None, last_error=None):
        import ctypes
        self.ctypes = ctypes
        self.Credential = credential_type()
        self.advapi = advapi or ctypes.WinDLL("advapi32", use_last_error=True)
        self.last_error = last_error or ctypes.get_last_error
        pointer = ctypes.POINTER(self.Credential)
        self.advapi.CredWriteW.argtypes = [pointer, ctypes.c_uint32]
        self.advapi.CredReadW.argtypes = [ctypes.c_wchar_p, ctypes.c_uint32, ctypes.c_uint32, ctypes.POINTER(pointer)]
        self.advapi.CredDeleteW.argtypes = [ctypes.c_wchar_p, ctypes.c_uint32, ctypes.c_uint32]
        self.advapi.CredFree.argtypes = [ctypes.c_void_p]
        for name in ("CredWriteW", "CredReadW", "CredDeleteW"):
            getattr(self.advapi, name).restype = ctypes.c_int

    def put(self, name, value, user, comment):
        ctypes = self.ctypes
        blob = value.encode("utf-16-le")
        buffer = (ctypes.c_ubyte * len(blob)).from_buffer_copy(blob)
        credential = self.Credential()
        credential.Type = CRED_TYPE_GENERIC
        credential.TargetName = name
        credential.Comment = comment
        credential.CredentialBlobSize = len(blob)
        credential.CredentialBlob = ctypes.cast(buffer, ctypes.POINTER(ctypes.c_ubyte))
        credential.Persist = CRED_PERSIST_LOCAL_MACHINE
        credential.UserName = user
        if not self.advapi.CredWriteW(ctypes.byref(credential), 0):
            raise TokensError(f"Credential Manager refused the copy (Windows error {self.last_error()})")

    def get(self, name):
        """The value kept under the name, or None where none is."""
        ctypes = self.ctypes
        found = ctypes.POINTER(self.Credential)()
        if not self.advapi.CredReadW(name, CRED_TYPE_GENERIC, 0, ctypes.byref(found)):
            error = self.last_error()
            if error == ERROR_NOT_FOUND:
                return None
            raise TokensError(f"Credential Manager could not read {name} (Windows error {error})")
        try:
            credential = found.contents
            blob = ctypes.string_at(credential.CredentialBlob, credential.CredentialBlobSize)
        finally:
            self.advapi.CredFree(found)
        try:
            return blob.decode("utf-16-le")
        except UnicodeDecodeError:
            raise TokensError(f"the copy under {name} is not UTF-16 text; delete it in Credential Manager") from None

    def drop(self, name):
        """True where a copy was there to delete."""
        if self.advapi.CredDeleteW(name, CRED_TYPE_GENERIC, 0):
            return True
        error = self.last_error()
        if error == ERROR_NOT_FOUND:
            return False
        raise TokensError(f"Credential Manager could not delete {name} (Windows error {error})")


def local_copies():
    """Where this machine keeps copies: Credential Manager on Windows, and
    nowhere elsewhere, where the token is only on GitHub."""
    if sys.platform != "win32":
        return None
    try:
        return CredentialManager()
    except (OSError, AttributeError):
        return None


def machine_copies(copies):
    """None asks this machine; False, --no-copy, keeps none."""
    return local_copies() if copies is None else (copies or None)


def keep_copy(copies, target, name, token):
    """The token, once on GitHub, kept on this machine too, so show can
    give it back. A copy that fails costs only the copy: the secret is
    stored already."""
    if copies is None:
        return
    try:
        copies.put(copy_name(target, name), token, name,
                   f"Claude token {name} of {place_key(target)}; python3 scripts/tokens.py show reads it")
        print(f"  {name}: a copy kept in Windows Credential Manager on this machine (show reads it back)")
    except (TokensError, OSError) as error:
        print(f"  {name}: no copy kept on this machine ({error}); the secret on GitHub is stored all the same")


def drop_copy(copies, target, name):
    if copies is None:
        return
    try:
        if copies.drop(copy_name(target, name)):
            print(f"  {name}: its copy on this machine deleted")
    except (TokensError, OSError) as error:
        print(f"  {name}: its copy on this machine could not be deleted ({error}); delete {copy_name(target, name)} in Credential Manager")


def has_copy(copies, target, name):
    try:
        return copies is not None and copies.get(copy_name(target, name)) is not None
    except (TokensError, OSError):
        return False


# What is there.

def slot_secret(slot):
    return f"{SLOT_PREFIX}{slot}"


def slot_label(slot):
    return f"{SLOT_PREFIX}{slot}{LABEL_SUFFIX}"


def list_words(text):
    """The words of a list of slots, split as the bench splits them."""
    return [word for word in LIST_SEPARATORS.split((text or "").strip()) if word]


def slot_number(text):
    """The slot a word names, or None. The word is looked up among the words
    the bench takes, never read by int(), which would take 08 and refuse a
    digit such as a superscript two that isdigit() takes, or more than 4,300
    digits."""
    return SLOT_WORDS.get(text)


def parse_slots(text):
    """The slots a list names, each once, in the list's own order, which is
    the order the bench deals its jobs to them."""
    return list(dict.fromkeys(slot for slot in map(slot_number, list_words(text)) if slot is not None))


def holds_other_words(text):
    """Whether a list holds a word that is not a slot, which stops the
    bench's real runs whatever else it names."""
    return any(slot_number(word) is None for word in list_words(text))


def slots_text(text):
    """A list as status shows it: the slots it names, never the text as
    written, since a list may be made by hand."""
    shown = ",".join(str(slot) for slot in parse_slots(text)) or "no slot"
    if holds_other_words(text):
        shown += ", and a word that is not a slot, which stops the bench's real runs"
    elif len(parse_slots(text)) != len(list_words(text)):
        shown += ", with a slot named twice, which gives it more of the bench's jobs; the next add, remove or order rewrites it"
    return shown


class State:
    """The secrets and variables at the target, by name, with the
    repositories each selected one reaches, read when asked."""

    def __init__(self, gh, target):
        self.gh = gh
        self.target = target
        fields = "name,updatedAt,visibility,numSelectedRepos"
        self.secrets = {item["name"]: item for item in gh.json(["secret", "list"] + target.flags + ["--json", fields])}
        self.variables = {item["name"]: item for item in gh.json(["variable", "list"] + target.flags + ["--json", fields + ",value"])}
        self.reached = {}

    def reach(self, kind, name):
        """The repositories a selected organization secret or variable
        reaches, listed once for the state: a write is always followed by a
        new state before anything is read from it again."""
        if (kind, name) not in self.reached:
            found, text = self.gh.read(["api", f"{self.target.api}/{kind}/{name}/repositories", "--paginate", "--jq", ".repositories[].name"])
            self.reached[(kind, name)] = sorted(text.split()) if found else []
        return self.reached[(kind, name)]

    def slots_set(self):
        return [slot for slot in SLOTS if slot_secret(slot) in self.secrets]

    def labels(self):
        """Each slot's label as status may print it. A label may be made by
        hand on the variables page, so one that is not a word of the label's
        shape is never printed as written."""
        found = {}
        for slot in SLOTS:
            item = self.variables.get(slot_label(slot))
            if item is not None:
                value = item.get("value", "")
                found[slot] = value if LABEL_SHAPE.fullmatch(value) else "a label not shown, not a word of a label's shape"
        return found


def when(item):
    """The moment GitHub last saw a value stored, or None."""
    text = (item or {}).get("updatedAt") or (item or {}).get("createdAt") or ""
    try:
        return datetime.fromisoformat(text.replace("Z", "+00:00"))
    except ValueError:
        return None


def age_line(item, now=None):
    """set 2026-10-16, renew by 2027-10-16; a token lasts about a year."""
    stored = when(item)
    if stored is None:
        return "set"
    now = now or datetime.now(timezone.utc)
    renew = stored + TOKEN_LIFE
    note = ""
    if renew <= now:
        note = ", renew now: it is older than a year"
    elif renew - now <= RENEW_SOON:
        note = ", renew soon"
    return f"set {stored.date().isoformat()}, renew by {renew.date().isoformat()}{note}"


def print_status(state, now=None):
    target = state.target
    print(f"Claude tokens on {target.describe()}")
    labels = state.labels()
    for slot in SLOTS:
        name = slot_secret(slot)
        if name not in state.secrets:
            continue
        label = f" ({labels[slot]})" if labels.get(slot) else ""
        print(f"  slot {slot}{label}: {name} {age_line(state.secrets[name], now)}{reach_note(state, 'secrets', name)}")
    free = [slot for slot in SLOTS if slot_secret(slot) not in state.secrets]
    print("  free slots: " + (", ".join(str(slot) for slot in free) if free else "none"))
    have = state.variables.get(ACCOUNTS_VARIABLE)
    if have is None:
        print(f"  {ACCOUNTS_VARIABLE}: not set")
    else:
        value = have.get("value", "")
        rewrite = "; the next add or remove rewrites it" if holds_other_words(value) else ""
        print(f"  {ACCOUNTS_VARIABLE}: {slots_text(value)}{rewrite}")
    if MAIN_SECRET in state.secrets:
        print(f"  {MAIN_SECRET}: {age_line(state.secrets[MAIN_SECRET], now)}{reach_note(state, 'secrets', MAIN_SECRET)} (the workflows' token)")
    else:
        print(f"  {MAIN_SECRET}: not set here, so @claude, the workers in Actions and the coordinator workflow do nothing")
    stale = [slot for slot in labels if slot_secret(slot) not in state.secrets]
    for slot in stale:
        print(f"  the label of slot {slot} ({labels[slot]}) outlives its secret; remove {slot} clears it")
    if target.level == "org":
        print_repository_own(state)


def reach_note(state, kind, name):
    item = state.secrets.get(name) if kind == "secrets" else state.variables.get(name)
    if state.target.level != "org" or (item or {}).get("visibility") != "selected":
        return ""
    try:
        reach = state.reach(kind, name)
    except TokensError as error:
        return f", seen by repositories this login could not list ({error})"
    return ", seen by " + (", ".join(reach) if reach else "no repository")


def print_repository_own(state):
    """A repository's own secret or variable of the same name outranks the
    organization's, so status lists them for every repository that reads
    what is kept here: a list of slots of a repository's own, made by hand
    before the wizard, would hide from its workflows the list the wizard
    keeps."""
    for full, names in repositories_reading(state).items():
        print_own_secrets(state, full, names)
        print_own_accounts(state, full, names)


def repositories_reading(state):
    """Each repository that reads what is kept here, with the names it reads
    from the organization. A repository named reads any of them (None),
    since the wizard writes each with it; another, found on the listing of
    a name of selected visibility, such as the bench's on the list's, reads
    the names whose listing shows it and the ones every repository sees."""
    reading = {full: None for full in state.target.repositories()}
    known = {full.lower(): full for full in reading}
    names = [MAIN_SECRET] + [slot_secret(slot) for slot in SLOTS]
    kept = [("variables", ACCOUNTS_VARIABLE, state.variables.get(ACCOUNTS_VARIABLE))]
    kept += [("secrets", name, state.secrets[name]) for name in names if name in state.secrets]
    everywhere = {name for _, name, item in kept if (item or {}).get("visibility") in ("all", "private")}
    for kind, name, item in kept:
        if (item or {}).get("visibility") != "selected":
            continue
        try:
            reached = state.reach(kind, name)
        except TokensError as error:
            if kind == "variables":
                print(f"  the repositories that see {name} could not be listed ({error})")
            continue
        for repo in reached:
            full = known.setdefault(f"{state.target.owner}/{repo}".lower(), f"{state.target.owner}/{repo}")
            if full not in reading:
                reading[full] = set(everywhere)
            if reading[full] is not None:
                reading[full].add(name)
    return reading


def print_own_secrets(state, full, names=None):
    """A repository's own token secrets that outrank the organization's, of
    the names it reads from the organization (any, for None)."""
    try:
        own = {item["name"].upper(): item for item in state.gh.json(["secret", "list", "--repo", full, "--json", "name,updatedAt"])}
    except TokensError as error:
        print(f"  {full}: its own secrets could not be read ({error})")
        return
    mine = [name for name in own if (name == MAIN_SECRET or name.startswith(SLOT_PREFIX)) and (names is None or name in names)]
    for name in mine:
        print(f"  {full} holds its own {name} ({age_line(own[name])}), which outranks the organization's")


def print_own_accounts(state, full, names=None):
    """A repository's own list of slots, when it would hide the
    organization's: a repository the organization's list does not reach
    keeps its own, which is its only list."""
    if names is not None and ACCOUNTS_VARIABLE not in names:
        return
    try:
        own = {item["name"].upper(): item for item in state.gh.json(["variable", "list", "--repo", full, "--json", "name,value"])}
    except TokensError as error:
        print(f"  {full}: its own variables could not be read ({error})")
        return
    if ACCOUNTS_VARIABLE in own:
        print(f"  {full} holds its own {ACCOUNTS_VARIABLE} ({slots_text(own[ACCOUNTS_VARIABLE].get('value', ''))}), which "
              f"outranks the organization's, so its workflows read that list and not the one kept here: delete it on "
              f"https://github.com/{full}/settings/variables/actions (the delete icon beside it)")


# Writes.

def scope_args(state, kind, name):
    """The visibility a write at the organization sends: what the secret or
    variable has, its list grown by the target's, or the target's for a new
    one. gh's default would be "private", which a write with no visibility
    would quietly set on an existing one."""
    target = state.target
    if target.level != "org":
        return []
    have = state.secrets.get(name) if kind == "secrets" else state.variables.get(name)
    visibility = (have or {}).get("visibility") or target.visibility
    if visibility != "selected":
        return ["--visibility", visibility]
    repos = sorted(set(target.repos) | set(state.reach(kind, name) if have else []))
    if not repos:
        raise TokensError(f"{name} would be seen by no repository: name them with --repos")
    return ["--visibility", "selected", "--repos", ",".join(repos)]


def store_secret(state, name, token, scope):
    """The value goes to gh on its standard input, where gh reads it without
    echo or a command line; the arguments carry the name and the place, the
    scope found before the token was."""
    args = ["secret", "set", name] + state.target.flags + scope
    state.gh.must(args, stdin=token)
    print(f"  {name}: stored on {state.target.describe()}")


def set_variable(state, name, value, scope=None):
    scope = scope_args(state, "variables", name) if scope is None else scope
    args = ["variable", "set", name] + state.target.flags + scope + ["--body", value]
    state.gh.must(args)
    print(f"  {name}: {value}")


def delete_secret(state, name):
    state.gh.must(["secret", "delete", name] + state.target.flags)
    print(f"  {name}: deleted")


def delete_variable(state, name):
    state.gh.must(["variable", "delete", name] + state.target.flags)
    print(f"  {name}: deleted")


def verify(state, names):
    """GitHub lists each stored name, by name alone, with a date later than
    the one it had: GitHub's own dates on both sides, so a clock that is off
    on this machine proves nothing."""
    fresh = State(state.gh, state.target)
    for name in names:
        before, after = when(state.secrets.get(name)), when(fresh.secrets.get(name))
        if after is None:
            raise TokensError(f"GitHub does not list {name} after storing it; open the secrets page and look")
        if before is not None and after <= before:
            raise TokensError(f"GitHub lists {name} as updated {after.isoformat()}, as before the store; it did not take")
        print(f"  {name}: GitHub lists it, updated {after.strftime('%Y-%m-%d %H:%M')} UTC")
    return fresh


def ordered_slots(state, order=None):
    """Every slot whose secret exists, in the order given, else in the
    list's own order, a slot new to it after the rest."""
    have = state.variables.get(ACCOUNTS_VARIABLE)
    listed = order if order is not None else parse_slots(have.get("value", "") if have else "")
    present = state.slots_set()
    kept = [slot for slot in listed if slot in present]
    return kept + [slot for slot in present if slot not in kept]


def update_accounts(state, scope=None):
    """The variable names every slot whose secret exists, as the bench reads
    it, and is written only when that changed: a list with any other word,
    which the bench refuses, is written again."""
    words = [str(slot) for slot in ordered_slots(state)]
    wanted = ",".join(words)
    have = state.variables.get(ACCOUNTS_VARIABLE)
    if have is not None and list_words(have.get("value", "")) == words:
        print(f"  {ACCOUNTS_VARIABLE}: {wanted or 'empty'} (unchanged)")
        return
    if not wanted:
        if have is not None:
            delete_variable(state, ACCOUNTS_VARIABLE)
        return
    set_variable(state, ACCOUNTS_VARIABLE, wanted, scope)


def refuse_crowd(gh, target, project=None):
    """Each repository named is one this login can see, and none holds a
    Claude token for being in crowd mode (docs/setup/crowdsourcing.md):
    this clone's mode line and each repository's variable say so."""
    try:
        text = (PROJECT if project is None else project).read_text(encoding="utf-8")
    except OSError:
        text = ""
    if CROWD_ON.search(text):
        raise TokensError("crowd mode is on in this clone's docs/project.md: a crowd repository holds no Claude token, so nothing was stored")
    for full in target.repositories():
        found, _ = gh.read(["api", f"repos/{full}", "--jq", ".full_name"])
        if not found:
            raise TokensError(f"GitHub lists no repository {full} for this login: check the name, or --repos")
        found, value = gh.read(["api", f"repos/{full}/actions/variables/{CROWD_VARIABLE}", "--jq", ".value"])
        if found and value.strip().lower() == "true":
            raise TokensError(f"{full} is in crowd mode (its variable {CROWD_VARIABLE} is true): it holds no Claude token, so nothing was stored")


def auth_status(gh):
    """gh auth status for the active github.com account alone: the plain
    command fails on a stale account of any host it knows, and lists that
    host first; a gh too old for --active takes the host alone."""
    done = gh.run(["auth", "status", "--hostname", "github.com", "--active"])
    if done.returncode != 0 and "unknown flag" in (done.stderr or ""):
        done = gh.run(["auth", "status", "--hostname", "github.com"])
    return done


def check_signed_in(gh):
    """What gh says of the github.com account, once it says one is signed in."""
    done = auth_status(gh)
    if done.returncode != 0:
        raise TokensError("gh is not signed in to github.com: run gh auth login (https://cli.github.com/manual/gh_auth_login), then run this again")
    return (done.stdout or "") + (done.stderr or "")


def check_scopes(status, target):
    """An organization's secrets need gh's token to carry admin:org, which
    gh auth login does not grant by default; better found here than after
    a sign-in. gh quotes each scope of an OAuth token; a token of another
    kind (a fine-grained token in GH_TOKEN) lists none, and gh says itself
    what such a token lacks when a call is refused."""
    if target.level != "org":
        return
    line = re.search(r"Token scopes: (.*)", status)
    scopes = re.findall(r"'([^']+)'", line.group(1)) if line else []
    if scopes and "admin:org" not in scopes:
        raise TokensError("gh's token lacks the admin:org scope an organization's secrets need: run gh auth refresh -h github.com -s admin:org, "
                          "approve it in the browser, then run this again")


def check_label(label):
    if label and not LABEL_SHAPE.fullmatch(label):
        raise TokensError("a label is a word of your own, up to 40 letters, digits, spaces, dots or dashes, and never an e-mail address")


def first_free(state):
    for slot in SLOTS:
        if slot_secret(slot) not in state.secrets:
            return slot
    raise TokensError(f"every slot from 1 to {SLOTS[-1]} is set; remove one, or name one to replace with --slot")


# The commands.

def sign_in_note(slot, opener=None):
    if opener is not None:
        print(f"Slot {slot}: the sign-in opens in {opener.browser.name}'s profile {opener.profile}. The first time, sign in to claude.ai")
        print("  there as the account this slot is for; the browser keeps it, so later the sign-in is one Approve.")
        print("  The page that opens names the account it will approve; approve it there and the token is stored.")
        return
    print(f"Slot {slot}: in your browser, sign in to claude.ai as the account this slot is for.")
    print("  If another account is signed in there, sign out first, or open the sign-in link in a private window.")
    print("  The page that opens names the account it will approve; approve it there and the token is stored.")


def add(gh, target, slot=None, label="", also_main=False, paste=False, write=None, stdin_fd=None, ask=None, noted=False, opener=None,
        copies=None):
    """One account into a slot. Every read and every scope is settled before
    the sign-in: a token minted and then not stored is a year's token lost
    with the process, so what can fail fails first."""
    refuse_crowd(gh, target)
    copies = machine_copies(copies)
    check_label(label)
    state = State(gh, target)
    slot = slot or first_free(state)
    name = slot_secret(slot)
    names = [name] + ([MAIN_SECRET] if also_main else [])
    scopes = {each: scope_args(state, "secrets", each) for each in names}
    label_scope = scope_args(state, "variables", slot_label(slot)) if label else None
    accounts_scope = scope_args(state, "variables", ACCOUNTS_VARIABLE)
    if name in state.secrets:
        print(f"Slot {slot} is set ({age_line(state.secrets[name])}); the new token replaces it.")
    if not noted:
        sign_in_note(slot, opener)
    token = get_token(argparse.Namespace(paste=paste), write=write, stdin_fd=stdin_fd, ask=ask, opener=opener)
    for each in names:
        store_secret(state, each, token, scopes[each])
        keep_copy(copies, target, each, token)
    token = None
    try:
        if label:
            set_variable(state, slot_label(slot), label, label_scope)
        state = verify(state, names)
        update_accounts(state, accounts_scope)
    except TokensError as error:
        raise TokensError(f"{error}; the token itself was sent to GitHub, so run status before storing another, and set the "
                          "label or the list of slots by hand on the variables page if status shows them missing") from None
    return state


def main_token(gh, target, paste=False, write=None, stdin_fd=None, ask=None, copies=None):
    refuse_crowd(gh, target)
    copies = machine_copies(copies)
    state = State(gh, target)
    scope = scope_args(state, "secrets", MAIN_SECRET)
    if MAIN_SECRET in state.secrets:
        print(f"{MAIN_SECRET} is set ({age_line(state.secrets[MAIN_SECRET])}); the new token replaces it.")
    print(f"{MAIN_SECRET}: in your browser, sign in to claude.ai as the account whose plan @claude, the workers in Actions and the coordinator workflow should draw on.")
    token = get_token(argparse.Namespace(paste=paste), write=write, stdin_fd=stdin_fd, ask=ask)
    store_secret(state, MAIN_SECRET, token, scope)
    keep_copy(copies, target, MAIN_SECRET, token)
    token = None
    return verify(state, [MAIN_SECRET])


def profile_opener(profile, browser=None):
    """add --profile: the profile, in the browser this machine has."""
    if not profile:
        return None
    check_profile(profile)
    browser = browser or find_browser()
    if browser is None:
        raise TokensError("--profile needs Chrome, Edge or Brave, and none was found where its installer puts it: "
                          "set SYMPHONY_BROWSER to the browser's program")
    return Opener(browser, profile)


def reorder(gh, target, text):
    """The list of slots in the order given: every slot that is set, each
    once. The bench deals its jobs to the slots in this order, one job to
    each in turn, so the first gets the first job and, when the jobs do
    not divide evenly, one more than the last."""
    state = State(gh, target)
    order = parse_slots(text)
    present = state.slots_set()
    if not present:
        raise TokensError("no slot is set, so there is nothing to order")
    if holds_other_words(text) or len(order) != len(list_words(text)) or sorted(order) != present:
        raise TokensError("the order names every slot that is set, each once: "
                          + (",".join(str(slot) for slot in present) or "none is set") + " in the order you want")
    scope = scope_args(state, "variables", ACCOUNTS_VARIABLE)
    set_variable(state, ACCOUNTS_VARIABLE, ",".join(str(slot) for slot in order), scope)
    print_status(State(gh, target))


def print_remembered(target, path=None, copies=None):
    """What this machine remembers for each slot of this place, and which
    tokens it keeps a copy of."""
    path = path or settings_path()
    copies = machine_copies(copies)
    settings = read_settings(path)
    print(f"Remembered on this machine, in {path}:")
    shown = False
    for slot in SLOTS:
        known = remembered(settings, target, slot)
        copied = has_copy(copies, target, slot_secret(slot))
        if known["profile"] or known["cswap"] or copied:
            shown = True
            parts = [f"browser profile {known['profile']}" if known["profile"] else "the default browser"]
            if known["cswap"]:
                parts.append(f"cswap account {known['cswap']}")
            if copied:
                parts.append("a copy of its token (show " + str(slot) + ")")
            print(f"  slot {slot}: " + ", ".join(parts))
    if has_copy(copies, target, MAIN_SECRET):
        shown = True
        print(f"  {MAIN_SECRET}: a copy of its token (show main)")
    if not shown:
        print("  nothing")


def forget(target, slot, path=None):
    """The slot's browser profile and cswap account, forgotten here; the
    browser's profile itself is the browser's to delete."""
    path = path or settings_path()
    settings = read_settings(path)
    remember(path, settings, target, slot, None, None)
    print(f"Slot {slot}: forgotten on this machine. Delete the browser profile itself in the browser, if you want it gone.")


def remove(gh, target, slot, copies=None):
    state = State(gh, target)
    name = slot_secret(slot)
    if name in state.secrets:
        delete_secret(state, name)
    else:
        print(f"  {name}: not set")
    drop_copy(machine_copies(copies), target, name)
    if slot_label(slot) in state.variables:
        delete_variable(state, slot_label(slot))
    update_accounts(State(gh, target))
    forget(target, slot)


def show(target, which, copies=None, clipboard=False):
    """A token back from this machine's copy: printed, or put on the
    clipboard, where it never reaches the screen or the terminal's
    scrollback."""
    slot = slot_number(which)
    if which != "main" and slot is None:
        raise UsageError(f"show takes a slot from 1 to {SLOTS[-1]}, or main")
    name = MAIN_SECRET if slot is None else slot_secret(slot)
    copies = machine_copies(copies)
    if copies is None:
        raise TokensError("copies are kept in Windows Credential Manager, and this machine is not Windows: "
                          "GitHub never gives a secret back, so mint a new token with add or main")
    token = copies.get(copy_name(target, name))
    if token is None:
        raise TokensError(f"this machine keeps no copy of {name} for {target.describe()}: a token stored before copies were kept "
                          "(2026-10-17), on another machine or with --no-copy has none; add --slot N, or main, mints one that keeps a copy")
    if clipboard:
        copy_to_clipboard(token)
        print(f"{name}: on the clipboard; paste it where you need it")
    else:
        print(token)
    token = None


def copy_to_clipboard(text):
    """Windows' clip reads what it puts on the clipboard from its standard input."""
    try:
        done = subprocess.run(["clip"], input=text, text=True, capture_output=True, timeout=30)
    except (OSError, subprocess.TimeoutExpired) as error:
        raise TokensError(f"clip could not be run ({error})") from None
    if done.returncode != 0:
        raise TokensError("clip did not take the token; run show without --copy to print it")


def share(gh, target, repos):
    """The slots' secrets, their labels and the list of slots seen by more
    repositories, such as a benchmark's, with no new sign-in: each is added
    to the selected list it already has. CLAUDE_CODE_OAUTH_TOKEN stays
    where it is, since no benchmark reads it."""
    if target.level != "org":
        raise TokensError("share adds repositories to an organization's secrets; a repository's own secrets are seen by it alone")
    refuse_crowd(gh, Target("org", target.owner, repos))
    state = State(gh, target)
    present = state.slots_set()
    if not present:
        raise TokensError("no slot is set, so there is nothing to share")
    ids = {}
    for repo in repos:
        found, text = gh.read(["api", f"repos/{target.owner}/{repo}", "--jq", ".id"])
        if not found or not re.fullmatch(r"[0-9]{1,20}", text.strip()):
            raise TokensError(f"GitHub lists no repository {target.owner}/{repo} for this login")
        ids[repo] = text.strip()
    names = [("secrets", slot_secret(slot), state.secrets[slot_secret(slot)]) for slot in present]
    names += [("variables", slot_label(slot), state.variables[slot_label(slot)]) for slot in present if slot_label(slot) in state.variables]
    if ACCOUNTS_VARIABLE in state.variables:
        names.append(("variables", ACCOUNTS_VARIABLE, state.variables[ACCOUNTS_VARIABLE]))
    for kind, name, item in names:
        if item.get("visibility") != "selected":
            print(f"  {name}: seen by every {'private ' if item.get('visibility') == 'private' else ''}repository already")
            continue
        for repo in repos:
            gh.must(["api", "-X", "PUT", f"{target.api}/{kind}/{name}/repositories/{ids[repo]}"])
        print(f"  {name}: now seen by {', '.join(repos)} too")
    print_status(State(gh, target))


def pool_module():
    """tools/pool.py, which keeps the readings' shape for the workflows too."""
    sys.path.insert(0, str(ROOT / "tools"))
    import pool
    return pool


def cswap_rows(run=None):
    """The accounts `cswap list --json` reports, each with its usage."""
    run = run or subprocess.run
    try:
        done = run(["cswap", "list", "--json"], capture_output=True, text=True, encoding="utf-8", timeout=120)
    except FileNotFoundError:
        raise TokensError("cswap is not installed here: usage reads its accounts' usage, which only cswap can") from None
    except subprocess.TimeoutExpired:
        raise TokensError("cswap list gave no answer in 120 seconds") from None
    try:
        data = json.loads(done.stdout or "")
    except json.JSONDecodeError:
        raise TokensError(f"cswap list --json printed no JSON (exit {done.returncode}); run it by hand to see why") from None
    rows = data.get("accounts") if isinstance(data, dict) else None
    if not isinstance(rows, list):
        raise TokensError("cswap list --json gave no list of accounts")
    return [row for row in rows if isinstance(row, dict)]


def usage_readings(target, rows, settings, pool):
    """One reading per slot whose cswap account the wizard remembers."""
    by_number = {row.get("number"): row for row in rows if type(row.get("number")) is int}
    readings, missing = [], []
    for slot in SLOTS:
        number = remembered(settings, target, slot)["cswap"]
        if number is None:
            continue
        record = pool.from_cswap(by_number.get(number, {}), str(slot))
        (readings if record else missing).append(record or (slot, number))
    return readings, missing


def send_usage(gh, repository, readings, pool, attempts=3):
    """The readings onto the branch pool through GitHub's git API, one
    commit on the branch's head, which fails and is tried again when a
    workflow wrote the branch in between."""
    files = [{"path": pool.file_name(record), "mode": "100644", "type": "blob",
              "content": json.dumps(record, indent=2, sort_keys=True) + "\n"} for record in readings]
    base = f"repos/{repository}/git"
    for _ in range(attempts):
        found, head = gh.read(["api", f"{base}/ref/heads/{pool.BRANCH}", "--jq", ".object.sha"])
        head = head.strip() if found else ""
        tree_body = {"tree": files}
        if head:
            tree_body["base_tree"] = gh.must(["api", f"{base}/commits/{head}", "--jq", ".tree.sha"]).strip()
        tree = gh.must(["api", "-X", "POST", f"{base}/trees", "--input", "-", "--jq", ".sha"], stdin=json.dumps(tree_body)).strip()
        commit_body = {"message": "Usage of the pool's slots, from cswap", "tree": tree, "parents": [head] if head else []}
        commit = gh.must(["api", "-X", "POST", f"{base}/commits", "--input", "-", "--jq", ".sha"], stdin=json.dumps(commit_body)).strip()
        if head:
            done = gh.run(["api", "-X", "PATCH", f"{base}/refs/heads/{pool.BRANCH}", "-f", f"sha={commit}", "-F", "force=false"])
        else:
            done = gh.run(["api", "-X", "POST", f"{base}/refs", "-f", f"ref=refs/heads/{pool.BRANCH}", "-f", f"sha={commit}"])
        if done.returncode == 0:
            return commit
        if "(HTTP 422)" not in (done.stderr or ""):
            raise TokensError(f"GitHub refused the write to the branch {pool.BRANCH}: {first_line(done.stderr) or 'no reason given'}")
    raise TokensError(f"the branch {pool.BRANCH} moved under each of {attempts} writes; run usage again")


def usage(gh, target, root=ROOT, path=None, run=None, also=()):
    """cswap's readings of each remembered account, sent to the branch pool
    of this clone's repository, and of each repository of the same owner
    named in also, such as coordbench, whose benchmark's auto reads its own
    branch; there the control panel and the auto mode read them; nothing
    but the slot, the windows and the per-model limits."""
    pool = pool_module()
    repository = origin_repository(root)
    if repository is None:
        raise TokensError("this clone has no GitHub origin, whose branch pool takes the readings")
    settings = read_settings(path or settings_path())
    readings, missing = usage_readings(target, cswap_rows(run), settings, pool)
    for slot, number in missing:
        print(f"  slot {slot}: cswap account {number} gave no usage; run cswap list to see why")
    if not readings:
        raise TokensError("no slot has a cswap account with usage: the wizard asks for each slot's cswap account when it adds one")
    for record in readings:
        parts = [f"{name.replace('_', ' ')} {round(each['utilization'] * 100)}%" for name, each in record["windows"].items()]
        parts += [f"{each['name']} {round(each['utilization'] * 100)}%" for each in record.get("scoped", [])]
        print(f"  slot {record['slot']}: " + (", ".join(parts) or "no windows"))
    for name in [repository[1]] + [each for each in also if each != repository[1]]:
        commit = send_usage(gh, f"{repository[0]}/{name}", readings, pool)
        print(f"Sent to the branch {pool.BRANCH} of {repository[0]}/{name} ({commit[:7]}).")


def chosen_slot(state, ask):
    """The slot to fill: the first free one by default, a set one to
    replace by its number, none to stop."""
    free = [slot for slot in SLOTS if slot_secret(slot) not in state.secrets]
    if free:
        answer = (ask(f"Which slot gets the next account? [{free[0]}; a set slot's number replaces it; n to stop] ") or "").strip().lower()
    else:
        answer = (ask("Every slot is set. Which one to replace? [its number, or Enter to stop] ") or "").strip().lower()
    if not answer:
        return free[0] if free else None
    if answer in ("n", "no", "q", "quit", "stop"):
        return None
    slot = slot_number(answer)
    if slot is not None:
        return slot
    raise TokensError(f"a slot is a number from 1 to {SLOTS[-1]}")


def yes(answer, default):
    text = (answer or "").strip().lower()
    return default if not text else text in ("y", "yes")


def chosen_label(account, ask):
    """The slot's label: a word typed, else the cswap alias where it is one."""
    offered = label_from(account)
    if offered:
        answer = (ask(f"A word of your own for this slot, so status says which account it is; anyone who reads the repository's "
                      f"workflows sees it (never an e-mail; y for cswap's {offered}; Enter for none): ") or "").strip()
        label = offered if answer.lower() in ("y", "yes") else answer
    else:
        label = (ask("A word of your own for this slot, so status says which account it is (never an e-mail; Enter for none): ") or "").strip()
    check_label(label)
    return label


def wizard(gh, target, ask=input, write=None, stdin_fd=None, local=None):
    """One account at a time: the state, then a slot, the cswap account, a
    label, the browser profile, the sign-in, the store, the check, until
    the person says no more."""
    refuse_crowd(gh, target)
    local = Local.discover() if local is None else local
    settings = read_settings(local.settings)
    state = State(gh, target)
    print_status(state)
    added = False
    while True:
        slot = chosen_slot(state, ask)
        if slot is None:
            break
        known = remembered(settings, target, slot)
        account = chosen_account(local.accounts, known["cswap"], ask)
        label = chosen_label(account, ask)
        profile = chosen_profile(local.browser, slot, known["profile"], ask)
        opener = Opener(local.browser, profile) if profile else None
        remember(local.settings, settings, target, slot, profile, account.number if account else None)
        sign_in_note(slot, opener)
        if opener is None:
            ask("Press Enter once the browser is signed in to that account ")
        also_main = MAIN_SECRET not in state.secrets and yes(
            ask(f"Store the same token as {MAIN_SECRET} too, the one @claude, the workers in Actions and the coordinator workflow run on? [y/N] "), False)
        state = add(gh, target, slot, label, also_main, write=write, stdin_fd=stdin_fd, noted=True, opener=opener)
        added = True
        if not yes(ask("Another account? [y/N] "), False):
            break
    # The list of slots is written after the check's read, so the closing
    # status reads everything again.
    print_status(State(gh, target) if added else state)


def add_place_options(parser, default):
    """Where the secrets go, accepted before the command and after it: the
    options after the command set nothing when absent, so the main parser's
    values stand."""
    parser.add_argument("--org", default=default, help="the organization the secrets go on, seen by --repos")
    parser.add_argument("--repos", default=default, help="repository names (without the owner, comma-separated) that see the organization's secrets")
    parser.add_argument("--visibility", default=default, choices=VISIBILITIES, help="for a new organization secret: selected (the default, --repos), all, or private")
    parser.add_argument("--repo", default=default, help="one repository, OWNER/NAME, that holds the secrets itself")


def build_parser():
    parser = argparse.ArgumentParser(
        prog="python3 scripts/tokens.py",
        description="Store your Claude accounts' subscription tokens as Actions secrets without seeing one; with no command, the wizard.",
    )
    add_place_options(parser, None)
    commands = parser.add_subparsers(dest="command")
    status = commands.add_parser("status", help="which slots are set, since when, and renew by when")
    adding = commands.add_parser("add", help="one account into a slot")
    adding.add_argument("--slot", type=int, choices=SLOTS, help="the slot, 1 to 8; by default the first free one")
    adding.add_argument("--label", default="", help="a word of your own for the slot, never an e-mail address")
    adding.add_argument("--also-main", action="store_true", help=f"store the same token as {MAIN_SECRET} too")
    adding.add_argument("--paste", action="store_true", help="take the token from a hidden prompt instead of running claude setup-token")
    adding.add_argument("--profile", help="open the sign-in in this profile of Chrome, Edge or Brave (its folder name, such as \"Profile 2\")")
    main_command = commands.add_parser("main", help=f"{MAIN_SECRET} alone, the token the Claude workflows run on")
    main_command.add_argument("--paste", action="store_true", help="take the token from a hidden prompt instead of running claude setup-token")
    for command in (adding, main_command):
        command.add_argument("--no-copy", action="store_true", help="keep no copy in Windows Credential Manager on this machine")
    removing = commands.add_parser("remove", help="a slot's secret and its label, and what this machine remembers of it")
    removing.add_argument("slot", type=int, choices=SLOTS)
    ordering = commands.add_parser("order", help="the order the bench deals its jobs to the slots, such as 3,1,2")
    ordering.add_argument("slots")
    forgetting = commands.add_parser("forget", help="a slot's browser profile and cswap account, on this machine only")
    forgetting.add_argument("slot", type=int, choices=SLOTS)
    profiles = commands.add_parser("profiles", help="what this machine remembers for each slot")
    showing = commands.add_parser("show", help="a token back from this machine's copy in Windows Credential Manager")
    showing.add_argument("which", help=f"a slot from 1 to {SLOTS[-1]}, or main")
    showing.add_argument("--copy", action="store_true", help="put it on the clipboard instead of printing it")
    using = commands.add_parser("usage", help="cswap's usage of each slot's account, sent to the branch pool")
    using.add_argument("--also", default="", help="more repositories of the same owner to send it to, comma-separated, such as coordbench")
    sharing = commands.add_parser("share", help="the slots' secrets seen by more repositories, such as coordbench")
    sharing.add_argument("repos", help="repository names without the owner, comma-separated")
    for command in (status, adding, main_command, removing, ordering, forgetting, profiles, showing, sharing, using):
        add_place_options(command, argparse.SUPPRESS)
    return parser


def show_target(args, gh):
    """The place a copy is kept under, with nothing asked of GitHub that the
    options already say: a copy is named by the organization alone, so show
    needs neither --repos nor a signed-in gh for it."""
    if args.org:
        return Target("org", args.org)
    if args.repo:
        owner, name = args.repo.split("/", 1)
        return Target("repo", owner, name=name)
    return target_from(args, gh)


def main_show(args):
    try:
        show(show_target(args, Gh()), args.which, clipboard=args.copy)
    except UsageError as error:
        print(f"tokens: {error}", file=sys.stderr)
        return 2
    except TokensError as error:
        print(f"tokens: {error}", file=sys.stderr)
        return 1
    return 0


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)
    if args.repo and not FULL_NAME.match(args.repo):
        parser.error("--repo takes OWNER/NAME")
    if args.repo and args.org:
        parser.error("--org or --repo, not both")
    if args.command == "show":
        return main_show(args)
    try:
        check_options(args)
    except UsageError as error:
        parser.error(str(error))
    gh = Gh()
    try:
        status = check_signed_in(gh)
        target = target_from(args, gh)
        if args.command not in ("forget", "profiles", "show", "usage"):
            check_scopes(status, target)
        if args.command == "status":
            print_status(State(gh, target))
        elif args.command == "add":
            add(gh, target, args.slot, args.label, args.also_main, args.paste, opener=profile_opener(args.profile),
                copies=False if args.no_copy else None)
        elif args.command == "main":
            main_token(gh, target, args.paste, copies=False if args.no_copy else None)
        elif args.command == "remove":
            remove(gh, target, args.slot)
        elif args.command == "order":
            reorder(gh, target, args.slots)
        elif args.command == "forget":
            forget(target, args.slot)
        elif args.command == "profiles":
            print_remembered(target)
        elif args.command == "usage":
            usage(gh, target, also=parse_repos(args.also) if args.also else ())
        elif args.command == "share":
            share(gh, target, parse_repos(args.repos))
        else:
            wizard(gh, target)
    except TokensError as error:
        print(f"tokens: {error}", file=sys.stderr)
        return 1
    except (EOFError, KeyboardInterrupt):
        print("tokens: stopped", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
