"""The protocol over GitHub's REST API, through urllib.

Four parts of the package read GitHub, and they read it the same way, so
the calls live once, here. The board collects a Snapshot, posts the
comments board.comments_due names and renames the branch of each claim
board.released_keys counts as released; the coordinator, in the same run,
releases the claims coordinator.releases_due names and tells the open work
that meets other work in a file; the team clerk calls the few methods
team.Client lists; and the merge gate asks GitHubWorld, gate.World over one
repository, how a branch stands, whether its required check passed, and
for the merge itself, which goes through GitHub's merge API and only when
the policy says the gate merges.

The token is an argument or GITHUB_TOKEN: a workflow's own token and an
App installation token are both sent as a bearer token. It travels once
per call, in the Authorization header, and never in a message: an error
names the method, the path and the status, and keeps GitHub's own words
apart, on the exception, for a caller that wants to say why.
"""

import base64
import collections
import dataclasses
import http.client
import json
import os
import time
import urllib.error
import urllib.parse
import urllib.request
from datetime import datetime, timedelta, timezone
from typing import NamedTuple

from .. import board, claims, coordinator, gate, labels, team
from ..gate import GateFailed, MainMoved, MergeRefused, MergeResult, NotAllowed

API = "https://api.github.com"
# Where an owner changes a ruleset, for the board's comment on a rename a
# ruleset refused.
WEB = "https://github.com"
# The answers GitHub gives a rename it refuses for good: a rule on the
# branch, or a token without the right. A 404 can mean only that the branch
# went away since it was read, which the next run sees for itself.
REFUSED_FOR_GOOD = {403, 422}
# What GitHub answers for an issue someone deleted.
GONE = 410
API_VERSION = "2022-11-28"
USER_AGENT = "sdlc-board"
TIMEOUT = 60
PER_PAGE = 100
# The branch creations one read of the repository's activity lists, the
# most GitHub gives on a page (began_times).
ACTIVITY_PAGE = 100
# The calls the board keeps in reserve for the writes a run makes after
# its reads, which the Dependencies section's reads never eat into: one
# read per open epic and goal, one per sub-issue with a blocker, and the
# epic bodies written after the page, against a token's thousand an hour
# (within_budget). A section left out this run is drawn the next.
DEPENDENCIES_RESERVE = 300
# What the Overlaps do with the work of a branch GitHub refused to date
# (began_refused_lines).
UNDATED_PULL = "date its pull request by when it was opened"
UNDATED_CLAIM = "put its claim last, as work nothing dates"
UNDATED_STRANGER = ("date a stranger's pull request from it by when it was opened and put its claim last, as work "
                    "nothing dates")
# GitHub computes mergeability after the first ask; one more ask after a
# pause is enough, and "unknown" is the honest word when it still has not.
MERGEABLE_RETRY_SECONDS = 3
# The repository's activity may not list a push the moment GitHub takes it:
# one more read after a pause before a claim is dated by its committer
# date (head_since, F8).
PUSH_RETRY_SECONDS = 3
NO_PUSH_LISTED = "the repository's activity lists no push that set the branch's head"
# In the run a branch's push started, an entry that set its head more than
# PUSH_SINCE_SLACK before the run began dates nothing: an older push of
# that head, the run's own not listed yet, or the run's own where the run
# began that much after it, queued behind other runs; either way the run's
# start dates the head, and the line says only what is so of both
# (head_since).
OLDER_PUSH_LISTED = "the newest push the repository's activity lists of the branch's head came before this run began"
# An entry of a branch's activity that names another head dates the claim
# only when it came after the run began, a push after the run read the
# branch; this much earlier still counts, for a runner's clock ahead of
# GitHub's. Dating by it can only hold a claim, never release one. In the
# run a branch's push started, an entry that set the head counts as that
# push only when it came as late as this too (OLDER_PUSH_LISTED).
PUSH_SINCE_SLACK = timedelta(minutes=5)
RECENT_MERGED = board.RECENT_MERGED
# The headers GitHub answers a REST call with that say what the token's
# rate limit has left: the calls left, the calls an hour allows and when
# the hour ends, in seconds since the epoch. A workflow's GITHUB_TOKEN is
# allowed 1,000 calls an hour per repository, shared by every workflow of
# the repository that uses it, the board's among them (calls_line).
LIMIT_HEADERS = ("x-ratelimit-remaining", "x-ratelimit-limit", "x-ratelimit-reset")
# The most digits a header of the rate limit is read with: a count or a
# time GitHub writes never comes near it, a time of ten digits is before
# the year 2287, which datetime writes, and int() of a longer string
# spends time no answer is worth (count_in).
LIMIT_DIGITS = 10
# GitHub refuses a call for a rate limit with a 403 or a 429, as its own
# documentation says ("Rate limits for the REST API", "Exceeding the rate
# limit"): its primary limit with x-ratelimit-remaining at 0, to be asked
# again only after x-ratelimit-reset; a secondary limit with a message that
# names it, to be asked again after retry-after's seconds where the answer
# gives them, else after at least a minute and longer each time, giving up
# after a few tries. The client waits as GitHub asks only while the wait
# is short (rate_wait): a run that waited longer would hold up the runs
# queued behind it, and the next run reads everything again anyway, so a
# longer wait fails the call as any other refusal does. A client's waits
# are bounded in all too (RATE_WAIT_BUDGET), at three minutes, so a job
# with a limit of five minutes, the shortest a consumer gives a job that
# runs the client, ends with GitHub's own words and not by its limit.
RATE_STATUSES = (403, 429)
RATE_WORDS = "rate limit"
RATE_RETRIES = 2
RATE_FIRST_WAIT = 60
RATE_LONGEST_WAIT = 120
RATE_WAIT_BUDGET = 180

# GitHub's mergeable_state words that say only that something besides the
# branch stops the merge, a rule, a review or the draft itself: how the
# branch stands to its base is then read from mergeable and a comparison.
OPAQUE_STATES = {"blocked", "draft"}
# The permissions whose holders' release comments count when GitHub's
# author_association does not vouch for them: write and above, in the old
# permission word and in the role name.
WRITE_PERMISSIONS = {"admin", "write"}
WRITE_ROLES = {"admin", "maintain", "write"}
# A required check run that finished with one of these let the merge go;
# any other finish (failure, timed_out, cancelled, action_required, stale)
# did not.
PASSING = {"success", "neutral", "skipped"}


class GitHubError(team.ClientError):
    """A call GitHub refused or did not answer. The message names the
    method, the path and the status, never the token; status is the HTTP
    status (None when there was no answer), detail GitHub's own message
    from the answer's body, "" when it gave none, and headers the
    answer's headers, None when there was no answer, for the wait a rate
    limit asks (rate_wait); after_rate_limit says the client made the
    call again after GitHub refused an earlier try of it for a rate limit
    (GitHub.request), so the refusal is the last of several answers and
    may follow a try GitHub acted on (wrote_nothing)."""

    def __init__(self, message, status=None, detail="", headers=None):
        super().__init__(message)
        self.status = status
        self.detail = detail
        self.headers = headers
        self.after_rate_limit = False


def next_link(headers):
    """The URL of the next page named by a Link header, or None."""
    link = header(headers, "link") or ""
    for part in link.split(","):
        pieces = part.split(";")
        if len(pieces) >= 2 and pieces[1].strip() == 'rel="next"':
            return pieces[0].strip().strip("<>")
    return None


def error_detail(problem):
    """GitHub's message from a refused call's JSON body, or "" when the
    body is empty, unreadable or not GitHub's shape. A body cut short
    raises http.client's IncompleteRead, which is no OSError, when it is
    read, and GitHub.request reads it inside its handler of the refusal,
    whose neighbours never catch what a handler raises: only this except
    keeps it from ending the run."""
    try:
        body = json.loads(problem.read() or b"null")
    except (AttributeError, OSError, TypeError, ValueError, http.client.HTTPException):
        return ""
    return str(body.get("message") or "") if isinstance(body, dict) else ""


def count_in(text):
    """The whole number a header of the rate limit says, or None for one
    missing or written otherwise: ASCII digits, at most LIMIT_DIGITS of
    them, around which GitHub may leave spaces."""
    text = text.strip() if isinstance(text, str) else ""
    if not 0 < len(text) <= LIMIT_DIGITS or not (text.isascii() and text.isdigit()):
        return None
    return int(text)


def header(headers, name):
    """One header of an answer by its lowercase name, read in any case as
    HTTP reads names, or None where the answer had none: urllib gives an
    HTTPError's headers as None when it has no answer to read them from,
    and a stand-in may give them as a plain mapping, whose names keep the
    case they were written in."""
    if headers is None:
        return None
    return next((value for key, value in headers.items() if isinstance(key, str) and key.lower() == name), None)


def rate_limited(status, headers, detail):
    """Whether GitHub refused a call for a rate limit, a refusal it asks to
    have waited out (rate_wait) and that says nothing against the call
    itself: a 429, Too Many Requests; or a 403 that says it is one, by a
    retry-after header or by a message that names a rate limit, as
    GitHub's refusals for its primary limit ("API rate limit exceeded")
    and its secondary limits do. That retry-after tells one whatever the
    message says is this project's assumption (validation/DESIGN.md,
    3.2): GitHub's documentation names it for a secondary limit's refusal
    and never says that no other refusal sends it. A 403 that says
    neither is GitHub's answer to a token without the right or to a rule,
    whatever its x-ratelimit-remaining says: every answer counts against
    the hour, so the hour's last call is answered with none left,
    whatever GitHub answers it with."""
    if status not in RATE_STATUSES:
        return False
    return status == 429 or header(headers, "retry-after") is not None or RATE_WORDS in (detail or "").lower()


def rate_wait(status, headers, detail, tries, now):
    """The seconds to wait before a call GitHub refused is made again, the
    call having been made again tries times already: for a refusal of a
    rate limit's (rate_limited), the longest of the waits GitHub asks
    for, the seconds retry-after names and, where the token has no call
    left, until the hour's end, since each holds whatever the other says
    ("Exceeding the rate limit"), or, where it asks for neither, a minute
    and then two; None for any other refusal, for a spent hour whose end
    the answer does not say, once the tries are spent (RATE_RETRIES), or
    where GitHub asks for a longer wait than a run waits
    (RATE_LONGEST_WAIT). now is the clock's time in seconds since the
    epoch, which x-ratelimit-reset is written in."""
    if tries >= RATE_RETRIES or not rate_limited(status, headers, detail):
        return None
    asked = count_in(header(headers, "retry-after"))
    left = count_in(header(headers, "x-ratelimit-remaining"))
    reset = count_in(header(headers, "x-ratelimit-reset"))
    if left == 0 and reset is None:
        return None
    waits = [] if asked is None else [asked]
    if left == 0:
        waits.append(max(reset - int(now), 0) + 1)
    wait = max(waits) if waits else RATE_FIRST_WAIT * 2 ** tries
    return wait if wait <= RATE_LONGEST_WAIT else None


def many(count, word):
    return f"{count} {word}{'' if count == 1 else 's'}"


def calls_line(api):
    """The run's last log line: how many of the run's calls to GitHub's
    REST API GitHub answered, its reads and its writes, and what the
    token's rate limit had left by the latest answer that said it, so
    whoever reads the run's log sees what a board run costs and how near
    the hour's limit it ran (LIMIT_HEADERS), and, where it waited for a
    rate limit, how many times and how long in all, each wait one call
    made again, while a refusal it did not wait out is a call like any
    other; a stand-in client that keeps no waits, as one written before
    the client waited, waited for none."""
    reads = api.made["GET"]
    writes = sum(api.made.values()) - reads
    line = f"GitHub's REST API: {many(reads + writes, 'call')} in this run, {many(reads, 'read')} and " \
           f"{many(writes, 'write')}"
    waited = getattr(api, "waited", None) or ()
    if waited:
        line += f"; it waited for GitHub's rate limits {many(len(waited), 'time')}, {sum(waited)} s in all"
    if api.limit is None:
        return line
    left, allowed, reset = api.limit
    until = datetime.fromtimestamp(reset, timezone.utc).strftime("%H:%M UTC")
    return f"{line}; GitHub last said the token had {left} of its {allowed} calls an hour left, the hour " \
           f"ending at {until}"


def quoted(name):
    """A branch or label name made safe for a URL path, slashes kept, since
    GitHub reads a ref's slashes as part of its name."""
    return urllib.parse.quote(name, safe="/")


def unsent(call, token):
    """Why a call was refused before it was sent, in words that quote no
    header, since one holds the token. http.client writes the request line
    in ASCII and refuses any other path, and every path here quotes what a
    person names (quoted), so such a path is a call built wrong; it writes a
    header in Latin-1 with no bare line break, and the token is the one
    header whose value comes from outside the adapter, so a token pasted
    with a curly quote or a newline is said as the token. Anything else
    urllib refuses to send, a proxy setting it cannot read for one, is said
    in no more words than that, since urllib's own may quote a password."""
    if not call.selector.isascii():
        return "the path holds a character outside ASCII, which http.client does not send"
    if not (token.isascii() and token.isprintable()):
        return ("the token holds a character http.client does not send in a header, such as a curly quote or a "
                "line break")
    return "urllib would not send the call"


class GitHub:
    """The REST calls the package makes, through urllib, and the methods
    team.Client lists.

    The opener, the sleep and the clock are arguments so a test can hand in
    fakes; the token is the argument, else GITHUB_TOKEN, and with neither
    there is no client, since every call would be refused.
    """

    def __init__(self, token, repo, opener=None, sleep=time.sleep, clock=time.time):
        self.token = token or os.environ.get("GITHUB_TOKEN")
        if not self.token:
            raise GitHubError("no token: pass one, or set GITHUB_TOKEN")
        self.repo = repo
        self.opener = opener or urllib.request.build_opener()
        self.sleep = sleep
        self.clock = clock
        # The seconds of each wait the client made for a rate limit, in
        # order, for the run's log (calls_line) and the run's bound
        # (RATE_WAIT_BUDGET).
        self.waited = []
        self.asked = permissions(self)
        # The calls GitHub answered this client, by method, and the rate
        # limit as the latest answer that said it left it: the calls left,
        # the calls an hour allows and when the hour ends (LIMIT_HEADERS),
        # for the run's log (calls_line).
        self.made = collections.Counter()
        self.limit = None

    def headers(self):
        return {
            "Authorization": f"Bearer {self.token}",
            "Accept": "application/vnd.github+json",
            "X-GitHub-Api-Version": API_VERSION,
            "User-Agent": USER_AGENT,
        }

    def url(self, path, params=None):
        return API + path + ("?" + urllib.parse.urlencode(params) if params else "")

    def request(self, method, url, body=None, missing_ok=False):
        """One call, made again after the wait GitHub asks where it refuses
        the call for a rate limit and the wait is one a run takes
        (wait_before_again); the decoded JSON and the answer's headers, or
        the last refusal's GitHubError. A write is made again as a read is,
        on an assumption: that GitHub refuses a rate limit's call before it
        acts on it, which its documentation does not say, and which the
        validation harness's fake assumes until a canary check settles it
        (validation/DESIGN.md, 3.2 and section 11, C4). The labels a
        release puts back after its comment failed do not rest on it
        (wrote_nothing), so the refusal raised says whether an earlier try
        was refused for a rate limit (GitHubError.after_rate_limit), since
        only the last try's answer is raised, and a 422 after a rate
        limit's refusal says nothing of what that earlier try did (F23's
        second review)."""
        tries = 0
        while True:
            try:
                return self.send(method, url, body, missing_ok)
            except GitHubError as problem:
                wait = self.wait_before_again(problem, tries)
                if wait is None:
                    problem.after_rate_limit = tries > 0
                    raise
            self.waited.append(wait)
            self.sleep(wait)
            tries += 1

    def wait_before_again(self, problem, tries):
        """The wait rate_wait asks before the call is made again, or None
        where it asks none or the run has waited as long as it waits in all
        (RATE_WAIT_BUDGET)."""
        wait = rate_wait(problem.status, problem.headers, problem.detail, tries, self.clock())
        if wait is None or sum(self.waited) + wait > RATE_WAIT_BUDGET:
            return None
        return wait

    def send(self, method, url, body=None, missing_ok=False):
        """One call, made once; the decoded JSON and the headers. A 404 is
        (None, headers) when missing_ok says absence is an answer. Every
        failure of the call, whether GitHub, the connection or http.client
        refused it, is a GitHubError, whichever part of urllib raised it;
        the URL is always the adapter's own or the next page a Link header
        names, which has the scheme urllib asks of it before trying
        anything. An answer cut short (IncompleteRead) and one that is no
        HTTP (BadStatusLine) are http.client's own exceptions, no OSError,
        which urllib passes on as they are, past every caller that turns a
        refusal into a line. A call http.client will not send is said by
        what holds the character it refuses (unsent), and the answer is read
        as JSON only once the call is done, so a call never sent is never
        said to have had an answer that was not JSON."""
        data = json.dumps(body).encode("utf-8") if body is not None else None
        headers = self.headers()
        if data is not None:
            headers["Content-Type"] = "application/json"
        call = urllib.request.Request(url, data=data, method=method, headers=headers)
        where = f"{method} {urllib.parse.urlsplit(url).path}"
        try:
            with self.opener.open(call, timeout=TIMEOUT) as answer:
                self.count_answer(method, answer.headers)
                raw, answered = answer.read(), answer.headers
        except urllib.error.HTTPError as problem:
            self.count_answer(method, problem.headers)
            if problem.code == 404 and missing_ok:
                return None, problem.headers
            raise GitHubError(f"{where}: HTTP {problem.code}", problem.code, error_detail(problem),
                              problem.headers) from None
        except OSError as problem:
            raise GitHubError(f"{where}: {getattr(problem, 'reason', problem)}") from None
        except http.client.HTTPException as problem:
            raise GitHubError(f"{where}: {problem!r}") from None
        except ValueError:
            raise GitHubError(f"{where}: {unsent(call, self.token)}") from None
        try:
            return (json.loads(raw) if raw else None), answered
        except ValueError:
            raise GitHubError(f"{where}: the answer was not JSON") from None

    def count_answer(self, method, headers):
        """Count a call GitHub answered, whatever it answered, and keep the
        rate limit its headers say, where they say all of it, for the run's
        log (calls_line). A call never sent, or one the connection lost
        before any answer, is not counted, as GitHub said nothing of it; a
        redirect urllib follows, which GitHub answers a renamed repository's
        old name with, is two calls counted as one. The headers are read
        in any case (header), as rate_wait reads them."""
        self.made[method] += 1
        limit = tuple(count_in(header(headers, name)) for name in LIMIT_HEADERS)
        if None not in limit:
            self.limit = limit

    def get(self, path, params=None, missing_ok=False):
        data, _ = self.request("GET", self.url(path, params), missing_ok=missing_ok)
        return data

    def pages(self, path, params=None):
        """A list endpoint page by page, following the Link header, so a
        caller that has seen enough can stop asking."""
        params = dict(params or {})
        params["per_page"] = PER_PAGE
        url = self.url(path, params)
        while url:
            data, headers = self.request("GET", url)
            yield data or []
            url = next_link(headers)

    def get_all(self, path, params=None):
        """Every item of a list endpoint."""
        return [item for page in self.pages(path, params) for item in page]

    def post(self, path, body):
        data, _ = self.request("POST", self.url(path), body=body)
        return data

    def put(self, path, body):
        data, _ = self.request("PUT", self.url(path), body=body)
        return data

    def patch(self, path, body):
        data, _ = self.request("PATCH", self.url(path), body=body)
        return data

    def delete(self, path, missing_ok=False):
        data, _ = self.request("DELETE", self.url(path), missing_ok=missing_ok)
        return data

    # --- team.Client ---------------------------------------------------------

    def issue(self, number):
        return self.get(f"/repos/{self.repo}/issues/{number}")

    def open_issues(self):
        """Every open issue, with the pull requests GitHub mixes into the
        list left in, as the API gives them; the clerk and the board each
        leave them out."""
        return self.get_all(f"/repos/{self.repo}/issues", {"state": "open"})

    def branch_names(self):
        return [branch["name"] for branch in self.get_all(f"/repos/{self.repo}/branches")]

    def open_pull(self, branch):
        """The open pull request from a branch of this repository, as
        pull_detail gives it, with its standing where GitHub's word says
        only that something else stops the merge (standing_of), so the
        clerk deals a claim back for a branch behind its base, as the board
        tells one behind the default branch "main moved"; or None. A draft
        is not compared, since the clerk carries a draft whatever its
        standing. No default branch is named, so a pull request stacked on
        another branch is compared with that base, which the clerk reads as
        a base it has fallen behind and the board leaves to GitHub's word;
        the template's rulesets hold only the default branch, so GitHub
        never calls such a pull request blocked there."""
        detail = self.pull_detail(branch)
        if detail and not detail.get("draft"):
            standing = standing_of(self, detail)
            if standing is not None:
                detail["standing"] = standing
        return detail

    def pull_detail(self, branch):
        """The open pull request from a branch of this repository, in the
        detail that holds mergeable_state and the head's sha, as GitHub
        gives it, or None. The list's head filter names the branch by its
        owner, which a fork in the same organization shares (github/docs,
        content/pull-requests/reference/forks.md:113), so a pull request
        from such a fork is left out by the head and base the list gives
        (board.from_fork; F17's review); and so, from a claim's branch, is
        one that someone who may not write opened from the branch here,
        which is no claim's either (claims.opened_by_writer, may_open;
        F19), a filter for claim branches alone, as mark_strangers's is,
        since another branch's pull request is no claim's whoever opened
        it. The gate's world reads a branch's pull request through this too
        (GitHubWorld.pull_of, merge), so a stranger's from a claim's branch
        has no standing and is never merged."""
        owner = self.repo.split("/")[0]
        found = self.get(f"/repos/{self.repo}/pulls", {"state": "open", "head": f"{owner}:{branch}"}) or []
        own = [pull for pull in found
               if not board.from_fork(pull)
               and (not claims.is_claim_branch(branch) or claims.opened_by_writer(pull, may_open(self.asked)))]
        return self.get(f"/repos/{self.repo}/pulls/{own[0]['number']}") if own else None

    def check_runs(self, sha):
        """Every check run on one commit, page by page: the endpoint wraps
        each page's runs in an object, and a busy commit has more than one
        page of them."""
        runs = []
        for page in self.pages(f"/repos/{self.repo}/commits/{sha}/check-runs"):
            if isinstance(page, dict):
                runs.extend(page.get("check_runs") or [])
        return runs

    def claim_text(self, branch):
        """The claim file on a claude/<issue> branch, as the board reads it
        (claim_text below): the clerk reads whose a claim is from it (F18)."""
        return claim_text(self, branch)

    def released(self, branch):
        """Whether a release comment counts for the claim on a claude/<issue>
        branch (claim_released below): the clerk carries a released claim
        for nobody (F18's review)."""
        return claim_released(self, branch)

    def create_label(self, name, colour):
        """Make the label; GitHub answers 422 for one that exists, which is
        the label being there, not a failure."""
        try:
            self.post(f"/repos/{self.repo}/labels", {"name": name, "color": colour})
        except GitHubError as problem:
            if problem.status != 422:
                raise

    def create_issue(self, title, body, labels):
        return self.post(f"/repos/{self.repo}/issues", {"title": title, "body": body, "labels": labels})

    def add_sub_issue(self, parent, child_id):
        self.post(f"/repos/{self.repo}/issues/{parent}/sub_issues", {"sub_issue_id": child_id})

    def add_blocked_by(self, number, blocker_id):
        """Record GitHub's blocked-by relation: the issue with the id
        blocks the issue with the number, as the issue's page shows it
        under Relationships and its dependency summary counts it."""
        self.post(f"/repos/{self.repo}/issues/{number}/dependencies/blocked_by", {"issue_id": blocker_id})

    def set_body(self, number, body):
        self.patch(f"/repos/{self.repo}/issues/{number}", {"body": body})

    def remove_label(self, number, name):
        """Take a label off an issue; GitHub's 404 for one the issue no
        longer carries is the label being off, not a failure."""
        self.delete(f"/repos/{self.repo}/issues/{number}/labels/{urllib.parse.quote(name, safe='')}",
                    missing_ok=True)


# --- the board ----------------------------------------------------------------

def settled(api, detail):
    """A pull request's detail with mergeable_state settled: asked once
    more after a pause when GitHub had not computed it, and "unknown" when
    it still has not."""
    if detail.get("mergeable") is None:
        api.sleep(MERGEABLE_RETRY_SECONDS)
        detail = api.get(f"/repos/{api.repo}/pulls/{detail['number']}")
    if detail.get("mergeable") is None:
        detail["mergeable_state"] = "unknown"
    return detail


def pull_with_detail(api, number, default_branch=None):
    """One pull request with its mergeable_state settled, its standing read
    where that word is opaque, and its file paths."""
    path = f"/repos/{api.repo}/pulls/{number}"
    detail = settled(api, api.get(path))
    standing = standing_of(api, detail, default_branch)
    if standing is not None:
        detail["standing"] = standing
    detail["files"] = [f["filename"] for f in api.get_all(path + "/files")]
    return detail


def standing_of(api, detail, default_branch=None):
    """How the pull request's branch stands to its base, in the board's
    words, where GitHub's mergeable_state says only that something else
    stops the merge (a rule, a review, the draft itself): mergeable tells
    a conflict apart, and the comparison tells behind from clean, as the
    gate's world reads such a pull request (F9: until this the board read
    blocked as blocked and draft as unknown, and neither heard "main
    moved"). None where GitHub's word is one the board reads as it is:
    behind and dirty say how the branch stands, unknown is GitHub still
    computing, and clean, unstable and has_hooks say the repository's
    rules do not hold the pull request for being behind, a draft's as any
    other's, since GitHub gives a draft the word its branch has (C0), which
    they say of a branch that is behind only under a rule that does not
    require an up-to-date branch, or under none, so the spec owes it
    nothing however far the default branch has moved (F16). None too for a pull request
    stacked on a branch other than the default branch the caller names,
    since "main moved" speaks of the default branch. A
    comparison GitHub refuses, a 5xx, a 422 or a rate limit's 403, reads
    as unknown, which keeps GitHub's own word, as the coordinator's own
    comparison does (holding): one pull request GitHub will not compare
    must not fail the board's whole run."""
    if detail.get("mergeable_state") not in OPAQUE_STATES:
        return None
    base = (detail.get("base") or {}).get("ref")
    if default_branch and base and base != default_branch:
        return None
    if not detail.get("mergeable"):
        return standing_word(detail)
    try:
        behind = behind_by(api, detail, default_branch)
    except GitHubError:
        behind = None
    return standing_word(detail, behind)


def behind_by(api, pull, default_branch=None):
    """How many commits of the pull request's base its head lacks, from
    GitHub's comparison of the two; None when there is nothing to compare
    or GitHub has no such comparison (404). Any other refusal raises."""
    base = (pull.get("base") or {}).get("ref") or default_branch
    head = (pull.get("head") or {}).get("sha")
    if not head or not base:
        return None
    found = api.get(f"/repos/{api.repo}/compare/{quoted(base)}...{head}", missing_ok=True)
    return (found or {}).get("behind_by")


def claim_text(api, branch):
    """The claim file on a claude/<issue> branch, or None when there is none
    yet, or when what stands at its path is no file: GitHub answers a folder
    with the list of what it holds, and a submodule, or a link to anything
    but a file, with an object of that type, and reading a list as a file
    raised in every run while the branch stood. A file GitHub gives no text
    of, one over 1 MB, which it answers with an empty content and the
    encoding "none", is board.UNREAD_CLAIM, never an empty claim file. The
    key is quoted, as a hand-made branch may put anything after claude/:
    http.client refuses a path it cannot send as ASCII, and a # would begin
    the URL's fragment, which asks for the claims folder."""
    name = board.branch_name(branch)
    found = api.get(
        f"/repos/{api.repo}/contents/{claims.CLAIM_DIR}/{quoted(claims.key_of(name))}.md",
        {"ref": name}, missing_ok=True,
    )
    if not isinstance(found, dict) or found.get("type", "file") != "file":
        return None
    if found.get("encoding") == "none":
        return board.UNREAD_CLAIM
    return base64.b64decode(found.get("content") or "").decode("utf-8", "replace")


def unread_claim_lines(snapshot):
    """A line for each claim file GitHub gave no text of, for the workflow's
    log, since the page says so only where no label names the worker: where
    one does, as the protocol's Announce step makes it, the row shows that
    worker and no files, as for a claim file that lists none."""
    snap = board.as_snapshot(snapshot)
    return [f"{claims.branch_of(key)}: GitHub gave no text for {claims.claim_path(key)}, a file over 1 MB, so the "
            f"board reads no worker and no files from it, and a release quotes no plan"
            for key, text in sorted(snap.claims.items(), key=lambda item: claims.by_number(item[0]))
            if text is board.UNREAD_CLAIM]


def recent_merged(api):
    """The pull requests merged last, newest first, at most RECENT_MERGED,
    from the first page of closed pull requests by last update (a hundred).
    The list depends on the repository alone, so two runs in a row render
    the same one: a merge fires two runs (the pull request closing and the
    push to the default branch), and a list measured from the previous
    board's commit lost, in the second, the merge the first had reported."""
    for page in api.pages(f"/repos/{api.repo}/pulls", {"state": "closed", "sort": "updated", "direction": "desc"}):
        merged = [p for p in page if p.get("merged_at")]
        merged.sort(key=lambda p: p["merged_at"], reverse=True)
        return merged[:RECENT_MERGED]
    return []


def claim_issues(api, branches, open_issues):
    """The issue behind each claim branch, by its key: the state and the
    label names, and GitHub's state_reason for one reopened, which
    may_be_due reads. An open issue comes from the list already read; any other
    takes one call, so the board can say that a branch outlived its issue.
    A branch whose key is no number (claims.number_in) has none and costs
    no call: a key such as a superscript would put a path in the URL that
    http.client refuses, and that refusal stopped the whole run."""
    by_number = {str(i["number"]): i for i in open_issues}
    found = {}
    for branch in board.claim_branches(branches):
        key = claims.key_of(branch)
        if claims.number_of(branch) is None:
            continue
        issue = by_number.get(key) or issue_or_none(api, key)
        if issue is not None:
            found[key] = {"state": issue.get("state", ""), "labels": labels.label_names(issue)} \
                | ({"state_reason": board.REOPENED} if issue.get("state_reason") == board.REOPENED else {})
    return found


def issue_or_none(api, number):
    """The issue with this number, or None when GitHub has none by it: 404
    for a number no issue was given, and 410 for an issue someone deleted,
    whose branch the template's tidy workflow, which runs when an issue
    closes, leaves standing; the 410 raised, and stopped every run while
    the branch stood."""
    try:
        return api.get(f"/repos/{api.repo}/issues/{number}", missing_ok=True)
    except GitHubError as problem:
        if problem.status == GONE:
            return None
        raise


def kept_names(listed):
    """The branches a release left standing, released/<issue>-<short sha>,
    by name, from the branch list as the API gives it."""
    return [name for name in (board.branch_name(b) for b in listed) if name.startswith(claims.RELEASED_PREFIX)]


def ended_releases(api, issues, branches, kept, asked):
    """The reopened issues, by key, whose released branches a release of a
    claim before the issue's last close left: an open issue GitHub says
    was reopened that carries claimed or in-review, with no claim branch
    and a branch released/<issue>-<short sha> standing, whose events
    record a close, with no label claimed put on after it but by a claim
    made since (claimed_by_a_claim), and whose comments hold no release
    comment after it. Such a branch is no
    release being undone, so the coordinator frees the issue wholly
    (coordinator.claim_stands; F22).

    The label claimed put on after the last close is a person undoing the
    release, as the release's note on the issue asks (board.released_text),
    unless the newest entry the repository's activity lists for the claim
    branch claude/<issue>, of any kind, came at or after that close and no
    earlier than the label: a claim made since puts its own claimed on
    after its branch's making, and the branch, gone now, was deleted or
    pushed after that, while an undo puts claimed back with no push to the
    branch after it, whether no claim came between or one came and ended
    before the person put claimed back (claimed_by_a_claim; F22's third
    review and its review). An entry in the label's own second reads as
    the claim's. A claim the activity does not list, one GitHub no longer
    keeps the entries of or one whose entries a rename took, reads as an
    undo: an empty list is what the ordinary undo, with no claim between,
    shows, so reading it as a claim's would take every such person's
    claimed off, while reading the unlisted claim as an undo only leaves
    claimed on, with in-review off, for a person to take off. Claimed that
    came through the close is the ended claim's (F22's second review).
    Where GitHub refuses the activity, keeps none of the repository (a
    404), or writes no time for the label or the entry, the issue is left
    out this run. A release comment after the close counts unless GitHub
    says its author may not write (may_open, from asked, the run's answers
    to who may write), so a refused permission read leaves the issue out
    this run, as a refused read of the events or the comments does, rather
    than taking a writer's release for a stranger's and the labels off on
    a guess. One left out keeps its released branch as a release being
    undone, and the freeing takes in-review off alone, as before this was
    read, until a later run reads it.

    Only such an issue costs calls: every page of its events; one read of
    the claim branch's newest entry in the activity where claimed went on
    after the close; and every page of its comments where the events and the
    activity leave the question open, with a permission asked, once a
    run, for each author of a release comment after the close that
    GitHub's word does not vouch for. A freed issue carries neither label
    and costs nothing again; one whose release is being undone costs its
    events and the activity read in every board run until a person
    renames the branch back, as the note says, which nothing bounds
    without a record kept from run to run (validation/DESIGN.md, F22)."""
    found = []
    for issue in issues:
        key = str(issue.get("number", ""))
        if not reopened_with_released_branch(issue, key, branches, kept):
            continue
        try:
            events = api.get_all(f"/repos/{api.repo}/issues/{key}/events")
            closed = latest_close(events)
            if not closed or (claimed_since_close(events)
                              and not claimed_by_a_claim(api, key, closed, last_claimed(events))):
                continue
            comments = api.get_all(f"/repos/{api.repo}/issues/{key}/comments")
        except GitHubError:
            continue
        if not released_after(comments, closed, may_open(asked)):
            found.append(key)
    return found


def reopened_with_released_branch(issue, key, branches, kept):
    """Whether an issue is one ended_releases reads: GitHub says it was
    reopened, it carries claimed or in-review, no claim branch stands for
    it and a branch released/<issue>-<short sha> does."""
    names = set(labels.label_names(issue))
    return (issue.get("state_reason") == board.REOPENED and bool({labels.CLAIMED, labels.IN_REVIEW} & names)
            and claims.branch_of(key) not in branches
            and any(name.startswith(claims.released_branch(key, "")) for name in kept))


def latest_close(events):
    """When an issue was last closed, as GitHub wrote the time of its
    latest event closed, or "" where its events record none."""
    times = [entry_time({"timestamp": e.get("created_at")}) for e in events
             if isinstance(e, dict) and e.get("event") == "closed"]
    return max((at for at in times if at), key=claims.parse_time, default="")


def claimed_since_close(events):
    """Whether the label claimed was put on after the issue's last close,
    by the order GitHub lists its events in, oldest first, so a close and
    a label in the same second still read in the order they came."""
    last_close = last_claimed = -1
    for at, event in enumerate(e for e in events if isinstance(e, dict)):
        if event.get("event") == "closed":
            last_close = at
        elif event.get("event") == "labeled" and (event.get("label") or {}).get("name") == labels.CLAIMED:
            last_claimed = at
    return last_claimed > last_close >= 0


def last_claimed(events):
    """When the label claimed was last put on the issue, as GitHub wrote
    the time of the latest such event, by the order GitHub lists them in;
    "" where none is listed or its time reads as none."""
    put_on = [e for e in events if isinstance(e, dict) and e.get("event") == "labeled"
              and (e.get("label") or {}).get("name") == labels.CLAIMED]
    return entry_time({"timestamp": put_on[-1].get("created_at")}) if put_on else ""


def claimed_by_a_claim(api, key, closed, claimed_at):
    """Whether the label claimed, last put on at claimed_at, after the
    time closed, is a claim's own rather than a person's undo: the newest
    entry the repository's activity lists for claude/<issue>, of any kind
    (newest_entry: its making, a push or its delete), came at or after
    the close and no earlier than the label, one read (ended_releases;
    F22's third review). False where the activity lists no entry for the
    branch, which an undo shows, an unlisted claim reading so too. None
    where GitHub keeps no activity of the repository (a 404) or no time is
    written for the label or the entry, which says neither, so the issue
    is left out this run and no label comes off on a guess; any other
    refusal raises, which leaves it out too."""
    if not claimed_at:
        return None
    entry = newest_entry(api, claims.branch_of(key))
    if entry is None:
        return None
    if not entry:
        return False
    at = entry_time(entry)
    if not at:
        return None
    when = claims.parse_time(at)
    return when >= claims.parse_time(closed) and when >= claims.parse_time(claimed_at)


def released_after(comments, closed, trusted):
    """Whether a release comment among an issue's comments came after the
    time closed, its author counted by trusted (claims.is_release_comment)."""
    after = [c for c in comments if isinstance(c, dict) and (at := entry_time({"timestamp": c.get("created_at")}))
             and claims.parse_time(at) > claims.parse_time(closed)]
    return bool(claims.latest_release(after, trusted))


def claim_heads(listed):
    """The head commit of each claim branch, by branch, from the branch
    list as the API gives it; a branch given by name alone has none."""
    return {name: sha for name, sha in branch_heads(listed).items() if claims.is_claim_branch(name)}


def branch_heads(listed):
    """The head commit of every branch, by name, from the branch list as
    the API gives it; a branch given by name alone has none."""
    return {branch["name"]: (branch.get("commit") or {}).get("sha") or ""
            for branch in listed if isinstance(branch, dict) and branch.get("name")}


def may_be_due(issue):
    """Whether a release can be due on a claim on this issue, or can have
    happened: the issue is open and not labelled in-review, which a
    finished claim waiting for review is, whether or not the label claimed
    is still on it; or open and reopened, since the coordinator takes the
    in-review an ended claim left off it before it releases, so its claim
    branch, which stands, is released in the same run once it has gone a
    day without a commit (coordinator.freeing; F22)."""
    names = (issue or {}).get("labels", [])
    return issue is not None and issue.get("state") == "open" and (
        labels.IN_REVIEW not in names or issue.get("state_reason") == board.REOPENED)


def ready_branches(pulls):
    """The head branches of the open pull requests marked ready, as the API
    details them; a pull request whose draft state is not given counts as
    ready, as the board counts it. One from a fork, or one a stranger
    opened (mark_strangers), is left out, whatever its branch is called,
    since it is no claim's (board.claim_of)."""
    return {(p.get("head") or {}).get("ref") for p in pulls if not p.get("draft") and is_claims(p)}


def is_claims(pull):
    """Whether a pull request, as the API details it, may be a claim's own:
    not from a fork, and not one mark_strangers found a stranger opened
    (board.claim_of)."""
    return not board.from_fork(pull) and not board.from_stranger(pull)


def permissions(api):
    """A callable saying whether a login may write to the repository, from
    GitHub's collaborator permission, asked once per login: True where it
    may, False where GitHub says it may not or knows no such login, and
    None where GitHub refused the call, which may_write and may_open each
    read their own way."""
    answers = {}

    def ask(login):
        if login not in answers:
            try:
                found = api.get(f"/repos/{api.repo}/collaborators/{urllib.parse.quote(login, safe='')}/permission",
                                missing_ok=True) or {}
            except GitHubError:
                answers[login] = None
            else:
                answers[login] = (found.get("permission") in WRITE_PERMISSIONS
                                  or found.get("role_name") in WRITE_ROLES)
        return answers[login]
    return ask


def may_write(api, asked=None):
    """A callable saying whether a login may write to the repository, from
    GitHub's collaborator permission (permissions), asked once per login:
    through asked where the caller passes it, else the client's own cache
    (GitHub.asked), so a board run and a clerk's look ask about an author
    once however many readers name them, and only a client with no cache
    pays for one of its own (F19's review). It stands in for
    author_association where that word does not vouch for an author:
    GitHub may show a private member of the organization to a workflow's
    token as a contributor. A call GitHub refuses counts as no, so a
    release by an author nobody can vouch for frees nothing."""
    asked = asked or getattr(api, "asked", None) or permissions(api)
    return lambda login: asked(login) is True


def may_open(asked):
    """A callable saying whether a login may have opened a claim's pull
    request (claims.opened_by_writer), from the answers permissions gives:
    unless GitHub says it may not write. A call GitHub refuses counts as
    yes, so a lost answer never takes a claim's own pull request for a
    stranger's, which would free a finished claim; it costs the runs in
    which GitHub refuses the call, each releasing no claim such a pull
    request holds."""
    return lambda login: asked(login) is not False


def mark_strangers(pulls, asked):
    """Write the board's key stranger (board.from_stranger) on each pull
    request, as the API details it, from a branch of this repository named
    as a claim's that someone who may not write opened
    (claims.opened_by_writer, may_open): it is no claim's, whatever its
    branch is called (board.claim_of; F19). An author author_association
    or an App's login vouches for costs no call; any other costs one, once
    a run, through asked. The pull requests are returned."""
    for pull in pulls:
        branch = (pull.get("head") or {}).get("ref") or ""
        if (claims.is_claim_branch(branch) and not board.from_fork(pull)
                and not claims.opened_by_writer(pull, may_open(asked))):
            pull["stranger"] = True
    return pulls


def committed_at(api, sha):
    """When the commit was committed, as GitHub writes it, or "" when GitHub
    does not know the commit."""
    commit = api.get(f"/repos/{api.repo}/git/commits/{sha}", missing_ok=True) or {}
    return (commit.get("committer") or {}).get("date") or ""


def head_since(api, branch, sha, refused=None, now=None, started_by=""):
    """When the branch's last commit came, as (pushed_at, committed_at),
    the one read and the other "" (board.Release.last_commit): when GitHub
    received the push that set the head, from the repository's activity
    (read_push, pushed_by), so no worker's clock decides how long a claim has
    been idle (F8: until its fix the head's committer date did, which the
    pusher sets, so a commit dated ahead held a claim past its day and one
    dated back had it released at once). now is when the run began, the
    wall clock unless given, taken before the run read the branch. A first
    read that lists no push that set the head, or that GitHub gives no
    answer to that holds (a 5xx, or none at all, but not a rate limit's
    refusal), is made once more after a pause, as settled asks again for
    mergeable: the push a run was started by is the one the activity is
    likeliest not to list yet, and the committer date in its place would
    let a commit dated back have its claim released in the run its own
    push started. Where the second read dates no push either, or GitHub
    refuses the read for good, the head's committer date stands in, the
    rule before the fix, so a token that may not read the activity still
    sees claims released and no claim is held for good; but never earlier
    than the newest entry either read listed for the branch, since the
    push that set the head came no sooner than that, so an activity that
    lags still dates a claim no earlier than the push it lists last
    (latest). In the run that a push to this branch started, started_by,
    which the workflow's event names (cli.pushed_branch), a head no read
    dates is dated at the run's start instead, for that run alone, since
    that push came just before the run, whatever the two reads list: a
    fresh branch's activity that lags both reads is an empty list, with no
    floor, so the committer date would have a claim commit dated back
    released in the run its own push started (F8's third review). So is a
    head a read dates by an entry that came more than PUSH_SINCE_SLACK
    before the run began: an older push of the same head that the activity
    lists while it does not list the run's own yet, as a head pushed back
    to one the branch had before gives, which dated by that entry would be
    released in the run its own push started where the entry is more than
    a day old (OLDER_PUSH_LISTED; F8's follow-up's second review); or the
    run's own push, where the run began that much after it, queued behind
    other runs, which the run's start then dates a little late, holding
    the claim no longer than the queue held the run (its review). The next run
    reads the activity again, so no claim is held for good. refused
    records why, by branch, for the run's log (pushed_refused_lines),
    except for a repository GitHub keeps no activity of (a 404), whose
    committer dates were always the rule outside the run its push
    started, which dates it at its start with no line either, and the
    second read's 404 after a first that listed nothing, which says the
    same. A rate limit's refusal is not read again after the pause: the
    client has waited it out already, as long as GitHub asks and a run
    waits, and made the call again (GitHub.request), so it is a refusal
    like any other, after which the committer date stands in
    (may_answer_later). It costs one read where the first read dates the
    push; two where it does not and the read is not made again, a refusal
    for good or a rate limit's, or a 404, the second the commit; and up to
    three with a pause where it is made again; in the run the branch's
    push started, the commit is not read. Each read is one call, or up to
    three where GitHub refuses it for a rate limit and the client asks
    again after each wait."""
    now = now or datetime.now(timezone.utc)
    read = first = read_push(api, branch, sha, now)
    if first.again:
        api.sleep(PUSH_RETRY_SECONDS)
        read = read_push(api, branch, sha, now)
    if read.pushed and (branch != started_by or came_since(read.pushed, now)):
        return read.pushed, ""
    floor, why = latest(read.floor, first.floor), read.why or (first.why if first.floor else "")
    if branch == started_by:
        began = github_time(now)
        why = OLDER_PUSH_LISTED if read.pushed else why
        if why and refused is not None:
            refused[branch] = f"{why}, so its idleness is read from the start of this run, {began}, which its push started"
        return began, ""
    committed = committed_at(api, sha)
    if floor and latest(committed, floor) != committed:
        if refused is not None:
            refused[branch] = (f"{why}, so its idleness is read from its newest entry's time, {floor}, "
                               + ("later than its committer date" if committed else "GitHub not showing its commit"))
        return floor, ""
    if why and refused is not None:
        refused[branch] = f"{why}, so its idleness is read from its committer date"
    return "", committed


def github_time(moment):
    """A moment as GitHub writes a time, ISO 8601 in UTC to the second, so
    a run's start can stand where the activity's times do."""
    return moment.astimezone(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


class PushRead(NamedTuple):
    """One read of when the push that set a branch's head came: GitHub's
    time of that push where the activity dates it; else why the committer
    date would stand in, in words for the run's log ("" for a repository
    GitHub keeps no activity of, a 404), whether a second read may tell
    more, and the time of the newest entry listed for the branch, which
    the push of the head came no sooner than ("" where none was listed)."""

    pushed: str = ""
    why: str = ""
    again: bool = False
    floor: str = ""


def read_push(api, branch, sha, now=None):
    """One read of when the push that set the branch's head came
    (PushRead), now being when the run began: dated where the newest entry
    dates the head (pushed_by); else owed a second read where GitHub's
    answer may change (may_answer_later) or the list names no push of the
    head yet."""
    try:
        entry = newest_entry(api, branch)
    except GitHubError as problem:
        return PushRead(why=f"GitHub refused to say when the push that set the branch's head came "
                            f"({problem}{github_words(problem)})", again=may_answer_later(problem))
    if entry is None:
        return PushRead()
    pushed = pushed_by(entry, sha, now)
    return PushRead(pushed=pushed) if pushed else PushRead(why=NO_PUSH_LISTED, again=True, floor=entry_time(entry))


def may_answer_later(problem):
    """Whether GitHub may answer the same call after a pause: it gave no
    answer, or a 5xx. Any other refusal is its answer for the run: a 401,
    a 403 for a token that may not read, a 422, and a rate limit's 429 or
    403 too, which the client has waited out already, as long as GitHub
    asks and a run waits, and made again (GitHub.request), so a read after
    a shorter pause would ask again while GitHub still asks for no call,
    which "may result in the banning of your integration" (GitHub's "Rate
    limits for the REST API")."""
    return problem.status is None or problem.status >= 500


def head_pushed_at(api, branch, sha, now=None):
    """When GitHub received the push that set the branch's head, as the
    repository's activity records it (pushed_by, on the newest entry of the
    branch's activity); "" when that entry does not date it, or there is
    no entry at all; None when GitHub keeps no activity of the repository
    (a 404). Any other refusal raises."""
    entry = newest_entry(api, branch)
    return None if entry is None else pushed_by(entry, sha, now)


def newest_entry(api, branch):
    """The newest entry of the branch's activity, of any kind; {} where the
    list is empty or its first item is no entry of this branch; None where
    GitHub keeps no activity of the repository (a 404). Any other refusal
    raises."""
    ref = f"refs/heads/{branch}"
    listed = api.get(f"/repos/{api.repo}/activity", {"ref": ref, "direction": "desc", "per_page": 1},
                     missing_ok=True)
    if listed is None:
        return None
    newest = listed[0] if isinstance(listed, list) and listed else None
    return newest if isinstance(newest, dict) and newest.get("ref") == ref else {}


def pushed_by(entry, sha, now):
    """The time the branch's newest entry dates its head by: its own, when
    it set the head the run read, a push, a forced push, the branch's
    creation or a merge into it; or when it came no earlier than
    PUSH_SINCE_SLACK before now, when the run began, whatever heads it
    names, a push after the run read the branch, which says the claim had
    a push then. An older entry that names another head is stale, since a
    later push the activity does not list yet set the head the run read,
    and dating by it would read the claim as idle while its holder pushes;
    with now None only the first rule dates. "" where neither holds, or
    the entry's time is none GitHub writes."""
    at = entry_time(entry)
    if at and (entry.get("after") == sha or came_since(at, now)):
        return at
    return ""


def entry_time(entry):
    """The entry's time as GitHub writes it, ISO 8601; "" for one missing
    or none such, which dates nothing rather than failing the run."""
    at = entry.get("timestamp") if isinstance(entry, dict) else None
    if not isinstance(at, str) or not at:
        return ""
    try:
        claims.parse_time(at)
    except ValueError:
        return ""
    return at


def came_since(timestamp, now):
    """Whether GitHub's time of an entry is no earlier than PUSH_SINCE_SLACK
    before now; a time that is missing or not one is not."""
    if not timestamp or now is None or not isinstance(timestamp, str):
        return False
    try:
        return claims.parse_time(timestamp) >= now - PUSH_SINCE_SLACK
    except ValueError:
        return False


def latest(first, second):
    """The later of two times as GitHub writes them, the first where they
    are the same instant, and either one where the other is "" or none
    GitHub writes."""
    try:
        later = claims.parse_time(first) >= claims.parse_time(second)
    except (ValueError, AttributeError):
        return first if entry_time({"timestamp": first}) else second
    return first if later else second


def pushed_refused_lines(snapshot):
    """A line for each claim branch whose head the run could not date by
    the push that set it, saying why and what it read in its place, the
    committer date, the newest entry listed or the start of the run its
    push started, for the workflow's log, since the page says nothing of
    it: a token that may not read the repository's activity would
    otherwise have every claim's idleness read from its committer date on
    every run unseen."""
    snap = board.as_snapshot(snapshot)
    return [f"{branch}: {why}" for branch, why in sorted(snap.pushed_refused.items())]


def claim_released(api, branch, now=None):
    """Whether a release comment counts for the claim on a claude/<issue>
    branch, by the part of the release rule its issue's comments and the
    branch's head decide, read as releases reads them for the board: the
    latest release comment by an author the repository trusts came after
    the branch had gone more than a day without a commit, as GitHub dates
    the push of its head (head_since), and nothing was committed since, or
    the head is one a release already kept (board.release_counts). The
    labels and the pull request decide the rest (board.is_released), which
    the clerk reads from its own look (team.Snapshot.released). A branch
    with no head GitHub knows is no released claim, as the board counts
    none. It costs a call for the branch's head, one for the released
    branch that would keep it, what head_since costs and, once the branch
    has gone more than a day without a commit or its head is kept, the
    issue's comments and a call for each author whose permission must be
    asked, once a client (GitHub.asked), so the clerk's look asks about an
    author once however many claims it reads."""
    now = now or datetime.now(timezone.utc)
    name = board.branch_name(branch)
    sha = branch_sha(api, name)
    if not sha:
        return False
    key = claims.key_of(name)
    kept = branch_sha(api, claims.released_branch(key, sha)) == sha
    pushed, committed = head_since(api, name, sha, now=now)
    release = board.Release(sha=sha, committed_at=committed, pushed_at=pushed, kept=kept)
    if not kept and (not release.last_commit or not claims.is_idle(release.last_commit, now)):
        return False
    comments = api.get_all(f"/repos/{api.repo}/issues/{quoted(key)}/comments")
    release.released_at = claims.latest_release(comments, may_write(api))
    return board.release_counts(release)


def finished_draft(api, pull):
    """Whether a draft's author said it is finished (claims.said_finished),
    from the pull request's own comments; None, a pull request, costs no
    call."""
    if not pull:
        return False
    author = (pull.get("user") or {}).get("login")
    return claims.said_finished(api.get_all(f"/repos/{api.repo}/issues/{pull['number']}/comments"), author)


def releases(api, heads, found_issues, pulls, now=None, asked=None, every_head=None, refused=None, started_by=""):
    """What the release rule reads about each claim that has gone more than
    a day without a commit, or whose head a release kept, by key, as a
    board.Release's JSON shape: the branch's head commit with the time its
    last commit came, and the time of the latest release comment on its
    issue, "" when there is none. The coordinator releases such a claim
    unless a comment already did (coordinator.releases_due), and the board
    counts it released once the label claimed is off too
    (board.released_keys). A finished claim, its issue labelled in-review
    or its pull request marked ready, is left out, since the board counts
    no release of it. So is a draft whose author said it is finished
    (finished_draft), which the coordinator never releases, unless a person
    has released it, as the protocol lets anyone who may write release any
    claim: it is read only once its issue no longer carries claimed, and
    kept only when its release comment counts (board.release_counts), so
    the board counts the release and the coordinator, which finds it done,
    adds none. A claim whose issue is labelled in-review or is not open,
    or whose pull request is marked ready, costs no call; any other costs
    what dating its head costs (head_since: one read where the activity
    dates the push at once, up to three where it does not, each up to
    three calls where GitHub refuses it for a rate limit), and, when its
    last commit is more than a day old or its head is kept, one for its
    draft's comments, then, unless the draft is finished and its issue
    still carries claimed, one for its issue's comments and one per
    author whose permission must be asked, so the run pays only where
    a release can be due or can have happened. A draft from a fork, or one
    a stranger opened, is no claim's (is_claims), so what its author says
    holds no claim. asked is the run's answers to who may write
    (permissions), shared with mark_strangers. When the last commit came is
    when GitHub received the push of the head (head_since, F8), and
    every_head, each branch's head by name, says which heads a release
    already kept, a branch released/<issue>-<short sha> standing at the
    same commit (board.Release.kept), which is read whatever its time, so
    the board frees the name of a branch pushed again at such a head;
    refused gathers, for the run's log, why each claim whose head no read
    dated by its push was dated by its committer date, by the newest
    entry listed or by the run's start: GitHub's refusal, or an activity
    that lists no push of the head, or, in the run the branch's push
    started, only one from before the run began. now is when the run
    began, taken before it read the branches, so a push after that read is
    told from a stale entry (pushed_by); started_by is the branch whose
    push started the run, "" for none, whose head is dated at the run's
    start where no read dates its push, or a read dates it by an entry from
    before the run began (head_since)."""
    now = now or datetime.now(timezone.utc)
    ready, trusted = ready_branches(pulls), may_write(api, asked)
    drafts = {(p.get("head") or {}).get("ref"): p for p in pulls if p.get("draft") and is_claims(p)}
    every_head = every_head or {}
    found = {}
    for branch, sha in sorted(heads.items()):
        key = claims.key_of(branch)
        if not sha or branch in ready or not may_be_due(found_issues.get(key)):
            continue
        kept = every_head.get(claims.released_branch(key, sha)) == sha
        pushed, committed = head_since(api, branch, sha, refused, now, started_by)
        release = board.Release(sha=sha, committed_at=committed, pushed_at=pushed, kept=kept)
        if not kept and (not release.last_commit or not claims.is_idle(release.last_commit, now)):
            continue
        finished = finished_draft(api, drafts.get(branch))
        if finished and labels.CLAIMED in found_issues[key]["labels"]:
            continue
        release.released_at = claims.latest_release(api.get_all(f"/repos/{api.repo}/issues/{key}/comments"), trusted)
        if finished and not board.release_counts(release):
            continue
        found[key] = release.to_dict()
    return found


def collect(api, sha, default_branch, now=None, started_by=""):
    """Everything the board renders, in the shape a --dry-run fixture has,
    so a live run's answers can be saved and rendered again offline; now is
    the moment a claim's idleness is measured from, the wall clock unless
    given, taken before the branches are read, so a push after that read
    dates its claim (pushed_by); started_by is the branch whose push
    started the run, "" for none (releases)."""
    now = now or datetime.now(timezone.utc)
    repo = api.repo
    asked = getattr(api, "asked", None) or permissions(api)
    issues = [i for i in api.get_all(f"/repos/{repo}/issues", {"state": "open"}) if "pull_request" not in i]
    pulls = mark_strangers([pull_with_detail(api, p["number"], default_branch)
                            for p in api.get_all(f"/repos/{repo}/pulls", {"state": "open"})], asked)
    listed = api.get_all(f"/repos/{repo}/branches")
    branches = board.claim_branches(listed)
    found_issues = claim_issues(api, branches, issues)
    undated = {}
    found = {
        "sha": sha,
        "default_branch": default_branch,
        "issues": issues,
        "pulls": pulls,
        "branches": branches,
        "claims": {claims.key_of(b): claim_text(api, b) for b in branches},
        "claim_issues": found_issues,
        "merged": recent_merged(api),
        "releases": releases(api, claim_heads(listed), found_issues, pulls, now, asked, branch_heads(listed),
                             undated, started_by),
    }
    if undated:
        found["pushed_refused"] = undated
    kept = kept_names(listed)
    if kept:
        found["kept_names"] = kept
        ended = ended_releases(api, issues, branches, kept, asked)
        if ended:
            found["ended_releases"] = ended
    found["began"], refused = began_times(api, found)
    if refused:
        found["began_refused"] = refused
    dependencies(api, found)
    return found


# --- the dependencies ---------------------------------------------------------

def within_budget(api, calls):
    """Whether the token has more than the reserve left after this many
    more calls, by the rate limit the latest answer said; a client that
    keeps no limit, a stand-in or one before its first answer, is within
    it."""
    limit = getattr(api, "limit", None)
    return limit is None or limit[0] - calls > DEPENDENCIES_RESERVE


def reserve_words(api):
    """What the log says of a read the reserve stopped."""
    left, allowed, _ = api.limit
    return f"{left} of {allowed} calls left, the reserve is {DEPENDENCIES_RESERVE}"


def pages_or_empty(api, path):
    """Every item of a list endpoint, or none where GitHub has no such
    issue (404) or someone deleted it (410), as issue_or_none reads one."""
    try:
        return api.get_all(path)
    except GitHubError as problem:
        if problem.status in (404, GONE):
            return []
        raise


def in_repository(issue, repo):
    """Whether an issue object is this repository's: its repository_url
    names the repository, as every issue GitHub gives carries, or it
    carries none, as a fixture's does. A sub-issue or an issue blocking
    one may live in another repository, where its number means nothing
    here and would be read and drawn as this repository's issue of the
    same number, so such an item is left out of the reads."""
    url = issue.get("repository_url")
    return not isinstance(url, str) or url.lower().endswith("/repos/" + repo.lower())


def owes_no_read(parent):
    """Whether a planned issue's own summary says it has no sub-issue, so
    its read would answer nothing: the open list carries the summary on
    every issue, and an issue with none, an older answer, is read."""
    summary = parent.get("sub_issues_summary")
    # A bool is not a count, though Python counts it: a summary that says
    # no whole number is read, as the issue is read when it carries none.
    return isinstance(summary, dict) and type(summary.get("total")) is int and summary["total"] == 0


class DependencyReads:
    """What the Dependencies section's reads found, as collect writes it
    into the snapshot: each open epic's and goal's sub-issues, the issues
    blocking each of an epic's sub-issues, and the issues those name that
    the open list does not hold, as board.Issue dictionaries."""

    def __init__(self, api, found):
        self.api = api
        self.repo = api.repo
        self.listed = {issue["number"] for issue in found["issues"]}
        self.sub_issues, self.blocked_by, self.others = {}, {}, {}
        self.reading = ""

    def read(self, planned):
        """Read every planned issue's sub-issues and, for an epic's, the
        issues blocking each sub-issue, which the epic's chart draws; a
        goal's sub-issues are epics, whose chart draws no relation of the
        epic's own, so their blocked-by relations are not read. False where
        the reserve stopped a read, with the rest left unread, since a
        section half drawn would mislead."""
        for parent in sorted(planned, key=lambda issue: issue["number"]):
            self.reading = f"#{parent['number']}'s sub-issues"
            children = [child for child in pages_or_empty(self.api, self.issue_path(parent["number"], "sub_issues"))
                        if "pull_request" not in child and in_repository(child, self.repo)]
            if not within_budget(self.api, 0):
                return False
            if children:
                self.sub_issues[parent["number"]] = [child["number"] for child in children]
            for child in children:
                self.remember(child)
                if labels.EPIC in labels.label_names(parent) and not self.read_blockers(child):
                    return False
        return True

    def read_blockers(self, child):
        """Read the issues blocking a sub-issue, where its summary counts
        one and no earlier parent's read brought them; False where the
        reserve stopped the read."""
        if labels.blocked_by_counts(child)[1] == 0 or child["number"] in self.blocked_by:
            return True
        if not within_budget(self.api, 1):
            return False
        self.reading = f"the issues blocking #{child['number']}"
        listed = pages_or_empty(self.api, self.issue_path(child["number"], "dependencies/blocked_by"))
        blockers = [blocker for blocker in listed if in_repository(blocker, self.repo)]
        if not within_budget(self.api, 0):
            return False
        self.blocked_by[child["number"]] = sorted(blocker["number"] for blocker in blockers)
        for blocker in blockers:
            self.remember(blocker)
        return True

    def remember(self, issue):
        """Keep an issue the open list does not hold, a closed sub-issue
        or an issue blocking one from outside the epic, so the chart draws
        it."""
        if issue["number"] not in self.listed:
            self.others[issue["number"]] = board.Issue.from_dict(issue).to_dict()

    def issue_path(self, number, rest):
        return f"/repos/{self.repo}/issues/{number}/{rest}"

    def write_into(self, found):
        for key in ("sub_issues", "blocked_by", "others"):
            if getattr(self, key):
                found[key] = getattr(self, key)


def dependencies(api, found):
    """The Dependencies section's reads, last in collect, so the page's
    tables are read whatever happens here: for every open issue labelled
    epic or goal whose summary does not say it has no sub-issue, its
    sub-issues, and for every sub-issue of an epic whose summary counts
    an issue blocking it, those issues. The reads stop short of the
    reserve (within_budget), and a read GitHub fails drops the section
    with its own notice; either way the run goes on and the page says why
    the section is not drawn (board.NOT_DRAWN). The chart each listed
    body holds is kept whatever is read (board.listed_graphs), so the
    graph writer knows what stands without reading again. Nothing here
    raises."""
    graphs = board.listed_graphs(found["issues"])
    if graphs:
        found["graphs"] = graphs
    planned = [issue for issue in found["issues"]
               if labels.planned(labels.label_names(issue)) and not owes_no_read(issue)]
    if not planned:
        return
    reads = DependencyReads(api, found)
    if not within_budget(api, len(planned)):
        found["dependencies_skipped"], found["dependencies_notice"] = "reserve", reserve_words(api)
        return
    try:
        complete = reads.read(planned)
    except GitHubError as problem:
        found["dependencies_skipped"] = "failed"
        found["dependencies_notice"] = f"reading {reads.reading} failed: {problem}{github_words(problem)}"
        return
    if not complete:
        found["dependencies_skipped"], found["dependencies_notice"] = "reserve", reserve_words(api)
        return
    reads.write_into(found)


def dependencies_skipped_lines(snapshot):
    """The line for the workflow's log where the Dependencies section was
    left out, and why, or none."""
    snap = board.as_snapshot(snapshot)
    if not snap.dependencies_skipped:
        return []
    return [f"board: Dependencies left out: {snap.dependencies_notice or snap.dependencies_skipped}"]


def write_graphs(api, snapshot, blocks):
    """Write each open epic's flowchart into its body, between the graph
    markers, after the page is written, and take the chart out of an
    issue that owes none any more, its sub-issues taken off; nothing in a
    run whose reads for the section were left out, since it knows nothing
    of what the bodies owe. An epic whose listed body holds the block
    already costs no call; any other is read afresh, since the clerk or a
    person may have written it since the run began, and written only when
    its body would change, so a run that finds nothing moved writes
    nothing. A closed or deleted epic, a body whose markers are not one
    clean pair, which a person wrote and the board leaves alone, a write
    that would reach the reserve and one GitHub refuses are each a line
    and never stop the next; nothing here raises."""
    snap = board.as_snapshot(snapshot)
    if snap.dependencies_skipped:
        return []
    lines = []
    for number in sorted(set(blocks) | set(snap.graphs)):
        block = blocks.get(number)
        if block is not None and snap.graphs.get(number) == block:
            continue
        if not within_budget(api, 2):
            lines.append(f"epic #{number}: graph not written: {reserve_words(api)}")
            continue
        try:
            lines += write_graph(api, number, block)
        except GitHubError as problem:
            lines.append(f"epic #{number}: graph write failed: {problem}{github_words(problem)}")
    return lines


def write_graph(api, number, block):
    """One epic's write: its body with the block, or without the chart
    where block is None, compared with the body as read with Windows
    line endings folded, so a body GitHub stored with them is not
    rewritten every run."""
    issue = issue_or_none(api, number)
    if issue is None or issue.get("state") == "closed":
        return []
    body = issue.get("body") or ""
    wanted = board.body_with_graph(body, block) if block is not None else board.body_without_graph(body)
    if wanted is None:
        return [] if block is None else [f"epic #{number}: graph block not written: markers unpaired"]
    if wanted.replace("\r\n", "\n") == body.replace("\r\n", "\n"):
        return []
    api.set_body(number, wanted)
    return [f"epic #{number}: graph {'written' if block is not None else 'taken out'}"]


def began_times(api, found):
    """When each piece of open work that meets other work in a file began,
    by its branch: the time GitHub writes for the newest push that made a
    branch of that name, read from the repository's activity, which GitHub
    records as the pushes reach it, so no clock of a worker's own comes in
    (F13); and GitHub's answer for each branch it refused to date. Only the
    work in the Overlaps is asked about, since the landing order is all the
    time is for, and a pull request from a fork never, since its branch is
    none of this repository's (board.open_work). One call reads the newest
    branch creations of the whole repository (newest_creations), and each
    branch it does not name is asked about on its own (branch_made_at), so
    a run pays one call where the work is new, as it is while it overlaps,
    and a refused first call costs only the reads it would have saved.
    When GitHub has no such push to give, or no activity at all (a 404),
    or refuses a branch's own read, whatever it answers, that branch is
    left out and the board orders its work as work whose beginning nothing
    says (board.landing_key), a pull request by when it was opened: the
    order is the one thing the time is for, and one branch GitHub will not
    answer for must not fail the board's whole run, as a comparison GitHub
    refuses does not (standing_of)."""
    wanted = list(dict.fromkeys(work.branch for works in board.overlaps(found).values()
                                for work in works if work.branch))
    if not wanted:
        return {}, {}
    try:
        times = newest_creations(api, wanted)
    except GitHubError:
        times = {}
    refused = {}
    for branch in wanted:
        if branch in times:
            continue
        try:
            began = branch_made_at(api, branch)
        except GitHubError as problem:
            refused[branch] = f"{problem}{github_words(problem)}"
            continue
        if began:
            times[branch] = began
    return times, refused


def newest_creations(api, branches):
    """When the newest push that made each of these branches came, for those
    the repository's newest ACTIVITY_PAGE branch creations name: one call,
    the activity filtered to branch creations and listed newest first, so
    the first entry for a branch is its newest creation. A refusal other
    than a 404 raises."""
    wanted = {f"refs/heads/{branch}": branch for branch in branches}
    listed = api.get(f"/repos/{api.repo}/activity", {"activity_type": "branch_creation", "direction": "desc",
                                                     "per_page": ACTIVITY_PAGE}, missing_ok=True) or []
    times = {}
    for activity in listed if isinstance(listed, list) else []:
        made = made_at(activity)
        branch = wanted.get(activity.get("ref")) if made and isinstance(activity.get("ref"), str) else None
        if branch and branch not in times:
            times[branch] = made
    return times


def branch_made_at(api, branch):
    """When the newest push that made the branch came, as GitHub writes the
    time, or "" when GitHub records none: the activity list, filtered to the
    branch and to branch creations and listed newest first. A refusal other
    than a 404 raises."""
    listed = api.get(f"/repos/{api.repo}/activity", {"ref": f"refs/heads/{branch}", "activity_type": "branch_creation",
                                                     "direction": "desc", "per_page": 1}, missing_ok=True) or []
    for activity in listed if isinstance(listed, list) else []:
        if made_at(activity) and activity.get("ref") == f"refs/heads/{branch}":
            return made_at(activity)
    return ""


def made_at(activity):
    """The time of an entry of the activity that records a push making a
    branch, "" for any other entry."""
    if not isinstance(activity, dict) or activity.get("activity_type") != "branch_creation":
        return ""
    return activity.get("timestamp") or ""


def began_refused_lines(snapshot):
    """A line for each branch of overlapping work GitHub refused to date in
    this run, for the workflow's log, since the page says nothing of why
    that work stands where it does: a token that may not read the
    repository's activity would otherwise leave every claim undated on
    every run unseen. A pull request from such a branch is then dated by
    its opening, and a claim with none yet is work nothing dates, which
    lands last (board.pull_began, board.UNDATED)."""
    snap = board.as_snapshot(snapshot)
    opened = {pull.branch for pull in snap.pulls if not pull.fork}
    own = {pull.branch for pull in snap.pulls if board.claim_of(pull) is not None}
    return [f"{branch}: GitHub refused to say when the branch was made ({answer}), so the Overlaps "
            f"{undated_words(branch, opened, own)}"
            for branch, answer in sorted(snap.began_refused.items())]


def undated_words(branch, opened, own):
    """What the Overlaps do with the work of a branch GitHub refused to
    date: date its pull request by its opening where the branch has one of
    its own, or is no claim's; and for a claim branch with none of its own,
    put the claim last, and date a stranger's pull request from it, where
    there is one, by its opening, since that is a piece of its own (F19;
    the fix's review)."""
    if branch in own or (branch in opened and not claims.is_claim_branch(branch)):
        return UNDATED_PULL
    return UNDATED_STRANGER if branch in opened else UNDATED_CLAIM


def snapshot(api, sha, default_branch, now=None, started_by=""):
    """The board's Snapshot of the repository as GitHub knows it now;
    started_by is the branch whose push started the run, "" for none
    (collect)."""
    return board.Snapshot.from_dict(collect(api, sha, default_branch, now, started_by))


def has_marker(api, number, marker):
    """Whether a comment on the issue or pull request carries the marker."""
    comments = api.get_all(f"/repos/{api.repo}/issues/{number}/comments")
    return any(marker in (c.get("body") or "") for c in comments)


def comment_once(api, number, marker, text):
    """Post the comment unless one on the pull request carries the marker."""
    if has_marker(api, number, marker):
        return False
    api.post(f"/repos/{api.repo}/issues/{number}/comments", {"body": marker + "\n" + text})
    return True


def post_comments(api, due):
    """Post each comment board.comments_due names, once each; the lines say
    which were posted and which failed, one each, for the workflow's log.
    Nothing here raises: a comment GitHub refuses, or whose pull request's
    comments it would not give, is a line, and the next comment is still
    tried. The board is written before the comments, and the workflow
    writes it to its branch only after a run that succeeded, so a refusal
    that escaped here would keep the board off its branch, and one GitHub
    gives for good would keep it off on every run. A refused comment is
    due again, idle or "main moved" alike, since both depend on the state
    and not on the run's event (board.comments_due), and a later run,
    which finds no marker yet, posts it."""
    lines = []
    for number, marker, text in due:
        kind = "idle" if marker == board.IDLE_MARKER else "main-moved"
        try:
            posted = comment_once(api, number, marker, text)
        except GitHubError as problem:
            lines.append(f"PR #{number}: {kind} comment failed: {problem}{github_words(problem)}")
            continue
        if posted:
            lines.append(f"PR #{number}: {kind} comment posted")
    return lines


# --- the coordinator ----------------------------------------------------------

def coordinate(api, snapshot, now):
    """What the coordinator does in the board's run, before the board is
    rendered: free every reopened issue coordinator.reopened_due names,
    release every claim coordinator.releases_due names, then tell each
    piece of open work about other work it meets that it has not been told
    of. The snapshot is brought up to date with each change, so the board
    this run renders shows it and a released branch is renamed in the same
    run, as a claim whose in-review a freeing took off is released once
    due, its release read with the others' (may_be_due). Nothing here raises: a
    call GitHub refuses is a line, and the next issue is still freed and
    the next claim still released. The lines say what was done, one each,
    for the workflow's log."""
    snap = board.as_snapshot(snapshot)
    lines = [free_reopened(api, snap, number) for number in coordinator.reopened_due(snap)]
    lines += [release_one(api, snap, key, now) for key in coordinator.releases_due(snap, now)]
    return lines + tell_overlaps(api, snap)


def free_reopened(api, snap, number):
    """Free one reopened issue whose labels still speak of a claim that
    ended (coordinator.freeing; F22): read it again first, since a merge
    can close an issue after the board read it, and leave one that is
    closed by then; where the board read no branch for it, read its claim
    branch again too, and leave the issue to the next run once one
    stands, since a worker that claimed it after the board read the
    branches has just put claimed and its label on; and read its claim's
    open pull request again, which the freeing reads in place of the
    board's, so the issue is left as it is once that is marked ready,
    since in-review is then the finished claim's own; then
    put todo on where it goes and take the labels off in freeing's order,
    claimed last, reading the claim branch once more just before claimed
    comes off and leaving claimed on once one stands, so a worker that
    claimed the issue while it was freed keeps the label its claim put
    on, and a run that fails partway leaves the issue carrying claimed or
    in-review for the next run to free. Nothing is put back: every
    step brings the issue nearer to free. The snapshot takes what was
    done, so the board this run renders shows it. Nothing here raises;
    what was done, in a line."""
    branch = claims.branch_of(str(number))
    try:
        issue = api.issue(number) or {}
        if issue.get("state") != "open":
            return f"#{number}: not freed, closed since the board read it"
        pull = api.pull_detail(branch)
        plan = coordinator.freeing(with_claims_pull(snap, number, pull), number, labels.label_names(issue))
        if not plan.taken_off and marked_ready(pull):
            return f"#{number}: not freed, its claim's pull request was marked ready since the board read it"
        if plan.taken_off and not plan.branch_stands and branch_sha(api, branch):
            return f"#{number}: not freed, claimed since the board read the branches"
    except GitHubError as problem:
        return f"#{number}: not freed, reading it again failed: {problem}{github_words(problem)}"
    put_on, taken_off = [], []
    try:
        if plan.put_on:
            api.post(f"/repos/{api.repo}/issues/{number}/labels", {"labels": list(plan.put_on)})
            put_on += plan.put_on
        for name in plan.taken_off:
            if name == labels.CLAIMED and not plan.branch_stands and branch_sha(api, branch):
                coordinator.after_labels(snap, str(number), taken_off, put_on)
                return (f"#{number}: reopened, not wholly freed, claimed while it was freed, so claimed stays"
                        + done_words(put_on, taken_off))
            take_label_off(api, number, name)
            taken_off.append(name)
    except GitHubError as problem:
        coordinator.after_labels(snap, str(number), taken_off, put_on)
        return (f"#{number}: reopened, not wholly freed, the next run goes on: {problem}{github_words(problem)}"
                + done_words(put_on, taken_off))
    coordinator.after_labels(snap, str(number), taken_off, put_on)
    if not taken_off:
        return f"#{number}: reopened, nothing left to free"
    if plan.branch_stands:
        return f"#{number}: reopened with a branch standing for it" + done_words(put_on, taken_off)
    return f"#{number}: reopened with its claim ended, freed" + done_words(put_on, taken_off)


def with_claims_pull(snap, number, pull):
    """The snapshot with the claim's own open pull request on this issue as
    read again (pull_detail), or none where it reads none, in place of the
    one the board read, so a pull request marked ready since then finishes
    the claim for the freeing (coordinator.freeing)."""
    others = [p for p in snap.pulls if not board.of_claim(p, str(number))]
    return dataclasses.replace(snap, pulls=others + ([board.Pull.from_dict(pull)] if pull else []))


def marked_ready(pull):
    """Whether the claim's own open pull request, as pull_detail gives it,
    is marked ready; None, no such pull request, is not."""
    return pull is not None and not pull.get("draft")


def done_words(put_on, taken_off):
    """What a freeing did, as the end of its line."""
    said = f"; {', '.join(put_on)} put on" if put_on else ""
    return said + (f"; {', '.join(taken_off)} taken off" if taken_off else "")


def holding(api, default_branch, sha):
    """What a branch holds beside the default branch, as coordinator.Holding,
    from GitHub's comparison of the two; ahead is None when GitHub would not
    compare them."""
    try:
        found = api.get(f"/repos/{api.repo}/compare/{quoted(default_branch)}...{sha}", missing_ok=True)
    except GitHubError:
        found = None
    if not found:
        return coordinator.Holding()
    newest = (found.get("commits") or [{}])[-1]
    message = (newest.get("commit") or {}).get("message") or ""
    files = tuple(sorted(f["filename"] for f in found.get("files") or [] if f.get("filename")))
    return coordinator.Holding(found.get("ahead_by"), (message.splitlines() or [""])[0], files)


def release_one(api, snap, key, now):
    """Release one claim: the labels first, claimed and every worker's off
    and todo on where the issue carries none of todo, epic and blocker
    (coordinator.labels_put_on), then the comment, in the protocol's
    order, and the board counts a release only once both are done. A run
    that fails at the labels, taking one off or putting todo on, puts
    back the labels it took off, even one whose removal GitHub answered
    with an error other than a refusal that says it wrote nothing
    (take_labels_off; F27), and takes off the todo it put on, even one
    whose putting on GitHub answered with an error, since GitHub may have
    written either (put_labels_on, put_back), so the issue is as the run
    found it and the next run starts the release over. A run whose
    comment fails reads the issue's comments, and the comment there after
    all, one GitHub wrote and lost its answer to, makes the release
    whole. Where they list none and GitHub refused the comment with a
    4xx other than a rate limit's, which says it wrote nothing
    (wrote_nothing), the run puts the labels back as a run that fails at
    the labels does. Where they list none after any other answer, a rate
    limit's refusal, a server error, a redirect, an answer that was not
    JSON or no answer at all, none of which says GitHub wrote nothing,
    and GitHub may list a comment it wrote only later (not_listed; F23
    and its review), or where the run cannot read them, whatever GitHub
    answered the comment with (unsettled; F20), the run cannot tell
    whether a release comment stands, its own, an earlier run's or one
    someone else wrote, so it leaves the labels as the release left them
    and counts no release of its own, since labels put back beside a
    comment that stands would hold the claim for good: the comment tells
    every later coordinator the claim is released, and the label claimed
    tells the board it is held; the next run, reading the comments,
    settles it as it settles any claim whose labels are off, finding the
    release whole or releasing the claim again unless its holder pushed
    or a person put the label claimed back in the meantime, while on a
    head a release already kept, where a release comment that counts
    stands, whenever it came (board.release_counts), that earlier release
    frees the claim in the same run once the labels are off. Where GitHub
    refuses the putting back too, the labels are left off with no
    comment, as an unsettled release leaves them where GitHub wrote none,
    until a run releases the claim or a person puts them back. What was
    done, in a line."""
    release = snap.releases[key]
    names, put_on = coordinator.labels_released(snap, key), coordinator.labels_put_on(snap, key)
    taken_off, problem = take_labels_off(api, key, names)
    if problem:
        return (f"#{key}: not released, taking the labels off failed: {problem}{github_words(problem)}"
                + put_back(api, key, taken_off))
    refused = put_labels_on(api, key, put_on, taken_off)
    if refused:
        return refused
    text = coordinator.release_text(snap, key, holding(api, snap.default_branch, release.sha))
    lost = ""
    try:
        released_at = (api.post(f"/repos/{api.repo}/issues/{key}/comments", {"body": text}) or {}).get("created_at")
    except GitHubError as problem:
        try:
            released_at = landed_release(api, key, release)
        except GitHubError as unread:
            return unsettled(snap, key, names, put_on, problem, unread)
        if not released_at and not wrote_nothing(problem):
            return not_listed(snap, key, names, put_on, problem)
        if not released_at:
            return (f"#{key}: not released, the comment failed: {problem}{github_words(problem)}"
                    + put_back(api, key, taken_off, put_on))
        lost = f" (GitHub's answer to the comment was lost, {problem}, but the comment is there)"
    coordinator.after_release(snap, key, released_at or now.isoformat(), names, put_on)
    return f"#{key}: claim released, idle since {release.last_commit}{lost}"


def take_labels_off(api, key, names):
    """Take a release's labels off its issue, in order, stopping at the
    first GitHub fails: the labels to put back should a step fail, and
    the failure, or None. The label whose removal failed is among them
    unless GitHub refused it in a way that says it wrote nothing
    (wrote_nothing), since a server error or no answer at all may come
    after GitHub took it off, and putting back a label still on changes
    nothing, where leaving off one GitHub took leaves the issue as
    neither the release nor the run found it (F27, until whose fix only
    the labels whose removal GitHub answered went back). A 404 that
    answers a try after a rate limit's refusal counts as a label taken
    off, for the same reason (take_label_off; F27's review)."""
    taken_off = []
    for name in names:
        try:
            if take_label_off(api, key, name, retried_404_taken=True):
                taken_off.append(name)
        except GitHubError as problem:
            return taken_off + ([] if wrote_nothing(problem) else [name]), problem
    return taken_off, None


def take_label_off(api, key, name, retried_404_taken=False):
    """Take one label off a claim's issue: True when it came off, False when
    GitHub answered 404, for an issue that no longer carried it, which is
    the label being off too, so a release goes on and puts back, if a later
    step fails, only what it took off. With retried_404_taken, as a
    release takes its labels off (take_labels_off), a 404 that answers a
    try the client made after a rate limit's refusal of an earlier one
    counts as the label taken off: GitHub may have taken it off on the
    refused try, so the 404 says it is off now, not that the issue never
    carried it, and the release goes on and puts it back should a later
    step fail, where a plain 404 left it off a held claim's issue (F27's
    review). Taking a todo off again, or freeing a reopened issue, wants
    only the label off, which such a 404 says either way."""
    try:
        api.delete(f"/repos/{api.repo}/issues/{key}/labels/{urllib.parse.quote(name, safe='')}")
    except GitHubError as problem:
        if problem.status == 404:
            return retried_404_taken and problem.after_rate_limit
        raise
    return True


def put_labels_on(api, key, put_on, taken_off):
    """Put on a claim's issue the labels its release puts on
    (coordinator.labels_put_on): nothing when they went on, else the
    release's line saying it failed there, with the labels taken off put
    back and those put on taken off again, since GitHub may have written
    them and lost its answer."""
    if not put_on:
        return ""
    try:
        api.post(f"/repos/{api.repo}/issues/{key}/labels", {"labels": put_on})
    except GitHubError as problem:
        return (f"#{key}: not released, putting {', '.join(put_on)} on failed: {problem}{github_words(problem)}"
                + put_back(api, key, taken_off, put_on))
    return ""


def put_back(api, key, names, put_on=()):
    """Put the labels a release took off back on its issue once a later step
    failed, so the claim stays its holder's until a run releases it whole,
    then take off the labels it put on, so the issue is as the run found
    it; what was done, as the end of the release's line. A label put on
    whose answer was an error is taken off too, since GitHub may have
    written it and lost the answer; one that is not there counts as taken
    off, as in the release itself."""
    said = ""
    if names:
        listed = ", ".join(names)
        try:
            api.post(f"/repos/{api.repo}/issues/{key}/labels", {"labels": names})
            said += f"; {listed} put back"
        except GitHubError as problem:
            said += f"; putting {listed} back failed too, so they may still be off: {problem}{github_words(problem)}"
    for name in put_on:
        try:
            if take_label_off(api, key, name):
                said += f"; {name} taken off again"
        except GitHubError as problem:
            said += f"; taking {name} off again failed too, so it may still be on: {problem}{github_words(problem)}"
    return said


def unsettled(snap, key, names, put_on, problem, unread):
    """The line of a release whose comment failed and whose issue's
    comments could not be read to see whether a release comment stands:
    the run's own, which GitHub may have written unless its refusal says
    it wrote nothing (wrote_nothing), or, whatever GitHub answered, a
    refusal among them, one an earlier run's post or someone else wrote;
    so the labels stay as the release left them (labels_left; F20, until
    whose fix the labels went back beside a comment that might stand,
    which held the claim for good). Until F23's review the line asked
    whether GitHub wrote the run's own comment, which after a refusal
    that says it wrote nothing it surely did not."""
    return labels_left(snap, key, names, put_on, problem,
                       f"its issue's comments could not be read to see whether a release comment stands: "
                       f"{unread}{github_words(unread)}")


def not_listed(snap, key, names, put_on, problem):
    """The line of a release whose comment GitHub answered with anything
    but a refusal that says it wrote nothing (wrote_nothing), a rate
    limit's refusal, a server error, a redirect, an answer that was not
    JSON or no answer at all, and whose issue's comments, read after it,
    list none: none of those says GitHub wrote nothing, and it may list a
    comment it wrote only later, so the read does not tell whether the
    comment stands, and the labels stay as the release left them
    (labels_left; F23, until whose fix the labels went back, which held
    the claim for good where GitHub had written the comment, and until
    whose review they went back after a rate limit's refusal too). The
    line says the same of a rate limit's refusal as of a server error,
    since GitHub may have acted on the call before it refused it."""
    return labels_left(snap, key, names, put_on, problem,
                       "its issue's comments list none yet, though GitHub may have written it before that failure "
                       "and list it only later")


def labels_left(snap, key, names, put_on, problem, why):
    """The line of a release the run cannot tell is whole, why being the
    reason it cannot: the labels stay as the release left them, and the
    snapshot shows them so, with no release of this run's, since the run
    cannot tell whether one stands, so the page lists the claim as held
    and nothing is renamed this run; the next run, reading the comments,
    finds the comment and the release is whole, or finds none and
    releases the claim again, unless its holder pushed since or a person
    put the label claimed back beside the comment. A head a release
    already kept, where a release comment that counts stands, is the
    exception: that comment counts whenever it came
    (board.release_counts), so with the labels off the snapshot counts
    that earlier release, and the same run lists the claim as released
    and deletes the branch, whose work is kept, with no word of this
    release on the issue where GitHub wrote none (F20's review). So the
    line says what was done and nothing of what comes next."""
    coordinator.after_labels(snap, key, names, put_on)
    return (f"#{key}: release unsettled, the comment failed: {problem}{github_words(problem)}; {why}; "
            f"the labels stay as the release left them")


def wrote_nothing(problem):
    """Whether a write that failed surely changed nothing: GitHub refused
    it with a 4xx, which it answers in place of making the write, where a
    server error, a 502 among them, or no answer at all may come after
    GitHub made it, and GitHub may list what it made only later, so a read
    that finds nothing then tells nothing (F23); a redirect of the post, a
    307 or a 308, which urllib does not follow for a POST, and an answer
    that was not JSON, which has no status either, say nothing of whether
    GitHub made it, so neither is a refusal (F23's review). Nor is a rate
    limit's refusal, a 429 or a 403 that says it is one (rate_limited),
    though it is a 4xx: GitHub's documentation says how long to wait and
    when to stop asking, never that GitHub refuses such a call before it
    acts on it, which the client's making a write again assumes
    (GitHub.request) and canary C4 would settle and is not run, so the
    labels a release puts back rest on no such assumption (F23's review;
    validation/DESIGN.md, 3.2 and 11). Nor is any refusal of a call the
    client made again after a rate limit's refusal of an earlier try
    (GitHubError.after_rate_limit): the 4xx that ends the call answers its
    last try alone, and GitHub may have acted on the try it refused for
    the rate limit (F23's second review)."""
    if problem.after_rate_limit or rate_limited(problem.status, problem.headers, problem.detail):
        return False
    return problem.status is not None and 400 <= problem.status < 500


def landed_release(api, key, release):
    """When the comment that releases this claim was made, by the rule the
    next run reads it with (claims.was_released), or "" when none is
    listed: a call that GitHub answered with an error, a 502 or a timeout,
    may still have written the comment, and once it stands the release is
    done; only a comment newer than the latest the run held is this run's
    (newer). "" says no more than that none is listed yet, since GitHub
    may list a comment it wrote only later (release_one, wrote_nothing;
    F23). A read GitHub refuses raises GitHubError, since the run then
    cannot tell whether the comment stands (release_one)."""
    comments = api.get_all(f"/repos/{api.repo}/issues/{key}/comments")
    released_at = claims.latest_release(comments, may_write(api))
    if not released_at or not newer(released_at, release.released_at):
        return ""
    return released_at if board.release_counts(dataclasses.replace(release, released_at=released_at)) else ""


def newer(found, held):
    """Whether a release comment found after a failed post is newer than
    the latest the run held before it: one no newer is an earlier
    release's, which for a kept head counts whenever it came
    (board.release_counts), so it would pass for this run's comment while
    none was written."""
    return not held or claims.parse_time(found) > claims.parse_time(held)


def tell_overlaps(api, snap):
    """Post each overlap comment that is due, on the piece's identity, the
    issue for a claim's work and the pull request itself for one that is
    no claim's (board.identity), for the work that meets other work,
    counting only the markers of authors the repository trusts
    (coordinator.told_numbers) and only those that came no earlier than
    the claims they speak of (coordinator.heard_of). The markers on the
    identity are read first; a claim's work that they leave something to
    tell is read on its own open pull requests too
    (coordinator.told_places), where the coordinator told it before F2's
    fix, so the upgrade tells nothing again; the extra read is paid only
    where the issue's markers alone would have a comment go out, which for
    a claim told on its pull request before the fix, and never on its
    issue since, is every run until its overlap ends. Work whose comments
    GitHub would not give is told nothing this run, so it is never told
    twice. The lines say what was posted."""
    lines, told, trusted = [], {}, may_write(api)
    for number in coordinator.meeting(snap):
        markers = markers_on(api, number, trusted, lines, f"#{number}: overlaps not told, its comments were not read")
        if markers is not None:
            told[number] = markers
    for number, _, _ in coordinator.overlap_comments_due(snap, told):
        for place in coordinator.told_places(snap, number) if number in told else ():
            failed = f"#{number}: overlaps not told, the comments of its pull request #{place} were not read"
            markers = markers_on(api, place, trusted, lines, failed)
            if markers is None:
                told.pop(number)
                break
            told[number] = coordinator.merged_times(told[number], markers)
    for number, marker, text in coordinator.overlap_comments_due(snap, told):
        if number not in told:
            continue
        try:
            api.post(f"/repos/{api.repo}/issues/{number}/comments", {"body": marker + "\n" + text})
        except GitHubError as problem:
            lines.append(f"#{number}: overlap comment failed: {problem}{github_words(problem)}")
            continue
        lines.append(f"#{number}: overlap comment posted")
    return lines


def markers_on(api, number, trusted, lines, failed):
    """When the latest trusted overlap marker on an issue or a pull request
    named each number (coordinator.told_times), or None, with a line
    saying why, when GitHub would not give its comments."""
    try:
        return coordinator.told_times(api.get_all(f"/repos/{api.repo}/issues/{number}/comments"), trusted)
    except GitHubError as problem:
        lines.append(f"{failed}: {problem}{github_words(problem)}")
        return None


def branch_sha(api, branch):
    """The branch's commit, or "" when the repository has no such branch."""
    found = api.get(f"/repos/{api.repo}/git/ref/heads/{quoted(branch)}", missing_ok=True)
    return ((found or {}).get("object") or {}).get("sha") or ""


def release_claims(api, snapshot):
    """Rename the branch of every claim board.released_keys counts as
    released, so its name is free for the next claim, and say so once on
    its issue; the lines say what was done, one per claim, for the
    workflow's log. GitHub closes a pull request whose branch is renamed,
    and keeps the commits under the new name."""
    snap = board.as_snapshot(snapshot)
    return [release_claim(api, key, snap.releases[key])
            for key in sorted(board.released_keys(snap), key=claims.by_number)]


def release_claim(api, key, release):
    """Free one released claim's name and say so once on its issue. Nothing
    here raises: a call GitHub refuses (a rule on the branch, a server
    error) is said in the line and left for the next run, since the board
    is written before this runs, and a failure that stopped the run would
    stop every board after it, and the comments it still owes. A rename
    GitHub refused for good is also said once on the issue, with the
    rulesets in its way, since nobody reads the workflow's log."""
    old, new = claims.branch_of(key), claims.released_branch(key, release.sha)
    try:
        done = free_name(api, old, new, release.sha)
    except GitHubError as problem:
        return f"{old}: not renamed to {new}: {problem}{github_words(problem)}" + \
            explain_refusal(api, key, old, new, problem)
    if done is None:
        return f"{old}: no longer at {release.sha[:7]}, so it keeps its name"
    try:
        comment_once(api, key, board.RELEASED_MARKER.format(branch=new), board.released_text(old, new))
    except GitHubError as problem:
        return f"#{key}: {done}, but the comment on the issue failed: {problem}{github_words(problem)}"
    return f"#{key}: {done}"


def free_name(api, old, new, sha):
    """Rename the claim branch old to new, reading it again first: None,
    and nothing done, when it moved since the snapshot, since a commit
    after the release means its holder came back. When new already keeps
    the same commit, old was pushed again at the released head, and the
    rename is refused for good; deleting old then frees the name and loses
    nothing. What was done, in words."""
    if branch_sha(api, old) != sha:
        return None
    try:
        api.post(f"/repos/{api.repo}/branches/{urllib.parse.quote(old, safe='')}/rename", {"new_name": new})
        return f"claim released; {old} renamed to {new}"
    except GitHubError as problem:
        if problem.status != 422 or branch_sha(api, new) != sha:
            raise
    api.delete(f"/repos/{api.repo}/git/refs/heads/{quoted(old)}")
    return f"claim released; {new} already keeps {sha[:7]}, so {old} was deleted"


def explain_refusal(api, key, old, new, problem):
    """Say once on the issue why its released claim's name is not free,
    when GitHub refused the rename for good: a ruleset over the claim
    branches, or a token without the right. Any other answer, or none, is
    for the moment, and the next run tries again without a word: a 403 of
    a rate limit's among them (rate_limited), which the client waited out
    as long as a run waits, since its comment would send the owner after
    a rule or a right that is not missing, and its marker would keep a
    later refusal for good unsaid. The rulesets are read only for the
    first saying, so a rename refused run after run costs one more call
    each time. What was said, for the log line, or ""."""
    if problem.status not in REFUSED_FOR_GOOD or rate_limited(problem.status, problem.headers, problem.detail):
        return ""
    marker = board.RENAME_REFUSED_MARKER.format(branch=new)
    words = f"HTTP {problem.status}{github_words(problem)}"
    try:
        if has_marker(api, key, marker):
            return ""
        text = board.rename_refused_text(old, new, rulesets_on(api, (old, new)), words)
        api.post(f"/repos/{api.repo}/issues/{key}/comments", {"body": marker + "\n" + text})
    except GitHubError as failure:
        return f"; saying so on #{key} failed: {failure}{github_words(failure)}"
    return f"; said so on #{key}"


def rulesets_on(api, branches):
    """The rulesets whose rules apply to any of the branches, each once,
    with the types of the rules it applies, as board.Ruleset holds them;
    None when GitHub would not say, which the comment then says instead."""
    applied = {}
    try:
        for branch in branches:
            for rule in api.get_all(f"/repos/{api.repo}/rules/branches/{urllib.parse.quote(branch, safe='')}"):
                source = (rule.get("ruleset_source_type") or "", rule.get("ruleset_source") or "",
                          rule.get("ruleset_id"))
                applied.setdefault(source, set()).add(rule.get("type") or "unnamed")
    except GitHubError:
        return None
    order = sorted(applied.items(), key=lambda item: tuple(str(part) for part in item[0]))
    return [ruleset(api, source, types) for source, types in order]


def ruleset(api, source, types):
    """One ruleset as the comment names it: its own name where the token
    may read the ruleset, else its number, and the page that changes it."""
    source_type, owner, number = source
    try:
        found = api.get(f"/repos/{api.repo}/rulesets/{number}", missing_ok=True) or {}
    except GitHubError:
        found = {}
    return board.Ruleset(found.get("name") or f"ruleset {number}", ruleset_link(source_type, owner, number),
                         tuple(sorted(types)))


def ruleset_link(source_type, owner, number):
    """The settings page of a repository's or an organization's ruleset,
    or "" for a source the board cannot link, such as an enterprise's."""
    if source_type == "Repository" and owner:
        return f"{WEB}/{owner}/settings/rules/{number}"
    if source_type == "Organization" and owner:
        return f"{WEB}/organizations/{owner}/settings/rules/{number}"
    return ""


def github_words(problem):
    """GitHub's own words on a refused call, after a colon, or nothing."""
    return f": {problem.detail}" if problem.detail else ""


# --- the gate's world ---------------------------------------------------------

def standing_word(pull, behind=None):
    """How a pull request's branch stands to its base, in the gate's four
    words, from the pull's mergeable_state as the board reads it. Where
    GitHub's word says only that something else blocks the merge (a rule,
    a review, a draft), mergeable tells a conflict apart, and behind, the
    commits of the base the head lacks, tells behind from clean; with
    neither known the answer is unknown."""
    if pull.get("mergeable_state") not in OPAQUE_STATES:
        return board.pr_state(pull)
    if pull.get("mergeable") is False:
        return gate.CONFLICTING
    if pull.get("mergeable") is None or behind is None:
        return gate.UNKNOWN
    return gate.BEHIND if behind > 0 else gate.CLEAN


def check_word(runs, name):
    """The required check's state on a commit, in the gate's words: none
    when no run has its name, red when a run of it finished without
    passing, pending while one is still running, green when every run of
    it passed. Red wins over pending, since a failed run already means
    the branch needs work."""
    mine = [run for run in runs if run.get("name") == name]
    if not mine:
        return gate.NO_CHECK
    if any(run.get("status") == "completed" and run.get("conclusion") not in PASSING for run in mine):
        return gate.RED
    if any(run.get("status") != "completed" for run in mine):
        return gate.PENDING
    return gate.GREEN


class GitHubWorld:
    """gate.World over one repository on GitHub.

    The policy says which check is required, who merges and which branch
    is the default one, so the world reads the same branch the gate's
    answers name; default_branch overrides it for a caller that knows
    better. The clock is an argument so a test can fix it.
    """

    def __init__(self, api, policy, default_branch=None, clock=None):
        self.api = api
        self.policy = policy
        self.default_branch = default_branch or policy.default_branch
        self.clock = clock or (lambda: datetime.now(timezone.utc))

    def main_sha(self):
        return self.branch_sha(self.default_branch)

    def branch_sha(self, branch):
        """The branch's commit, or "" when the repository has no such branch."""
        return branch_sha(self.api, branch)

    def pull_of(self, branch):
        """The branch's open pull request, settled, or None."""
        found = self.api.pull_detail(branch)
        return settled(self.api, found) if found else None

    def behind_by(self, pull):
        """How many commits of the pull request's base its head lacks, from
        GitHub's comparison of the two; None when there is nothing to
        compare or GitHub has no such comparison (a 404). Any other refusal
        raises."""
        return behind_by(self.api, pull, self.default_branch)

    def standing(self, branch):
        """Unknown for a branch with no open pull request: GitHub computes
        mergeability only for one, and the gate asks only about those.
        Where GitHub's word is opaque the branch is read through mergeable
        and one comparison (standing_of), a comparison GitHub refuses
        reading as unknown, the same whether GitHub had computed the pull
        request's mergeability at the first read or only after the pause."""
        pull = self.pull_of(branch)
        if pull is None:
            return gate.UNKNOWN
        found = standing_of(self.api, pull)
        return found if found is not None else standing_word(pull)

    def standing_detail(self, branch):
        """The branch's gate.Standing, which never says a clean branch is
        level with the default branch: GitHub's word clean says only that
        it merges as it stands, behind or not under a rule that does not
        ask for an up-to-date branch (F16)."""
        return gate.Standing(self.standing(branch), level=False)

    def check(self, branch):
        sha = self.branch_sha(branch)
        return check_word(self.api.check_runs(sha), self.policy.required_check) if sha else gate.NO_CHECK

    def merge(self, branch):
        """Merge the branch's pull request through GitHub's merge API, with
        the head the gate looked at, so a push that lands meanwhile is
        refused rather than merged unseen. Only where the policy says the
        gate merges: elsewhere the owner or auto-merge does, and a merge
        here would go behind their back."""
        if self.policy.who_merges != gate.GATE:
            raise NotAllowed(f"{gate.MOVERS[self.policy.who_merges]} merges here, not the gate; "
                             f"{branch} waits for it")
        pull = self.api.pull_detail(branch)
        if pull is None:
            raise MergeRefused(f"{branch} has no open pull request to merge; open one first")
        try:
            done = self.api.put(f"/repos/{self.api.repo}/pulls/{pull['number']}/merge",
                                {"sha": pull["head"]["sha"], "merge_method": "merge"})
        except GitHubError as problem:
            raise self.refusal(branch, pull, problem) from None
        return MergeResult(done["sha"])

    def refusal(self, branch, pull, problem):
        """What a refused merge means for the gate. GitHub answers 405 for
        a pull request it will not merge and 409 when the head moved: main
        moved when the branch, read again, no longer fits the default
        branch, and otherwise the author has something to fix (a review, a
        rule, a push to look at), which is also the answer when GitHub will
        not compare the branch, since the standing is then unknown. Any
        other answer is the gate failing."""
        said = f": {problem.detail}" if problem.detail else ""
        if problem.status not in (405, 409):
            return GateFailed(f"GitHub did not merge #{pull['number']}: {problem}{said}")
        if self.standing(branch) in (gate.BEHIND, gate.CONFLICTING):
            return MainMoved(f"{branch} no longer merges cleanly into {self.default_branch}{said}")
        return MergeRefused(f"GitHub would not merge #{pull['number']}{said or ' as it stands'}")

    def now(self):
        return self.clock()
