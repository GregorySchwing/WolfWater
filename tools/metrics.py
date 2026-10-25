"""Field metrics: how the claim protocol runs in this repository, as numbers.

OpenSpineConsortium/coordbench measures ways of coordinating coding agents
on CooperBench's task pairs, this protocol among them; this tool measures
the same protocol where it is used for real, in aggregates only, so a
project can share how the protocol ran without sharing anything of its
own. docs/setup/metrics.md is the reader's page: every key, a sample of
each file, how each measure lines up with the benchmark's, and the switch.

  collect  run daily by .github/workflows/metrics.yml with GITHUB_TOKEN
           alone. It reads the window's pull requests, reviews, check
           runs, claim branches and comments from GitHub's REST API and
           writes METRICS.json: counts, durations in seconds and two rates,
           and the same pull requests counted again by work setting and in
           pairs (what else was at work beside each, and whose: one Claude
           Code session's, one person's several sessions', several
           people's), under the keys below and no other. No name, login,
           session link, title, branch name, file name, URL, issue number or
           text is ever written: the logins, session lines and changed files
           it reads to sort the pull requests stay in the run's memory, so
           the file says how the protocol ran and nothing about what was
           built or by whom. The workflow puts it on the branch metrics,
           where OpenSpineConsortium's collector reads it, and only where
           the repository installed the Symphony App. GITHUB_TOKEN's hour
           is the whole repository's, the board's runs among it, so a run
           stops before it would leave the board fewer requests than
           RESERVE, says so, and writes nothing. The breakdown reads only
           once every other key is measured, so it never costs one: when
           its reads would reach the reserve, or GitHub fails one of them,
           the record is written without it, and when only the changed
           files, read last, are left, with the counts that need them null.
           With --profile-out it writes PROFILE.json too, beside the
           breakdown and left out with it: the Symphony-Profile lines on the
           commits of the pull requests that ended in the window and the
           profiles of the Actions runs that name them, summed by work
           setting and alone against parallel (tools/profile.py, group and
           coordination). The workflow puts it on the branch metrics too;
           the collector does not read it, since its keys are the
           repository's own until the Apps' privacy statement lists them.
  judge    run by a maintainer on their own machine, with their own gh and
           their own Claude Code (claude -p) on their own subscription, never
           with ANTHROPIC_API_KEY. It reads each recent pull request's
           conversation, has the local model code it with CooperBench's App F
           taxonomy of communication errors, and keeps only how many
           conversations fell in each code, in all and by the work setting
           of their pull request, which never reaches the model. The
           conversation text goes to the maintainer's local claude process
           and nowhere else, which asks the model on their own account as
           any of their sessions does: JUDGE.json holds counts alone, and
           --push puts it on the branch metrics beside METRICS.json as an
           ordinary commit, never a forced one. It refuses while any
           variable that would bill an API or a cloud account is set, and
           --push only in a clone whose origin is the repository judged,
           since JUDGE.json names no repository, by design.
  push     put a JUDGE.json the judge already wrote on the branch metrics,
           after the same origin check, so a push that was refused, or a
           run kept without --push, never costs a second judgement.
  unchanged OLD NEW  exit 0 when two METRICS.json files, or two
           PROFILE.json files, differ in nothing but generated_at, so the
           workflow writes only when something moved.

The switch: sharing is on unless the file .github/symphony-metrics (or
.github/sdlc-metrics, its name before the rename, which a repository set up
before it carries) or the repository variable SDLC_METRICS (handed to the
workflow's environment) says otherwise. Off from any of them wins, and an
off run writes exactly {"schema": "sdlc-metrics/1", "sharing": "off"}, and
{"schema": "symphony-profile-sums/1", "sharing": "off"} as PROFILE.json,
without a single API call.

Standard library only, like every tool here. GitHub's REST calls go through
the vendored protocol library's client, the one the board uses, so there is
one way of calling GitHub; the judge reads through gh instead, because it
runs as the maintainer, signed in as themselves.
"""

import argparse
import http.client
import io
import json
import math
import os
import re
import shlex
import shutil
import subprocess
import sys
import tempfile
import time
import urllib.error
import urllib.parse
import urllib.request
import zipfile
from dataclasses import dataclass, field
from datetime import datetime, timedelta, timezone
from itertools import chain, combinations
from pathlib import Path
from typing import NamedTuple

TOOLS = Path(__file__).resolve().parent
sys.path.insert(0, str(TOOLS))

import economics  # noqa: E402
import usage  # noqa: E402
from sdlc_protocol import board, claims  # noqa: E402
from sdlc_protocol.adapters.github import TIMEOUT, USER_AGENT, GitHub, GitHubError, error_detail, settled  # noqa: E402

# tools/profile.py, loaded by its path, since the standard library has a
# module of the same name.
profile = usage.profiler()

SCHEMA = "sdlc-metrics/1"
JUDGE_SCHEMA = "sdlc-judge/1"
ON, OFF = "on", "off"
OFF_RECORD = {"schema": SCHEMA, "sharing": OFF}
SWITCH_FILE = ".github/symphony-metrics"
# The switch's name before the rename to Symphony, still read, since a
# repository set up before it carries that file.
OLD_SWITCH_FILE = ".github/sdlc-metrics"
SWITCH_VARIABLE = "SDLC_METRICS"
METRICS_BRANCH = "metrics"
METRICS_FILE = "METRICS.json"
JUDGE_FILE = "JUDGE.json"
# The profile's sums, beside METRICS.json on the branch metrics, which the
# collector does not read: its keys are the repository's own until the
# Apps' privacy statement lists them.
PROFILE_FILE = "PROFILE.json"
PROFILE_SCHEMA = "symphony-profile-sums/1"
PROFILE_OFF = {"schema": PROFILE_SCHEMA, "sharing": OFF}
VENDORED = TOOLS / "VENDORED"
WINDOW_DAYS = 30

# The keys of METRICS.json, in order; the file holds these and nothing else.
COUNT_KEYS = ("claims", "pull_requests_opened", "merged", "closed_unmerged", "reviews", "changes_requested",
              "check_runs_passed", "check_runs_failed", "conflicts_seen", "released_claims", "contributors")
DURATIONS = ("claim_to_first_pull_request", "open_to_merge", "first_review")
DURATION_KEYS = tuple(f"{name}_{kind}" for name in DURATIONS for kind in ("median", "p90"))
RATE_KEYS = ("merge_rate", "first_check_pass_rate")
# The work settings, in order. A person's pull request is in one of the
# first six, by what else was at work beside it and whose that work was; an
# app's or a bot's is in one of the last two, which are nobody's partner.
SETTINGS = ("solo", "one_session", "one_person_sessions", "several_people", "several_people_and_sessions",
            "undetermined", "workflow_agents", "automation")
PERSON_SETTINGS = SETTINGS[:6]
# What by_setting counts for each work setting.
SETTING_KEYS = ("pull_requests", "merged", "closed_unmerged", "open", "conflicted", "released", "changes_requested",
                "claim_branches", "partners", "work_seconds", "landings_under", "landings_before_opening",
                "file_overlap")
# How two people's pull requests at work together stand to each other, and
# what pairs counts for each.
RELATIONS = ("same_session", "one_person_sessions", "several_people", "undetermined")
PAIR_KEYS = ("pairs", "both_active", "file_overlap", "landings", "landings_file_overlap", "conflict_onsets",
             "conflict_onsets_file_overlap")
# The pair counts that need the changed files, null with them when the
# board's reserve left some of them unread.
FILE_PAIR_KEYS = ("file_overlap", "landings_file_overlap", "conflict_onsets_file_overlap")
# What parallel counts for the repository: who and what was at work.
PARALLEL_KEYS = ("people", "sessions", "sessions_per_person_max", "pull_requests_with_session",
                 "pull_requests_without_session", "pull_requests_multi_session", "sessions_across_people",
                 "peak_parallel", "peak_parallel_people", "peak_sessions_one_person", "conflict_marks",
                 "conflict_marks_unattributed", "open_conflicting", "open_conflicting_drafts",
                 "open_conflicting_dirty", "truncated", "files_unread")
# Every name in these tuples is a constant of the tool and never comes from
# data, so the file's words are the same in every repository's.
RECORD_KEYS = ("schema", "generated_at", "window_days", "sharing", "protocol_version", "counts",
               "durations_seconds", "rates", "parallel", "by_setting", "pairs", "economics", "coordination")
# The keys written with the breakdown, and left out with it.
BREAKDOWN_KEYS = RECORD_KEYS[-5:]
# The keys of PROFILE.json, in order; the last four are measured with the
# breakdown and left out with it.
PROFILE_KEYS = ("schema", "generated_at", "window_days", "sharing", "sessions", "actions_runs",
                "actions_runs_unread", "by_setting", "coordination")

# A check run counts as passed or failed by its conclusion; neutral,
# skipped, cancelled, stale and action_required say neither, and a run
# still going has no conclusion yet.
PASSED = {"success"}
FAILED = {"failure", "timed_out"}
# A commit's date is the committer's clock, not GitHub's: a minute's grace
# keeps a slightly fast clock from hiding the head a pull request opened with.
CLOCK_GRACE = timedelta(seconds=60)
# The coordinator routine starts every release comment with this line
# (docs/agents/routines/coordinator.md, step 2), so releases can be counted
# without reading what the comment says. It is the vendored protocol's, the
# line its board renames a released claim's branch on.
RELEASE_MARKER = claims.RELEASE_MARKER
# The board's "main moved" comment carries a marker and says "behind" or
# "conflicting with"; only the second is a conflict. Both pieces come from
# the vendored board, so a change of its words there changes them here.
MAIN_MOVED_PREFIX = board.MAIN_MOVED_MARKER.split("{sha}")[0]
CONFLICT_WORDS = board.main_moved_text("conflicting", "{base}", "{sha}").split("{base}")[0]
DIRTY = "dirty"
# The kind of a pull request a person opened. One an app or a bot opened is
# the repository's own Claude workflows at work when its branch starts with
# AGENT_BRANCH, the prefix the claim protocol and claude-code-action share,
# so the workers in Actions (workers.yml), the coordinator (coordinator.yml)
# and @claude among them; any other is automation, such as dependabot or
# the adoption App.
PERSON = "person"
AGENT_BRANCH = "claude/"
# An author date further back than this before the opening belongs to an
# old commit carried in, not to the start of the work.
WORK_HORIZON = timedelta(days=7)
# A Claude Code session is named by a whole line of a commit message or a
# description, never by a link inside a sentence: the trailer a cloud
# session ends its commits with (its key in any letter case) or its link
# alone on a line; the footer a cloud session ends its descriptions with;
# and the line a session with no claude.ai link, a local one or a claude -p
# worker, can write to name itself, whose name is kept behind TOKEN_PREFIX
# so it never meets a claude.ai id.
LINK_LINE = re.compile(r"(?im)^[ \t]*(?:claude-session:[ \t]*)?https://claude\.ai/code/(session_[A-Za-z0-9]{1,64})"
                       r"[ \t]*$")
FOOTER_LINE = re.compile(r"(?m)^[ \t]*_Generated by \[Claude Code\]\(https://claude\.ai/code/"
                         r"(session_[A-Za-z0-9]{1,64})\)_[ \t]*$")
TOKEN_LINE = re.compile(r"(?im)^[ \t]*symphony-session:[ \t]*([A-Za-z0-9][A-Za-z0-9_-]{7,63})[ \t]*$")
TOKEN_PREFIX = "token:"
# The board's mark of the default branch moving under a pull request, with
# the short sha it moved to, which names the merge that moved it. It is
# built from the vendored marker, so a change of its words in the protocol
# changes it here.
MARK_PREFIX, MARK_SUFFIX = board.MAIN_MOVED_MARKER.split("{sha}")
MARK = re.compile(re.escape(MARK_PREFIX) + r"([0-9a-f]{7,40}|unknown)" + re.escape(MARK_SUFFIX))
# GITHUB_TOKEN has 1,000 requests an hour for the whole repository, and
# the board's runs draw on the same hour. A run stops before it would leave
# fewer than this many, writes nothing, and the branch metrics keeps its
# last record until a run that fits; only the breakdown's reads, made last,
# give way instead, and the record is written without what they would have
# read. Each answer's headers say how many are left and when the hour turns.
HOURLY_BUDGET = 1000
RESERVE = 300
REMAINING_HEADER = "X-RateLimit-Remaining"
RESET_HEADER = "X-RateLimit-Reset"
# What Claude spent in GitHub Actions comes from the record each run
# uploads as an artifact named with this prefix (tools/usage.py). Each
# download costs a call of the hour, so a run reads at most USAGE_LIMIT new
# ones, newest first, and keeps what it read of each, the parts
# economics.kept_record keeps, in the file --usage-cache names, which the
# workflow carries from one run to the next, so the next run reads only
# what is new; an artifact or a file bigger than USAGE_MAX_BYTES is no
# record of ours and is left unread. The cache's version rises when the
# kept form keeps more, as version 2 keeps each run's profile. A record an
# older version kept is still read, without what that version dropped, and
# is marked old: its artifact is downloaded again once the new ones are
# read, within the same USAGE_LIMIT, while it is listed and unexpired, and
# the old record stays where that download fails or the artifact has gone,
# so no run drops out of the economics because the cache's form changed.
USAGE_PREFIX = "symphony-usage-"
USAGE_LIMIT = 200
USAGE_MAX_BYTES = 65536
USAGE_CACHE_VERSION = 2
# The versions whose records a run still reads (read_usage_cache).
USAGE_CACHE_VERSIONS = (1, 2)
REDIRECTS = (301, 302, 303, 307, 308)
# GitHub's revert button opens a pull request whose title starts so.
REVERT_TITLE = 'Revert "'


# --- time --------------------------------------------------------------------

def parse_time(text):
    """GitHub's ISO time as an aware datetime, or None for none."""
    if not text:
        return None
    try:
        moment = datetime.fromisoformat(str(text).replace("Z", "+00:00"))
    except ValueError:
        return None
    return moment if moment.tzinfo else moment.replace(tzinfo=timezone.utc)


def iso(moment):
    return moment.astimezone(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


class Window(NamedTuple):
    """The days a run looks back over, ending at its own moment."""

    start: datetime
    end: datetime
    days: int

    @classmethod
    def ending(cls, end, days):
        return cls(end - timedelta(days=days), end, days)

    def holds(self, moment):
        return moment is not None and self.start <= moment <= self.end


def look_back_start(window):
    """Where the look-back begins: the window's finished pull requests are
    read back one more window, so their life before it (what was at work
    beside them, the board's marks, the coordinator's releases) is whole.
    Nothing the existing keys read changes."""
    return window.start - timedelta(days=window.days)


def time_of(item, key):
    return parse_time((item or {}).get(key))


def seconds_between(earlier, later):
    """Whole seconds from one moment to a later one; None when either is
    missing or the order is wrong, which a skewed clock can make it."""
    if earlier is None or later is None or later < earlier:
        return None
    return int((later - earlier).total_seconds())


# --- the switch --------------------------------------------------------------

def first_word(text):
    words = (text or "").split()
    return words[0].lower() if words else ""


def sharing(root, environ):
    """'on' or 'off'. The variable turns sharing off with any word but on,
    and says nothing when empty; the file, under its name or its old one,
    says on only with the word on first, says nothing when missing, and is
    off when it says anything else or cannot be read. Off from any of them
    wins, and none saying anything is on: a repository set up before the
    rename may carry both files, and an off in either is the owner's word,
    which an on in the other never lifts. Each doubt falls to off, since
    off shares nothing."""
    said = first_word(environ.get(SWITCH_VARIABLE))
    if said and said != ON:
        return OFF
    for name in (SWITCH_FILE, OLD_SWITCH_FILE):
        try:
            text = (Path(root) / name).read_text(encoding="utf-8")
        except FileNotFoundError:
            continue
        except (OSError, UnicodeDecodeError):
            return OFF
        if first_word(text) != ON:
            return OFF
    return ON


def protocol_version(path=VENDORED):
    """The vendored protocol's release and short commit, "0.2.0 f530a66",
    from tools/VENDORED's first line; None when it cannot be read."""
    try:
        first = path.read_text(encoding="ascii").splitlines()[0].split()
    except (OSError, UnicodeDecodeError, IndexError):
        return None
    if len(first) != 3 or first[0] != "sdlc-protocol":
        return None
    return f"{first[1]} {first[2][:7]}"


# --- reading GitHub ----------------------------------------------------------

def header_number(headers, name):
    """A header's whole number, or None when it is missing or not one."""
    try:
        return int((headers or {}).get(name))
    except (TypeError, ValueError):
        return None


def refusal_detail(problem):
    """GitHub's words on a refused call, as error_detail reads them, or ""
    when the refusal's body was cut short: reading it raises http.client's
    IncompleteRead, which is no OSError, and a call GitHub refused is a
    refusal whatever became of its words."""
    try:
        return error_detail(problem)
    except http.client.HTTPException:
        return ""


class OutOfBudget(GitHubError):
    """The hour's requests are down to the board's reserve."""


class Unfollowed(urllib.request.HTTPRedirectHandler):
    """Leaves a redirect for the caller, so the address it gives is fetched
    without the token: urllib sends a redirected request every header of
    the first, Authorization among them, to whatever host it names. It
    reads nothing of the redirect, not even its address, which urllib
    parses before asking whether to follow it and raises ValueError on
    when it cannot, leaving the answer open: the caller reads it."""

    def http_error_302(self, req, fp, code, msg, headers):
        return None

    http_error_301 = http_error_303 = http_error_307 = http_error_308 = http_error_302


class Budgeted(GitHub):
    """The protocol library's client with the hour's budget watched: each
    answer says how many requests are left, and the next call is refused
    before it is made once that is down to the reserve. GitHub's own
    refusal of a spent hour, a 403, is read the same way, so a run that
    starts after the board spent the hour says so too. The calls it
    reports are every call GitHub answered: the client's own count (made),
    which holds each try it makes again after a rate limit's wait, and the
    redirects location asks for itself (the first rate-limit review's sixth
    finding: a count of this class's own requests left the tries out)."""

    def __init__(self, token, repo, reserve=RESERVE, **kwargs):
        super().__init__(token, repo, **kwargs)
        self.reserve, self.located, self.remaining, self.reset = reserve, 0, None, None
        self.unfollowed = kwargs.get("opener") or urllib.request.build_opener(Unfollowed)

    @property
    def calls(self):
        return sum(self.made.values()) + self.located

    def request(self, method, url, body=None, missing_ok=False):
        if self.remaining is not None and self.remaining <= self.reserve:
            raise OutOfBudget(self.refusal())
        try:
            data, headers = super().request(method, url, body, missing_ok)
        except GitHubError as problem:
            if problem.status == 403 and "api rate limit exceeded" in problem.detail.lower():
                self.remaining = 0
                raise OutOfBudget(f"GitHub refused the call ({problem.detail}); {self.refusal()}") from None
            raise
        self.counted(headers)
        return data, headers

    def counted(self, headers):
        if (remaining := header_number(headers, REMAINING_HEADER)) is not None:
            self.remaining, self.reset = remaining, header_number(headers, RESET_HEADER)

    def location(self, path):
        """Where GitHub sends the download at path, asked without following
        the redirect, since the file is served from another host, which is
        never sent the token; None when GitHub answers with no redirect.
        One call of the hour, counted like any other. It asks around the
        library's request, so it says each failure as a GitHubError itself,
        an answer that is not HTTP among them: http.client raises its own
        exceptions there, which are no OSError, and urllib passes them on
        as they are."""
        if self.remaining is not None and self.remaining <= self.reserve:
            raise OutOfBudget(self.refusal())
        request = urllib.request.Request(self.url(path), headers=self.headers())
        try:
            with self.unfollowed.open(request, timeout=TIMEOUT) as answer:
                headers, target = answer.headers, None
        except urllib.error.HTTPError as problem:
            if problem.code not in REDIRECTS:
                self.located += 1
                detail = refusal_detail(problem)
                if problem.code == 403 and "api rate limit exceeded" in detail.lower():
                    self.remaining = 0
                    raise OutOfBudget(f"GitHub refused the call ({detail}); {self.refusal()}") from None
                raise GitHubError(f"GET {path}: HTTP {problem.code}", problem.code, detail) from None
            headers, target = problem.headers, problem.headers.get("Location")
        except OSError as problem:
            raise GitHubError(f"GET {path}: {getattr(problem, 'reason', problem)}") from None
        except http.client.HTTPException as problem:
            raise GitHubError(f"GET {path}: {problem!r}") from None
        self.located += 1
        self.counted(headers)
        return target

    def download(self, url, limit):
        """The bytes at a storage address a redirect gave, fetched with no
        token at all, or None when they are more than limit; no call of the
        hour, since the storage host is not GitHub's API. An answer cut
        short is a GitHubError, as a refused one is, so the record counts as
        unread: http.client raises its own IncompleteRead, which is no
        OSError, for one sent in chunks and its BadStatusLine for one that
        is not HTTP, and gives what came of one whose length it was told,
        read up to a limit, without a word. The address is fetched only
        over https and asked, as location asks, without following a
        redirect, which would take it wherever the storage host said, plain
        http among the places; and an address urllib cannot read, where it
        raises ValueError, is a GitHubError too."""
        try:
            if urllib.parse.urlsplit(url).scheme != "https":
                raise GitHubError("a download address that is not https is not fetched")
            request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
            with self.unfollowed.open(request, timeout=TIMEOUT) as answer:
                data, told = answer.read(limit + 1), header_number(answer.headers, "Content-Length")
        except urllib.error.HTTPError as problem:
            problem.close()
            raise GitHubError(f"download: HTTP {problem.code}") from None
        except OSError as problem:
            raise GitHubError(f"download: {getattr(problem, 'reason', problem)}") from None
        except http.client.HTTPException as problem:
            raise GitHubError(f"download: {problem!r}") from None
        except ValueError as problem:
            raise GitHubError(f"download: an address urllib cannot read ({problem})") from None
        if told is not None and len(data) < min(told, limit + 1):
            raise GitHubError(f"download: the answer was cut short, {len(data)} of {told} bytes")
        return data if len(data) <= limit else None

    def refusal(self):
        back = f"at {iso(datetime.fromtimestamp(self.reset, timezone.utc))}" if self.reset else "within the hour"
        return (f"stopped after {self.calls} calls, with {self.remaining} of GITHUB_TOKEN's {HOURLY_BUDGET:,} "
                f"requests an hour left for the board, whose reserve is {self.reserve}: nothing is written, and "
                f"the branch {METRICS_BRANCH} keeps its last record. The hour's budget comes back {back}; run "
                "the workflow again once it has, at a quiet moment.")


def pull_listing(api):
    """Every pull request, most lately updated first, one at a time. A page
    is asked for only once the one before it is read through, so a reader
    that stops at the window asks for no page more than it needs, and can
    come back later to read on into the look-back."""
    for page in api.pages(f"/repos/{api.repo}/pulls", {"state": "all", "sort": "updated", "direction": "desc"}):
        yield from page


def updated_since(listing, since):
    """The pull requests the listing gives next whose last update is at or
    after since, or unknown, and the rest of the listing from the first one
    updated before it. GitHub lists them by last update, newest first, so
    that one ends them; a creation, a merge, a close or a review each moves
    the update time."""
    found = []
    for pull in listing:
        updated = time_of(pull, "updated_at")
        if updated is not None and updated < since:
            return found, chain([pull], listing)
        found.append(pull)
    return found, iter(())


def once_each(items, seen, number=lambda pull: pull["number"]):
    """The items whose pull request has not been listed before, each from
    its first listing, the numbers noted in seen. GitHub pages its listing
    by page number, so an update while it is read can push a pull request
    onto the next page a second time."""
    found = []
    for item in items:
        if number(item) not in seen:
            seen.add(number(item))
            found.append(item)
    return found


def listed_pulls(api, window):
    """Every pull request changed in the window (recent), and every one
    changed in the look-back before it (earlier), each once. The
    look-back's pull requests are read only to see what was at work beside
    the window's before it began."""
    recent, rest = updated_since(pull_listing(api), window.start)
    earlier, _ = updated_since(rest, look_back_start(window))
    seen = set()
    return once_each(recent, seen), once_each(earlier, seen)


@dataclass
class PullFacts:
    """What one pull request contributes: the pull as listed, its commits
    (oldest first), its submitted reviews, the check runs on each of its
    commits made in the window, by commit, the check runs on the head it
    was opened with (None when that head cannot be told) and, while it is
    open, GitHub's mergeable_state and mergeable."""

    pull: dict
    commits: list = field(default_factory=list)
    reviews: list = field(default_factory=list)
    runs: dict = field(default_factory=dict)
    first_runs: list = None
    mergeable_state: str = None
    mergeable: bool = None


def commit_time(commit, who):
    return time_of(((commit or {}).get("commit") or {}).get(who), "date")


def first_head_sha(pull, commits):
    """The head the pull request was opened with. GitHub keeps no record of
    it, so it is the newest of the pull request's commits committed by the
    time it opened; None when a rebase has replaced them all, and the pull
    request then stays out of the first-check rate."""
    cutoff = time_of(pull, "created_at")
    if cutoff is None:
        return None
    before = [c for c in commits if (when := commit_time(c, "committer")) and when <= cutoff + CLOCK_GRACE]
    return max(before, key=lambda c: commit_time(c, "committer"))["sha"] if before else None


def every_check_run(api, sha):
    """Every check run on one commit, earlier attempts included: GitHub's
    default listing, the latest attempt of each check, hides a failed run
    behind its re-run, and a failure that was run again is one
    check_runs_failed is for. The vendored client's check_runs asks with
    no filter, so the endpoint is paged here, the way it pages it."""
    runs = []
    for page in api.pages(f"/repos/{api.repo}/commits/{sha}/check-runs", {"filter": "all"}):
        if isinstance(page, dict):
            runs.extend(page.get("check_runs") or [])
    return runs


def pull_facts(api, pull, window):
    """The calls one pull request takes, so the budget can be reckoned
    (RESERVE): one for its commits, one for its reviews, one for the check
    runs of each commit made in the window, one more for the head it was
    opened with when that was made before the window, and while it is
    open one for its detail, two when GitHub has not computed its
    mergeability yet. The runs are read on every commit made in the
    window, not the head alone: the default branch takes only a pull
    request whose check passed, so a merged one's head is green by
    construction, and the failures its earlier heads had would never be
    seen."""
    path = f"/repos/{api.repo}/pulls/{pull['number']}"
    facts = PullFacts(pull, api.get_all(path + "/commits"),
                      [r for r in api.get_all(path + "/reviews") if r.get("state") != "PENDING"])
    made = [c["sha"] for c in facts.commits if c.get("sha") and window.holds(commit_time(c, "committer"))]
    facts.runs = {sha: every_check_run(api, sha) for sha in made}
    if window.holds(time_of(pull, "created_at")):
        first = first_head_sha(pull, facts.commits)
        if first is not None:
            facts.first_runs = facts.runs[first] if first in facts.runs else every_check_run(api, first)
    if pull.get("state") == "open":
        detail = settled(api, api.get(path))
        facts.mergeable_state = detail.get("mergeable_state")
        facts.mergeable = detail.get("mergeable")
    return facts


def claim_number(pull):
    return claims.number_of((pull.get("head") or {}).get("ref") or "")


def branch_claims(api, window, known):
    """The issue numbers of claude/<n> branches with no pull request among
    the window's whose head commit was made in the window. GitHub's REST API
    says when a commit was made, never when a branch was, so each such
    branch costs one call. A closed issue's branch is not skipped: tidy.yml
    deletes it when the issue closes, so nearly every branch here is an
    open issue's, and the listing of open issues that would tell them apart
    would cost calls and save none."""
    found = set()
    for branch in api.get_all(f"/repos/{api.repo}/branches"):
        number = claims.number_of(branch.get("name") or "")
        sha = (branch.get("commit") or {}).get("sha")
        if number is None or number in known or not sha:
            continue
        if window.holds(commit_time(api.get(f"/repos/{api.repo}/commits/{sha}", missing_ok=True), "committer")):
            found.add(number)
    return found


def window_comments(api, window):
    """The issue and pull request comments made in the window. GitHub's
    since is the last update, so an old comment edited lately is listed and
    then left out here."""
    listed = api.get_all(f"/repos/{api.repo}/issues/comments",
                         {"since": iso(window.start), "sort": "created", "direction": "asc"})
    return [c for c in listed if window.holds(time_of(c, "created_at"))]


def look_back_comments(api, window):
    """The comments made in the look-back, before the window. GitHub's
    since is the last update and the listing goes from the oldest made, so
    it starts with older comments edited lately, which are left out, and
    ends at the first comment made in the window, which window_comments
    already has."""
    start = look_back_start(window)
    found = []
    for page in api.pages(f"/repos/{api.repo}/issues/comments",
                          {"since": iso(start), "sort": "created", "direction": "asc"}):
        for c in page:
            made = time_of(c, "created_at")
            if made is not None and made >= window.start:
                return found
            if made is not None and made >= start:
                found.append(c)
    return found


def read_repository(api, window):
    """What the keys before the breakdown measure, read from GitHub in the
    calls they always took, in the same order: the window's pull requests
    with what each holds, the claims on branches without one and the
    window's comments; and the listing from where it stopped, at the first
    pull request changed before the window, for the breakdown to read on."""
    recent, rest = updated_since(pull_listing(api), window.start)
    facts = [pull_facts(api, pull, window) for pull in recent]
    known = {n for f in facts if (n := claim_number(f.pull)) is not None}
    return facts, branch_claims(api, window, known), window_comments(api, window), rest


# --- work settings -------------------------------------------------------------

@dataclass(eq=False)
class Unit:
    """One pull request as the breakdown sorts it, held in memory only: who
    opened it ("login:" and the lowercased login, "" for an app or a bot),
    its kind, the sessions it names, its opening, the start of its work, its
    merge and its end (None while open), its merge commit, whether its base
    is the default branch, its claim, its branch (the head's repository and
    name), its commits' author dates, whether it counts as changed in the
    window (window_unit), whether it is a draft, GitHub's mergeability while
    it is open, and the reviews by others. Only start changes once it is
    built, and units are told apart by identity."""

    number: int
    person: str
    kind: str
    sessions: frozenset
    created: datetime
    start: datetime
    merged: datetime
    end: datetime
    merge_sha: str
    into_default: bool
    claim: int
    branch: tuple
    authored: tuple
    recent: bool
    draft: bool
    mergeable: bool
    mergeable_state: str
    reviews: tuple


def well_formed(pull, recent):
    """Whether a pull request is a unit: opened, and ended (merged or
    closed), or open and changed in the window (recent, as window_unit
    tells it). One GitHub gives no end for is none, and so is one open and
    neither opened nor given a commit in the window, which is stale and at
    work beside nobody."""
    if time_of(pull, "created_at") is None:
        return False
    if time_of(pull, "merged_at") or time_of(pull, "closed_at"):
        return True
    return pull.get("state") == "open" and recent


def kind_of(pull):
    """PERSON for a person's pull request, and for an app's or a bot's the
    setting it counts under, by its branch."""
    if person(pull.get("user")):
        return PERSON
    return "workflow_agents" if ((pull.get("head") or {}).get("ref") or "").startswith(AGENT_BRANCH) else "automation"


def session_names(text):
    """The sessions a commit message or a description names, each on a line
    of its own. A line ended with CRLF, as a description saved from
    GitHub's web page is, counts as the same line ended with LF."""
    text = (text or "").replace("\r\n", "\n")
    names = {found.group(1) for line in (LINK_LINE, FOOTER_LINE) for found in line.finditer(text)}
    names.update(TOKEN_PREFIX + found.group(1) for found in TOKEN_LINE.finditer(text))
    return names


def sessions_of(pull, commits):
    """Every session the pull request's commit messages and description
    name. A commit's author never stands for a session, nor for a person:
    every Claude Code cloud session commits as the one account claude."""
    texts = [((c or {}).get("commit") or {}).get("message") for c in commits] + [pull.get("body")]
    return frozenset(name for text in texts for name in session_names(text))


def into_default_branch(pull):
    """Whether the pull request's base is the repository's default branch,
    which the listing says: only a merge there moves the branch the board
    marks the others against."""
    base = pull.get("base") or {}
    default = (base.get("repo") or {}).get("default_branch")
    return bool(default) and base.get("ref") == default


def work_start(created, authored):
    """When the work began: the earliest author date of its commits, or its
    opening when that is earlier. An author date more than WORK_HORIZON
    before the opening is an old commit carried in, and set aside."""
    return min([created] + [d for d in authored if d >= created - WORK_HORIZON])


def unit_of(pull, commits, recent, reviews=(), mergeable=None, mergeable_state=None):
    """The pull request as a unit, or None when it is not well formed."""
    if not well_formed(pull, recent):
        return None
    created, merged = time_of(pull, "created_at"), time_of(pull, "merged_at")
    authored = tuple(sorted(d for c in commits if (d := commit_time(c, "author"))))
    head = pull.get("head") or {}
    return Unit(number=pull["number"], person=person(pull.get("user")) or "", kind=kind_of(pull),
                sessions=sessions_of(pull, commits), created=created, start=work_start(created, authored),
                merged=merged, end=merged or time_of(pull, "closed_at"),
                merge_sha=pull.get("merge_commit_sha") or "", into_default=into_default_branch(pull),
                claim=claim_number(pull),
                branch=(((head.get("repo") or {}).get("full_name")) or "", head.get("ref") or ""),
                authored=authored, recent=recent, draft=bool(pull.get("draft")), mergeable=mergeable,
                mergeable_state=mergeable_state, reviews=tuple(reviews))


def floor_starts(units):
    """Raise each unit's start to the latest end of another pull request
    from the same branch that ended by its opening: a branch used again
    lists commits its earlier pull request already landed."""
    for unit in units:
        ends = [other.end for other in units if other is not unit and other.branch == unit.branch
                and other.end is not None and other.end <= unit.created]
        if ends:
            unit.start = max(unit.start, max(ends))


def stop(unit, now):
    """When the unit's work stopped: its end, or now while it is open."""
    return unit.end or now


def co_active(a, b, now):
    """Whether two people's pull requests were at work at once for more
    than a clock's grace; an app's or a bot's is at work beside nobody."""
    return (a is not b and a.kind == PERSON and b.kind == PERSON
            and min(stop(a, now), stop(b, now)) - max(a.start, b.start) > CLOCK_GRACE)


def partners_of(units, now):
    """Each unit's number, mapped to the units at work beside it."""
    return {unit.number: [other for other in units if co_active(unit, other, now)] for unit in units}


def relation(a, b):
    """How two people's pull requests at work together stand: another
    person's; the same person's from the same session or from another; or
    the same person's where one names no session, so it could be either."""
    if a.person != b.person:
        return "several_people"
    if not a.sessions or not b.sessions:
        return "undetermined"
    return "same_session" if a.sessions & b.sessions else "one_person_sessions"


def setting_of(unit, partners):
    """The unit's work setting, decided by its relations to its partners in
    this order: another person's work and the same person's other session's
    both beside it; the other session's; a session that cannot be told;
    another person's; and the same session's alone."""
    if unit.kind != PERSON:
        return unit.kind
    if not partners:
        return "solo"
    kinds = {relation(unit, other) for other in partners}
    if {"one_person_sessions", "several_people"} <= kinds:
        return "several_people_and_sessions"
    for name in ("one_person_sessions", "undetermined", "several_people"):
        if name in kinds:
            return name
    return "one_session"


def window_unit(pull, commits, window, reviews=(), mergeable=None, mergeable_state=None):
    """A pull request listed in the window as a unit, or None. An open one
    counts as changed in the window only when it was opened or given a
    commit in it, by its own work: the board comments on a pull request
    behind the default branch or in conflict with it at every push there,
    which renews GitHub's updated_at on one nobody has touched for months,
    and it would otherwise be at work beside everything since it opened."""
    changed = (pull.get("state") != "open" or window.holds(time_of(pull, "created_at"))
               or any(made_in(c, window) for c in commits))
    return unit_of(pull, commits, changed, reviews, mergeable, mergeable_state)


def earlier_units(api, earlier, units, window):
    """The look-back's pull requests as units, beside the window's (units).
    Commits cost a call each, so they are read only for a person's pull
    request that ended more than a clock's grace after the earliest start of
    the window's people's finished or open work, since only such a one can
    have been at work beside it. That start is taken after the branch floor,
    which the listing's ends and branches give without a call, so a branch
    used again costs no call for its pull request before. Without commits,
    a unit starts at its opening."""
    listed = [unit_of(pull, [], False) for pull in earlier]
    floor_starts(units + [unit for unit in listed if unit is not None])
    first = min((u.start for u in units if u.kind == PERSON and (u.end is None or window.holds(u.end))),
                default=None)
    found = []
    for pull, unit in zip(earlier, listed):
        if unit is None:
            continue
        if first is not None and unit.kind == PERSON and unit.end - first > CLOCK_GRACE:
            unit = unit_of(pull, api.get_all(f"/repos/{api.repo}/pulls/{pull['number']}/commits"), False)
        found.append(unit)
    return found


def work_units(api, facts, earlier, window):
    """collect's units: the window's pull requests, from what was read for
    them already, then the look-back's, then the branch floor. A pull
    request listed twice is one unit, from its first listing, while the
    keys before the breakdown take the listing as it came, as they always
    have."""
    seen = set()
    units = [unit for f in once_each(facts, seen, lambda f: f.pull["number"])
             if (unit := window_unit(f.pull, f.commits, window, others_reviews(f), f.mergeable,
                                     f.mergeable_state)) is not None]
    units += earlier_units(api, once_each(earlier, seen), units, window)
    floor_starts(units)
    return units


def judge_units(gh, recent, earlier, window):
    """The judge's units, by collect's rules: it reads nothing else of the
    pull requests, so each well-formed listed one's commits cost a call."""
    units = [unit for pull in recent if well_formed(pull, True)
             and (unit := window_unit(pull, gh.get_all(f"/repos/{gh.repo}/pulls/{pull['number']}/commits"),
                                      window)) is not None]
    units += earlier_units(gh, earlier, units, window)
    floor_starts(units)
    return units


def finished(units, window):
    """The units that ended, merged or closed, in the window, of any kind:
    each pull request counts once, whole, in the window where it ended."""
    return [unit for unit in units if window.holds(unit.end)]


def files_needed(units, window):
    """The pull requests whose changed files the breakdown needs, in number
    order: every finished person's one with a partner, and its partners."""
    partners = partners_of(units, window.end)
    needed = set()
    for unit in finished(units, window):
        if unit.kind == PERSON and partners[unit.number]:
            needed.update([unit.number, *(other.number for other in partners[unit.number])])
    return sorted(needed)


def changed_files(api, number):
    """The paths a pull request changed, renamed ones by both names, held
    in memory to see whether two pull requests' files meet."""
    names = set()
    for item in api.get_all(f"/repos/{api.repo}/pulls/{number}/files"):
        names.update(name for name in (item.get("filename"), item.get("previous_filename")) if name)
    return frozenset(names)


def read_files(api, numbers):
    """The changed files of each pull request named, read last and in this
    order, and how many were left unread: the files are the one read that
    stops at the board's reserve without failing the run, since only the
    counts that need them go without."""
    files = {}
    for position, number in enumerate(numbers):
        try:
            files[number] = changed_files(api, number)
        except OutOfBudget:
            return files, len(numbers) - position
    return files, 0


# --- measuring ---------------------------------------------------------------

# The median's one rule lives in economics, which takes the coordination
# block's medians by it too.
median_seconds = economics.median_seconds


def p90_seconds(values):
    """The 90th percentile by nearest rank: the smallest value at or above
    nine tenths of them; None for no values."""
    if not values:
        return None
    ordered = sorted(values)
    return ordered[math.ceil(0.9 * len(ordered)) - 1]


def ratio(part, whole):
    return round(part / whole, 4) if whole else None


def is_bot(user):
    login = (user or {}).get("login") or ""
    return (user or {}).get("type") == "Bot" or login.endswith("[bot]")


def person(user):
    """Who a GitHub user is, for counting alone: None for an app or a bot,
    which is no contributor."""
    if not isinstance(user, dict) or not user.get("login") or is_bot(user):
        return None
    return "login:" + user["login"].lower()


def commit_author(commit):
    """Who made a commit: its GitHub account when it has one, else the
    author's email, so an agent that commits under a name of its own counts
    once whichever account ran it."""
    account = commit.get("author")
    if isinstance(account, dict) and account.get("login"):
        return person(account)
    email = ((((commit.get("commit") or {}).get("author")) or {}).get("email") or "").lower()
    return "email:" + email if email and "[bot]" not in email else None


def made_in(commit, window):
    """Whether the commit was authored or committed in the window; a rebase
    keeps the author date and renews the committer date."""
    return any(window.holds(commit_time(commit, who)) for who in ("author", "committer"))


def contributors(facts, window):
    """How many people did something in the window on its pull requests:
    opened one, submitted a review on one, or made a commit on one. A pull
    request enters the window on any change, a board comment on a stale
    one among them, so each person is dated by their own act, not by the
    pull request's; the identities are counted here and never leave."""
    found = set()
    for f in facts:
        if window.holds(time_of(f.pull, "created_at")):
            found.add(person(f.pull.get("user")))
        found.update(person(r.get("user")) for r in f.reviews if window.holds(time_of(r, "submitted_at")))
        found.update(commit_author(c) for c in f.commits if made_in(c, window))
    return len(found - {None})


def issue_number(comment):
    tail = (comment.get("issue_url") or "").rstrip("/").rsplit("/", 1)[-1]
    return int(tail) if tail.isdecimal() else None


def conflicts(facts, comments):
    """The pull requests seen in conflict: open ones GitHub calls dirty
    now, and any the board's "main moved" comment called conflicting in
    the window. A pull request only behind is no conflict."""
    dirty = {f.pull["number"] for f in facts if f.mergeable_state == DIRTY}
    marked = {issue_number(c) for c in comments
              if MAIN_MOVED_PREFIX in (c.get("body") or "") and CONFLICT_WORDS in (c.get("body") or "")}
    return len((dirty | marked) - {None})


def released_claims(comments):
    """How many claims were released: the comments that open with the
    release marker (claims.starts_with_release_marker), where one on an
    issue that had another less than a day before is the same release said
    again and is not counted twice. An hourly coordinator said it again
    before it learnt to skip a claim it had released; a later claim on the
    issue cannot be released within a day of the last release, since it
    must go a day without a commit first. Whoever wrote the comment counts:
    the board asks GitHub whether an author may write, which a count need
    not, and the workflow's token may show the owner's account, a private
    member of the organization, as a mere contributor."""
    times = {}
    for c in comments:
        number, when = issue_number(c), time_of(c, "created_at")
        if claims.starts_with_release_marker(c) and number is not None and when is not None:
            times.setdefault(number, []).append(when)
    count = 0
    for found in times.values():
        found.sort()
        count += sum(1 for before, after in zip([None] + found, found)
                     if before is None or after - before > claims.IDLE_AFTER)
    return count


def window_runs(facts, window):
    """The check runs on the window's commits that finished in the window,
    each once: a commit in two pull requests brings its runs twice."""
    runs = {}
    for f in facts:
        for run in (r for found in f.runs.values() for r in found):
            if window.holds(time_of(run, "completed_at")):
                runs[run.get("id") or id(run)] = run
    return list(runs.values())


def claim_time(f):
    """When the claim was made: the earliest author date among the pull
    request's commits, which on a claude/<n> branch is the claim file's own
    commit; a rebase keeps author dates and changes committer dates."""
    dates = [d for c in f.commits if (d := commit_time(c, "author"))]
    return min(dates) if dates else None


def claim_waits(opened):
    """Seconds from each claim to its first pull request opened in the
    window, one per claimed issue."""
    first = {}
    for f in sorted(opened, key=lambda f: time_of(f.pull, "created_at")):
        number = claim_number(f.pull)
        if number is not None and number not in first:
            first[number] = f
    return [s for f in first.values()
            if (s := seconds_between(claim_time(f), time_of(f.pull, "created_at"))) is not None]


def others_reviews(f):
    """The reviews by anyone but the pull request's author: an author's
    answer in a review thread is filed by GitHub as a review of their own."""
    author = (f.pull.get("user") or {}).get("login")
    return [r for r in f.reviews if not author or (r.get("user") or {}).get("login") != author]


def review_waits(opened):
    """Seconds from opening to the first review by someone other than the
    author, for each pull request opened in the window that has one."""
    waits = []
    for f in opened:
        times = [t for r in others_reviews(f) if (t := time_of(r, "submitted_at"))]
        wait = seconds_between(time_of(f.pull, "created_at"), min(times)) if times else None
        if wait is not None:
            waits.append(wait)
    return waits


def merge_waits(merged):
    """Seconds from opening to merging, for each pull request merged in the
    window."""
    return [s for f in merged
            if (s := seconds_between(time_of(f.pull, "created_at"), time_of(f.pull, "merged_at"))) is not None]


def first_check_passed(f):
    """Whether every counted check run on the head the pull request opened
    with passed, every attempt among them, so a check that passed only when
    run again did not pass first; None when that head is unknown or none of
    its runs counted."""
    counted = [r for r in f.first_runs or [] if r.get("conclusion") in PASSED | FAILED]
    return all(r["conclusion"] in PASSED for r in counted) if counted else None


def summary(values):
    return median_seconds(values), p90_seconds(values)


def measures(facts, extra_claims, comments, window):
    """The counts, the durations and the rates over what was read."""
    opened = [f for f in facts if window.holds(time_of(f.pull, "created_at"))]
    merged = [f for f in facts if window.holds(time_of(f.pull, "merged_at"))]
    closed = [f for f in facts if not f.pull.get("merged_at") and window.holds(time_of(f.pull, "closed_at"))]
    reviews = [r for f in facts for r in others_reviews(f) if window.holds(time_of(r, "submitted_at"))]
    runs = window_runs(facts, window)
    claimed = {n for f in facts if (n := claim_number(f.pull)) is not None} | set(extra_claims)
    counts = dict(zip(COUNT_KEYS, (
        len(claimed), len(opened), len(merged), len(closed), len(reviews),
        sum(1 for r in reviews if r.get("state") == "CHANGES_REQUESTED"),
        sum(1 for r in runs if r.get("conclusion") in PASSED),
        sum(1 for r in runs if r.get("conclusion") in FAILED),
        conflicts(facts, comments),
        released_claims(comments),
        contributors(facts, window),
    )))
    durations = dict(zip(DURATION_KEYS, (*summary(claim_waits(opened)), *summary(merge_waits(merged)),
                                         *summary(review_waits(opened)))))
    verdicts = [v for f in opened if (v := first_check_passed(f)) is not None]
    rates = {"merge_rate": ratio(len(merged), len(merged) + len(closed)),
             "first_check_pass_rate": ratio(sum(verdicts), len(verdicts))}
    return {"counts": counts, "durations_seconds": durations, "rates": rates}


# --- the breakdown -------------------------------------------------------------

def conflict_marks(comments):
    """The board's conflicting marks, oldest first, as (pull request number,
    when, the short sha the default branch moved to). Only a mark an app or
    a bot posted counts: a person's paste of one marks nothing."""
    marks = []
    for c in comments:
        body = c.get("body") or ""
        found = MARK.search(body)
        number, when = issue_number(c), time_of(c, "created_at")
        if found and CONFLICT_WORDS in body and is_bot(c.get("user")) and number is not None and when is not None:
            marks.append((number, when, found.group(1)))
    return sorted(marks, key=lambda mark: mark[1])


def first_marks(marks):
    """The earliest mark of each pull request, by number, as (when, sha)."""
    first = {}
    for number, when, sha in marks:
        first.setdefault(number, (when, sha))
    return first


def names_landing(sha, merge_sha):
    """Whether a mark's short sha names this merge commit."""
    return sha != "unknown" and merge_sha != "" and merge_sha.startswith(sha)


def at_risk(unit, when, first):
    """Whether a landing at this moment could still put the unit in
    conflict: the board had not marked it conflicting yet."""
    mark = first.get(unit.number)
    return mark is None or when < mark[0]


def lands(unit):
    """Whether the unit's merge moved the default branch, the one the board
    marks the others against: a merge into another branch, a release branch
    say, moves nothing it marks."""
    return unit.merged is not None and unit.into_default


def release_times(comments):
    """(issue number, when) of each coordinator release comment, by any
    author, as counts.released_claims takes them: those that open with the
    marker (claims.starts_with_release_marker)."""
    return [(issue_number(c), when) for c in comments
            if claims.starts_with_release_marker(c) and (when := time_of(c, "created_at")) is not None]


def empty_settings():
    return {name: dict.fromkeys(SETTING_KEYS, 0) for name in SETTINGS}


def add_finished(row, unit, partners, first, releases, landings, files):
    """Count one finished unit into its setting's row. The landings it
    counts are other branches': a merge of its own branch, its earlier pull
    request's where the floor starts it, carries its own earlier work and
    puts it at no risk."""
    merges = [other.merged for other in landings if other.branch != unit.branch]
    mine = files.get(unit.number, frozenset())
    row["pull_requests"] += 1
    row["merged"] += unit.merged is not None
    row["closed_unmerged"] += unit.merged is None
    row["conflicted"] += unit.number in first
    row["released"] += unit.claim is not None and any(
        number == unit.claim and unit.start <= when <= unit.end for number, when in releases)
    row["changes_requested"] += any(r.get("state") == "CHANGES_REQUESTED" for r in unit.reviews)
    row["claim_branches"] += unit.claim is not None
    row["partners"] += len(partners)
    row["work_seconds"] += seconds_between(unit.start, unit.end) or 0
    row["landings_under"] += sum(1 for when in merges if unit.created <= when < unit.end and at_risk(unit, when, first))
    row["landings_before_opening"] += sum(1 for when in merges if unit.start <= when < unit.created)
    row["file_overlap"] += any(mine & files.get(other.number, frozenset()) for other in partners)


def add_pair(row, a, b, first, files):
    """Count one pair of people's pull requests at work together into its
    relation's row. Its landing is the merge of the first to end into the
    default branch, from another branch, while the other was open and
    before the board marked the other conflicting; its conflict onset is
    the other's first conflicting mark naming that merge."""
    early, late = sorted((a, b), key=lambda unit: (unit.end, unit.number))
    since, until = max(a.start, b.start), min(a.end, b.end)
    shared = bool(files.get(a.number, frozenset()) & files.get(b.number, frozenset()))
    landing = (lands(early) and early.branch != late.branch and late.created <= early.merged < late.end
               and at_risk(late, early.merged, first))
    onset = landing and late.number in first and names_landing(first[late.number][1], early.merge_sha)
    row["pairs"] += 1
    row["both_active"] += all(any(since <= d <= until for d in unit.authored) for unit in (a, b))
    row["file_overlap"] += shared
    row["landings"] += landing
    row["landings_file_overlap"] += landing and shared
    row["conflict_onsets"] += onset
    row["conflict_onsets_file_overlap"] += onset and shared


def session_groups(units):
    """How many sessions some units hold, units that share a session joined
    however many steps apart; a unit naming no session is left out, so this
    is a lower bound."""
    groups = []
    for unit in units:
        if not unit.sessions:
            continue
        joined, apart = set(unit.sessions), []
        for group in groups:
            if group & joined:
                joined |= group
            else:
                apart.append(group)
        groups = apart + [joined]
    return len(groups)


def peak_counts(person_units, window):
    """The most people's pull requests, people, and sessions of one person
    at work at one moment, over the window's start and each start of
    people's work the window holds."""
    moments = sorted({window.start} | {unit.start for unit in person_units if window.holds(unit.start)})
    peaks = dict.fromkeys(("peak_parallel", "peak_parallel_people", "peak_sessions_one_person"), 0)
    for moment in moments:
        live = [unit for unit in person_units if unit.start <= moment < stop(unit, window.end)]
        people = {unit.person for unit in live}
        sessions = max((session_groups([u for u in live if u.person == who]) for who in people), default=0)
        peaks["peak_parallel"] = max(peaks["peak_parallel"], len(live))
        peaks["peak_parallel_people"] = max(peaks["peak_parallel_people"], len(people))
        peaks["peak_sessions_one_person"] = max(peaks["peak_sessions_one_person"], sessions)
    return peaks


def parallel_counts(units, ended, open_units, marks, landings, window, files_unread):
    """The repository's counts of who and what was at work, over the
    people's pull requests that ended in the window or are open now and
    were opened or given a commit in it."""
    population = [unit for unit in ended + open_units if unit.kind == PERSON]
    sessions_by_person, people_by_session = {}, {}
    for unit in population:
        sessions_by_person.setdefault(unit.person, set()).update(unit.sessions)
        for name in unit.sessions:
            people_by_session.setdefault(name, set()).add(unit.person)
    in_window = [mark for mark in marks if window.holds(mark[1])]
    unmergeable = [unit for unit in open_units if unit.mergeable is False]
    return {
        "people": len(sessions_by_person),
        "sessions": len(people_by_session),
        "sessions_per_person_max": max(map(len, sessions_by_person.values()), default=0),
        "pull_requests_with_session": sum(1 for unit in population if unit.sessions),
        "pull_requests_without_session": sum(1 for unit in population if not unit.sessions),
        "pull_requests_multi_session": sum(1 for unit in population if len(unit.sessions) >= 2),
        "sessions_across_people": sum(1 for people in people_by_session.values() if len(people) >= 2),
        **peak_counts([unit for unit in units if unit.kind == PERSON], window),
        "conflict_marks": len(in_window),
        "conflict_marks_unattributed": sum(1 for _, _, sha in in_window
                                           if not any(names_landing(sha, q.merge_sha) for q in landings)),
        "open_conflicting": len(unmergeable),
        "open_conflicting_drafts": sum(1 for unit in unmergeable if unit.draft),
        "open_conflicting_dirty": sum(1 for unit in unmergeable if unit.mergeable_state == DIRTY),
        "truncated": sum(1 for unit in ended if unit.start < look_back_start(window)),
        "files_unread": files_unread,
    }


def null_file_counts(by_setting, pairs):
    """With some changed files left unread, the counts that need them are
    not known: null, never a zero that would read as no overlap. The apps'
    and bots' settings keep their 0, since they are nobody's partner."""
    for name in PERSON_SETTINGS:
        by_setting[name]["file_overlap"] = None
    for row in pairs.values():
        for key in FILE_PAIR_KEYS:
            row[key] = None


def breakdown(units, comments, files, files_unread, window):
    """parallel, by_setting and pairs: the same pull requests counted again
    by work setting and in pairs, so one person's several sessions are
    never pooled with several people's work."""
    partners = partners_of(units, window.end)
    marks = conflict_marks(comments)
    first = first_marks(marks)
    releases = release_times(comments)
    landings = [unit for unit in units if lands(unit)]
    ended = finished(units, window)
    open_units = [unit for unit in units if unit.end is None]
    by_setting = empty_settings()
    for unit in ended:
        add_finished(by_setting[setting_of(unit, partners[unit.number])], unit, partners[unit.number], first,
                     releases, landings, files)
    for unit in open_units:
        by_setting[setting_of(unit, partners[unit.number])]["open"] += 1
    pairs = {name: dict.fromkeys(PAIR_KEYS, 0) for name in RELATIONS}
    for a, b in combinations(sorted(units, key=lambda unit: unit.number), 2):
        if a.end is not None and b.end is not None and co_active(a, b, window.end) \
                and window.holds(max(a.end, b.end)):
            add_pair(pairs[relation(a, b)], a, b, first, files)
    if files_unread:
        null_file_counts(by_setting, pairs)
    return {"parallel": parallel_counts(units, ended, open_units, marks, landings, window, files_unread),
            "by_setting": by_setting, "pairs": pairs}


def began_before(units, window):
    """Whether a unit the breakdown counts, ended in the window or open,
    began its work before the window: only then can a comment made in the
    look-back, a mark on it or a release of its claim, bear on a count."""
    return any(unit.start < window.start for unit in units if unit.end is None or window.holds(unit.end))


def read_breakdown(api, facts, rest, comments, window):
    """The breakdown, how many changed-file lists it left unread, and the
    units it sorted, from reads made once every other key is measured, so
    it never costs one: the
    look-back's pull requests, read on from where the window's listing
    stopped; the commits of those that can have been at work beside the
    window's; the look-back's comments, when the window's work began before
    it; and last the changed files. Only the changed files give way to the
    board's reserve here; OutOfBudget, or any other GitHubError, from
    another read is the caller's, which then leaves the breakdown out."""
    earlier, _ = updated_since(rest, look_back_start(window))
    units = work_units(api, facts, earlier, window)
    if began_before(units, window):
        comments = look_back_comments(api, window) + comments
    files, unread = read_files(api, files_needed(units, window))
    return breakdown(units, comments, files, unread, window), unread, units


# --- economics -----------------------------------------------------------------

def usage_artifacts(api, since):
    """The usage records' artifacts made since since, with when each was
    made, newest first as GitHub lists them; expired ones are left out. The
    listing stops at the first page made wholly before since."""
    found = []
    for page in api.pages(f"/repos/{api.repo}/actions/artifacts"):
        listed = (page.get("artifacts") or []) if isinstance(page, dict) else []
        made = [(artifact, time_of(artifact, "created_at")) for artifact in listed]
        found += [(artifact, when) for artifact, when in made
                  if when is not None and when >= since and not artifact.get("expired")
                  and isinstance(artifact.get("id"), int) and str(artifact.get("name") or "").startswith(USAGE_PREFIX)]
        if not listed or all(when is not None and when < since for _, when in made):
            break
    return found


def artifact_record(api, artifact):
    """The usage record an artifact holds, in the form the cache keeps
    (economics.kept_record), or None when it holds none. Any workflow run
    of the repository can upload an artifact under this name, so its bytes
    are anyone's, and whatever zipfile, its decompressors and json raise on
    them is no record of ours: an entry marked encrypted or packed by a
    method zipfile does not know, a broken stream, JSON nested deeper than
    the interpreter's stack, among what was seen raised."""
    if (artifact.get("size_in_bytes") or 0) > USAGE_MAX_BYTES:
        return None
    target = api.location(f"/repos/{api.repo}/actions/artifacts/{artifact['id']}/zip")
    data = api.download(target, USAGE_MAX_BYTES) if target else None
    if data is None:
        return None
    try:
        with zipfile.ZipFile(io.BytesIO(data)) as archive:
            files = [info for info in archive.infolist() if info.filename.endswith(".json")]
            if len(files) != 1 or files[0].file_size > USAGE_MAX_BYTES:
                return None
            record = json.loads(archive.read(files[0]).decode("utf-8"))
    except Exception:
        return None
    return economics.kept_record(record)


def read_usage_cache(path):
    """The records an earlier run kept, by artifact id, each with when its
    artifact was made and in the form the cache keeps, which a cache an
    older run wrote, holding records as they were uploaded, is brought to;
    a record an older version of the cache kept is marked old and read
    with no profile, which that version dropped; an unreadable file is an
    empty one."""
    try:
        data = json.loads(Path(path).read_text(encoding="utf-8")) if path else None
    except (OSError, ValueError):
        return {}
    version = data.get("v") if isinstance(data, dict) else None
    entries = data.get("records") if version in USAGE_CACHE_VERSIONS else None
    kept = {}
    for key, entry in (entries if isinstance(entries, dict) else {}).items():
        record = economics.kept_record(entry.get("record")) if isinstance(entry, dict) else None
        if record is None or not isinstance(entry.get("at"), str) or parse_time(entry["at"]) is None:
            continue
        if version != USAGE_CACHE_VERSION or entry.get("old") is True:
            kept[key] = {"at": entry["at"], "record": {**record, "profile": None}, "old": True}
        else:
            kept[key] = {"at": entry["at"], "record": record}
    return kept


def write_usage_cache(path, kept):
    if path:
        write_json(path, {"v": USAGE_CACHE_VERSION, "records": kept})


def read_usage(api, since, cached):
    """The usage records of the Actions runs made since since, by artifact
    id, the cached ones taken as they are and new ones downloaded, then the
    old ones (read_usage_cache) downloaded again, at most USAGE_LIMIT
    downloads in all, newest first; and how many new ones were left unread,
    None when the artifacts could not be listed at all (a workflow without
    actions: read, or the board's reserve reached), when only the cached
    ones are kept. An old record whose download fails, or is left for a
    later run, is kept as it was."""
    kept = {key: entry for key, entry in cached.items() if parse_time(entry["at"]) >= since}
    try:
        listed = usage_artifacts(api, since)
    except GitHubError:
        return kept, None
    fresh = [(artifact, when) for artifact, when in listed if str(artifact["id"]) not in kept]
    old = [(artifact, when) for artifact, when in listed if kept.get(str(artifact["id"]), {}).get("old")]
    unread = fetched = 0
    for artifact, when in fresh:
        if fetched >= USAGE_LIMIT:
            unread += 1
            continue
        try:
            record = artifact_record(api, artifact)
        except OutOfBudget:
            fetched, unread = USAGE_LIMIT, unread + 1
            continue
        except GitHubError:
            unread += 1
            continue
        fetched += 1
        if record is not None:
            kept[str(artifact["id"])] = {"at": iso(when), "record": record}
    for artifact, when in old:
        if fetched >= USAGE_LIMIT:
            break
        try:
            record = artifact_record(api, artifact)
        except OutOfBudget:
            break
        except GitHubError:
            continue
        fetched += 1
        if record is not None:
            kept[str(artifact["id"])] = {"at": iso(when), "record": record}
    return kept, unread


def work_of(unit, facts):
    """A finished unit as the economics read it, with its commits and failed
    check runs from what collect read for it."""
    commits = facts.commits if facts else []
    failures = sum(1 for runs in (facts.runs.values() if facts else ()) for run in runs
                   if run.get("conclusion") in FAILED)
    return economics.Work(unit.number, claim=unit.claim, kind=unit.kind, start=unit.start, end=unit.end,
                          merged=unit.merged, opened=unit.created, commits=commits, check_failures=failures,
                          first_check=first_check_passed(facts) if facts else None)


def works_and_others(units, facts, window):
    """The units of people and of the repository's Claude workflows: those
    that ended in the window as works, and the rest, open or ended before
    it, as what can have been at work beside them."""
    by_number = {f.pull["number"]: f for f in facts}
    works, others = [], []
    for unit in units:
        if unit.kind not in economics.WORKING_KINDS:
            continue
        if window.holds(unit.end):
            works.append(work_of(unit, by_number.get(unit.number)))
        else:
            others.append(economics.Work(unit.number, kind=unit.kind, start=unit.start, end=unit.end))
    return works, others


def window_commits(facts, window):
    """The commits made in the window, of every pull request read, each
    once: they carry the trailers of what sessions spent in it."""
    seen, found = set(), []
    for f in facts:
        for commit in f.commits:
            sha = commit.get("sha")
            if sha not in seen and window.holds(commit_time(commit, "committer")):
                seen.add(sha)
                found.append(commit)
    return found


def reverts(facts, window):
    return sum(1 for f in facts if window.holds(time_of(f.pull, "merged_at"))
               and str(f.pull.get("title") or "").startswith(REVERT_TITLE))


class Records(NamedTuple):
    """The Actions runs' usage records since the look-back began, oldest
    first, those of them made in the window, how many were left unread
    (read_usage), and how many of those read are kept from an older cache
    without their profile, not yet downloaded again (read_usage_cache)."""

    all: list
    in_window: list
    unread: object
    unprofiled: int = 0


def usage_records(api, window, cache_path=None):
    """The usage records of the Actions runs since the look-back began, as
    Records; read after the breakdown, so they never cost a key before it."""
    kept, unread = read_usage(api, look_back_start(window), read_usage_cache(cache_path))
    write_usage_cache(cache_path, kept)
    entries = sorted(kept.values(), key=lambda entry: entry["at"])
    records = [economics.actions_record(entry["record"]) for entry in entries]
    in_window = [record for entry, record in zip(entries, records) if window.holds(parse_time(entry["at"]))]
    return Records(records, in_window, unread, sum(1 for entry in entries if entry.get("old")))


def economic_keys(facts, units, window, records):
    """economics and coordination, from the units the breakdown sorted, the
    commits read for them and the usage records (Records)."""
    works, others = works_and_others(units, facts, window)
    return {"economics": economics.economics(works, window_commits(facts, window), records.in_window, records.all,
                                             reverts(facts, window), records.unread),
            "coordination": economics.coordination(works, records.all, window.end, others, records.in_window)}


def profile_keys(facts, units, window, records):
    """PROFILE.json's measures: each pull request's profile, the
    Symphony-Profile lines on its own commits and the profiles of the
    Actions runs that name it (profile.pull_profile), summed per work
    setting over the pull requests that ended in the window, each counted
    once in its setting as the breakdown counts it, and in the coordination
    block's two sides (profile.coordination); how many sessions' lines and
    runs' profiles the sums hold, the coordinator's runs that count on the
    parallel side among them; and how many runs' profiles are not read
    yet, those left unread and those an older cache kept without one, None
    where the runs could not be listed."""
    by_number = {f.pull["number"]: f for f in facts}
    partners = partners_of(units, window.end)
    ended = finished(units, window)
    works = [work_of(unit, by_number.get(unit.number)) for unit in ended]
    profiles = {work.number: profile.pull_profile(work, records.all) for work in works}
    settings = {unit.number: setting_of(unit, partners[unit.number]) for unit in ended}
    working, others = works_and_others(units, facts, window)
    overhead = {id(run) for run in economics.overhead(working, records.in_window)}
    summed = [record for record in records.all if record["profile"] and (
        id(record) in overhead or any(economics.names(record, w) for w in works))]
    lines = [profile.commit_lines(commit) for work in works for commit in work.commits]
    sessions = {session for found in lines for session, _ in found}
    unread = None if records.unread is None else records.unread + records.unprofiled
    return {"sessions": len(sessions - {""}), "actions_runs": len(summed), "actions_runs_unread": unread,
            "by_setting": {name: profile.group([w for w in works if settings[w.number] == name], profiles)
                           for name in SETTINGS},
            "coordination": profile.coordination(working, profiles, window.end, others, records.in_window)}


def metrics_record(measured, window, version):
    return {"schema": SCHEMA, "generated_at": iso(window.end), "window_days": window.days, "sharing": ON,
            "protocol_version": version, **measured}


def profile_record(measured, window):
    return {"schema": PROFILE_SCHEMA, "generated_at": iso(window.end), "window_days": window.days, "sharing": ON,
            **measured}


def write_json(path, record):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(record, indent=2) + "\n", encoding="ascii")


def same_apart_from_time(old_text, new_text):
    """Whether two METRICS.json texts, or two PROFILE.json texts, hold the
    same record once generated_at is set aside."""
    try:
        old, new = json.loads(old_text), json.loads(new_text)
    except ValueError:
        return False
    if not isinstance(old, dict) or not isinstance(new, dict):
        return False
    old.pop("generated_at", None)
    new.pop("generated_at", None)
    return old == new


# --- the judge: App F, on the maintainer's own Claude Code --------------------

# App F's codes, their names and their three groups, as
# OpenSpineConsortium/coordbench's coordbench/judge.py has them (its
# APPF_7, ERROR_NAMES and ERROR_GROUPS), which follow CooperBench's paper.
APPF_CODES = ("C1a", "C1b", "C2", "C3a", "C3b", "C4a", "C4b")
APPF_NAMES = {"C1a": "Unanswered direct question (no reply)", "C1b": "Unanswered direct question (ignored)",
              "C2": "Non-answer or vague answer", "C3a": "Incorrect claim (uncorrected)",
              "C3b": "Incorrect claim (corrected)", "C4a": "Spammy repetition (repeats same information)",
              "C4b": "Spammy repetition (near-duplicate status blocks)"}
GROUPS = ("unresponsiveness", "hallucination", "repetition")
APPF_GROUPS = {"C1a": "unresponsiveness", "C1b": "unresponsiveness", "C2": "unresponsiveness",
               "C3a": "hallucination", "C3b": "hallucination", "C4a": "repetition", "C4b": "repetition"}
PROMPT_VARIANT = "appF-7"
DEFAULT_MODEL = "claude-sonnet-5-5"
# JUDGE.json's keys as the judge wrote them before the breakdown, which push
# still takes, and as it writes them now: the window and the limit the run
# covered, the conversations with two or more people, and the counts again
# by work setting, each setting's in JUDGE_SETTING_KEYS.
LEGACY_JUDGE_KEYS = ("schema", "generated_at", "model", "prompt_variant", "conversations", "counts", "groups")
JUDGE_KEYS = LEGACY_JUDGE_KEYS + ("window_days", "limit", "multi_speaker", "by_setting")
JUDGE_SETTING_KEYS = ("pull_requests", "conversations", "multi_speaker", "counts")
# The model is the one free-form string JUDGE.json carries, and the branch
# metrics is as public as the repository: a name is kept to the shape
# OpenSpineConsortium's collector keeps, so a cloud model ARN, with an
# account number in it, is refused before it can be written.
MODEL_NAME = re.compile(r"[A-Za-z0-9][A-Za-z0-9._-]{0,63}(\[[a-z0-9]{1,8}\])?")
# Each of these makes claude bill an API or a cloud account instead of the
# maintainer's subscription, and send the conversations there; the judge
# refuses while one is set, saying where they would have gone. It never
# passes one on, nor ANTHROPIC_BASE_URL, which alone would send them to a
# gateway: that one is dropped, and the dry run says so.
BILLING_VARIABLES = {
    "ANTHROPIC_API_KEY": "bill the API with it",
    "ANTHROPIC_AUTH_TOKEN": "bill the API, or the gateway ANTHROPIC_BASE_URL names, with it",
    "CLAUDE_CODE_USE_BEDROCK": "send every conversation to Amazon Bedrock, billed to that account",
    "CLAUDE_CODE_USE_VERTEX": "send every conversation to Google Vertex AI, billed to that account",
    "CLAUDE_CODE_USE_FOUNDRY": "send every conversation to Microsoft Foundry, billed to that account",
}
DROPPED_VARIABLES = (*BILLING_VARIABLES, "ANTHROPIC_BASE_URL")
TOOL_ROLES = {"gh": "reads the pull requests as you", "claude": "is the model that codes them"}
PUSH_COMMAND = "python3 tools/metrics.py push"
# App F's text word for word, as coordbench's prompts/commerrors-appF-7.txt
# holds it, up to its Output paragraph. That paragraph asked for quotes and
# a description; here the answer is the codes alone, so no text of the
# conversation comes back even into the reply this tool reads and drops.
APPF_PROMPT = "\n".join([
    "You are a careful reviewer of two agent collaboration conversations. This is a precision-first detector of "
    "bad conversation patterns. Prefer returning no issue unless the evidence is strong and explicit.",
    "Important exclusion. Do not label state mismatch or visibility confusion itself as an error (e.g., agents on "
    "separate branches unable to see each other's changes). Bad conversation patterns around these topics should "
    "still be labeled.",
    "Taxonomy. Label at most one category per conversation.",
    *(f"- {code} {APPF_NAMES[code]}" for code in APPF_CODES),
    "Evidence requirements. Include at least two exact quotes that make the issue undeniable. C1a/C1b require the "
    "question plus demonstration of missing or irrelevant response. C3a requires the incorrect claim and later "
    "contradiction. C4a/C4b require two quotes showing the repetition.",
    'Output. The quotes are for your own check and stay out of the answer: return only the JSON object '
    '{"codes": ["<category id>"]} naming the one category the evidence makes undeniable, or {"codes": []} when '
    'it is weak. No other key and no other text.',
])
CODE = re.compile(r"\bC([1-4])([ab])?\b", re.IGNORECASE)
FENCE = re.compile(r"^```(?:json)?\s*|\s*```$", re.IGNORECASE)
MENTION = re.compile(r"@([A-Za-z0-9](?:[A-Za-z0-9-]*[A-Za-z0-9])?)")
# coordbench's bounds: a conversation above 100,000 characters keeps its
# first and last 50,000.
CONVERSATION_LIMIT = 100_000
CONVERSATION_KEEP = 50_000
EMPTY_MCP = '{"mcpServers": {}}'
# The system prompt goes to claude as a file in the folder it runs in, not
# as an argument: a shell in the way (cmd.exe runs npm's claude.cmd on
# Windows) would split a multi-line argument, and a process listing shows
# every argument.
PROMPT_FILE = "system-prompt.txt"
CALL_TIMEOUT = 600
API_VERSION = "2022-11-28"
PER_PAGE = 100
REVIEW_WORDS = {"APPROVED": "review: approved", "CHANGES_REQUESTED": "review: changes requested",
                "COMMENTED": "review: commented", "DISMISSED": "review: dismissed"}
HEADER = ("Conversation of a pull request: {count} messages in time order, each as [message id] time sender -> "
          "recipient (all for a broadcast): content. The people are named A, B, C and so on in the order they "
          "first speak; bot is any app or bot account.")
# The protocol's own comments are posted by rule: the board's "main moved"
# on every move of the default branch and its "idle" on every quiet day,
# again and again, and the coordinator's release. They are nobody's words,
# and App F would read them as spammy repetition, so they stay out.
PROTOCOL_MARKERS = (MAIN_MOVED_PREFIX, board.IDLE_MARKER, RELEASE_MARKER)


class JudgeError(Exception):
    """A judge step that failed outright; the message says which."""


class Gh:
    """GitHub's REST API through the maintainer's own gh, signed in as
    them; the shape of the protocol library's client, so the window's
    pull requests are read by the same code collect uses."""

    def __init__(self, repo=None, run=subprocess.run, binary="gh"):
        self.repo, self.run, self.binary = repo, run, binary

    def call(self, *args):
        try:
            done = self.run([self.binary, *args], capture_output=True, text=True, encoding="utf-8")
        except OSError as problem:
            raise JudgeError(f"gh did not run: {problem}") from None
        if done.returncode != 0:
            raise JudgeError(f"gh {args[0]} {args[-1]}: {(done.stderr or done.stdout or '').strip()[:300]}")
        return done.stdout

    def get(self, path, params=None):
        query = "?" + urllib.parse.urlencode(params) if params else ""
        text = self.call("api", "-H", "Accept: application/vnd.github+json",
                         "-H", f"X-GitHub-Api-Version: {API_VERSION}", path.lstrip("/") + query)
        try:
            return json.loads(text) if text.strip() else None
        except ValueError:
            raise JudgeError(f"gh api {path}: the answer was not JSON") from None

    def pages(self, path, params=None):
        """A list endpoint page by page, until a page comes back short."""
        page = 1
        while True:
            items = self.get(path, {**(params or {}), "per_page": PER_PAGE, "page": page}) or []
            yield items
            if len(items) < PER_PAGE:
                return
            page += 1

    def get_all(self, path, params=None):
        return [item for page in self.pages(path, params) for item in page]

    def default_repo(self):
        return self.call("repo", "view", "--json", "nameWithOwner", "--jq", ".nameWithOwner").strip()


@dataclass
class Message:
    when: datetime
    user: dict
    note: str
    text: str


def is_protocol_comment(text):
    return any(marker in text for marker in PROTOCOL_MARKERS)


def pull_messages(pull, comments, review_comments, reviews):
    """A pull request's conversation in time order: its description, its
    comments, its review comments and its reviews, without the protocol's
    own automated comments. A review with no words still says something
    when it approves or asks for changes."""
    messages = []
    if (pull.get("body") or "").strip():
        messages.append(Message(time_of(pull, "created_at"), pull.get("user"), "description", pull["body"]))
    messages += [Message(time_of(c, "created_at"), c.get("user"), None, c.get("body") or "") for c in comments]
    messages += [Message(time_of(c, "created_at"), c.get("user"), "review comment", c.get("body") or "")
                 for c in review_comments]
    for review in reviews:
        state = review.get("state")
        if state in REVIEW_WORDS and ((review.get("body") or "").strip() or state != "COMMENTED"):
            messages.append(Message(time_of(review, "submitted_at"), review.get("user"), REVIEW_WORDS[state],
                                    review.get("body") or ""))
    epoch = datetime.min.replace(tzinfo=timezone.utc)
    said = (m for m in messages if (m.text.strip() or m.note) and not is_protocol_comment(m.text))
    return sorted(said, key=lambda m: m.when or epoch)


def letters(index):
    """A, B, ... Z, then AA, AB and on, for the index-th speaker."""
    name = ""
    index += 1
    while index:
        index, rest = divmod(index - 1, 26)
        name = chr(ord("A") + rest) + name
    return name


def speaker_names(messages):
    """Each person's letter, by login, in the order they first speak."""
    names = {}
    for m in messages:
        login = ((m.user or {}).get("login") or "ghost").lower()
        if not is_bot(m.user) and login not in names:
            names[login] = letters(len(names))
    return names


def speaker(user, names):
    return "bot" if is_bot(user) else names[((user or {}).get("login") or "ghost").lower()]


def anonymised(text, names):
    """The text with every mention of a speaker made their letter."""
    return MENTION.sub(lambda m: "@" + names.get(m.group(1).lower(), m.group(1)), text)


def bounded(text):
    """coordbench's cut: the text whole, or its two ends with a line saying
    how much is left out between them."""
    if len(text) <= CONVERSATION_LIMIT:
        return text
    left_out = len(text) - 2 * CONVERSATION_KEEP
    note = (f"\n[this conversation is {len(text):,} characters long; the {left_out:,} characters after the "
            f"first {CONVERSATION_KEEP:,} are left out here, and the last {CONVERSATION_KEEP:,} follow.]\n")
    return text[:CONVERSATION_KEEP] + note + text[-CONVERSATION_KEEP:]


def render(messages):
    """The conversation as the local model reads it, or "" when there is
    none, which no call is made for."""
    if not messages:
        return ""
    names = speaker_names(messages)
    lines = [HEADER.format(count=len(messages)), ""]
    for number, m in enumerate(messages, start=1):
        note = f"({m.note}) " if m.note else ""
        when = iso(m.when) if m.when else "-"
        lines.append(f"[{number}] {when} {speaker(m.user, names)} -> all: {note}{anonymised(m.text, names)}")
    return bounded("\n".join(lines) + "\n")


def conversation(gh, pull):
    """One pull request's conversation, read with gh and rendered, and how
    many people speak in it: App F's unanswered-question codes need two."""
    path = f"/repos/{gh.repo}"
    number = pull["number"]
    messages = pull_messages(pull, gh.get_all(f"{path}/issues/{number}/comments"),
                             gh.get_all(f"{path}/pulls/{number}/comments"),
                             gh.get_all(f"{path}/pulls/{number}/reviews"))
    return render(messages), len(speaker_names(messages))


def claude_command(binary, model, prompt_file):
    """claude -p with App F as the whole system prompt, read from a file,
    one turn, no tool, no MCP server, no settings (so no hook), no skill
    and no session kept; it runs in the folder that holds the prompt file
    and nothing else, so no project's CLAUDE.md is read either."""
    return [binary, "-p", "--output-format", "json", "--model", model, "--max-turns", "1",
            "--system-prompt-file", str(prompt_file), "--tools", "", "--strict-mcp-config", "--mcp-config",
            EMPTY_MCP, "--setting-sources", "", "--disable-slash-commands", "--no-session-persistence"]


class LocalClaude:
    """The maintainer's own Claude Code, signed in to their own account.
    The conversation goes to this one process on stdin and nowhere else.
    The binary is the path the PATH check found, so the process starts
    on a system that resolves a command name differently from a path."""

    def __init__(self, model, run=subprocess.run, binary="claude"):
        self.model, self.run, self.binary = model, run, binary

    def command(self, prompt_file):
        return claude_command(self.binary, self.model, prompt_file)

    def ask(self, conversation):
        with tempfile.TemporaryDirectory(prefix="sdlc-judge-") as folder:
            prompt_file = Path(folder) / PROMPT_FILE
            prompt_file.write_text(APPF_PROMPT, encoding="utf-8")
            try:
                done = self.run(self.command(prompt_file), input=conversation, capture_output=True, text=True,
                                encoding="utf-8", cwd=folder, timeout=CALL_TIMEOUT,
                                env=judge_environment(os.environ))
            except (OSError, subprocess.SubprocessError) as problem:
                raise JudgeError(f"claude did not run: {problem}") from None
        return claude_result(done)


def judge_environment(environ):
    """The maintainer's environment, which holds their sign-in, without
    the variables that would bill or redirect the call; the refusal has
    already kept the billing ones out, and the gateway URL goes here."""
    return {name: value for name, value in environ.items() if name not in DROPPED_VARIABLES}


def claude_result(done):
    """The reply's text from claude's JSON result; when claude failed, the
    error says what claude said, which reaches the maintainer's terminal
    and nowhere else."""
    try:
        result = json.loads(done.stdout)
    except (TypeError, ValueError):
        raise JudgeError(f"claude exited {done.returncode} without a JSON result: "
                         f"{(done.stderr or '').strip()[-300:]}") from None
    if done.returncode != 0 or not isinstance(result, dict) or result.get("is_error") \
            or not isinstance(result.get("result"), str):
        raise JudgeError(f"claude failed ({claude_failure(result, done.stderr)}); nothing is written")
    return result["result"]


def claude_failure(result, stderr):
    """Claude Code's own account of a failure: the subtype, then the end of
    its result text, where a sign-in or usage-limit message goes under the
    subtype success, its errors and stderr, whichever it filled."""
    if not isinstance(result, dict):
        return "no result object"
    errors = result.get("errors") if isinstance(result.get("errors"), list) else []
    said = [text.strip()[-300:] for text in (result.get("result"), *map(str, errors), stderr)
            if isinstance(text, str) and text.strip()]
    return "; ".join([str(result.get("subtype")), *said])


def reply_object(text):
    """The reply's JSON object, from the whole text, a fenced block or the
    span between its outer braces; None when there is none."""
    text = FENCE.sub("", (text or "").strip())
    for candidate in (text, text[text.find("{"): text.rfind("}") + 1]):
        try:
            value = json.loads(candidate)
        except ValueError:
            continue
        if isinstance(value, dict):
            return value
    return None


def reply_code(value):
    """The one App F code a reply names, or None. App F labels at most one
    category per conversation, so the first code of the taxonomy counts;
    anything else in the reply is dropped here and never kept."""
    codes = value.get("codes") if isinstance(value, dict) else None
    for item in codes if isinstance(codes, list) else []:
        match = CODE.search(str(item))
        code = match and f"C{match.group(1)}{(match.group(2) or '').lower()}"
        if code in APPF_CODES:
            return code
    return None


def judged_code(claude, conversation):
    """(readable, code) for one conversation: one call, and one more when
    the reply is not JSON, as coordbench's judge allows."""
    for _ in range(2):
        value = reply_object(claude.ask(conversation))
        if value is not None:
            return True, reply_code(value)
    return False, None


def empty_judge_settings():
    return {name: {"pull_requests": 0, "conversations": 0, "multi_speaker": 0, "counts": dict.fromkeys(APPF_CODES, 0)}
            for name in SETTINGS}


def judge_record(model, conversations, counts, now, window_days, limit, multi_speaker, by_setting):
    return {"schema": JUDGE_SCHEMA, "generated_at": iso(now), "model": model, "prompt_variant": PROMPT_VARIANT,
            "conversations": conversations, "counts": {code: counts.get(code, 0) for code in APPF_CODES},
            "groups": {group: sum(counts.get(code, 0) for code in APPF_CODES if APPF_GROUPS[code] == group)
                       for group in GROUPS},
            "window_days": window_days, "limit": limit, "multi_speaker": multi_speaker, "by_setting": by_setting}


def judge_settings(gh, recent, earlier, window):
    """The work setting of each listed pull request that ended in the
    window, by number, sorted as collect sorts them; it never reaches the
    model."""
    units = judge_units(gh, recent, earlier, window)
    partners = partners_of(units, window.end)
    return {unit.number: setting_of(unit, partners[unit.number]) for unit in finished(units, window)}


def judged_by_setting(pulls, settings):
    """JUDGE.json's by_setting with each judged pull request that ended in
    the window counted under its setting, its conversation still to come."""
    by_setting = empty_judge_settings()
    for pull in pulls:
        if pull["number"] in settings:
            by_setting[settings[pull["number"]]]["pull_requests"] += 1
    return by_setting


def add_conversation(by_setting, settings, pull, speakers, code=None):
    """Count one conversation under its pull request's setting, when the
    pull request ended in the window."""
    name = settings.get(pull["number"])
    if name is None:
        return
    row = by_setting[name]
    row["conversations"] += 1
    row["multi_speaker"] += speakers >= 2
    if code:
        row["counts"][code] += 1


def judge_refusal(environ, which):
    """Why the judge will not run here, or None."""
    for name, outcome in BILLING_VARIABLES.items():
        if environ.get(name):
            return (f"judge: {name} is set, and claude would {outcome}, which this project never does. Unset it "
                    f"(`unset {name}`) and run again: claude then runs on your own subscription, signed in with "
                    "`claude` as usual.")
    for tool, why in TOOL_ROLES.items():
        if not which(tool):
            return f"judge: {tool} is not on PATH; it {why}. Install it, sign in, and run again."
    return None


def repository_of(url):
    """OWNER/NAME from a remote's URL in any form git takes (https, ssh,
    scp-like, a path; with or without .git), or None when it ends in no
    such pair."""
    parts = re.split(r"[/:]", re.sub(r"(\.git)?/?$", "", url.strip(), count=1))
    return "/".join(parts[-2:]) if len(parts) >= 2 and all(parts[-2:]) else None


def origin_refusal(git, repo):
    """Why a push would land on the wrong repository, or None: the push
    goes to origin, and JUDGE.json names no repository, by design, so
    origin must be the repository judged."""
    try:
        url = git("remote", "get-url", REMOTE).strip()
    except JudgeError:
        url = ""
    found = repository_of(url)
    if found and found.lower() == repo.lower():
        return None
    where = f"is {url}" if url else "has no URL"
    return (f"the repository judged is {repo}, but this clone's {REMOTE} {where}, so the push would put its counts "
            f"on another repository's branch {METRICS_BRANCH}. Run this in a clone of {repo}.")


def show_dry_run(claude, conversations, pulls, record, environ, out):
    """What a run would send and write, and nothing sent or written."""
    print(f"judge: a dry run; nothing is sent, written or pushed. {len(conversations)} of {pulls} pull requests "
          "have a conversation, and a run would send each to your local claude, one call each.", file=out)
    tally = {name: record["by_setting"][name]["pull_requests"] for name in SETTINGS}
    print(f"judge: {sum(tally.values())} of them ended in the window, by work setting: "
          + ", ".join(f"{name} {count}" for name, count in tally.items()), file=out)
    if environ.get("ANTHROPIC_BASE_URL"):
        print("judge: ANTHROPIC_BASE_URL is set here and is dropped from claude's environment, so the "
              "conversations go where your sign-in does, not to that URL.", file=out)
    print(f"--- the command, run in a temporary folder holding {PROMPT_FILE} and nothing else ---", file=out)
    print(shlex.join(claude.command(f"<that folder>/{PROMPT_FILE}")), file=out)
    print(f"--- the system prompt, written to {PROMPT_FILE} ---", file=out)
    print(APPF_PROMPT, file=out)
    print("--- stdin: the first conversation ---", file=out)
    print(conversations[0] if conversations else "(no conversation in the window)", file=out)
    print("--- the JUDGE.json a run would write (its counts are the model's to fill) ---", file=out)
    print(json.dumps(record, indent=2), file=out)
    return 0


def judge(args, environ, gh, claude, which=shutil.which, git=None, out=sys.stdout):
    """The judge command over injected gh, claude and git, for the tests."""
    refusal = judge_refusal(environ, which)
    if refusal:
        print(refusal, file=sys.stderr)
        return 1
    gh.binary, claude.binary = which("gh"), which("claude")
    gh.repo = args.repo or gh.default_repo()
    # A mismatch is found before a pull request is read or a conversation
    # sent, so it costs nothing.
    refusal = args.push and origin_refusal(git or run_git, gh.repo)
    if refusal:
        print(f"judge: --push refused: {refusal}", file=sys.stderr)
        return 1
    now = parse_time(args.now) if args.now else datetime.now(timezone.utc)
    window = Window.ending(now, args.window_days)
    recent, earlier = listed_pulls(gh, window)
    pulls = recent[:args.limit or None]
    # Analyses leave out a record whose limit is set, as one that judged the
    # newest pull requests alone, so a limit that cut nothing is none.
    limit = args.limit if args.limit is not None and len(recent) > args.limit else None
    found = [(pull, text, speakers) for pull in pulls for text, speakers in [conversation(gh, pull)] if text]
    # Every listed pull request is sorted, whatever --limit says, since its
    # setting turns on what else was at work beside it.
    settings = judge_settings(gh, recent, earlier, window)
    by_setting = judged_by_setting(pulls, settings)
    if args.dry_run:
        for pull, _, speakers in found:
            add_conversation(by_setting, settings, pull, speakers)
        multi = sum(1 for _, _, speakers in found if speakers >= 2)
        record = judge_record(claude.model, len(found), {}, now, args.window_days, limit, multi, by_setting)
        return show_dry_run(claude, [text for _, text, _ in found], len(pulls), record, environ, out)
    counts, judged, unreadable, multi = {}, 0, 0, 0
    for number, (pull, text, speakers) in enumerate(found, start=1):
        readable, code = judged_code(claude, text)
        judged += readable
        unreadable += not readable
        multi += readable and speakers >= 2
        if code:
            counts[code] = counts.get(code, 0) + 1
        if readable:
            add_conversation(by_setting, settings, pull, speakers, code)
        print(f"judge: {number} of {len(found)} coded", file=sys.stderr)
    if unreadable:
        print(f"judge: {unreadable} replies were not JSON after one retry and are left out", file=sys.stderr)
    write_json(args.out, judge_record(claude.model, judged, counts, now, args.window_days, limit, multi,
                                      by_setting))
    print(f"judge: {args.out} written: {judged} conversations, counts only", file=out)
    if args.push:
        print(f"judge: pushed to {METRICS_BRANCH} as {push_judge(Path(args.out), git or run_git)}", file=out)
    return 0


# --- the judge's push: an ordinary commit on the branch metrics ----------------

REMOTE = "origin"
TRACKING = f"refs/remotes/{REMOTE}/{METRICS_BRANCH}"


def run_git(*args, input=None):
    """One git command in the current clone; its output, or JudgeError."""
    done = subprocess.run(["git", *args], input=input, capture_output=True, text=True, encoding="utf-8")
    if done.returncode != 0:
        raise JudgeError(f"git {args[0]}: {(done.stderr or '').strip()[:300]}")
    return done.stdout


def fetched_metrics(git):
    """The branch metrics' commit on origin, or None when it has none yet.
    The tracking ref is dropped first, so a branch deleted since the last
    fetch is not mistaken for one still there."""
    for step in (("update-ref", "-d", TRACKING), ("fetch", "--quiet", "--no-tags", REMOTE,
                                                   f"+refs/heads/{METRICS_BRANCH}:{TRACKING}")):
        try:
            git(*step)
        except JudgeError:
            if step[0] == "fetch":
                return None
    return git("rev-parse", "--verify", "--quiet", TRACKING + "^{commit}").strip() or None


def judge_commit(git, path, parent):
    """A commit whose tree is the parent's with JUDGE.json put in, so the
    METRICS.json and PROFILE.json the workflow wrote stay exactly as they
    are."""
    blob = git("hash-object", "-w", str(path)).strip()
    kept = [line for line in (git("ls-tree", parent).splitlines() if parent else [])
            if not line.endswith("\t" + JUDGE_FILE)]
    tree = git("mktree", input="\n".join(kept + [f"100644 blob {blob}\t{JUDGE_FILE}"]) + "\n").strip()
    # The push is made with a person's credentials, so [skip ci] keeps any
    # workflow a repository starts on every push off a tree of three files.
    message = f"Judge counts as of {iso(datetime.now(timezone.utc))} [skip ci]"
    return git("commit-tree", tree, *(("-p", parent) if parent else ()), "-m", message).strip()


def push_arguments(commit):
    """The push: a plain one, which only a fast-forward passes, so the push
    guard (tools/guard_push.py) lets it through and the workflow's lease
    sees it."""
    return ["push", "--quiet", REMOTE, f"{commit}:refs/heads/{METRICS_BRANCH}"]


def push_judge(path, git):
    """Push JUDGE.json onto the branch metrics; when the workflow rewrote
    the branch in between, read it again and try once more."""
    refused = None
    for _ in range(2):
        commit = judge_commit(git, path, fetched_metrics(git))
        try:
            git(*push_arguments(commit))
            return commit
        except JudgeError as problem:
            refused = problem
    raise JudgeError(f"the push to {METRICS_BRANCH} was refused twice, the last time with: {refused}; "
                     f"{path} is written, so push it alone with `{PUSH_COMMAND} --out {path}`")


def judge_file(path):
    """The JUDGE.json at path, checked to be one the judge wrote: the
    branch metrics takes the judge's aggregates and nothing else."""
    try:
        record = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeDecodeError, ValueError) as problem:
        raise JudgeError(f"{path} is not a JUDGE.json the judge wrote: {problem}") from None
    # A file judged before the breakdown holds the legacy keys alone, and
    # still pushes: it was the judge's then.
    if not isinstance(record, dict) or tuple(record) not in (JUDGE_KEYS, LEGACY_JUDGE_KEYS) \
            or record.get("schema") != JUDGE_SCHEMA or not numbers_alone(record):
        raise JudgeError(f"{path} is not a JUDGE.json the judge wrote; run the judge first")
    return record


def whole(value):
    """Whether a value is a whole number as JSON writes one: true and false
    are not, though Python counts them among its ints."""
    return type(value) is int


def numbers_alone(record):
    """Whether the keys the breakdown added hold the judge's numbers in the
    judge's shape and nothing else: the branch metrics is as public as the
    repository, so a file carrying words under them is no judge's."""
    if tuple(record) == LEGACY_JUDGE_KEYS:
        return True
    rows = record["by_setting"]
    return (whole(record["window_days"]) and whole(record["multi_speaker"])
            and (record["limit"] is None or whole(record["limit"]))
            and isinstance(rows, dict) and tuple(rows) == SETTINGS and all(map(judged_row, rows.values())))


def judged_row(row):
    """Whether one setting's row is the judge's: its keys in order, each a
    whole number but counts, which holds App F's codes in order, each a
    whole number."""
    return (isinstance(row, dict) and tuple(row) == JUDGE_SETTING_KEYS
            and all(whole(row[key]) for key in JUDGE_SETTING_KEYS if key != "counts")
            and isinstance(row["counts"], dict) and tuple(row["counts"]) == APPF_CODES
            and all(map(whole, row["counts"].values())))


def push(args, gh, git=None, out=sys.stdout):
    """The push command: a JUDGE.json the judge already wrote, onto the
    branch metrics, after the same origin check as judge --push, so a
    push that was refused, or a run kept without --push, is never judged
    again to be pushed once. gh is asked only for the default repository."""
    path = Path(args.out)
    judge_file(path)
    gh.repo = args.repo or gh.default_repo()
    refusal = origin_refusal(git or run_git, gh.repo)
    if refusal:
        print(f"push: refused: {refusal}", file=sys.stderr)
        return 1
    print(f"push: {path} pushed to {METRICS_BRANCH} as {push_judge(path, git or run_git)}", file=out)
    return 0


# --- the command line -----------------------------------------------------------

def collect(args, environ, opener=None, sleep=time.sleep):
    if sharing(args.root, environ) == OFF:
        write_json(args.out, OFF_RECORD)
        if args.profile_out:
            write_json(args.profile_out, PROFILE_OFF)
        print(f"metrics: sharing is off ({SWITCH_FILE}, {OLD_SWITCH_FILE} or {SWITCH_VARIABLE}); {args.out} says so"
              " and nothing else")
        return 0
    now = parse_time(args.now) if args.now else datetime.now(timezone.utc)
    window = Window.ending(now, args.window_days)
    api = Budgeted(environ.get("GITHUB_TOKEN"), args.repo, opener=opener, sleep=sleep)
    facts, extra_claims, comments, rest = read_repository(api, window)
    measured, records = measures(facts, extra_claims, comments, window), None
    # The breakdown never costs a key measured before it: when its reads
    # would reach the board's reserve, or GitHub fails one of them, the
    # record is written without it, the record a copy of the template from
    # before the breakdown writes.
    try:
        broken_out, unread, units = read_breakdown(api, facts, rest, comments, window)
    except OutOfBudget:
        left = ("; the breakdown was left out to keep the board's reserve, so the record has no parallel, "
                "by_setting, pairs, economics or coordination")
    except GitHubError as problem:
        answer = f"HTTP {problem.status}" if problem.status else "no answer"
        left = (f"; the breakdown was left out because GitHub failed one of its reads ({answer}), so the record "
                "has no parallel, by_setting, pairs, economics or coordination")
    else:
        records = usage_records(api, window, args.usage_cache)
        measured.update(broken_out)
        measured.update(economic_keys(facts, units, window, records))
        left = (f"; the changed files of {unread} pull requests were left unread to keep the board's reserve, so "
                "the counts that need them are null") if unread else ""
        if records.unread is None:
            left += "; the Actions runs' usage records could not be listed, so only those read before count"
        elif records.unread:
            left += f"; {records.unread} Actions runs' usage records were left for the next run"
    write_json(args.out, metrics_record(measured, window, protocol_version()))
    written = args.out
    if args.profile_out:
        # PROFILE.json is summed only once METRICS.json is written, from
        # what collect already read, so that no fault in its sums can stop
        # the record the collector reads; a run that writes no PROFILE.json
        # leaves the branch's own in place (metrics.yml).
        try:
            profiled = profile_keys(facts, units, window, records) if records is not None else {}
        except Exception as fault:
            # A run that exits 0 would leave the branch's older PROFILE.json
            # in place unnoticed, so the fault is an annotation on the run.
            left += f"; {PROFILE_FILE} was not written: {type(fault).__name__}: {fault}"
            print(f"::warning::{PROFILE_FILE} was not written, so the branch keeps its older one: "
                  f"{type(fault).__name__}")
        else:
            write_json(args.profile_out, profile_record(profiled, window))
            written = f"{args.out} and {args.profile_out}"
    print(f"metrics: {written} written over {args.window_days} days, in {api.calls} calls{left}")
    return 0


def unchanged(args):
    try:
        old, new = Path(args.old).read_text(encoding="utf-8"), Path(args.new).read_text(encoding="utf-8")
    except (OSError, UnicodeDecodeError):
        return 1
    return 0 if same_apart_from_time(old, new) else 1


def positive(text):
    value = int(text)
    if value < 1:
        raise argparse.ArgumentTypeError("must be 1 or more")
    return value


def model_name(text):
    if not MODEL_NAME.fullmatch(text):
        raise argparse.ArgumentTypeError("a model name starts with a letter or digit, goes on with letters, digits, "
                                         "dots, dashes and underscores, at most 64 in all, and may end in [a tag]; "
                                         f"{text!r} is written to JUDGE.json as it is, so it is refused")
    return text


def parse_args(argv):
    parser = argparse.ArgumentParser(description="Field metrics of the claim protocol (docs/setup/metrics.md).")
    commands = parser.add_subparsers(dest="command", required=True)
    gather = commands.add_parser("collect", help="write METRICS.json from GitHub's REST API (GITHUB_TOKEN)")
    gather.add_argument("--repo", required=True, help="OWNER/NAME")
    gather.add_argument("--out", required=True, help="where METRICS.json goes")
    gather.add_argument("--root", default=".", help=f"the checkout holding {SWITCH_FILE} (default: here)")
    gather.add_argument("--usage-cache", help="a file of the Actions runs' usage records read before, read and "
                                              "written back, so a run downloads only the new ones")
    gather.add_argument("--profile-out", help=f"where {PROFILE_FILE} goes, the profile's sums (default: not written)")
    code = commands.add_parser("judge", help="code recent pull request conversations with your own claude -p")
    code.add_argument("--repo", help="OWNER/NAME (default: the repository gh sees here)")
    code.add_argument("--limit", type=positive,
                      help="judge at most this many pull requests, newest first; every listed pull request's "
                           "commits are still read to sort them into work settings")
    code.add_argument("--dry-run", action="store_true", help="show what would be sent and written; send nothing")
    code.add_argument("--model", default=DEFAULT_MODEL, type=model_name,
                      help=f"the model claude runs (default: {DEFAULT_MODEL})")
    code.add_argument("--out", default=JUDGE_FILE, help=f"where JUDGE.json goes (default: ./{JUDGE_FILE})")
    code.add_argument("--push", action="store_true",
                      help=f"commit JUDGE.json onto the branch {METRICS_BRANCH} of origin, the repository judged")
    for each in (gather, code):
        each.add_argument("--window-days", type=positive, default=WINDOW_DAYS, help="days to look back (default: 30)")
        each.add_argument("--now", help="the window's end, an ISO time (default: now; for tests)")
    send = commands.add_parser("push", help=f"push a JUDGE.json the judge wrote onto the branch {METRICS_BRANCH}")
    send.add_argument("--repo", help="OWNER/NAME the file was judged for (default: the repository gh sees here)")
    send.add_argument("--out", default=JUDGE_FILE, help=f"the file, as judge --out wrote it (default: ./{JUDGE_FILE})")
    same = commands.add_parser("unchanged", help=f"exit 0 when two {METRICS_FILE}, or two {PROFILE_FILE}, differ "
                                                 "only in generated_at")
    same.add_argument("old")
    same.add_argument("new")
    return parser.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    try:
        if args.command == "collect":
            return collect(args, os.environ)
        if args.command == "judge":
            return judge(args, os.environ, Gh(), LocalClaude(args.model))
        if args.command == "push":
            return push(args, Gh())
        return unchanged(args)
    except (GitHubError, JudgeError) as problem:
        print(f"metrics: {problem}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
