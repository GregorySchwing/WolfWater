"""The board: what a repository's open work looks like, as one page.

The board is the file STATE.md on the branch "board": the open claims, the
open pull requests and how each stands to the default branch, the files
that more than one piece of open work changes and the order that work
should land in, the queue of unclaimed issues, the issues given to
workers, the blockers and the pull requests merged last, so every session
and every worker reads the same page instead of asking the API.
This module is the pure part. A Snapshot holds what GitHub knows, filled by
an adapter or read back from JSON; render turns it into the page, with no
wall-clock value and every row sorted, so the same state gives the same
bytes and a workflow writes the branch only when something moved; and
comments_due says which pull requests the board would speak to, and what
it would say, for the adapter to post once each. A claim the coordinator
released (released_keys) leaves the page as soon as the snapshot shows
the release, its issue back in the queue, and the adapter renames its
branch so the name is free for the next claim. Nothing here calls GitHub.
"""

import json
import unicodedata
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import NamedTuple

from . import labels
from .claims import (
    CLAIM_DIR, IDLE_AFTER, PREFIX, RELEASED_PREFIX, branch_of, by_number, claim_files, claim_worker,
    is_claim_branch, is_idle, key_of, number_in, number_of, parse_time, was_released,
)

RECENT_MERGED = 5
# Where a workflow writes the board when nothing else is named: a file in
# the current directory, since the workflow checks out the branch "board".
STATE_PATH = Path("STATE.md")
NO_CLAIM_FILE = "no claim file yet"
# A claim file GitHub gave no text of: its contents API gives the text of a
# file of up to 1 MB, and answers a bigger one with an empty content and
# the encoding "none". The snapshot holds it as False, where a text is the
# file's and None is no file yet, and the Worker cell says so when no label
# names the worker either.
UNREAD_CLAIM = False
CLAIM_FILE_NOT_READ = "claim file over 1 MB, not read"
# A claim branch that outlived its issue: the row stays, and says what to do.
ISSUE_CLOSED = "issue closed"
DELETE_BRANCH = "delete the branch"
NO_PULL_YET = "(no pull request yet)"
# The issue cell of a branch under claude/ that names no issue, as
# claude/x does, or claude/03, which reads as a claim on #3 and is none.
NOT_A_CLAIM = "(not a claim)"

MAIN_MOVED_MARKER = "<!-- sdlc-board: main-moved {sha} -->"
IDLE_MARKER = "<!-- sdlc-board: idle -->"
# The coordinator releases a draft's claim in the same run the board would
# call it idle, so the comment reaches a pull request marked ready, which
# waits however long its review takes, or one from a branch that is no claim.
IDLE_TEXT = (
    "idle: no change for more than a day. A pull request marked ready waits for its review and merge, "
    "however long that takes; a draft whose branch goes a day without a commit is released by the coordinator."
)
# What the board says on an issue once it has renamed a released claim's
# branch; the marker names the new branch, so each release is told once.
RELEASED_MARKER = "<!-- sdlc-board: released {branch} -->"
# What the board says on an issue when GitHub refused for good to rename
# its released claim's branch; the marker names the branch the rename was
# to make, so each refusal is told once, and a newer head's again.
RENAME_REFUSED_MARKER = "<!-- sdlc-board: rename-refused {branch} -->"

CLAIM_HEADING = "## Claims (open claude/<issue> branches and docs/agents/claims/<issue>.md)"
CLAIM_COLUMNS = ("Issue", "Branch", "Worker", "PR", "Files it touches")
PULL_HEADING = "## Open pull requests"
PULL_COLUMNS = ("PR", "Title", "Branch", "State (clean, behind, conflicting, blocked, unknown)", "Last update")
OVERLAP_HEADING = ("## Overlaps (files that more than one piece of open work changes; land the work in this "
                   "order: a blocker first, then work marked ready, then the oldest)")
OVERLAP_COLUMNS = ("File", "Open work, in the order it should land")
QUEUE_HEADING = "## Queue (open issues labelled todo, unclaimed and given to no worker)"
ISSUE_COLUMNS = ("Issue", "Title", "Labels")
GIVEN_HEADING = "## Given to workers (open issues a worker's label names, not yet claimed)"
GIVEN_COLUMNS = ("Issue", "Title", "Worker")
BLOCKER_HEADING = "## Blockers (open issues labelled blocker)"
# The last section, drawn only where an open epic or goal has sub-issues,
# or where the run left it out: each open epic's tasks and their blocked-by
# relations as a Mermaid flowchart, and the goals above their epics.
DEPENDENCIES_HEADING = "## Dependencies"
DEPENDENCIES_PREFACE = (
    "Each open epic's tasks as GitHub's sub-issues and blocked-by relations hold them: an arrow runs from the "
    "issue that blocks to the issue it blocks, dotted once the blocker is closed; a node's shade is its state, "
    "and an epic's heading counts its closed tasks."
)
GOALS_HEADING = "### Goals"
# Why the section was left out, by the word the adapter writes in the
# snapshot (dependencies_skipped): its reads would have reached the calls
# the board keeps in reserve, or GitHub failed one of them.
NOT_DRAWN = {
    "reserve": ("Not drawn this run: its reads would reach the calls the board keeps in reserve; "
                "the next run tries again."),
    "failed": "Not drawn this run: GitHub did not answer one of its reads; the next run tries again.",
}
# A node's class is its state, one of these, in this order of precedence
# once closed is settled: review, claimed, blocked, todo; an issue in no
# state is drawn plain.
NODE_CLASSES = (
    "classDef todo fill:#ffffff,stroke:#0e8a16,color:#000000",
    "classDef claimed fill:#fff8d6,stroke:#fbca04,color:#000000",
    "classDef review fill:#e3f0ff,stroke:#1d76db,color:#000000",
    "classDef blocked fill:#ffe7e7,stroke:#e99695,color:#000000",
    "classDef closed fill:#eeeeee,stroke:#999999,color:#666666",
)
# A node's label is cut to this many characters, so a wide title keeps
# the chart readable.
NODE_LABEL_WIDTH = 40
# The same block, written into each open epic's body between these
# markers, so the epic's own page shows its graph; the plan section the
# clerk writes (team.PLAN_MARKER) stands beside it.
GRAPH_MARKER = "<!-- sdlc-board: graph -->"
GRAPH_END_MARKER = "<!-- sdlc-board: /graph -->"
MERGED_HEADING = f"## Recently merged (the last {RECENT_MERGED})"
MERGED_COLUMNS = ("PR", "Title", "Branch", "Merged")

# The five words the board uses for how a pull request stands to the
# default branch, and GitHub's mergeable_state words that mean each.
STATE_WORDS = ("clean", "behind", "conflicting", "blocked", "unknown")
GITHUB_STATES = {"has_hooks": "clean", "unstable": "clean", "dirty": "conflicting"}
# GitHub's state_reason for an open issue that was closed and opened again.
REOPENED = "reopened"


# --- the snapshot -----------------------------------------------------------

@dataclass
class Issue:
    """An open issue as the queue and the blockers show it, and whether
    GitHub says it was reopened (its state_reason), which the coordinator
    reads to free one whose labels still speak of a claim that ended with
    its close (coordinator.reopened_due; F22); with, for the Dependencies
    section, its state (a closed sub-issue or an issue blocking one is
    drawn too), how many issues GitHub's dependencies say block it, open
    and in all, from the summary GitHub carries on every issue
    (labels.blocked_by_counts), the relation and never the label blocker,
    and its sub-issues' progress, (completed, total), from the other
    summary."""

    number: int
    title: str = ""
    labels: list = field(default_factory=list)
    updated_at: str = ""
    reopened: bool = False
    state: str = "open"
    blocked_by_open: int = 0
    blocked_by_total: int = 0
    progress: tuple | None = None

    @classmethod
    def from_dict(cls, data):
        blocked_by_open, blocked_by_total = labels.blocked_by_counts(data)
        return cls(
            number=int(data["number"]),
            title=data.get("title") or "",
            labels=labels.label_names(data),
            updated_at=data.get("updated_at") or "",
            reopened=data.get("state_reason") == REOPENED,
            state="closed" if data.get("state") == "closed" else "open",
            blocked_by_open=blocked_by_open,
            blocked_by_total=blocked_by_total,
            progress=progress_of(data.get("sub_issues_summary")),
        )

    def to_dict(self):
        """Its JSON shape; state_reason only for an issue reopened, so a
        snapshot saved before F22's fix is written back as it was, and the
        state, the dependency summary and the sub-issue summary only when
        they say something."""
        return {"number": self.number, "title": self.title, "labels": list(self.labels), "updated_at": self.updated_at} \
            | ({"state_reason": REOPENED} if self.reopened else {}) \
            | ({"state": "closed"} if self.closed else {}) \
            | ({"issue_dependencies_summary": {"blocked_by": self.blocked_by_open,
                                               "total_blocked_by": self.blocked_by_total}}
               if self.blocked_by_total else {}) \
            | ({"sub_issues_summary": {"total": self.progress[1], "completed": self.progress[0]}}
               if self.progress else {})

    @property
    def closed(self):
        return self.state == "closed"


def progress_of(summary):
    """(completed, total) from GitHub's sub_issues_summary, or None where
    there is no summary or it counts no sub-issue."""
    if not isinstance(summary, dict):
        return None
    total, completed = labels.whole_count(summary.get("total")), labels.whole_count(summary.get("completed"))
    return (completed, total) if total else None


@dataclass
class Pull:
    """An open pull request with its head branch, GitHub's word for how it
    stands to the default branch (mergeable_state, as the API gives it),
    the paths it changes, and whether it is still a draft: one marked ready
    waits for review however long it takes, so its claim is never released,
    and a pull request whose draft state is not given counts as ready. The
    dict shape is the API's, head and all, so a pull request detail from
    GitHub is a fixture as it is; standing is the adapter's own key, how
    the branch stands where GitHub's word says only that something else
    stops the merge (adapters.github.standing_of). created_at, when the
    pull request was opened, is when its work began where GitHub says
    nothing of its branch (open_work); fork, whether its branch stands in
    another repository, or in none, a fork since deleted, whose branch a
    branch of the same name here is not (from_fork), and stranger, whether
    someone who may not write to the repository opened it (from_stranger),
    so that either is never a claim's (claim_of)."""

    number: int
    title: str = ""
    branch: str = ""
    head_sha: str = ""
    mergeable_state: str | None = None
    updated_at: str = ""
    files: list = field(default_factory=list)
    draft: bool = False
    standing: str | None = None
    created_at: str = ""
    fork: bool = False
    stranger: bool = False

    @property
    def state(self):
        """The board's word for the pull request: GitHub's mergeable_state
        as state_word reads it, except where that word says only that
        something else stops the merge (blocked, or the draft itself) and
        the adapter read how the branch stands through mergeable and a
        comparison (F9): behind and conflicting win, since the main-moved
        rule speaks to them, and a draft whose branch is clean reads clean
        rather than unknown; a blocked branch that is clean stays blocked,
        which says what still stops it."""
        word = state_word(self.mergeable_state)
        if self.standing in ("behind", "conflicting"):
            return self.standing
        if word == "unknown" and self.standing == "clean":
            return "clean"
        return word

    @classmethod
    def from_dict(cls, data):
        head = data.get("head") or {}
        return cls(
            number=int(data["number"]),
            title=data.get("title") or "",
            branch=head.get("ref") or "",
            head_sha=head.get("sha") or "",
            mergeable_state=data.get("mergeable_state"),
            updated_at=data.get("updated_at") or "",
            files=[file_name(f) for f in data.get("files") or []],
            draft=bool(data.get("draft")),
            standing=data.get("standing"),
            created_at=str(data.get("created_at") or ""),
            fork=from_fork(data),
            stranger=from_stranger(data),
        )

    def to_dict(self):
        data = {
            "number": self.number,
            "title": self.title,
            "head": {"ref": self.branch, "sha": self.head_sha},
            "mergeable_state": self.mergeable_state,
            "updated_at": self.updated_at,
            "files": list(self.files),
            "draft": self.draft,
        }
        if self.standing is not None:
            data["standing"] = self.standing
        if self.created_at:
            data["created_at"] = self.created_at
        if self.fork:
            data["fork"] = True
        if self.stranger:
            data["stranger"] = True
        return data


def from_fork(data):
    """Whether a pull request's branch stands in another repository than the
    one it would merge into, or in none, a fork since deleted, as the API's
    head and base name them; a dict that does not name both by their full
    names, as a fixture or a snapshot written by hand, which may hold
    anything there, says so by the board's own key fork, which
    Pull.to_dict writes."""
    head, base = mapping(data.get("head")), mapping(data.get("base"))
    base_name = mapping(base.get("repo")).get("full_name")
    if "repo" not in head or not isinstance(base_name, str) or not base_name:
        return bool(data.get("fork"))
    if head["repo"] is None:
        return True
    head_name = mapping(head["repo"]).get("full_name")
    if not isinstance(head_name, str):
        return bool(data.get("fork"))
    return head_name.lower() != base_name.lower()


def from_stranger(data):
    """Whether someone who may not write to the repository opened a pull
    request, by the board's own key stranger, which the adapter writes on a
    pull request from a branch here named as a claim's
    (adapters.github.mark_strangers) and Pull.to_dict writes back: who may
    write is an answer GitHub gives the adapter, not one the pull request's
    own fields always hold (claims.opened_by_writer)."""
    return bool(data.get("stranger"))


def mapping(value):
    """The value when it is a dict, else an empty one."""
    return value if isinstance(value, dict) else {}


@dataclass
class Merged:
    """A pull request that landed: enough for the recently merged table."""

    number: int
    title: str = ""
    branch: str = ""
    merged_at: str = ""

    @classmethod
    def from_dict(cls, data):
        return cls(
            number=int(data["number"]),
            title=data.get("title") or "",
            branch=(data.get("head") or {}).get("ref") or "",
            merged_at=data.get("merged_at") or "",
        )

    def to_dict(self):
        return {"number": self.number, "title": self.title, "head": {"ref": self.branch}, "merged_at": self.merged_at}


@dataclass
class ClaimIssue:
    """The issue behind a claim branch: its state and its label names, so
    the board can say that a branch outlived its issue, and can name the
    account from the issue's label once the claim file is gone."""

    state: str = ""
    labels: list = field(default_factory=list)

    @classmethod
    def from_dict(cls, data):
        return cls(state=data.get("state") or "", labels=labels.label_names(data))

    def to_dict(self):
        return {"state": self.state, "labels": list(self.labels)}


@dataclass
class Release:
    """What the release rule reads about one claim: when the latest release
    comment on its issue was made (claims.latest_release), "" when none
    stands; the head of its branch and when its last commit came
    (last_commit), which is when GitHub received the push that set the
    head, as the repository's activity records it (pushed_at), so no
    worker's clock decides it (F8), or what the adapter dates the head by
    in that push's place (adapters.github.head_since): an entry of the
    branch that came after the run began, whatever heads it names, the
    newest entry the reads listed where it is later than the committer
    date, or, in the run the branch's own push started, the run's start,
    each in pushed_at too, since each dates the last commit as a push
    would, or, where none of those holds, when the head was committed, as
    its committer date says (committed_at), one of the two read and the
    other ""; and whether a
    release already kept that head, a branch released/<issue>-<short sha>
    standing at it (kept), so the branch was pushed again with nothing a
    release has not kept. An adapter fills it for every claim a release can
    be due on or can have freed (adapters.github.releases): one gone more
    than a day without a commit, or whose head a release kept, whose issue
    is open and not labelled in-review, and whose pull request, if it has
    one, is a draft its author has not called finished, or one its author
    called finished that a person has released. The coordinator releases
    those that no release comment has released yet
    (coordinator.releases_due), and the board counts a release once the
    label claimed is off too (released_keys)."""

    released_at: str = ""
    sha: str = ""
    committed_at: str = ""
    pushed_at: str = ""
    kept: bool = False

    @property
    def last_commit(self):
        """When the branch's last commit came: GitHub's time of the push
        that set the head, or what the adapter dates the head by in its
        place, an entry after the run began, the newest entry listed or
        the start of the run the branch's push started (pushed_at), else
        the head's committer date, "" when neither is known."""
        return self.pushed_at or self.committed_at

    @classmethod
    def from_dict(cls, data):
        return cls(released_at=data.get("released_at") or "", sha=data.get("sha") or "",
                   committed_at=data.get("committed_at") or "", pushed_at=data.get("pushed_at") or "",
                   kept=data.get("kept") is True)

    def to_dict(self):
        """Its JSON shape; pushed_at and kept only when set, so a snapshot
        read before F8 is written back as it was."""
        return {"released_at": self.released_at, "sha": self.sha, "committed_at": self.committed_at} \
            | ({"pushed_at": self.pushed_at} if self.pushed_at else {}) | ({"kept": True} if self.kept else {})


@dataclass
class Snapshot:
    """Everything the board renders, in the shape a --dry-run fixture has:
    the default branch's commit and name, the open issues (pull requests
    left out), the open pull requests with their standing and files, the
    branch names, the claim files by issue key (None for a branch with no
    file yet, UNREAD_CLAIM for one GitHub gave no text of), the issue
    behind each claim branch, what a release reads about a claim by issue
    key, the merges last, and when the work that meets other work began,
    by its branch: when GitHub says the push that made the branch came, as
    the adapter reads it for the overlaps' landing order (began_at), with
    GitHub's answer for each branch it refused to date, which the run's log
    names (adapters.github.began_refused_lines); and why each claim
    branch whose head no read dated by its push was dated otherwise, by
    its committer date, by the newest entry listed or by the start of the
    run its push started, GitHub's refusal or an activity that lists no
    push of the head, with what was read in its place
    (adapters.github.pushed_refused_lines); and the branches a release
    left standing, released/<issue>-<short sha>, by name, whose work a
    person may put back, so the coordinator frees no reopened issue whose
    release is being undone (coordinator.claim_stands; F22), with the
    reopened issues, by key, whose released branches an earlier claim's
    release left, its latest release comment and the label claimed put on
    no later than the issue's last close, or claimed put on by a claim
    made since, so they undo nothing
    (adapters.github.ended_releases)."""

    sha: str = ""
    default_branch: str = "main"
    issues: list = field(default_factory=list)
    pulls: list = field(default_factory=list)
    branches: list = field(default_factory=list)
    claims: dict = field(default_factory=dict)
    claim_issues: dict = field(default_factory=dict)
    merged: list = field(default_factory=list)
    releases: dict = field(default_factory=dict)
    began: dict = field(default_factory=dict)
    began_refused: dict = field(default_factory=dict)
    pushed_refused: dict = field(default_factory=dict)
    kept_names: list = field(default_factory=list)
    ended_releases: list = field(default_factory=list)
    # The Dependencies section's reads (adapters.github.dependencies): each
    # open epic's and goal's sub-issues, by number; the issues that block
    # each sub-issue, by number, ascending; the issues those name that the
    # open list does not hold (closed sub-issues, blockers outside), by
    # number; the word for why the section was left out (NOT_DRAWN), ""
    # when it was not, and what the run's log says of it.
    sub_issues: dict = field(default_factory=dict)
    blocked_by: dict = field(default_factory=dict)
    others: dict = field(default_factory=dict)
    dependencies_skipped: str = ""
    dependencies_notice: str = ""
    # The chart each open issue's body held between the graph markers when
    # the run listed the issues, by number, folded to LF (graph_block):
    # what the adapter's graph writer compares its block with before it
    # reads the epic again, and the issues whose chart it takes out once
    # they owe none.
    graphs: dict = field(default_factory=dict)

    @classmethod
    def from_dict(cls, data):
        """A snapshot from its JSON shape, which is also what the API gives:
        labels may be objects or names, branches objects or names."""
        return cls(
            sha=data.get("sha") or "",
            default_branch=data.get("default_branch") or "main",
            issues=[Issue.from_dict(i) for i in data.get("issues") or []],
            pulls=[Pull.from_dict(p) for p in data.get("pulls") or []],
            branches=[branch_name(b) for b in data.get("branches") or []],
            claims={str(k): v for k, v in (data.get("claims") or {}).items()},
            claim_issues={str(k): ClaimIssue.from_dict(v or {}) for k, v in (data.get("claim_issues") or {}).items()},
            merged=[Merged.from_dict(m) for m in data.get("merged") or []],
            releases={str(k): Release.from_dict(v or {}) for k, v in (data.get("releases") or {}).items()},
            began={str(k): str(v) for k, v in (data.get("began") or {}).items() if v},
            began_refused={str(k): str(v) for k, v in (data.get("began_refused") or {}).items()},
            pushed_refused={str(k): str(v) for k, v in (data.get("pushed_refused") or {}).items()},
            kept_names=[branch_name(b) for b in data.get("kept_names") or []],
            ended_releases=[str(k) for k in data.get("ended_releases") or []],
            sub_issues={int(k): [int(n) for n in v] for k, v in (data.get("sub_issues") or {}).items()},
            blocked_by={int(k): sorted(int(n) for n in v) for k, v in (data.get("blocked_by") or {}).items()},
            others={int(k): Issue.from_dict(v) for k, v in (data.get("others") or {}).items()},
            dependencies_skipped=str(data.get("dependencies_skipped") or ""),
            dependencies_notice=str(data.get("dependencies_notice") or ""),
            graphs={int(k): str(v) for k, v in (data.get("graphs") or {}).items()},
        )

    def to_dict(self):
        return {
            "sha": self.sha,
            "default_branch": self.default_branch,
            "issues": [i.to_dict() for i in self.issues],
            "pulls": [p.to_dict() for p in self.pulls],
            "branches": list(self.branches),
            "claims": dict(self.claims),
            "claim_issues": {k: v.to_dict() for k, v in self.claim_issues.items()},
            "merged": [m.to_dict() for m in self.merged],
            "releases": {k: v.to_dict() for k, v in self.releases.items()},
            "began": dict(self.began),
        } | ({"began_refused": dict(self.began_refused)} if self.began_refused else {}) \
            | ({"pushed_refused": dict(self.pushed_refused)} if self.pushed_refused else {}) \
            | ({"kept_names": list(self.kept_names)} if self.kept_names else {}) \
            | ({"ended_releases": list(self.ended_releases)} if self.ended_releases else {}) \
            | ({"sub_issues": {str(k): list(v) for k, v in self.sub_issues.items()}} if self.sub_issues else {}) \
            | ({"blocked_by": {str(k): list(v) for k, v in self.blocked_by.items()}} if self.blocked_by else {}) \
            | ({"others": {str(k): v.to_dict() for k, v in self.others.items()}} if self.others else {}) \
            | ({"dependencies_skipped": self.dependencies_skipped} if self.dependencies_skipped else {}) \
            | ({"dependencies_notice": self.dependencies_notice} if self.dependencies_notice else {}) \
            | ({"graphs": {str(k): v for k, v in self.graphs.items()}} if self.graphs else {})

    def issue_by_number(self, number):
        """The issue with this number, among the open issues or the others
        the Dependencies reads brought, or a bare one when neither holds
        it, so a node is drawn for every number an edge names."""
        for issue in self.issues:
            if issue.number == number:
                return issue
        return self.others.get(number) or Issue(number)

    @classmethod
    def from_json(cls, text):
        return cls.from_dict(json.loads(text))

    def to_json(self):
        return json.dumps(self.to_dict(), indent=2) + "\n"


def branch_name(branch):
    """A branch as a name, whether the API's object or a fixture's string."""
    return branch if isinstance(branch, str) else branch["name"]


def file_name(changed):
    """A changed file as a path, whether the API's object from the pull
    request's files or a fixture's string."""
    return changed if isinstance(changed, str) else changed["filename"]


def as_snapshot(value):
    """A Snapshot from one, or from the dict shape a fixture or an API
    adapter holds, so render takes either."""
    return value if isinstance(value, Snapshot) else Snapshot.from_dict(value)


def as_pull(value):
    return value if isinstance(value, Pull) else Pull.from_dict(value)


# --- the rules --------------------------------------------------------------

def state_word(mergeable_state):
    """One of the five words the board uses for a pull request, from
    GitHub's mergeable_state. GitHub's has_hooks and unstable are mergeable
    states (the checks are what they are), so they read as clean; dirty is
    conflicting; blocked is a computed state of its own; unknown is kept for
    a state GitHub has not computed yet. The board's own words pass through,
    so a snapshot that already speaks them renders the same."""
    if mergeable_state in STATE_WORDS:
        return mergeable_state
    return GITHUB_STATES.get(mergeable_state, "unknown")


def pr_state(pull):
    """The board's word for a pull request, given as a Pull or as the dict
    the API returns, with the adapter's standing where it read one."""
    if isinstance(pull, Pull):
        return pull.state
    return Pull(0, mergeable_state=pull.get("mergeable_state"), standing=pull.get("standing")).state


def short_sha(sha):
    return (sha or "")[:7] or "unknown"


def issue_ref(key):
    """The issue cell of a claim's row: #<number>, or NOT_A_CLAIM for a key
    that names no issue, which shown as itself could read as one (03, #3,
    or a digit of another script that the page folds into an ASCII one)."""
    return f"#{key}" if number_in(key) is not None else NOT_A_CLAIM


def escaped(name):
    """A name with each character outside ASCII written as its escape, so
    the page, which folds such characters into ASCII ones (cell), never
    shows a branch such as claude/ and a fullwidth 3 as claude/3."""
    return "".join(ch if ch.isascii() else escape(ch) for ch in name)


def escape(ch):
    """One character as Python writes it in a string: \\uXXXX, or \\UXXXXXXXX
    for one past U+FFFF."""
    code = ord(ch)
    return f"\\u{code:04x}" if code <= 0xFFFF else f"\\U{code:08x}"


def claim_branches(branches):
    return [branch_name(b) for b in branches if is_claim_branch(branch_name(b))]


def claimed_numbers(branches, released=()):
    """The issue numbers that have a claim branch, leaving out the claims
    named in released, by key, since a released claim holds nothing; a
    branch whose key is no number (claims.number_in) claims none."""
    numbers = (number_of(b) for b in claim_branches(branches) if key_of(b) not in released)
    return {number for number in numbers if number is not None}


def released_keys(snapshot, idle_after=IDLE_AFTER):
    """The claims released, by the coordinator or by a person who may
    write, by issue key: a claim branch whose issue is open and no longer
    labelled claimed, whose latest release comment came after the branch
    had gone without a commit for more than idle_after, and which nothing
    has been committed to since (claims.was_released), or whose head a
    release already kept, whenever its comment came (release_counts); the
    last commit is when GitHub received the push that set the head, or
    what the adapter dates the head by in its place, else its committer
    date (Release.last_commit, F8). It takes both acts
    of a release, the comment and the label claimed taken off, so a stray
    comment frees nothing; a commit after the release, its holder's or a
    new claim's, ends it; and putting the label claimed back undoes it. A
    finished claim is never released, however long it waits for review:
    one whose issue is labelled in-review or whose pull request is marked
    ready; a draft its author called finished only by a person, since the
    coordinator never releases one. The board leaves a released claim off
    the page and puts its issue back in the queue, and the GitHub adapter
    renames its branch, so the name is free for the next claim; the page
    reads the same before the rename and after it, so one release rewrites
    the board once."""
    snap = as_snapshot(snapshot)
    return {key_of(b) for b in claim_branches(snap.branches) if is_released(snap, key_of(b), idle_after)}


def is_released(snap, key, idle_after=IDLE_AFTER):
    """Whether the claim with this key was released, by the coordinator
    or by a person who may write, by released_keys' rule; a claim with no
    head commit known is not, since its branch could not be renamed."""
    release, issue = snap.releases.get(key), snap.claim_issues.get(key)
    if release is None or issue is None or not release.sha:
        return False
    if issue.state != "open" or labels.CLAIMED in issue.labels or labels.IN_REVIEW in issue.labels:
        return False
    if any(of_claim(p, key) and not p.draft for p in snap.pulls):
        return False
    return release_counts(release, idle_after)


def release_counts(release, idle_after=IDLE_AFTER):
    """Whether the latest release comment that counts released the claim
    as its branch stands: it came after the branch had gone more than
    idle_after without a commit, with none since (was_released), or the
    head is one a release already kept (Release.kept), since a branch
    pushed again at it brings nothing new, and its time, when GitHub
    received that push, says nothing of the work, which the release that
    kept the head already judged idle."""
    if not release.released_at:
        return False
    return release.kept or was_released(release.released_at, release.last_commit, idle_after)


def claim_of(pull):
    """The key of the claim whose pull request this is: the key its branch
    names (key_of), or None for a branch that is no claim, for a pull
    request from a fork, whatever its branch is called, and for one that
    someone who may not write to the repository opened, from the claim's
    own branch here. A claim branch lives in this repository, and a fork's
    branch of the same name is another branch; and a claim's pull request
    is opened by someone who may write, or by an App the repository
    installed, while GitHub's documentation says anyone who may read the
    repository may open one from any of its branches. So such a pull
    request, marked ready or said finished, finishes no claim, and the
    board lists it and tells it what it tells any other (F17: until its fix
    a ready pull request from a fork's claude/<issue> kept that claim from
    ever being released, and the board took it for the claim's own on the
    page and in the Overlaps; F19, the same for a stranger's pull request
    from the claim's own branch)."""
    return None if pull.fork or pull.stranger else key_of(pull.branch)


def of_claim(pull, key):
    """Whether the pull request is the claim's with this key (claim_of)."""
    return claim_of(pull) == key


def closed_by_release(pull, released):
    """Whether the rename of a released claim's branch closes the pull
    request: every open one from a branch released_branches names, here,
    the claim's own and a stranger's alike, since GitHub closes each pull
    request whose head is a branch it renames (renaming-a-branch.md:19);
    never one from a fork, whose branch of the same name stands in the
    fork (F17, F19)."""
    return not pull.fork and pull.branch in released


def released_branches(snap):
    """The claim branches of the released claims, whose pull requests the
    board no longer lists or speaks to: renaming the branch closes them."""
    return {branch_of(key) for key in released_keys(snap)}


def released_text(old, new):
    """What the board says on an issue once it has renamed the branch of
    its released claim."""
    return (
        f"released: the claim on this issue was released, so its branch {old} is now {new}, "
        f"which keeps its work; GitHub closes a pull request from the branch with the rename. "
        f"The issue is back in the queue, and the next claim pushes {old} anew from the default branch. "
        f"To undo the release, put the label claimed back on the issue first, then rename the branch "
        f"back and reopen its pull request."
    )


class Ruleset(NamedTuple):
    """A ruleset whose rules apply to a branch the board could not rename:
    its name, the page where an owner changes it, "" when there is none to
    link, and the types of the rules it applies to the branch."""

    name: str
    link: str
    types: tuple


def rename_refused_text(old, new, rulesets, words):
    """What the board says on an issue when GitHub refused for good to
    rename the branch of its released claim, in GitHub's words: the name
    stays taken, so nobody can claim the issue, until a board run renames
    it. A ruleset over the claim branches is what makes a rename an
    administrator's, so each one that applies to either name is named with
    its rules; rulesets is None when the board could not read them, and
    empty when none applies, which leaves the token's own right."""
    said = " ".join(str(words).split())
    lines = [f"rename refused: the claim on this issue was released, but GitHub refused to rename its "
             f"branch {old} to {new} ({said}), so the name {old} is still taken and nobody can claim the "
             f"issue yet. Every board run tries the rename again."]
    if rulesets:
        lines += ["", "These rulesets apply to the two names:"]
        lines += [f"- {f'[{r.name}]({r.link})' if r.link else r.name}: {', '.join(r.types)}" for r in rulesets]
        lines += ["", f"Take `{PREFIX}**` and `{RELEASED_PREFIX}**` out of their targets, so that they hold the "
                      f"default branch alone, and the next board run renames the branch."]
    elif rulesets is None:
        lines += ["", f"The board could not read which rulesets apply to the two names. A ruleset whose target "
                      f"takes in `{PREFIX}**` or `{RELEASED_PREFIX}**` makes the rename an administrator's: take "
                      f"those names out of its target, and the next board run renames the branch."]
    else:
        lines += ["", "No ruleset applies to the two names, so the board's token may lack the right to rename "
                      "a branch: its workflow needs contents: write."]
    return "\n".join(lines)


def pulls_behind_main(pulls):
    """The pull requests the main-moved rule speaks to: behind or
    conflicting, and only those, a draft or a blocked one included once
    the adapter has read how its branch stands (F9)."""
    return [p for p in pulls if pr_state(p) in ("behind", "conflicting")]


def idle_pulls(pulls, now, idle_after=IDLE_AFTER):
    """The pull requests with no change for more than a day."""
    return [p for p in pulls if pull_is_idle(as_pull(p), now, idle_after)]


def pull_is_idle(pull, now, idle_after=IDLE_AFTER):
    """Whether one pull request has waited too long. One with no update
    time is not idle rather than failing the run: there is nothing to
    measure, and GitHub always gives one."""
    return bool(pull.updated_at) and is_idle(pull.updated_at, now, idle_after)


def main_moved_text(state, default_branch, sha7):
    how = "behind" if state == "behind" else "conflicting with"
    return (
        f"main moved: this branch is {how} {default_branch} at {sha7}; "
        f"rebase or merge {default_branch} before the next push."
    )


class DueComment(NamedTuple):
    """A comment the board owes a pull request: posted only when no comment
    on it carries the marker, which is the adapter's job."""

    number: int
    marker: str
    text: str

    @property
    def kind(self):
        return "idle" if self.marker == IDLE_MARKER else "main-moved"


def comments_due(snapshot, event_name, ref_name, now, idle_after=IDLE_AFTER):
    """What the board says after this run, as (pull number, marker, text).
    Every pull request behind the default branch or in conflict with it, as
    its state word says (Pull.state), hears so, once per commit of the
    default branch (the marker carries the short sha), in whichever run
    first finds it so; behind is GitHub's word, or, where GitHub says only
    that something else stops the merge, the comparison's (F9), so one
    GitHub calls clean, which the repository's rules do not hold for being
    behind, hears nothing (F16): the word used to
    be owed only by the run the push to the default branch started, and
    that run is lost when the workflow's concurrency group replaces it
    behind a later event's (F12), so the word depends on the state alone,
    and event_name and ref_name, which the callers still pass, decide
    nothing. A pull request nobody touched for a day is asked once whether
    it is parked. Main-moved comments come first, then idle ones, each in
    the order the snapshot lists the pull requests, as the template's
    board posted them. A released claim's pull request hears nothing:
    renaming its branch closes it; a pull request from a fork, whatever
    its branch is called, hears what any other does (closed_by_release)."""
    snap = as_snapshot(snapshot)
    released = released_branches(snap)
    pulls = [p for p in snap.pulls if not closed_by_release(p, released)]
    due = []
    if snap.sha:
        sha7 = short_sha(snap.sha)
        marker = MAIN_MOVED_MARKER.format(sha=sha7)
        for pull in pulls_behind_main(pulls):
            due.append(DueComment(pull.number, marker, main_moved_text(pull.state, snap.default_branch, sha7)))
    for pull in idle_pulls(pulls, now, idle_after):
        due.append(DueComment(pull.number, IDLE_MARKER, IDLE_TEXT))
    return due


# --- the page ---------------------------------------------------------------

def ascii_fold(text):
    """Text folded to ASCII, accented letters to their plain ones and the
    rest dropped, so the page passes a repository's own ASCII check."""
    return unicodedata.normalize("NFKD", text).encode("ascii", "ignore").decode()


def cell(value):
    """A table cell: one line, ASCII, pipes escaped, so the board stays a
    valid table and passes a repository's own ASCII check."""
    text = str(value).replace("\r", " ").replace("\n", " ").replace("|", "\\|")
    return ascii_fold(text).strip()


def table(columns, rows):
    """A Markdown table; an empty one keeps its header and says (none), so
    the layout is the same whatever the repository holds."""
    lines = ["| " + " | ".join(columns) + " |", "|" + " --- |" * len(columns)]
    for row in rows or [["(none)"] + [""] * (len(columns) - 1)]:
        lines.append("| " + " | ".join(cell(c) for c in row) + " |")
    return lines


def claim_row(key, branch, text, pull, issue):
    """One claim as the board shows it. The worker is the claim file's, or
    the issue's label once the finishing pull request has deleted the file
    or when GitHub gave no text of it (UNREAD_CLAIM).
    A branch whose issue is closed is a leftover: the row says so in place
    of the pull request and asks for the branch to go in place of the files.
    A branch whose key names no issue claims nothing, and the row says so in
    place of the issue, with the branch's name escaped (escaped)."""
    closed = issue is not None and issue.state == "closed"
    if pull:
        files = sorted(pull.files)
    else:
        files = claim_files(text) if text else []
    if pull:
        pull_cell = f"#{pull.number}"
    elif closed:
        pull_cell = ISSUE_CLOSED
    else:
        pull_cell = NO_PULL_YET
    if files:
        files_cell = ", ".join(files)
    elif closed and not pull:
        files_cell = DELETE_BRANCH
    else:
        files_cell = "-"
    if text:
        worker = claim_worker(text)
    else:
        missing = CLAIM_FILE_NOT_READ if text is UNREAD_CLAIM else NO_CLAIM_FILE
        worker = labels.worker_of({"labels": issue.labels} if issue else None) or missing
    return [issue_ref(key), escaped(branch), worker, pull_cell, files_cell]


def claim_rows(snap):
    """A row per claim branch, in issue order, but a released claim's; the
    pull request a row names is the claim's own, never one from a fork
    (claim_of)."""
    pulls_by_branch = {p.branch: p for p in snap.pulls if claim_of(p) is not None}
    released = released_keys(snap)
    rows = []
    for branch in sorted(claim_branches(snap.branches), key=lambda b: by_number(key_of(b))):
        key = key_of(branch)
        if key in released:
            continue
        rows.append(claim_row(key, branch, snap.claims.get(key), pulls_by_branch.get(branch), snap.claim_issues.get(key)))
    return rows


def pull_rows(snap):
    """A row per open pull request, in number order, but one from a
    released claim's branch, which the rename of that branch closes; one
    from a fork stays, whatever its branch is called (closed_by_release)."""
    released = released_branches(snap)
    return [
        [f"#{p.number}", p.title, p.branch, p.state, p.updated_at]
        for p in sorted(snap.pulls, key=lambda p: p.number) if not closed_by_release(p, released)
    ]


def blocked_note(snap, issue):
    """What a row says after an issue's title when GitHub's dependencies
    say an open issue blocks it: the open blockers by number where the
    Dependencies reads brought them, else how many there are, since the
    summary every listed issue carries counts them and costs no read. ""
    for an issue blocked by nothing."""
    if issue.blocked_by_open <= 0:
        return ""
    blocking = [n for n in snap.blocked_by.get(issue.number, []) if not snap.issue_by_number(n).closed]
    if blocking:
        return " (blocked by " + ", ".join(f"#{n}" for n in blocking) + ")"
    count = issue.blocked_by_open
    return f" (blocked by {count} open issue{'s' if count != 1 else ''})"


def in_pick_order(issues):
    """The issues in the order Pick takes them within a place: by priority
    label, p1 first and the issues with none last (labels.priority), then
    by number, so the board lists what a worker will take first at the
    top. Without priority labels this is number order, as the board always
    listed them."""
    return sorted(issues, key=lambda i: (labels.priority(i.labels), i.number))


def issue_rows(snap, label, unclaimed_only):
    """The issues carrying a label, in Pick's order (in_pick_order); the
    queue leaves out an issue with the label claimed, a claim branch that
    is not released, an epic or a goal, or a worker's label, which gives
    the issue to that worker. An issue an open issue blocks stays listed
    and says so beside its title (blocked_note), since the blocked-by
    relation keeps it from being taken, not from being seen."""
    claimed = claimed_numbers(snap.branches, released_keys(snap))
    rows = []
    for issue in in_pick_order(snap.issues):
        names = sorted(issue.labels)
        if label not in names:
            continue
        if unclaimed_only and (labels.CLAIMED in names or issue.number in claimed or labels.reserved(names)):
            continue
        rows.append([f"#{issue.number}", issue.title + blocked_note(snap, issue), ", ".join(names)])
    return rows


def given_rows(snap):
    """The open issues a worker's label gives to that worker and that
    nobody has claimed yet, in Pick's order (in_pick_order), whatever else
    they are labelled: the spec says an issue a worker's label names is
    that worker's alone, left out of the Queue and shown here until the
    worker claims it, so the queue leaves them out and this is where the
    worker and everyone else see them. Until F14's fix the rows needed
    todo and left epics out, so an issue labelled worker:<name> alone was
    shown nowhere."""
    claimed = claimed_numbers(snap.branches, released_keys(snap))
    rows = []
    for issue in in_pick_order(snap.issues):
        names = sorted(issue.labels)
        workers = labels.holder_names(names)
        if labels.CLAIMED in names or not workers or issue.number in claimed:
            continue
        rows.append([f"#{issue.number}", issue.title + blocked_note(snap, issue), ", ".join(workers)])
    return rows


class Work(NamedTuple):
    """One piece of open work, as the overlaps read it: its number, the pull
    request's or, for a claim with no pull request yet, the issue's; the
    files it changes or, before it changes any, the files its claim says it
    will touch; whether it is marked ready and whether its issue is a
    blocker; whether it is a pull request at all; its branch, "" for a
    pull request from a fork; when it began, as GitHub writes a time, ""
    when nothing says (open_work); and the issue it is a claim's work for,
    None for a pull request that is no claim's (identity)."""

    number: int
    files: tuple
    ready: bool = False
    blocker: bool = False
    is_pull: bool = True
    branch: str = ""
    began: str = ""
    issue: int | None = None


def identity(work):
    """The number that names a piece of open work for as long as it is
    open, in the overlap markers and as the place it is told: the issue,
    for a claim and the pull request later opened from its branch, which
    are one piece, so opening the pull request tells nothing again; the
    pull request's own number for one that is no claim's, from a fork,
    from a branch that is no claim or opened by someone who may not write
    (spec, Overlaps; F2: until the fix the work was named by its pull
    request's number as soon as it had one, and told on it, so a claim
    that opened its pull request was told again of the work it was told
    of, and that work told again of it). An issue's number names each
    claim on it in turn, which the coordinator tells apart by when the
    claim's branch was made (coordinator.heard_of)."""
    return work.number if work.issue is None else work.issue


# Where work whose beginning nothing says stands in the landing order: with
# the newest. Only a claim with no pull request yet can be such work, since
# a pull request began at its opening where GitHub says nothing of its
# branch, and the protocol has a claim open its draft at once (Work in the
# open), so such a claim is all but always the newest piece; counted as the
# oldest, a claim GitHub would not date in the run that tells an overlap
# was told as landing before work days older (F13's review). The number
# orders such work among itself.
UNDATED = datetime.max.replace(tzinfo=timezone.utc)


def landing_key(work):
    """The order open work that meets in a file should land in: a blocker
    first, since it stops everyone; then work marked ready, which is done
    and waits only for its merge; then the oldest, the work that began
    first, by the time GitHub gives (began_at), and not by its number,
    which for a claim is its issue's and for a pull request the pull
    request's, two numbers that tell nothing of which began first (F13);
    the number breaks a tie."""
    return (not work.blocker, not work.ready, began_at(work.began), work.number)


def began_at(text):
    """When a piece of work began, as an aware datetime, from the time
    GitHub wrote; UNDATED when there is none or it does not read as a
    time, since a fixture or a snapshot saved by hand may hold anything,
    a number or a list included, and the board never stops for one
    value."""
    try:
        return parse_time(text) if text else UNDATED
    except (AttributeError, TypeError, ValueError):
        return UNDATED


def own_files(files):
    """The files a piece of work changes, without claim files: each claim
    has its own, so they never meet."""
    return tuple(sorted(f for f in files if not f.startswith(CLAIM_DIR + "/")))


def is_blocker(snap, key):
    issue = snap.claim_issues.get(key)
    return bool(issue) and labels.BLOCKER in issue.labels


def open_work(snapshot):
    """Every piece of open work, in number order: each open pull request but
    a released claim's, with the files it changes, or the files its claim
    file lists while it changes nothing but that file; and each claim with
    no pull request yet, with the files its claim file lists. Each began
    when the push that made its branch came, as the snapshot holds it by
    branch, and a pull request where pull_began says; a pull request from a
    fork carries no branch, since its branch is none of this
    repository's, and is a piece of its own, whatever its branch is called,
    since it is no claim's (claim_of), as is a stranger's pull request from
    a claim's branch here, which a release's rename closes with the claim's
    own (closed_by_release). A claim's own pull request carries the issue
    it is for, which names the piece (identity)."""
    snap = as_snapshot(snapshot)
    released = released_keys(snap)
    renamed = released_branches(snap)
    pulls = {p.branch: p for p in snap.pulls if claim_of(p) is not None}
    work = []
    for pull in snap.pulls:
        key = claim_of(pull)
        if closed_by_release(pull, renamed):
            continue
        files = own_files(pull.files) or (own_files(claim_files(snap.claims.get(key) or "")) if key else ())
        work.append(Work(pull.number, files, not pull.draft, bool(key) and is_blocker(snap, key), True,
                         "" if pull.fork else pull.branch, pull_began(snap, pull),
                         number_of(pull.branch) if key else None))
    for branch in claim_branches(snap.branches):
        key, number = key_of(branch), number_of(branch)
        if branch in pulls or key in released or number is None:
            continue
        files = own_files(claim_files(snap.claims.get(key) or ""))
        work.append(Work(number, files, False, is_blocker(snap, key), False, branch, snap.began.get(branch, ""),
                         number))
    return sorted(work, key=lambda w: w.number)


def pull_began(snap, pull):
    """When a pull request's work began: when the push that made its branch
    came, as the snapshot holds it, so a claim and the pull request later
    opened from its branch are one piece that began at the claim (F13);
    else, for a pull request from a branch the snapshot says nothing of,
    or nothing that reads as a time, or from a fork, whose branch a branch
    of the same name here is not, when it was opened."""
    made = "" if pull.fork else snap.began.get(pull.branch, "")
    return made if began_at(made) != UNDATED else pull.created_at


def overlaps(snapshot):
    """The files that more than one piece of open work changes, each with
    that work in the order it should land (landing_key), in file order; a
    claim's two pull requests, which GitHub allows to two bases, are one
    piece (identity), so a file they alone change meets nothing."""
    meeting = {}
    for work in open_work(snapshot):
        for path in work.files:
            meeting.setdefault(path, []).append(work)
    return {path: sorted(works, key=landing_key) for path, works in sorted(meeting.items())
            if len({identity(work) for work in works}) > 1}


def overlap_rows(snap):
    return [[path, ", ".join(f"#{w.number}" for w in works)] for path, works in overlaps(snap).items()]


def merged_rows(snap):
    return [
        [f"#{m.number}", m.title, m.branch, m.merged_at]
        for m in sorted(snap.merged, key=lambda m: m.merged_at, reverse=True)
    ]


# --- the dependencies -------------------------------------------------------

def node_id(number):
    return f"i{number}"


def node_label(issue):
    """A node's text: the number, the title folded to ASCII and cut to
    NODE_LABEL_WIDTH, and the sub-issues' progress where there is any.
    Mermaid reads a quoted label's entities, so the characters that would
    end the label or the node are written as entities, "#" first since
    the entities carry one, and so are "%", which ends nothing but would
    let a title open a %%{ directive that restyles or breaks the chart,
    and "$", since Mermaid typesets $$...$$ as math and a formula it
    cannot parse stops the whole chart from drawing; both go after the
    fold, which turns the small and the fullwidth percent and dollar
    signs into them. The control panel writes the same entities in the
    same order, so its charts and the board's are the same bytes."""
    title = " ".join(ascii_fold(issue.title).split())
    if len(title) > NODE_LABEL_WIDTH:
        title = title[:NODE_LABEL_WIDTH - 3].rstrip() + "..."
    for plain, entity in (("#", "#35;"), ("%", "#37;"), ("$", "#36;"), ('"', "#quot;"), ("|", "#124;"),
                          ("<", "#lt;"), (">", "#gt;")):
        title = title.replace(plain, entity)
    label = f"#{issue.number}" + (f" {title}" if title else "")
    if issue.progress and issue.progress[1] > 0:
        label += f" ({issue.progress[0]}/{issue.progress[1]})"
    return label


def node_class(issue, claimed):
    """A node's state as its class: closed first, then in review, then
    claimed (the label, or an unreleased claude/<issue> branch), then
    labelled blocked, then todo; "" for an issue in none of them."""
    names = issue.labels
    if issue.closed:
        return "closed"
    if labels.IN_REVIEW in names:
        return "review"
    if labels.CLAIMED in names or issue.number in claimed:
        return "claimed"
    if labels.BLOCKED in names:
        return "blocked"
    return "todo" if labels.TODO in names else ""


def node_line(issue, claimed):
    """One node of a flowchart: its id, its label in quotes, and its class
    where it has one."""
    suffix = f":::{node_class(issue, claimed)}" if node_class(issue, claimed) else ""
    return f'{node_id(issue.number)}["{node_label(issue)}"]{suffix}'


def mermaid_block(nodes, edges, claimed):
    """A fenced Mermaid flowchart: the nodes in number order, the edges
    sorted and each once, and the five classes, so the same state draws
    the same bytes."""
    lines = ["```mermaid", "flowchart TD"]
    lines += ["  " + node_line(issue, claimed) for issue in sorted(nodes, key=lambda i: i.number)]
    lines += ["  " + edge for edge in sorted(set(edges))]
    lines += ["  " + line for line in NODE_CLASSES]
    lines.append("```")
    return "\n".join(lines)


def blocked_by_edge(snap, blocker, blocked):
    """An arrow from the issue that blocks to the issue it blocks, dotted
    once the blocker is closed."""
    arrow = "-.->" if snap.issue_by_number(blocker).closed else "-->"
    return f"{node_id(blocker)} {arrow} {node_id(blocked)}"


def epic_block(snap, epic_number):
    """One open epic's flowchart: its sub-issues and the issues outside
    the epic that block one of them, with every blocked-by relation among
    them as an edge."""
    numbers = set(snap.sub_issues.get(epic_number, []))
    edges = []
    for blocked in sorted(numbers):
        for blocker in snap.blocked_by.get(blocked, []):
            numbers.add(blocker)
            edges.append(blocked_by_edge(snap, blocker, blocked))
    nodes = [snap.issue_by_number(n) for n in numbers]
    return mermaid_block(nodes, edges, claimed_numbers(snap.branches, released_keys(snap)))


def goals_block(snap, goals):
    """The goals above their epics: each open goal with sub-issues and the
    issues under it, joined by thick arrows and no blocked-by edge, since
    an epic's own block draws those."""
    numbers, edges = set(), []
    for goal in goals:
        numbers.add(goal.number)
        for epic in snap.sub_issues.get(goal.number, []):
            numbers.add(epic)
            edges.append(f"{node_id(goal.number)} ==> {node_id(epic)}")
    nodes = [snap.issue_by_number(n) for n in numbers]
    return mermaid_block(nodes, edges, claimed_numbers(snap.branches, released_keys(snap)))


def planned_with_sub_issues(snap, label):
    """The open issues carrying the label that the reads found sub-issues
    for, in number order."""
    return [issue for issue in sorted(snap.issues, key=lambda i: i.number)
            if label in issue.labels and snap.sub_issues.get(issue.number)]


def epic_heading(issue):
    """An epic's heading in the section, counting its closed tasks where
    GitHub's summary gives the count."""
    heading = f"### Epic #{issue.number}" + (f": {cell(issue.title)}" if cell(issue.title) else "")
    if issue.progress and issue.progress[1] > 0:
        heading += f" ({issue.progress[0]}/{issue.progress[1]})"
    return heading


def epic_blocks(snap):
    """The flowchart of every open epic with sub-issues, by epic number:
    what the section draws, and what the adapter writes into each epic's
    body (write_graphs)."""
    return {issue.number: epic_block(snap, issue.number) for issue in planned_with_sub_issues(snap, labels.EPIC)}


def dependencies_lines(snap):
    """The Dependencies section's lines, or none where there is nothing to
    draw and the reads were not left out, so a repository with no epic
    reads the page it always did."""
    if snap.dependencies_skipped:
        return ["", DEPENDENCIES_HEADING, "", NOT_DRAWN.get(snap.dependencies_skipped, NOT_DRAWN["failed"])]
    goals = planned_with_sub_issues(snap, labels.GOAL)
    blocks = epic_blocks(snap)
    if not goals and not blocks:
        return []
    lines = ["", DEPENDENCIES_HEADING, "", DEPENDENCIES_PREFACE]
    if goals:
        lines += ["", GOALS_HEADING, "", goals_block(snap, goals)]
    for number, block in blocks.items():
        lines += ["", epic_heading(snap.issue_by_number(number)), "", block]
    return lines


def graph_section(body):
    """Where an epic's body holds its graph: the span from the one
    GRAPH_MARKER through the GRAPH_END_MARKER after it; None where the
    body holds neither, or anything but one clean pair, a marker without
    the other, an end before the start or a second pair, which a person
    wrote and the board leaves alone, since writing into it could take
    out what they meant to keep."""
    start = body.find(GRAPH_MARKER)
    if start < 0 or body.find(GRAPH_END_MARKER) < start:
        return None
    end = body.find(GRAPH_END_MARKER, start + len(GRAPH_MARKER))
    if end < 0 or body.find(GRAPH_MARKER, start + len(GRAPH_MARKER)) >= 0:
        return None
    return (start, end + len(GRAPH_END_MARKER))


def line_ending(body):
    """The line ending a body is written with, Windows's where it holds
    one, so what is put into it is written the same way."""
    return "\r\n" if "\r\n" in body else "\n"


def graph_block(body):
    """The chart a body holds between its graph markers, folded to LF, or
    None where it holds no clean pair of them (graph_section)."""
    span = graph_section(body)
    if span is None:
        return None
    inside = body[span[0] + len(GRAPH_MARKER):span[1] - len(GRAPH_END_MARKER)].replace("\r\n", "\n")
    return inside.removeprefix("\n").removesuffix("\n")


def listed_graphs(issues):
    """The chart each listed issue's body holds between the graph markers,
    by number, as the API's issue objects carry the body (Snapshot.graphs);
    an issue with no clean pair of markers is left out."""
    found = {}
    for issue in issues:
        block = graph_block(issue.get("body") or "")
        if block is not None:
            found[int(issue["number"])] = block
    return found


def body_with_graph(body, block):
    """An epic's body with its graph: the block between the markers, put
    where the body holds them already or appended after a blank line, in
    the body's own line ending, so nothing outside the markers changes;
    None where the body holds markers that are not one clean pair, which
    a person wrote and the board leaves alone."""
    newline = line_ending(body)
    graph = newline.join((GRAPH_MARKER, block.replace("\n", newline), GRAPH_END_MARKER))
    span = graph_section(body)
    if span is None:
        if GRAPH_MARKER in body or GRAPH_END_MARKER in body:
            return None
        return body.rstrip("\r\n") + (newline * 2 if body.strip() else "") + graph + newline
    return body[:span[0]] + graph + body[span[1]:]


def body_without_graph(body):
    """An epic's body with its graph taken out, the blank line before it
    too, for an epic that owes no chart any more, its sub-issues taken
    off; None where the body holds no clean pair of markers."""
    span = graph_section(body)
    if span is None:
        return None
    newline = line_ending(body)
    before, after = body[:span[0]].rstrip("\r\n"), body[span[1]:].lstrip("\r\n")
    if before and after:
        return before + newline * 2 + after
    return before + newline if before else after


# --- the page ---------------------------------------------------------------

def render(snapshot, default_branch=None):
    """STATE.md for a snapshot: the layout the protocol's STATE.md describes.
    The default branch named here wins over the snapshot's own, for a
    workflow that knows it from the event; None takes the snapshot's."""
    snap = as_snapshot(snapshot)
    branch = cell(default_branch or snap.default_branch or "main")
    lines = [
        "# Board",
        "",
        f"As of commit {cell(short_sha(snap.sha))} on {branch}. "
        "Rewritten by .github/workflows/board.yml; do not edit by hand.",
    ]
    sections = (
        (CLAIM_HEADING, CLAIM_COLUMNS, claim_rows(snap)),
        (PULL_HEADING, PULL_COLUMNS, pull_rows(snap)),
        (OVERLAP_HEADING, OVERLAP_COLUMNS, overlap_rows(snap)),
        (QUEUE_HEADING, ISSUE_COLUMNS, issue_rows(snap, labels.TODO, True)),
        (GIVEN_HEADING, GIVEN_COLUMNS, given_rows(snap)),
        (BLOCKER_HEADING, ISSUE_COLUMNS, issue_rows(snap, labels.BLOCKER, False)),
        (MERGED_HEADING, MERGED_COLUMNS, merged_rows(snap)),
    )
    for heading, columns, rows in sections:
        lines += ["", heading, ""] + table(columns, rows)
    lines += dependencies_lines(snap)
    return "\n".join(lines) + "\n"


def write_if_changed(path, text):
    """Write the board only when its bytes differ, so a caller comparing
    files sees no change when nothing moved; True when it was written. The
    one thing in this module that touches a file, kept here because the
    rendered page is ASCII by construction and this is where that is known."""
    path = Path(path)
    data = text.encode("ascii")
    if path.is_file() and path.read_bytes() == data:
        return False
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(data)
    return True
