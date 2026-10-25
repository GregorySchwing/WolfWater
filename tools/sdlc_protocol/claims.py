"""The claim: a branch named for an issue and the file that describes it.

A claim is the first push of the branch claude/<issue number>, the number
as GitHub writes it, in ASCII digits with no leading zero. GitHub refuses
the second push of the same name, so the branch is the lock, and
its first commit adds docs/agents/claims/<issue number>.md, which names
the worker, the plan and the files the work expects to touch. The board
reads that file until the pull request has files of its own, and the
coordinator releases a claim whose branch has had no commit for a day.
A release is a comment on the issue that starts with RELEASE_MARKER;
once the board sees one, it renames the branch to
released/<issue>-<short sha>, which keeps the work and frees the name, so
the next claim is again a first push of claude/<issue number>.
This module holds those conventions and the parsing, and nothing that
talks to git or to GitHub.
"""

import re
import unicodedata
from datetime import datetime, timedelta, timezone

PREFIX = "claude/"
CLAIM_DIR = "docs/agents/claims"
# A claim with no commit for this long is released by the coordinator,
# and the board asks a pull request idle this long whether it is parked.
IDLE_AFTER = timedelta(hours=24)
# The first line of every release comment the coordinator posts, so the
# board and the field metrics tell a release from what people write
# without reading what the comment says.
RELEASE_MARKER = "<!-- sdlc-coordinator: claim-released -->"
# Whose release comments, and overlap markers, count by GitHub's
# author_association word alone: the repository's owner, the members of its
# organization and its collaborators. Anyone else's count only when an
# adapter finds that they may write to the repository (is_trusted_author).
RELEASERS = frozenset({"OWNER", "MEMBER", "COLLABORATOR"})
# The bot a repository's own workflows write as, with their GITHUB_TOKEN:
# the coordinator posts its releases and overlap comments as this bot.
# Only a workflow of the repository can write as it, and a workflow only
# runs once someone who may write has put it there, so what it posts
# counts like theirs. No person's login can end in "[bot]", which GitHub
# keeps for its Apps.
RELEASING_BOTS = frozenset({"github-actions[bot]"})
# How every GitHub App's bot login ends, and no person's can, since a login
# is letters, digits and hyphens. An App opens a pull request only in a
# repository it is installed on, which only someone who administers the
# repository or its organization does, so a pull request one opened is one
# the repository asked for, as the workers in Actions open theirs through
# the Claude App or the repository's own workflows (opened_by_writer).
APP_LOGIN_SUFFIX = "[bot]"
# Where a released claim's branch goes: out of the claude/ names, so the
# issue can be claimed again, and under a name of its own, so the work it
# holds stays for whoever carries it on.
RELEASED_PREFIX = "released/"
# What a worker that may not mark its pull request ready, a cloud session
# whose gh pr calls are refused, comments on its draft when the work is
# done, for a maintainer to mark: a draft whose author's latest comment
# starts with it is finished, and waits for its review, which the
# coordinator never cuts short with a release; a person who may write
# can still release it by hand, and that release counts.
FINISHED_SIGNAL = "finished; ready to mark"

# The line that names the claim's worker: "- Worker: name", or the older
# "- Account: name" that claim files carried before workers had one name.
WORKER_LINE = re.compile(r"^\W*(?:worker|account)\W*:\s*(.+?)\s*$", re.IGNORECASE)
ACCOUNT_LINE = WORKER_LINE
# A note in parentheses after the name, like "(label worker:nightly)",
# is for the reader, not the board.
TRAILING_NOTE = re.compile(r"\s*\(.*\)\s*$")
FILES_HEADING = re.compile(r"^#+\s*files to touch\s*$", re.IGNORECASE)
PLAN_HEADING = re.compile(r"^#+\s*plan\s*$", re.IGNORECASE)
# Where an issue lists the files its work will touch: the Task form's
# "Files likely touched", or "Files to touch" in a task the coordinator
# planned out of an epic.
ISSUE_FILES_HEADING = re.compile(r"^#+\s*(?:files likely touched|files to touch)\s*$", re.IGNORECASE)
NO_RESPONSE = "_No response_"
LIST_ITEM = re.compile(r"^\s*[-*]\s+(.+?)\s*$")
UNNAMED = "unnamed"
# An issue number as a claim branch's key must spell it: as branch_of
# writes it, in ASCII digits with no leading zero, one to nine of them.
# GitHub refuses only a second push of the same name, so claude/03 beside
# claude/3 would be a second claim on #3. str.isdigit would also take a
# superscript, which int() refuses, and str.isdecimal another script's
# digits, which int() reads as a number; a run of digits longer than int()
# reads would raise where an answer belongs. Nine digits is far more than
# any queue has issues, as the gate reads an issue too (gate.ISSUE_TEXT).
NUMBER = re.compile(r"0|[1-9][0-9]{0,8}")


def branch_of(number, prefix=PREFIX):
    """The claim branch for an issue: claude/<number>. The prefix is an
    argument so a gate policy that names another one uses the same rule."""
    return f"{prefix}{number}"


def is_claim_branch(branch, prefix=PREFIX):
    return branch.startswith(prefix)


def key_of(branch, prefix=PREFIX):
    """What follows claude/ in a claim branch: the issue number by protocol,
    kept as text because a hand-made branch may put something else there
    and the board still shows it. None for a branch that is no claim."""
    return branch[len(prefix):] if is_claim_branch(branch, prefix) else None


def number_in(key):
    """The issue number a claim branch's key spells, or None when it spells
    none (NUMBER): every place that reads a branch name as an issue number
    reads it here, so no branch name can raise there and only one name
    holds an issue."""
    text = str(key)
    return int(text) if NUMBER.fullmatch(text) else None


def number_of(branch, prefix=PREFIX):
    """The issue number a claim branch names, or None when the branch is no
    claim or what follows claude/ is no number as branch_of writes one
    (number_in): a number comes back only when branch_of(number) is the
    branch itself, so claude/03 names no issue, and a branch name never
    raises."""
    key = key_of(branch, prefix)
    return None if key is None else number_in(key)


def reads_as_claim(branch, prefix=PREFIX):
    """Whether a branch reads as a claim and claims nothing: digits after
    claude/, of any script, that are not an issue's number as branch_of
    writes it (number_in), as in claude/03, a name its author meant for a
    claim, which the protocol check and the push guard refuse."""
    key = key_of(branch, prefix)
    return bool(key) and key.isdigit() and number_in(key) is None


def meant_number(key):
    """The issue number a key of digits seems to mean, read in ASCII with
    no leading zero, as #3 for claude/03, or None when it has more digits
    than an issue's number has, or is not all digits."""
    if not key or not key.isdigit():
        return None
    text = "".join(str(unicodedata.digit(digit)) for digit in key).lstrip("0") or "0"
    return number_in(text)


def claim_path(number):
    """The claim file's path inside the repository, as text with slashes."""
    return f"{CLAIM_DIR}/{number}.md"


def by_number(value):
    """A sort key that puts numbers in order and anything else after them,
    so claude/9 comes before claude/100 and claude/x after both; a key that
    spells no number (number_in) sorts by its text, and never raises."""
    number = number_in(value)
    return (0, number, "") if number is not None else (1, 0, str(value))


def claim_worker(text):
    """The worker a claim file names, from the first line like "Worker:
    main", or the older "Account: main", in any of the Markdown spellings a
    hand might give it, without a note in parentheses after it; "unnamed"
    when no line does. A line whose name is empty once the note and the
    markup are gone, as "- Worker: " or "- Worker: (to be filled)", names
    nobody: until F18's review it named "", which the board showed as the
    worker, a release named as whose the claim was, and the clerk took for
    the holder, so nobody carried the claim."""
    for line in text.splitlines():
        match = WORKER_LINE.match(line)
        if match:
            name = TRAILING_NOTE.sub("", match.group(1)).strip("`*_ ")
            if name:
                return name
    return UNNAMED


# The name the board read the same line under while it named an account.
claim_account = claim_worker


def section_lines(text, heading):
    """The lines of a claim file under the heading the pattern matches, up
    to the next heading."""
    lines, inside = [], False
    for line in text.splitlines():
        if heading.match(line):
            inside = True
            continue
        if inside and line.startswith("#"):
            break
        if inside:
            lines.append(line)
    return lines


def claim_files(text):
    """The paths a claim file lists under its "Files to touch" heading, one
    list line each, up to the next heading, sorted: what the board shows
    for a claim whose pull request does not exist yet."""
    matches = (LIST_ITEM.match(line) for line in section_lines(text, FILES_HEADING))
    return sorted(match.group(1).strip("`") for match in matches if match)


def issue_files(body):
    """The paths an issue's body lists under "Files likely touched" or
    "Files to touch": the first block of lines under the heading, one path
    a line, without list marks or code spans, sorted and each once; what
    the clerk deals work to a worker's area by. The prose after a blank
    line is not a path, nor is a line in parentheses or GitHub's "_No
    response_"."""
    paths, started = set(), False
    for line in section_lines(body or "", ISSUE_FILES_HEADING):
        if not line.strip():
            if started:
                break
            continue
        started = True
        text = line.strip().lstrip("-*").strip().strip("`").strip()
        if text and text != NO_RESPONSE and not text.startswith("("):
            paths.add(text)
    return sorted(paths)


def claim_plan(text):
    """The plan a claim file states under its "Plan" heading, as one line
    of prose, or "" when it states none: what a release leaves behind for
    the next worker."""
    return " ".join(" ".join(section_lines(text or "", PLAN_HEADING)).split())


def parse_time(value):
    """An ISO 8601 time as GitHub or git writes it, as an aware datetime; a
    datetime given is returned as it is. A time without a zone is UTC."""
    if isinstance(value, datetime):
        moment = value
    else:
        text = value.strip()
        if text.endswith("Z"):
            text = text[:-1] + "+00:00"
        moment = datetime.fromisoformat(text)
    if moment.tzinfo is None:
        moment = moment.replace(tzinfo=timezone.utc)
    return moment


def is_idle(last_change, now, idle_after=IDLE_AFTER):
    """Whether more than idle_after has passed since the last change: the
    one rule behind the board's idle comment, the coordinator's release
    and the board's reading of a release (was_released)."""
    return now - parse_time(last_change) > idle_after


def release_due(claims, now, idle_after=IDLE_AFTER, released=None):
    """The claims the coordinator releases: those whose branch has had no
    commit for more than idle_after. The claims are a mapping, or pairs,
    from the claim's key (its branch, or the issue number) to the time of
    the branch's last commit, as a datetime or as GitHub writes it; the
    answer is the keys due, in issue order, so the comments the coordinator
    posts come out the same for the same state. released maps the same
    keys to the time of the latest release comment on each claim's issue;
    a claim that comment already released (was_released) is not due
    again, so a coordinator that runs every hour releases it once."""
    pairs = claims.items() if hasattr(claims, "items") else claims
    released = released or {}
    due = [key for key, last_commit in pairs
           if is_idle(last_commit, now, idle_after) and not was_released(released.get(key, ""), last_commit, idle_after)]
    return sorted(due, key=issue_order)


def issue_order(key):
    """A sort key for a claim given as its branch or as its issue number."""
    text = str(key)
    return by_number(key_of(text) if is_claim_branch(text) else text)


def starts_with_release_marker(comment):
    """Whether an issue comment, as the API gives it, has RELEASE_MARKER as
    its first line. A comment that quotes the marker further down, a reply
    for instance, talks about a release and is none."""
    lines = ((comment or {}).get("body") or "").strip().splitlines()
    return bool(lines) and lines[0].strip() == RELEASE_MARKER


def is_release_comment(comment, trusted=None):
    """Whether an issue comment releases a claim: its first line is the
    marker and the repository trusts its author (is_trusted_author). A
    comment by anyone else cannot take a claim away from its holder."""
    return starts_with_release_marker(comment) and is_trusted_author(comment, trusted)


def is_trusted_author(comment, trusted=None):
    """Whether the repository trusts the author of a comment, as the API
    gives it, to speak for the coordinator: by GitHub's author_association
    (RELEASERS), as the repository's own workflows (RELEASING_BOTS), or, for
    an author neither vouches for, by trusted, a callable given the author's
    login that an adapter answers from the author's permission on the
    repository: GitHub may show a private member of the organization to a
    workflow's token as a mere contributor."""
    if (comment or {}).get("author_association") in RELEASERS:
        return True
    login = ((comment or {}).get("user") or {}).get("login")
    if login in RELEASING_BOTS:
        return True
    return bool(trusted and login and trusted(login))


def opened_by_writer(pull, trusted=None):
    """Whether someone who may write to the repository opened a pull
    request, as the API gives it, so that it may be a claim's own
    (board.claim_of): by GitHub's author_association (RELEASERS), as a
    GitHub App (APP_LOGIN_SUFFIX), or, for an author neither vouches for, by
    trusted, a callable given the author's login that an adapter answers
    from the author's permission on the repository. GitHub's documentation
    says anyone who may read a repository may open a pull request in it,
    from one of its own branches as from a fork, so a stranger's pull
    request from a claim's branch is no claim's (F19). A pull request whose
    author the API does not name, as a fixture written by hand, counts as
    a writer's, as one does whose author's permission GitHub will not say
    (adapters.github.may_open), so a missing answer never frees a finished
    claim."""
    if (pull or {}).get("author_association") in RELEASERS:
        return True
    login = ((pull or {}).get("user") or {}).get("login")
    if not isinstance(login, str) or not login or login.endswith(APP_LOGIN_SUFFIX):
        return True
    return bool(trusted and trusted(login))


def said_finished(comments, author):
    """Whether the author's latest comment among a pull request's comments
    starts with FINISHED_SIGNAL, so the draft is finished, and the
    coordinator never releases its claim however long it waits; a
    person's release of it still counts."""
    theirs = [c for c in comments if (c.get("user") or {}).get("login") == author and author]
    if not theirs:
        return False
    latest = max(theirs, key=lambda c: parse_time(c.get("created_at") or "1970-01-01T00:00:00Z"))
    return (latest.get("body") or "").strip().lower().startswith(FINISHED_SIGNAL)


def latest_release(comments, trusted=None):
    """When the latest release comment among an issue's comments was made,
    as GitHub wrote the time, or "" when none of them releases a claim."""
    times = [c["created_at"] for c in comments if c.get("created_at") and is_release_comment(c, trusted)]
    return max(times, key=parse_time) if times else ""


def was_released(released_at, last_commit, idle_after=IDLE_AFTER):
    """Whether a release comment made at released_at released the claim
    whose branch was last committed to at last_commit: the claim had gone
    without a commit for more than idle_after when the comment came, as a
    release is due only then, and nothing has been committed since, which
    would mean that its holder came back or another account claimed the
    issue anew. A comment made sooner than that releases nothing, and with
    either time unknown no release can be told."""
    if not released_at or not last_commit:
        return False
    return is_idle(last_commit, parse_time(released_at), idle_after)


def released_branch(key, sha):
    """The name a released claim's branch takes: released/<issue>-<the
    first seven characters of its head commit>, so an issue released twice
    keeps both branches apart."""
    return f"{RELEASED_PREFIX}{key}-{sha[:7]}"
