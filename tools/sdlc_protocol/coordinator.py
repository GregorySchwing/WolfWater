"""The coordinator's rules: which claims to release and what to say when it
does, and what to tell open work that meets other work in the same files.

The coordinator is a step of the board's own run, on the repository's own
token, so it needs no Claude account and no routine: before the board is
rendered, it releases every claim the release rule finds due, and after
the release it tells each piece of open work that changes a file another
piece changes which one should land first. What needs judgement, splitting
a request into issues, planning an epic, settling a blocked task, tidying
the shared memory, is Claude's, in the repository's coordinator workflow
or a coordinator routine, and is not here.

A release is the protocol's (spec/PROTOCOL.md, Release): a claim whose
branch has gone more than a day without a commit, whose issue is open and
not labelled in-review, and whose pull request, if it has one, is still a
draft. The coordinator takes the label claimed and the holder's label off
the issue, and puts todo on an issue that carries none of todo, epic and
blocker, so it is listed again (labels_put_on); then it posts one comment
on it, whose first line is claims.RELEASE_MARKER, saying why, whose the
claim was, what its branch holds and what its plan was, so the next worker
can carry the work on or start over, and the board renames the branch in
the same run. The board counts a release only once both are done. A run
that fails at the labels puts back the labels it took off, the one whose
taking off GitHub answered with an error among them unless GitHub refused
it in place of taking it off (F27), and takes off the todo it put on, so
the issue is as it was and the next run starts the release over. A run
whose comment fails reads the issue's comments: a release comment there
makes the release whole; where they list none after GitHub refused the
comment, a 4xx other than a rate limit's on a try that followed none
refused for a rate limit, which it answers in place of writing anything,
the run puts the labels back as at the labels; and where
it cannot read them, whatever GitHub answered (F20), or they list none
after any other answer, a rate limit's refusal, on the last try or on one
the run made the comment again after, a server error or no answer among
them, none of which says GitHub wrote nothing,
and a read may list a comment GitHub wrote only later (F23), it cannot
tell whether a release comment stands, its own, an earlier run's or one
someone else wrote, so it leaves the labels as the release left them for
the next run to settle (adapters/github.py, release_one).

Before it releases, the coordinator frees each issue that was closed and
opened again while its labels still speak of a claim that ended with the
close: claimed, in-review and its worker's label would keep it out of the
Queue and from every worker, and in-review would keep a claim branch that
still stands from ever being released (reopened_due, freeing; F22). Like
the board, this module only decides and words; adapters/github.py makes
the calls.
"""

import dataclasses
import re
from typing import NamedTuple

from . import board, claims, labels

# The first line of the comment that tells a piece of open work about the
# other work it meets; the marker names every number it told of, so a later
# run tells again only when work it has not told of starts to meet it.
OVERLAP_MARKER = "<!-- sdlc-coordinator: overlap {numbers} -->"
OVERLAP_LINE = re.compile(r"<!-- sdlc-coordinator: overlap ([0-9 ]+) -->")
# An issue carrying one of these labels gets no todo from a release: todo
# and blocker each list it on the board with no claim, and the Queue never
# lists an epic.
NEEDS_NO_TODO = frozenset({labels.TODO, labels.EPIC, labels.BLOCKER})
# How many commit subjects and files a release note names before it counts
# the rest, so a long-lived branch does not bury the comment.
FILES_NAMED = 20
# How much of the claim's plan, and of each name (the worker, the newest
# commit's subject, a file), a release note quotes before it cuts the rest:
# GitHub refuses a comment of more than 65,536 characters, a refused comment
# releases nothing, and the claim file, the branch and its history keep the
# whole. With these the comment stays far under that limit, whatever the
# claim file and the branch hold.
PLAN_QUOTED = 4000
NAME_QUOTED = 200


# --- releases -----------------------------------------------------------------

def is_candidate(snap, key):
    """Whether a release can be due on the claim with this key: its branch
    exists, its issue is open and not labelled in-review, no pull request
    from its branch is marked ready, one from a fork whose branch bears the
    same name being no claim's (board.of_claim), and its head commit and
    when its last commit came (board.Release.last_commit) are known."""
    issue, release = snap.claim_issues.get(key), snap.releases.get(key)
    if issue is None or release is None or not release.sha or not release.last_commit:
        return False
    if claims.branch_of(key) not in snap.branches or issue.state != "open" or labels.IN_REVIEW in issue.labels:
        return False
    return not any(board.of_claim(p, key) and not p.draft for p in snap.pulls)


def releases_due(snapshot, now, idle_after=claims.IDLE_AFTER):
    """The keys of the claims the coordinator releases now, in issue order:
    claims.release_due over every candidate, so a claim already released
    (claims.was_released) is not released twice, and one released too soon
    is released again. A draft its author called finished is a candidate
    only once a person's release of it counts (adapters.github.releases),
    so it is never released here."""
    snap = board.as_snapshot(snapshot)
    keys = [key for key in snap.releases if is_candidate(snap, key) and not kept_released(snap, key)]
    last_commits = {key: snap.releases[key].last_commit for key in keys}
    released = {key: snap.releases[key].released_at for key in keys}
    return claims.release_due(last_commits, now, idle_after, released)


def kept_released(snap, key):
    """Whether the board counts released a claim whose head a release
    kept (board.is_released): such a branch was pushed again with nothing
    new, and the board deletes it, so a release of it would be a second
    one for the same work. With no release comment that counts, or with
    the label claimed put back, which undoes a release, it is released as
    any idle claim is, its push dating it, so a person's undoing never
    holds it for good."""
    release = snap.releases[key]
    return release.kept and board.is_released(snap, key)


class Holding(NamedTuple):
    """What a released branch holds, as a comparison with the default branch
    found it: how many commits it is ahead, the subject of its newest one
    and the files it changes; ahead is None when GitHub would not compare."""

    ahead: int | None = None
    last_subject: str = ""
    files: tuple = ()


def holder(snap, key):
    """Whose a claim was: the worker its claim file names, else the
    worker its issue's label names, else nobody known."""
    text = snap.claims.get(key)
    if text:
        named = claims.claim_worker(text)
        if named != claims.UNNAMED:
            return named
    issue = snap.claim_issues.get(key)
    return labels.worker_of({"labels": issue.labels} if issue else None) or "a worker nobody named"


def cut(text, limit):
    """The text, or its first limit characters and "..." when it is longer."""
    return text if len(text) <= limit else text[:limit].rstrip() + "..."


def named_files(files):
    """Up to FILES_NAMED files, in code spans, each cut to NAME_QUOTED
    characters, and how many more there are."""
    shown = ", ".join(f"`{cut(path, NAME_QUOTED)}`" for path in files[:FILES_NAMED])
    rest = len(files) - FILES_NAMED
    return shown + (f" and {rest} more" if rest > 0 else "")


def holding_text(holding, default_branch):
    """What the branch holds, in one sentence."""
    files = board.own_files(holding.files)
    if holding.ahead is None:
        held = "GitHub would not compare the branch with the default branch, so look at the branch itself"
    else:
        held = f"{holding.ahead} commit{'s' if holding.ahead != 1 else ''} ahead of {default_branch}"
        if holding.last_subject:
            held += f", the newest \"{cut(holding.last_subject, NAME_QUOTED)}\""
    if files:
        return f"{held}; it changes {named_files(files)}."
    if holding.ahead is not None:
        return f"{held}; it changes nothing but its claim file."
    return held + "."


def release_text(snapshot, key, holding, idle_after=claims.IDLE_AFTER):
    """The comment that releases the claim with this key: the marker as
    its first line, then why, whose, what the branch holds, what the plan
    was, and how to carry the work on or undo the release. What it quotes
    is cut to PLAN_QUOTED and NAME_QUOTED, so GitHub never refuses it for
    its length."""
    snap = board.as_snapshot(snapshot)
    release = snap.releases[key]
    old, new = claims.branch_of(key), claims.released_branch(key, release.sha)
    hours = int(idle_after.total_seconds() // 3600)
    pull = next((p for p in snap.pulls if board.of_claim(p, key)), None)
    pull_clause = f"its pull request #{pull.number} is still a draft" if pull else "it has no pull request"
    return "\n".join([
        claims.RELEASE_MARKER,
        f"released: the claim on this issue, held by {cut(holder(snap, key), NAME_QUOTED)}, has had no commit on "
        f"{old} since {release.last_commit}, more than {hours} hours, and {pull_clause}, so the coordinator "
        f"released it. The board renames the branch to {new}, which keeps the work, and the issue goes back to "
        f"the queue.",
        "",
        f"What the branch holds: {holding_text(holding, snap.default_branch)}",
        plan_line(snap.claims.get(key), key, new),
        "",
        f"To carry the work on, claim the issue anew, a new {old} from {snap.default_branch}, and bring the "
        f"commits in from {new} after the claim commit; or start over. To undo the release, put the label "
        f"claimed back on the issue.",
    ]) + "\n"


def plan_line(text, key, new):
    """The release note's line on the claim's plan, cut to PLAN_QUOTED with
    a word on where the rest is, the claim file the branch keeps; for a
    claim file GitHub gave no text of (board.UNREAD_CLAIM), where it is."""
    if text is board.UNREAD_CLAIM:
        return (f"The claim file is over 1 MB, more than GitHub gives the text of, so its plan is not quoted; "
                f"{claims.claim_path(key)} on {new} has it.")
    plan = claims.claim_plan(text or "")
    if not plan:
        return "The claim file states no plan."
    if len(plan) <= PLAN_QUOTED:
        return f"The claim's plan: {plan}"
    return f"The claim's plan, cut short ({claims.claim_path(key)} on {new} has it whole): {cut(plan, PLAN_QUOTED)}"


def labels_released(snapshot, key):
    """The labels a release takes off the issue: claimed and every label
    that names a worker, so the issue is nobody's, but a bare worker:,
    which names nobody (labels.is_holder_label); labels_put_on says what
    goes on, so it is back in the open queue."""
    issue = board.as_snapshot(snapshot).claim_issues.get(key)
    names = issue.labels if issue else []
    return [name for name in names if name == labels.CLAIMED or labels.is_holder_label(name)]


def labels_put_on(snapshot, key):
    """The labels a release puts on the issue: todo, when it carries none
    of todo, epic and blocker, since the Queue lists only an issue labelled
    todo and the spec has a released claim's issue listed there again,
    while a task a worker's label alone gave carries none (F15: until its
    fix such an issue was shown nowhere once its claim was released). An
    epic gets none, since the Queue never lists one and its tasks are the
    work; a blocker none, since the board lists every open blocker under
    Blockers and the clerk deals one first without todo, so todo would
    only list it twice; and nor does an issue the snapshot knows nothing
    of, which no release is due on (is_candidate)."""
    issue = board.as_snapshot(snapshot).claim_issues.get(key)
    if issue is None or NEEDS_NO_TODO & set(issue.labels):
        return []
    return [labels.TODO]


def after_release(snapshot, key, released_at, taken_off, put_on=()):
    """The snapshot as GitHub knows it once the claim with this key is
    released: the comment made at released_at, the labels taken off its
    issue and put on it in both places the board reads them, so the same
    run renders the release and renames the branch without asking GitHub
    again. A Snapshot given is changed in place, and returned."""
    snap = after_labels(snapshot, key, taken_off, put_on)
    snap.releases[key] = dataclasses.replace(snap.releases[key], released_at=released_at)
    return snap


def after_labels(snapshot, key, taken_off, put_on=()):
    """The snapshot as GitHub knows it once a release has taken taken_off
    off the issue with this key and put put_on on it, in both places the
    board reads them, with no new release comment known: after_release's
    labels, all of an unsettled release (adapters.github.labels_left),
    whose comment the run cannot tell stands, and what a run freeing a
    reopened issue changed (adapters.github.free_reopened), where the issue
    may have no claim branch and so no place among the claims' issues. The
    release comment the snapshot held stays, so a head a release already
    kept, whose latest release comment counts whenever it came
    (board.release_counts), is released by it once the labels are off. A
    Snapshot given is changed in place, and returned."""
    snap = board.as_snapshot(snapshot)
    issue = snap.claim_issues.get(key)
    if issue is not None:
        issue.labels = relabelled(issue.labels, taken_off, put_on)
    for listed in snap.issues:
        if str(listed.number) == key:
            listed.labels = relabelled(listed.labels, taken_off, put_on)
    return snap


def relabelled(names, taken_off, put_on):
    """An issue's labels once a release took taken_off off and put put_on
    on: those kept, then those put on that it did not carry."""
    kept = [name for name in names if name not in taken_off]
    return kept + [name for name in put_on if name not in kept]


# --- reopened issues --------------------------------------------------------

class Freeing(NamedTuple):
    """What the coordinator does to a reopened issue, in that order: the
    labels it puts on, then those it takes off; and whether a branch
    stands for the issue (claim_stands), where it takes in-review off
    alone."""

    put_on: tuple = ()
    taken_off: tuple = ()
    branch_stands: bool = False


def reopened_due(snapshot):
    """The open issues GitHub says were reopened whose labels still speak of
    a claim that ended with the close, in number order: those freeing
    takes a label off. The issues are the snapshot's, as the board read
    them; the adapter reads each again before it frees it, since a merge
    can close an issue after the board read it."""
    snap = board.as_snapshot(snapshot)
    return [issue.number for issue in sorted(snap.issues, key=lambda i: i.number)
            if issue.reopened and freeing(snap, issue.number, issue.labels).taken_off]


def freeing(snapshot, number, names):
    """What freeing the reopened issue with this number, which carries the
    labels names, does. An open pull request marked ready from its claim
    branch finishes the claim, and in-review may be its own, so nothing.
    Where a branch stands for the issue (claim_stands), in-review alone
    comes off: a claim branch is a claim whatever the labels say, and once
    in-review is off the release frees it when its branch has gone a day
    without a commit, in the same run where it has (adapters.github.
    may_be_due), while a released branch beside the label claimed is a
    release a person is undoing, unless that release came before the
    issue's last close and claimed came through it, or was put on by a
    claim made since (claim_stands). Where
    none stands and the issue carries claimed or in-review, the claim
    ended with the close: todo goes on as
    a release puts it on (labels_put_on), then every worker's label,
    in-review and claimed come off, claimed last, so a run that fails
    partway leaves the issue carrying claimed or in-review, which the next
    run frees."""
    snap, key = board.as_snapshot(snapshot), str(number)
    if any(board.of_claim(p, key) and not p.draft for p in snap.pulls):
        return Freeing()
    if claim_stands(snap, key):
        return Freeing(taken_off=(labels.IN_REVIEW,) if labels.IN_REVIEW in names else (), branch_stands=True)
    if labels.CLAIMED not in names and labels.IN_REVIEW not in names:
        return Freeing()
    put_on = () if NEEDS_NO_TODO & set(names) else (labels.TODO,)
    holders = tuple(name for name in names if labels.is_holder_label(name))
    return Freeing(put_on, holders + tuple(name for name in (labels.IN_REVIEW, labels.CLAIMED) if name in names))


def claim_stands(snap, key):
    """Whether a branch stands for the issue with this key: its claim
    branch, or a branch a release left (board.Snapshot.kept_names) where
    that release may be being undone. A release whose latest comment came
    no later than the issue's last close, with the label claimed put on no
    later than the close either, or put on by a claim made since, no later
    than the newest entry the activity lists for its branch since the close
    (board.Snapshot.ended_releases), belongs to a
    claim before the close, which a person undoing a release on the open
    issue cannot mean, since an undo puts claimed back after the close with
    no claim branch made, the claimed that came through it is the ended
    claim's and a claim made since puts its own: counting its
    branch would take in-review off alone and leave claimed and the
    worker's label on an issue with no claim branch for good, listed in
    no section of the board."""
    if claims.branch_of(key) in snap.branches:
        return True
    kept = claims.released_branch(key, "")
    return key not in snap.ended_releases and any(name.startswith(kept) for name in snap.kept_names)


# --- overlaps -----------------------------------------------------------------

def meeting(snapshot):
    """For each piece of open work that meets other work in a file, by its
    identity (board.identity: the issue for a claim's work, the pull
    request's number for one that is no claim's): the work it meets with
    the files they share, each file once however many of the piece's pull
    requests change it, in landing order."""
    found = {}
    for path, works in board.overlaps(snapshot).items():
        for work in works:
            for other in works:
                if board.identity(other) != board.identity(work):
                    found.setdefault(board.identity(work), {}).setdefault(other, set()).add(path)
    return {identity: sorted(((other, tuple(sorted(paths))) for other, paths in others.items()),
                             key=lambda pair: board.landing_key(pair[0]))
            for identity, others in sorted(found.items())}


def identities(snapshot):
    """What each number a piece of open work is shown by names now: the
    piece's identity (board.identity), so a marker that names a claim's
    pull request, as the markers written before F2's fix do, is read as
    naming the claim."""
    return {work.number: board.identity(work) for work in board.open_work(snapshot)}


def told_places(snapshot, identity):
    """The other places a piece of open work has been told on, beside its
    identity: the claim's own open pull requests, where the coordinator
    told a claim's work before F2's fix and whose markers count as the
    issue's (spec, Overlaps), so the upgrade tells nothing again; none for
    a pull request that is no claim's, which is told on itself."""
    return [work.number for work in board.open_work(snapshot) if work.issue == identity and work.is_pull]


def claims_made(snapshot):
    """When the push that made each claim's branch came, by the claim's
    issue, for the claims among the open work GitHub dates (snapshot
    began): an issue's number names each claim on it in turn, so a marker
    tells a claim's work, or tells of it, only when it came no earlier than
    that push; a claim GitHub does not date takes every marker."""
    snap = board.as_snapshot(snapshot)
    made = {}
    for work in board.open_work(snap):
        if work.issue is None or not work.branch:
            continue
        when = board.began_at(snap.began.get(work.branch, ""))
        if when != board.UNDATED:
            made[work.issue] = when
    return made


def told_numbers(comments, trusted=None):
    """Every number an overlap comment among these has told of, read from
    the marker on its first line: a comment that quotes a marker further
    down talks about an overlap and told nothing. A marker counts only from
    an author the repository trusts (claims.is_trusted_author; trusted as
    there), as the coordinator writes as the repository's own workflows, so
    nobody else can keep a piece of work from being told, and only when
    every number in it is an issue's or a pull request's (claims.number_in):
    anyone may comment, and a number of thousands of digits would raise,
    stopping every board run."""
    told = set()
    for comment in comments:
        lines = ((comment or {}).get("body") or "").strip().splitlines()
        match = OVERLAP_LINE.fullmatch(lines[0].strip()) if lines else None
        numbers = [claims.number_in(number) for number in match.group(1).split()] if match else []
        if numbers and None not in numbers and claims.is_trusted_author(comment, trusted):
            told.update(numbers)
    return told


def told_times(comments, trusted=None):
    """For each number the overlap markers among these name, when the latest
    of them came, as GitHub writes the time: each comment is read alone by
    told_numbers, so a marker counts by the same rules, and its time is
    the comment's created_at, "" when it has none."""
    times = {}
    for comment in comments:
        for number in told_numbers([comment], trusted):
            times[number] = latest(times.get(number), (comment or {}).get("created_at") or "")
    return times


def latest(first, second):
    """The later of two times GitHub wrote, a time that does not read as one
    counting as the latest of all (board.began_at), so a marker GitHub
    dates oddly still counts, and None for no time yet."""
    if first is None:
        return second
    return first if board.began_at(first) >= board.began_at(second) else second


def merged_times(first, second):
    """Two readings of told_times as one, each number at its latest marker."""
    times = dict(first)
    for number, at in second.items():
        times[number] = latest(times.get(number), at)
    return times


def heard_of(markers, named, made, identity):
    """The identities a piece has been told of, from the markers on it:
    markers maps each number they name to when the latest of them came
    (told_times), or is a set of numbers, which carries no times; each
    number is read as the identity of the work shown by it now (named),
    and counts only when its marker came no earlier than the push that
    made the claim branch of the piece and of the work named, where GitHub
    dates that push (claims_made): a claim released and claimed anew is
    new work, told and told of again."""
    times = markers if isinstance(markers, dict) else dict.fromkeys(markers, "")
    heard = set()
    for number, at in times.items():
        other = named.get(number, number)
        when = board.began_at(at)
        if all(when >= made[i] for i in (identity, other) if i in made):
            heard.add(other)
    return heard


def describe(work):
    """A piece of work in a few words, for an overlap comment."""
    if not work.is_pull:
        kind = "a claim with no pull request yet"
    elif work.ready:
        kind = "a pull request marked ready"
    else:
        kind = "a draft pull request"
    return kind + (", for a blocker" if work.blocker else "")


def this_work(works):
    """What the comment calls the piece it is posted on: the claim, its pull
    request or its pull requests, or the pull request that is no claim's."""
    pulls = [f"#{work.number}" for work in works if work.is_pull]
    if not pulls:
        return "claim"
    if works[0].issue is None:
        return "pull request"
    if len(pulls) == 1:
        return f"claim's pull request {pulls[0]}"
    return f"claim's pull requests {', '.join(pulls[:-1])} and {pulls[-1]}"


def overlap_text(works, others):
    """The comment that tells a piece of open work, shown by works, about the
    other work it meets: each other piece with the files they share, and
    the order the whole group should land in."""
    order = sorted(list(works) + [other for other, _ in others], key=board.landing_key)
    lines = [f"overlap: other open work changes the same files as this {this_work(works)}:"]
    lines += [f"- #{other.number}, {describe(other)}: {named_files(list(paths))}" for other, paths in others]
    lines += ["", f"Land them in this order: {', '.join(f'#{w.number}' for w in order)}: a blocker first, then "
                  f"work marked ready, then the oldest. Each one after the first brings the default branch in "
                  f"once the one before it has landed, and the board's Overlaps section shows the order as it "
                  f"changes."]
    return "\n".join(lines) + "\n"


def overlap_comments_due(snapshot, told):
    """The overlap comments to post, as board.DueComment, each on the
    piece's identity, the issue for a claim's work (board.identity): one
    for each piece of open work that meets work it has not been told of.
    told maps a piece's identity to what the overlap markers on it name,
    read by heard_of: each number as the identity of the piece it names
    now (identities), since a marker written before F2's fix named a
    claim's pull request, and only from a marker that came no earlier than
    the claims it speaks of (claims_made). The marker names the identities
    of the work met, each once, and the comment names all of it, told
    before or not, so it reads whole on its own."""
    snap = board.as_snapshot(snapshot)
    works = {}
    for work in board.open_work(snap):
        works.setdefault(board.identity(work), []).append(work)
    named, made = identities(snap), claims_made(snap)
    due = []
    for identity, others in meeting(snap).items():
        met = sorted({board.identity(other) for other, _ in others})
        if set(met) <= heard_of(told.get(identity, set()), named, made, identity):
            continue
        marker = OVERLAP_MARKER.format(numbers=" ".join(str(n) for n in met))
        due.append(board.DueComment(identity, marker, overlap_text(works[identity], others)))
    return due
