"""The clerk: the clerical work of the coordinator and the workers, done
by code so no Claude usage goes to it.

The coordinator plans an epic the way CooperBench's best-scoring
configuration does: Claude turns it into a plan and tasks given to
workers, the workers carry them in parallel on claim branches, and the
coordinator integrates at the end. A worker is any worker of the
registry (workers.py): a person's sessions, a routine or a worker in
GitHub Actions. Claude does the thinking; this module does the rest, three
ways, each run by a workflow step through the command line (cli.py's
"team" subcommand):

  seed   read the coordinator's plan (the JSON claude-code-action returns
         under --json-schema), write it into the epic, open one sub-issue
         per task labelled todo and worker:<owner>, or todo alone where no
         worker is registered, and link each under the epic.
  route  after an issue event, say what the coordinator does next:
         "plan <epic>", "integrate <epic>", "unblock <task>" or "none".
  next   for one worker, say what it does next: "carry <issue>" (its own
         claim, one whose claim file names it or, where no file names a
         worker, whose issue's label names it, needs work), "take <issue>"
         (the task it was given, else its share of the unclaimed blockers
         and todos that fit its area) or "none".
  deal   for every worker of a run, from one look at the repository, the
         same answer next gives each, as JSON by worker: the workers'
         workflow deals once, in its first job, and hands each leg its
         move, so no two legs looking for themselves at different moments
         are dealt one issue (F4).

Every call to GitHub goes through a Client, the small interface below that
adapters/github.py implements over REST and a test implements in memory,
so the rules here can be pinned without a network. The answers are the
one-line strings the workflows read; what the clerk cannot do is a
TeamError, whose message never holds a token.
"""

import json
import re
from typing import Protocol, runtime_checkable

from .claims import PREFIX as CLAIM_PREFIX, UNNAMED, claim_worker, issue_files, number_of
# The label names and the prefixes are spelled once, in labels.py: a
# worker names itself with worker:<name>, and the older account:<name> and
# agent:<name> are read as the same label.
from .labels import (ACCOUNT_PREFIX, AGENT_PREFIX, BLOCKED as BLOCKED_LABEL,
                     BLOCKER as BLOCKER_LABEL, CLAIMED as CLAIMED_LABEL, EPIC as EPIC_LABEL,
                     IN_REVIEW as IN_REVIEW_LABEL, TODO as TODO_LABEL, WORKER_COLOUR, agent_label as agent_of,
                     blocked_by_counts, holders, label_names, planned, priority, worker_label, worker_of)
from .workers import fits

# The one check the ruleset requires and a member can fix; the board's own
# job and Claude's runs also leave check runs, which no member can mend.
REQUIRED_CHECK = "check"
PLAN_MARKER = "<!-- team:plan -->"
TASKS_LINE = re.compile(r"<!-- team:tasks ([0-9 ]*) -->")
EPIC_LINE = re.compile(r"^Epic: #(\d+)\s*$", re.MULTILINE)
# A ready pull request with one of these states has fallen behind its
# base or conflicts with it; any other state waits on checks or a merge.
NEEDS_REBASE = {"behind", "dirty"}
# The same, read where GitHub's word says only that something else stops
# the merge: the standing the adapter reads through mergeable and a
# comparison (adapters.github.standing_of), as the board does (F9).
STANDING_NEEDS_REBASE = {"behind", "conflicting"}
FAILED = {"failure", "timed_out"}
# The colour a worker's label is made in; kept under team mode's name too.
AGENT_COLOUR = WORKER_COLOUR


class TeamError(Exception):
    """What the clerk could not do; the message never holds the token."""


class ClientError(Exception):
    """A call GitHub refused or did not answer. An adapter raises this, or a
    subclass of it, so the clerk can tell the refusals it tolerates (a
    sub-issue link that failed) from a bug in its own code; the message
    never holds the token."""


@runtime_checkable
class Client(Protocol):
    """The few calls the clerk makes, as one object handed in by the caller.

    Issues, pull requests and check runs are the dictionaries GitHub's REST
    API returns, so the clerk reads the same fields whoever made them: an
    issue has "number", "id", "title", "body", "state" and "labels" (each a
    dictionary with a "name"); a pull request has "draft", "mergeable_state"
    and "head" with its "sha"; a check run has "name" and "conclusion".
    A call GitHub refuses raises ClientError.
    """

    def issue(self, number: int) -> dict:
        """One issue, open or closed, by number."""

    def open_issues(self) -> list[dict]:
        """Every open issue of the repository."""

    def branch_names(self) -> list[str]:
        """The names of the repository's branches."""

    def open_pull(self, branch: str) -> dict | None:
        """The open pull request from the branch, in the detail that holds
        "mergeable_state" and the head's "sha", or None when there is none;
        where that word says only that something else stops the merge, a
        "standing" may say how the branch stands to its base."""

    def check_runs(self, sha: str) -> list[dict]:
        """The check runs on one commit."""

    def claim_text(self, branch: str) -> str | bool | None:
        """The claim file on a claude/<issue> branch, at its head, or None
        when there is none; a file GitHub gives no text of, one over 1 MB,
        is False (board.UNREAD_CLAIM), never an empty claim file."""

    def released(self, branch: str) -> bool:
        """Whether a release comment counts for the claim on a claude/<issue>
        branch, by its issue's comments and the branch's head
        (board.release_counts): the latest by an author the repository
        trusts came after the branch had gone more than a day without a
        commit, as GitHub dates the push of its head (F8), and nothing was
        committed since, or the head is one a release already kept. The
        clerk asks it only of a claim whose labels and pull request leave a
        release possible (Snapshot.released)."""

    def create_label(self, name: str, colour: str) -> None:
        """Make the label; one that exists already is left as it is, since
        GitHub answers a second creation with 422 and that is not a failure."""

    def create_issue(self, title: str, body: str, labels: list[str]) -> dict:
        """Open an issue and return it, with its "number" and "id"."""

    def add_sub_issue(self, parent: int, child_id: int) -> None:
        """Hang the issue with that id under the parent issue's number."""

    def add_blocked_by(self, number: int, blocker_id: int) -> None:
        """Record that the issue with that id blocks the issue with this
        number, as GitHub's blocked-by relation, the one an issue's page
        shows under Relationships."""

    def set_body(self, number: int, body: str) -> None:
        """Replace an issue's body."""

    def remove_label(self, number: int, name: str) -> None:
        """Take a label off an issue; one the issue does not carry is fine."""


# --- label rules --------------------------------------------------------------
# label_names, agent_of (labels.agent_label) and holders come from labels.py
# and are re-exported here, so a reader of the clerk finds them by the names
# the clerk uses them under.

def is_mine(issue, agent):
    """An issue the labels give to this worker alone, so no peer's move can
    take it from the worker. One that several workers' labels give is
    each of theirs (given_to), and whose a claim is, its claim file says,
    whatever the labels say (Snapshot.holder; F18)."""
    return holders(issue) == {agent}


def given_to(issue, agent):
    """Whether the issue's labels give it to this worker, alone or beside
    others: the board shows it to each of them under Given to workers, and
    Pick has each take it, the claim deciding whose it is (F18)."""
    return agent in holders(issue)


def epic_of(issue):
    """The epic a task belongs to, from the "Epic: #N" line seed writes."""
    found = EPIC_LINE.search(issue.get("body") or "")
    return int(found.group(1)) if found else None


def task_numbers(epic):
    """The task issues seed opened for an epic, from its machine line."""
    found = TASKS_LINE.search(epic.get("body") or "")
    return [int(n) for n in found.group(1).split()] if found else []


def open_issues(client):
    """Every open issue (pull requests, which GitHub's issue list mixes in
    and an adapter may pass on, left out)."""
    return [issue for issue in client.open_issues() if "pull_request" not in issue]


def claim_branches(client):
    """The issue numbers that have a claude/<number> branch on the remote."""
    numbers = (number_of(name) for name in client.branch_names())
    return {number for number in numbers if number is not None}


# --- route ------------------------------------------------------------------

def route(client, event, number, label=None):
    """What the lead does after an issue event: plan a newly labelled epic,
    help a member that labelled its task blocked, integrate an epic whose
    last task just closed, or nothing. A dispatch with the epic's own
    number and the event closed integrates it too, which is the way to
    recover an integration whose event was lost."""
    issue = client.issue(number)
    is_open = issue.get("state") == "open"
    if event == "labeled" and label == EPIC_LABEL:
        planned = PLAN_MARKER in (issue.get("body") or "")
        return f"plan {number}" if is_open and not planned else "none"
    if event == "labeled" and label == BLOCKED_LABEL:
        waiting = is_open and BLOCKED_LABEL in label_names(issue) and epic_of(issue)
        return f"unblock {number}" if waiting else "none"
    if event == "closed":
        return integration_due(client, epic_of(issue) or (number if task_numbers(issue) else None))
    return "none"


def integration_due(client, epic_number):
    """"integrate <epic>" when every task of an open epic is closed and the
    lead has not yet claimed the integration, or its claim needs work: two
    tasks closing together must not start two integrations."""
    if epic_number is None:
        return "none"
    epic = client.issue(epic_number)
    tasks = task_numbers(epic)
    if epic.get("state") != "open" or not tasks:
        return "none"
    if any(client.issue(n).get("state") != "closed" for n in tasks):
        return "none"
    if epic_number in claim_branches(client) and not claim_needs_work(client, epic_number):
        return "none"
    return f"integrate {epic_number}"


# --- next -------------------------------------------------------------------

def claim_needs_work(client, number):
    """Whether an open claim has something for its holder to do: no pull
    request yet, a draft, a base it has fallen behind or conflicts with, or
    a failed required check. One that is ready and only waiting for its
    checks or its merge needs nothing, and a run spent on it would be usage
    wasted. A branch behind its base counts whether GitHub says so or the
    pull request's standing does, where GitHub's word says only that a
    check or a review holds it, since the board tells such a pull request
    "main moved" (F9) and its holder should hear it from the clerk too. So
    a ready claim whose check is still running is carried at once when its
    branch is behind, where until F9 it waited for the check and, under a
    rule that asks for an up-to-date branch, for GitHub to say behind; a
    busy default branch can then take a worker in Actions round after
    round, each push starting the check again, up to its limit of rounds.
    integration_due reads the integration's claim the same way."""
    return needs_work(client, client.open_pull(f"{CLAIM_PREFIX}{number}"))


def needs_work(client, pull):
    """claim_needs_work over the claim's open pull request already read, None
    for a claim with none: the Snapshot reads each claim's pull request once
    a look and shares it between this rule and Snapshot.released, where
    until F18's second review a look read it twice."""
    if not pull:
        return True
    if pull.get("draft") or pull.get("mergeable_state") in NEEDS_REBASE:
        return True
    if pull.get("standing") in STANDING_NEEDS_REBASE:
        return True
    runs = client.check_runs(pull["head"]["sha"]) or []
    return any(run.get("name") == REQUIRED_CHECK and run.get("conclusion") in FAILED
               for run in runs)


def in_pick_order(issue):
    """The key that sorts the issues of one of Pick's places into the order
    a worker takes them in: by priority label (labels.priority), then by
    number."""
    return (priority(label_names(issue)), issue["number"])


class Snapshot:
    """One look at the open issues and the claim branches, shared by every
    member's answer so the members of one run divide the work alike."""

    def __init__(self, client):
        self.client = client
        # In Pick's order within each place: p1 first, then p2, p3 and the
        # issues with no priority, the lower number first among equals; a
        # repository whose issues carry no priority label sees number order.
        self.issues = sorted(open_issues(client), key=in_pick_order)
        self.claimed = claim_branches(client)
        # What the look has read so far, by issue number, each read once
        # the first time a member's answer asks for it.
        self.needs_work, self.pulls, self.names, self.releases = {}, {}, {}, {}

    def open_pull(self, number):
        """The claim's open pull request, or None, read once a look and
        shared by claim_needs_work and released, which until F18's second
        review each read it for themselves, so one look paid twice for one
        pull request."""
        if number not in self.pulls:
            self.pulls[number] = self.client.open_pull(f"{CLAIM_PREFIX}{number}")
        return self.pulls[number]

    def claim_needs_work(self, number):
        if number not in self.needs_work:
            self.needs_work[number] = needs_work(self.client, self.open_pull(number))
        return self.needs_work[number]

    def holds(self, issue, agent):
        """Whether the member holds the claim on this issue: the claim names
        the member (named) and does not count as released (released), which
        is nobody's until the board renames its branch, as the board shows
        its issue unclaimed again, though its claim file still names its
        holder, so the clerk sends nobody back to a claim a release freed
        (F18's review). The name is read first, so a peer's claim costs the
        member's look one read of its claim file and nothing more, and
        whether a claim is released, which costs its pull request and what
        Client.released costs, is asked only of a claim that names the
        member (F18's second review). Until F18's fix the clerk read the
        labels alone, so a worker whose label an issue still carried was
        told to carry the claim another worker had pushed and not yet
        announced, and the claim's holder was told nothing."""
        return self.named(issue) == agent and not self.released(issue)

    def holder(self, issue):
        """Whose the claim on this issue is (named), or None for a claim
        that counts as released (released) or that names nobody; holds is
        the question a member's answer asks, which costs a peer's claim
        less."""
        named = self.named(issue)
        return None if named is None or self.released(issue) else named

    def named(self, issue):
        """The worker a claim names, as a release names whose it was
        (coordinator.holder): the one its claim file names, else, where the
        finishing pull request has deleted the file, no line names one or
        GitHub gave no text of it, the one the issue's label names; None
        when neither does. Each claim's file is read once a look."""
        number = issue["number"]
        if number not in self.names:
            text = self.client.claim_text(f"{CLAIM_PREFIX}{number}")
            named = claim_worker(text) if text else UNNAMED
            self.names[number] = named if named != UNNAMED else worker_of(issue)
        return self.names[number]

    def released(self, issue):
        """Whether the claim on this issue counts as released, by the
        board's rule (board.released_keys): its issue is labelled neither
        claimed nor in-review, its own pull request, if any, is still a
        draft, and a release comment that counts came after its branch had
        gone more than a day without a commit, with none since, or its head
        is one a release already kept (Client.released). A claim whose
        issue is labelled claimed or in-review costs no call, and any other
        one its pull request, read once a look (open_pull), and what
        Client.released costs; the answer is kept for the look."""
        number = issue["number"]
        if number not in self.releases:
            self.releases[number] = self.is_released(issue)
        return self.releases[number]

    def is_released(self, issue):
        names = label_names(issue)
        if CLAIMED_LABEL in names or IN_REVIEW_LABEL in names:
            return False
        pull = self.open_pull(issue["number"])
        if pull and not pull.get("draft"):
            return False
        return self.client.released(f"{CLAIM_PREFIX}{issue['number']}")

    def free(self, issue):
        """Whether an issue may be taken: no claim branch holds it, and it
        is labelled neither blocked nor claimed. The board shows an issue
        labelled claimed in neither its Queue nor Given to workers, branch
        or no branch, as when a release is undone before its branch is
        renamed back, so a worker reading the board never takes it, and the
        clerk never deals it: a worker in Actions that leaves such an issue
        alone would otherwise be dealt it again on every run, its share of
        the queue held back behind it. An issue reopened after its merged
        branch was deleted is such an issue only until the coordinator's
        next run frees it (coordinator.freeing; F22). Nor is an issue free
        that GitHub's dependencies say an open issue blocks, the relation
        the coordinator writes from a plan's depends_on and a person draws
        on the issue's page: it waits until every open issue blocking it
        is closed, and the summary GitHub carries on every listed issue
        says how many are open (labels.blocked_by_counts), so the clerk
        reads no more for it."""
        labels = label_names(issue)
        return (issue["number"] not in self.claimed and BLOCKED_LABEL not in labels
                and CLAIMED_LABEL not in labels and blocked_by_counts(issue)[0] == 0)

    def carry(self, agent):
        """"carry <issue>" for the first claim of the member's own that needs
        work, one the member holds (holds), whatever the issue's labels say,
        and that is not labelled blocked, whose holder waits on a peer or the
        coordinator and so has nothing else to do (spec/PROTOCOL.md,
        Planning); None when it has none."""
        for issue in self.issues:
            number = issue["number"]
            if (number in self.claimed and BLOCKED_LABEL not in label_names(issue)
                    and self.holds(issue, agent) and self.claim_needs_work(number)):
                return f"carry {number}"
        return None

    def given(self, agent, taken=(), alone=False):
        """The first task the labels give the member, beside other workers
        or, with alone, to the member alone (is_mine), that is free and no
        epic, which the coordinator plans into tasks, nor goal, which a
        person breaks into epics, since nobody takes either whole, and that
        no member before it in the team took (taken), whether or not it also
        carries todo, since the board shows it under Given to workers either
        way (F14) and the spec has the worker take it first; None when there
        is none."""
        mine = is_mine if alone else given_to
        for issue in self.issues:
            if (issue["number"] not in taken and mine(issue, agent) and self.free(issue)
                    and not planned(label_names(issue))):
                return issue
        return None

    def own_move(self, agent):
        """A member's move that no peer's can change, which next_for answers
        from the member's own look: a claim of its own that needs work, else
        a task the labels give it alone. A task several workers' labels give
        goes to the first of them in the team's order with nothing else to
        do (moves_from), so it is no member's own move. None otherwise."""
        carried = self.carry(agent)
        if carried:
            return carried
        issue = self.given(agent, alone=True)
        return f"take {issue['number']}" if issue is not None else None

    def pool(self):
        """The shared work: unclaimed issues labelled blocker or todo that
        no worker was given and that are neither epics, which the
        coordinator plans, nor goals, which a person breaks into epics,
        blockers first, since a blocker stops everyone (spec/PROTOCOL.md),
        then each in Pick's order (in_pick_order)."""
        def shared(issue):
            labels = label_names(issue)
            return (self.free(issue) and not holders(issue) and not planned(labels)
                    and (TODO_LABEL in labels or BLOCKER_LABEL in labels))
        return sorted((i for i in self.issues if shared(i)),
                      key=lambda i: (BLOCKER_LABEL not in label_names(i), in_pick_order(i)))


def deal(members, pool, areas=None):
    """The shared work dealt out, by member: each member in turn takes the
    first issue no member before it took whose files fit its area
    (workers.fits; an issue that names no files fits every area), and a
    member left with none is not in the answer. With no areas this is the
    first issue to the first member, the second to the second, and so on."""
    areas = areas or {}
    dealt, taken = {}, set()
    for member in members:
        for issue in pool:
            if issue["number"] not in taken and fits(areas.get(member, ()), issue_files(issue.get("body"))):
                dealt[member] = issue
                taken.add(issue["number"])
                break
    return dealt


def deals(client, team, areas=None):
    """Every worker's next move, by name, from one look at the repository:
    its own claims first, the ones whose claim files name it or, where no
    file names a worker, whose issues' labels name it, then the task it was
    given, one given to it alone before one that several workers'
    labels give, which goes to the first of them in the order given with
    nothing else to do, then its share of the shared work, which is dealt
    out (deal) among the workers with no work of their own, in the order
    given, each taking the next issue that fits its area; a worker with no
    issue left for it gets "none" and waits, rather than start a session
    that would race a peer for the same claim. areas maps a worker's name
    to its area; a worker not in it takes any issue.

    One look for the whole run is the point: a worker that looked for
    itself could see a peer's ready pull request clean, and be dealt the
    issue after the peer's, while a later look saw the same pull request
    behind the default branch, the peer carrying its own claim, and dealt
    the peer's issue again (F4). So the workers' workflow asks once, in its
    first job, and hands each leg its move; next_for is one worker's answer
    from the same rule."""
    return moves_from(Snapshot(client), team, areas)


def moves_from(snap, team, areas=None):
    """deals over a look already taken: each member's own move, its claim
    that needs work or else a task the labels give it alone; then the tasks
    several workers' labels give, each member left with nothing, in the
    team's order, taking the first given to it that no member before it
    took, so such a task goes to the first of them with nothing else to do
    (F18: until its fix it went to none of them); then the shared work
    dealt among the members with none. The team is every worker the deal
    is for, whether or not a run starts it (deal_moves), so a worker
    outside the run keeps its share and two runs dealt from one state
    never give one issue to two workers. A name repeated in team is one
    member, dealt once, so no issue goes to a repetition and to nobody."""
    team = list(dict.fromkeys(team))
    moves = {member: snap.own_move(member) for member in team}
    taken = set()
    for member in team:
        issue = None if moves[member] else snap.given(member, taken)
        if issue is not None:
            moves[member] = f"take {issue['number']}"
            taken.add(issue["number"])
    idle = [member for member in team if not moves[member]]
    dealt = deal(idle, snap.pool(), areas)
    for member in idle:
        moves[member] = f"take {dealt[member]['number']}" if member in dealt else "none"
    return moves


def next_for(client, agent, team=None, areas=None):
    """The worker's next move, as deals gives it: the same answer the run's
    one look would give this worker, from the worker's own look. A move no
    peer's can change (Snapshot.own_move) is answered first, before any
    peer's pull request is read: deals gives a member such a move whatever
    the others do; a peer's claim costs the look one read of its claim
    file (Snapshot.holds), and its pull request, its head and its comments
    cost nothing, since whether a claim is released is asked only of one
    that names the worker."""
    team = list(team or [agent])
    if agent not in team:
        team.append(agent)
    snap = Snapshot(client)
    return snap.own_move(agent) or moves_from(snap, team, areas)[agent]


def next_move(client, agent, team=None, label=None, areas=None):
    """The worker's next move as a round of the workers' workflow asks for
    it: next_for, and when the move is a take and the worker named its
    label, the label is made first, so the session that follows can put it
    on the issue (a label gh cannot find makes `gh issue edit --add-label`
    add none)."""
    answer = next_for(client, agent, team, areas)
    if label and answer.startswith("take "):
        ensure_label(client, label)
    return answer


def deal_moves(client, team, areas=None, labels=None, run=None):
    """The moves of the workers of one run, as the workers' workflow asks
    for them in its first job: deals over the whole team, since a worker not
    in this run is still dealt its share and so left it, as next_for leaves
    it; the answer holds the workers of the run (run, default the whole
    team; a name not in the team is dealt after it), and each take's worker
    has its label made first, as next_move does, labels giving the label by
    name."""
    team = list(dict.fromkeys(team))
    run = list(dict.fromkeys(team if run is None else run))
    moves = deals(client, team + [name for name in run if name not in team], areas)
    for name in run:
        if moves[name].startswith("take ") and (labels or {}).get(name):
            ensure_label(client, labels[name])
    return {name: moves[name] for name in run}


# --- seed -------------------------------------------------------------------

def plan_schema(owners):
    """The JSON schema the coordinator's plan answers to, as Claude Code's
    --json-schema takes it, and as seed reads the plan: one or two
    paragraphs, then the tasks, each with a title, the worker it goes to,
    the files it owns, its body and, where one task must close before
    another can be taken, the positions of the tasks it depends on, which
    seed writes as GitHub's blocked-by relations. With workers registered,
    a task's owner is one of them; with none, any text, which seed leaves
    out, so every task goes to the open queue."""
    owner = {"type": "string", "enum": list(owners)} if owners else {"type": "string"}
    # No minimum and no uniqueItems: the API's structured output takes
    # neither keyword for sure, and dependency_problems and positions hold
    # the plan to both whatever wrote it.
    depends_on = {"type": "array", "items": {"type": "integer"},
                  "description": ("positions of the tasks of this plan that must close before this one "
                                  "can be taken, counted from 1 in the order the tasks are listed")}
    task = {"type": "object",
            "properties": {"title": {"type": "string"}, "owner": owner,
                           "files": {"type": "array", "items": {"type": "string"}}, "body": {"type": "string"},
                           "depends_on": depends_on},
            "required": ["title", "owner", "files", "body"]}
    return {"type": "object",
            "properties": {"plan": {"type": "string"}, "tasks": {"type": "array", "minItems": 1, "items": task}},
            "required": ["plan", "tasks"]}


def owners_for(tasks, agents):
    """Each task's owner: the one the plan named when it is a worker here;
    a task whose owner is not goes first to a worker no task names yet, so
    no worker sits idle while another holds two, then to the workers in
    turn; with no workers registered, no task has an owner."""
    if not agents:
        return [None] * len(tasks)
    named = [task.get("owner") if task.get("owner") in agents else None for task in tasks]
    spare = [agent for agent in agents if agent not in named]
    owners, turn = [], 0
    for owner in named:
        if owner is None and spare:
            owner = spare.pop(0)
        elif owner is None:
            owner, turn = agents[turn % len(agents)], turn + 1
        owners.append(owner)
    return owners


def depends_on(task):
    """What a task's depends_on lists, as the plan wrote it, or none; a
    plan written by hand may list anything, so dependency_problems reads
    this and the rest read positions."""
    return task.get("depends_on") or []


def positions(task):
    """The positions a task depends on, each once: a position listed
    twice, which the schema asks against and a hand-written plan may
    carry, is one relation, not two writes of it. Only called once
    dependency_problems has found the list to hold positions."""
    return list(dict.fromkeys(depends_on(task)))


def is_position(value):
    """Whether a value is a position: a whole number from 1, and not a
    bool, which Python counts as one."""
    return isinstance(value, int) and not isinstance(value, bool) and value >= 1


def dependency_problems(tasks):
    """What is wrong with the plan's dependencies, each in one line, in
    the order a reader finds them: a task whose depends_on lists something
    other than positions, one that depends on a position the plan has no
    task of, one that depends on itself, and a cycle, which no order of
    work could ever close. The first problem refuses the plan before any
    task is opened, since a plan half written is harder to mend than one
    refused whole."""
    problems = []
    for position, task in enumerate(tasks, 1):
        listed = depends_on(task)
        # A depends_on that is no list, a bare number for one dependency,
        # is the same problem as a list of the wrong things.
        if not isinstance(listed, list) or not all(is_position(value) for value in listed):
            problems.append(f"plan: task {position}'s depends_on must list positions, whole numbers from 1")
            continue
        for value in listed:
            if value > len(tasks):
                problems.append(f"plan: task {position} depends on {value}, which the plan has no task of")
            elif value == position:
                problems.append(f"plan: task {position} depends on itself")
    if not problems:
        cycle = dependency_cycle(tasks)
        if cycle:
            problems.append("plan: tasks " + ", ".join(str(p) for p in cycle) + " depend on each other in a cycle")
    return problems


def dependency_cycle(tasks):
    """The positions of one cycle in the plan's dependencies, ascending,
    or none: a depth-first walk from each position along what it depends
    on, which meets a position still on its path exactly when there is a
    cycle. Only called once every position listed is a task's."""
    edges = {position: positions(task) for position, task in enumerate(tasks, 1)}
    done = set()
    for start in edges:
        stack = [(start, iter(edges[start]))]
        while stack:
            position, following = stack[-1]
            following_position = next(following, None)
            path = [p for p, _ in stack]
            if following_position is None:
                stack.pop()
                done.add(position)
            elif following_position in path:
                return sorted(path[path.index(following_position):])
            elif following_position not in done:
                stack.append((following_position, iter(edges[following_position])))
    return []


def task_body(epic_number, task):
    files = "\n".join(f"- {path}" for path in task.get("files", [])) or "- (none named)"
    return (f"Epic: #{epic_number}\n\n{task.get('body', '').strip()}\n\n"
            f"## Files to touch\n\n{files}\n\n"
            "Stay inside these files; other workers own the rest of the epic "
            "(docs/agents/WORKERS.md).\n")


def ensure_label(client, name):
    """Make a worker's label; one that exists already is fine, and the
    client treats it so. A label gh cannot find makes `gh issue edit
    --add-label` add none of the labels it was given, so the label is made
    before it is needed."""
    client.create_label(name, WORKER_COLOUR)


def open_task(client, epic_number, task, owner):
    """Open one task issue, given to its owner when it has one, and hang it
    under the epic as a sub-issue; the link is a convenience for people, so
    its failure is not one. The answer is the issue as GitHub returned it,
    with its "number" and its "id", which a blocked-by relation names."""
    names = [TODO_LABEL]
    if owner:
        names.append(worker_label(owner))
        ensure_label(client, names[-1])
    made = client.create_issue(task["title"], task_body(epic_number, task), names)
    try:
        client.add_sub_issue(epic_number, made["id"])
    except ClientError:
        pass
    return made


def refusal_words(problem):
    """How a refused call is told in the plan and the clerk's answer: the
    status GitHub answered, which a GitHubError carries, or that it did not
    answer."""
    status = getattr(problem, "status", None)
    return f"GitHub answered {status}" if status else "GitHub did not answer"


def write_relations(client, tasks, made):
    """Write every blocked-by relation the plan's depends_on names, once
    every task is open and has an id, and say which GitHub refused: each
    relation alone, so one refusal stops no other, and none retried, since
    the relation is a convenience the plan's text keeps too. The answer is
    the relations written, by the blocked task's number, and the refusals
    as (blocked, blocker, words)."""
    numbers = [issue["number"] for issue in made]
    blockers, refused = {}, []
    for position, task in enumerate(tasks, 1):
        for depended in positions(task):
            blocked, blocker = numbers[position - 1], numbers[depended - 1]
            try:
                client.add_blocked_by(blocked, made[depended - 1]["id"])
                blockers.setdefault(blocked, []).append(blocker)
            except ClientError as problem:
                refused.append((blocked, blocker, refusal_words(problem)))
    return blockers, refused


def task_line(number, title, owner, blockers=(), refused=()):
    """One task's line in the plan section: its number and title, its
    worker, the tasks that block it, and the relations GitHub refused to
    record, so a reader sees the plan's order even where GitHub does not."""
    line = f"- [ ] #{number} {title}" + (f" ({worker_label(owner)})" if owner else "")
    if blockers:
        line += ", blocked by " + ", ".join(f"#{blocker}" for blocker in blockers)
    if refused:
        line += "; not written: " + ", ".join(f"#{blocker} ({words})" for blocker, words in refused)
    return line


def plan_section(plan, made, blockers=None, refused=()):
    """What seed appends to the epic: the coordinator's plan, the task
    list, with the tasks that block each task where the plan has any, and the
    machine line route reads."""
    blockers = blockers or {}
    lines = [task_line(number, title, owner, blockers.get(number, ()),
                       [(blocker, words) for blocked, blocker, words in refused if blocked == number])
             for number, title, owner in made]
    numbers = " ".join(str(number) for number, _, _ in made)
    return (f"\n\n{PLAN_MARKER}\n## Plan, by the coordinator\n\n{plan.strip()}\n\n"
            "## Tasks\n\n" + "\n".join(lines) + f"\n\n<!-- team:tasks {numbers} -->\n")


def body_as_it_stands(client, epic_number, first):
    """The epic's body right before the plan is written: read again, since
    the tasks took a while to open and a person, or the board's graph, may
    have changed it meanwhile; where GitHub refuses that read the body
    read when seed began stands in, since the tasks and their relations
    exist by now and a plan not written would have the next run open them
    all again, which costs more than an edit lost."""
    try:
        return client.issue(epic_number).get("body") or ""
    except ClientError:
        return first


def write_plan(client, epic_number, first, section):
    """Write the plan section onto the epic's body as it stands, then read
    the body back once and write the section again where the marker is
    gone. The board's graph writer reads and writes the whole body too,
    and one that read the epic before this write and wrote after it would
    have put the plan out; the plan is written once and route reads it
    ever after, so losing it would open the tasks twice, while the graph
    is drawn again every run and costs nothing lost. A read-back GitHub
    refuses is let go, since the write itself went through."""
    client.set_body(epic_number, body_as_it_stands(client, epic_number, first) + section)
    try:
        written = client.issue(epic_number).get("body") or ""
    except ClientError:
        return
    if PLAN_MARKER not in written:
        client.set_body(epic_number, written + section)


def seed(client, epic_number, plan, agents):
    """Open the epic's tasks, write the dependencies the plan names as
    GitHub's blocked-by relations, and write the plan into the epic, once.
    A plan whose dependencies do not hold together is refused before any
    task is opened (dependency_problems); the plan is written onto the
    epic as it stands once the tasks are open and read back once, so a
    body written over it meanwhile does not lose it (write_plan)."""
    epic = client.issue(epic_number)
    first = epic.get("body") or ""
    if PLAN_MARKER in first:
        return "already planned"
    tasks = plan.get("tasks") or []
    if not tasks:
        raise TeamError("the plan has no tasks")
    problems = dependency_problems(tasks)
    if problems:
        raise ValueError(problems[0])
    owners = owners_for(tasks, agents)
    opened = [open_task(client, epic_number, task, owner) for task, owner in zip(tasks, owners)]
    made = [(issue["number"], task["title"], owner) for issue, task, owner in zip(opened, tasks, owners)]
    blockers, refused = write_relations(client, tasks, opened)
    write_plan(client, epic_number, first, plan_section(plan.get("plan", ""), made, blockers, refused))
    # An epic opened with the Task form carries todo, which would put it on
    # the board's queue for any worker to take whole; its tasks are the work.
    client.remove_label(epic_number, TODO_LABEL)
    answer = "seeded " + " ".join(str(number) for number, _, _ in made)
    for blocked, blocker, words in refused:
        answer += f"; #{blocked} blocked by #{blocker} not written: {words}"
    return answer


# --- the team variable ----------------------------------------------------------

def team_members(text):
    """The members from the CLAUDE_TEAM variable's JSON list."""
    try:
        members = json.loads(text or "[]")
    except ValueError:
        raise TeamError("CLAUDE_TEAM is not a JSON list") from None
    if not isinstance(members, list) or not all(isinstance(m, str) and m for m in members):
        raise TeamError("CLAUDE_TEAM must be a JSON list of names, as in [\"w1\",\"w2\"]")
    return members
