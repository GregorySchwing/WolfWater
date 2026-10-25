"""The command line: python3 -m sdlc_protocol <board|check|team|workers|guard|verify-vendored>.

Each subcommand is a thin shell over the modules and takes exactly the
arguments the template's tool for the same job took, so a workflow line
changes only from "python3 tools/board.py ..." to "python3 -m
sdlc_protocol board ...":

  board  release the claims gone a day without a commit and tell open
         work that meets other work in a file (the coordinator), render
         STATE.md from GitHub (--default-branch, --sha, --out,
         --render-only), rename the branches of released claims and
         comment on the pull requests that wait, or render a JSON
         snapshot with --dry-run FIXTURE and no network
  check  hold the pull request in GITHUB_EVENT_PATH to the claim protocol
  team   [--repo OWNER/NAME] seed, route, next, deal or schema, the clerk;
         seed, next, deal and schema read the workers from the registry
         with --registry DIR
  workers check, list, describe, get or from-form: the registry of
         workers, .github/symphony/workers/<name>.json, read offline
  guard  the PreToolUse hook: the tool call's JSON on stdin, and exit 2
         with the reason on stderr to refuse a push
  locks  status, plan, acquire, release, sweep, verify or overlap: the locks on
         files, symbols and named resources in LOCKS.json on the branch
         locks; verify obeys the switch SYMPHONY_LOCKS, off by default
  verify-vendored VENDOR_DIR
         hold a vendored copy, VENDOR_DIR/sdlc_protocol, to the record
         VENDOR_DIR/VENDORED beside it (vendoring.py, which a copy also
         runs as a script of its own)

The exit codes are the old tools' too: 0 when done, 1 for a failure told
in one line on stderr (or, for check, the problems on stdout), 2 for a
push the guard refuses or arguments that cannot be read.

Only board and team import the GitHub adapter, and only when they run:
the guard answers before every Bash call of a session, and the start-up
time of urllib is no use to it.
"""

import argparse
import json
import os
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

from . import board, claims, guard, protocol_check, team, vendoring, workers

PROG = "python3 -m sdlc_protocol"
USAGE = (f"usage: {PROG} {{board,check,team,workers,guard,locks,verify-vendored}} [arguments]; "
         f"{PROG} <command> -h says more")


class CommandError(Exception):
    """What a subcommand cannot start without, in one line, never a token."""


def github_adapter():
    """The GitHub adapter, imported when a subcommand first calls GitHub."""
    from .adapters import github
    return github


def load_json(path):
    return json.loads(Path(path).read_text(encoding="utf-8"))


def environment(given):
    """The environment a subcommand reads: the one a test hands in, else
    the process's own."""
    return os.environ if given is None else given


# --- board --------------------------------------------------------------------

def board_arguments(argv):
    parser = argparse.ArgumentParser(
        prog=f"{PROG} board",
        description="Release idle claims and tell overlapping work (the coordinator), render the board "
                    "(STATE.md), rename released claims' branches and comment on waiting pull requests.")
    parser.add_argument("--dry-run", metavar="FIXTURE", help="render this JSON snapshot instead of asking GitHub")
    parser.add_argument("--render-only", action="store_true",
                        help="write the board, and release no claim, rename no branch and post no comment")
    parser.add_argument("--out", help=f"where to write the board (default: {board.STATE_PATH} in the current "
                                      "directory; stdout with --dry-run)")
    parser.add_argument("--default-branch", default="main", help="the default branch's name (default: main)")
    parser.add_argument("--sha", help="the default branch's head, which the board describes and names in its "
                        "main-moved comments on every event (default: GITHUB_SHA where that is the default "
                        "branch's head, which it is not on a pull_request event or a push to another branch; "
                        "then no main-moved comment is made)")
    return parser.parse_args(argv)


def board_dry_run(args):
    """The fixture's board, as it stands: the fixture names its own default
    branch, as the template's dry run read it."""
    text = board.render(load_json(args.dry_run))
    if args.out:
        board.write_if_changed(Path(args.out), text)
    else:
        sys.stdout.write(text)
    return 0


def plan_board(api, sha, default_branch, now, render_only=False, say=print, started_by=""):
    """The first half of a board run: read the repository as GitHub holds
    it, release the claims due and tell overlapping work (the coordinator),
    and render the page; with render_only, the command's --render-only,
    the coordinator is skipped and the half only reads and renders. The
    snapshot as the coordinator left it, which the second half acts on,
    and the page. The run is in two halves so that the validation harness
    can run them on its own clock with other work between them, as the
    half-minute of a workflow run allows; say takes each log line as it is
    made, and started_by is the branch whose push started the run
    (pushed_branch), "" for none, whose claim no read of the activity
    dates, or one dates only by an entry from before the run began, is
    dated at the run's start (adapters.github.head_since)."""
    github = github_adapter()
    snap = github.snapshot(api, sha, default_branch, now, started_by)
    for line in (github.unread_claim_lines(snap) + github.pushed_refused_lines(snap)
                 + github.began_refused_lines(snap) + github.dependencies_skipped_lines(snap)):
        say(line)
    if not render_only:
        for line in github.coordinate(api, snap, now):
            say(line)
    return snap, board.render(snap)


def act_on_board(api, snap, event_name, ref_name, now, say=print):
    """The second half, once the page is written: rename the released
    claims' branches, post the comments the pull requests are owed, and
    write each open epic's flowchart into its body, all from the snapshot
    the first half read; the event and the ref are passed on, though
    since F12's fix no comment depends on them."""
    github = github_adapter()
    for line in github.release_claims(api, snap):
        say(line)
    for line in github.post_comments(api, board.comments_due(snap, event_name, ref_name, now)):
        say(line)
    for line in github.write_graphs(api, snap, board.epic_blocks(snap)):
        say(line)


def default_head_sha(environ, default_branch):
    """GITHUB_SHA where it is the default branch's head: on a push to that
    branch, and on every event but a pull request's and a push to another
    branch, which carry the pull request's merge commit or that branch's
    head, no commit of the default branch at all; there it is "", so the
    board names no commit and tells nobody that main moved (board.
    comments_due), since the word now goes out on every event and would
    otherwise carry a sha that is not main's, once per such event. A
    workflow that passes --sha from a checkout of the default branch, as
    the template's does, never comes here."""
    event, ref = environ.get("GITHUB_EVENT_NAME", ""), environ.get("GITHUB_REF_NAME", "")
    if event.startswith("pull_request") or (event == "push" and ref != default_branch):
        return ""
    return environ.get("GITHUB_SHA", "")


def pushed_branch(environ):
    """The branch whose push started the run, as the workflow's event names
    it, "" where none did: on a push or a branch's creation the ref pushed
    (GITHUB_REF), when it is a branch's and no tag's; on a pull request's
    synchronize, the head branch the push moved, read from the event's
    payload, when the head is this repository's, since a fork's branch is
    no claim here; and nothing on any other event, nor on a pull request
    event whose action is not synchronize, whose payload is read all the
    same, since the action is in it, nor on a payload that cannot be read
    or where the environment does not name the repository the run is in
    (GITHUB_REPOSITORY). A claim on that branch whose push no read of the
    activity dates, or one dates only by an entry from before the run
    began, is dated at the run's start, for that run alone
    (spec/PROTOCOL.md, Release; adapters.github.head_since): the push
    that started the run is the one the activity is likeliest not to list
    yet. A synchronize that GitHub fires for a change of the pull request's
    base, with no push, holds such a claim one run too, which the next run
    reads afresh."""
    event = environ.get("GITHUB_EVENT_NAME", "")
    if event in ("push", "create"):
        ref = environ.get("GITHUB_REF", "")
        return ref[len("refs/heads/"):] if ref.startswith("refs/heads/") else ""
    if event != "pull_request":
        return ""
    payload, _ = read_event(environ)
    if payload is None or payload.get("action") != "synchronize":
        return ""
    head = (payload.get("pull_request") or {}).get("head") or {}
    repo = environ.get("GITHUB_REPOSITORY")
    if not repo or (head.get("repo") or {}).get("full_name") != repo:
        return ""
    ref = head.get("ref")
    return ref if isinstance(ref, str) else ""


def say_what_it_cost(api):
    """The board run's last line (adapters.github.calls_line), which never
    takes the place of the run's own outcome: a client that keeps no count,
    or keeps one in another shape, as a stand-in may, or a time no clock
    writes, says nothing."""
    try:
        print(github_adapter().calls_line(api))
    except (AttributeError, LookupError, TypeError, ValueError, OverflowError, OSError):
        return


def board_live_run(args, environ):
    """A board run on GitHub: the two halves, the page written between
    them, and last, whether the run ended or GitHub refused it on the way,
    the calls GitHub answered and what the token had left (say_what_it_cost),
    the one figure of what a board run costs taken in the repository
    itself."""
    token, repo = environ.get("GITHUB_TOKEN"), environ.get("GITHUB_REPOSITORY")
    if not token or not repo:
        raise CommandError("GITHUB_TOKEN and GITHUB_REPOSITORY must be set")
    api = github_adapter().GitHub(token, repo)
    now = datetime.now(timezone.utc)
    try:
        snap, page = plan_board(api, args.sha or default_head_sha(environ, args.default_branch), args.default_branch,
                                now, args.render_only, started_by=pushed_branch(environ))
        out = Path(args.out) if args.out else board.STATE_PATH
        changed = board.write_if_changed(out, page)
        print(f"{out.as_posix()}: {'rewritten' if changed else 'unchanged'}")
        if not args.render_only:
            act_on_board(api, snap, environ.get("GITHUB_EVENT_NAME", ""), environ.get("GITHUB_REF_NAME", ""), now)
    finally:
        say_what_it_cost(api)
    return 0


def run_board(argv, environ=None):
    args = board_arguments(argv)
    try:
        return board_dry_run(args) if args.dry_run else board_live_run(args, environment(environ))
    except (CommandError, team.ClientError) as problem:
        print(f"board: {problem}", file=sys.stderr)
        return 1


# --- check --------------------------------------------------------------------

def read_event(environ):
    """The event GitHub hands the run, and what to say when there is none
    to check: no event, or one that is not a JSON object."""
    where = environ.get("GITHUB_EVENT_PATH")
    if not where or not Path(where).is_file():
        return None, "no pull request event to check"
    try:
        event = load_json(where)
    except ValueError:
        event = None
    if not isinstance(event, dict):
        return None, "the event is not readable; nothing checked"
    return event, ""


# How long one git call of the check may take: a fetch is of a few
# commits' objects, or of the history down to a merge base.
GIT_SECONDS = 120
# What git is run with whatever the environment or the checkout's own
# configuration say: no automatic housekeeping, and pathspecs read with
# their magic, never literally (F31's second review).
GIT_SETTINGS = ("--no-literal-pathspecs", "-c", "gc.auto=0", "-c", "maintenance.auto=false",
                "-c", "diff.ignoreSubmodules=none")
UNREAD = "the claim files this pull request changes could not be read: "
# How far a shallow clone is deepened, in turn, to find where a pull
# request's branch left its base; the last step takes the whole history.
DEEPENINGS = (("--deepen=64",), ("--deepen=512",), ("--unshallow",))
# How much of git's last words a reason quotes.
STDERR_TAIL = 200


class Git:
    """git run in one checkout, with every GIT_* variable of the
    environment left out, since an earlier step of the job can set any of
    them through $GITHUB_ENV (GIT_LITERAL_PATHSPECS, which empties the
    diff, GIT_DIR, which reads another repository, GIT_CONFIG_*), and the
    last failure's words kept for the reason the check gives."""

    def __init__(self, root):
        self.root = root
        self.environ = {k: v for k, v in os.environ.items() if not k.upper().startswith("GIT_")}
        self.why = ""

    def __call__(self, *args):
        """git's standard output, or None, with why saying what went
        wrong: git missing, a timeout, or the tail of what git wrote to
        its standard error."""
        try:
            done = subprocess.run(["git", *GIT_SETTINGS, *args], cwd=self.root, env=self.environ,
                                  capture_output=True, timeout=GIT_SECONDS)
        except subprocess.TimeoutExpired:
            self.why = f"git {args[0]} timed out after {GIT_SECONDS} s"
            return None
        except (OSError, subprocess.SubprocessError, ValueError) as error:
            self.why = f"git could not be run ({error.__class__.__name__})"
            return None
        if done.returncode != 0:
            self.why = f"git {args[0]} failed: {last_words(done.stderr)}"
            return None
        return done.stdout

    def said(self, what):
        """A reason, with git's own last words where it gave any."""
        return f"{what} ({self.why})" if self.why else what


def last_words(stderr):
    """The last line git wrote to its standard error, as plain ASCII of
    bounded length, or how it ended where it wrote nothing."""
    lines = [line.strip() for line in stderr.decode("utf-8", "replace").splitlines() if line.strip()]
    text = (lines[-1] if lines else "no message").encode("ascii", "replace").decode("ascii")
    return text[:STDERR_TAIL]


def event_sha(side):
    """A commit an event's head or base names, or None where it names
    none git could take as one."""
    sha = (side or {}).get("sha")
    return sha if isinstance(sha, str) and 7 <= len(sha) <= 64 and sha.isascii() and sha.isalnum() else None


def parents_of(git, commit):
    """The parents a commit's own object names, which a clone of depth 1
    still holds where commit^1 is not there; None where git cannot read
    it."""
    text = git("cat-file", "commit", commit)
    if text is None:
        return None
    return [line.split()[1] for line in text.decode("ascii", "replace").splitlines()
            if line.startswith("parent ") and len(line.split()) == 2]


def fetched(git, sha):
    """Whether the commit is in the clone, fetched one deep from origin
    where it was not, as actions/checkout's default depth of 1 leaves it
    out; origin is the base repository, which holds a fork's head too."""
    if git("cat-file", "-e", f"{sha}^{{commit}}") is not None:
        return True
    if git("fetch", "--quiet", "--no-tags", "--depth=1", "origin", sha) is None:
        return False
    return git("cat-file", "-e", f"{sha}^{{commit}}") is not None


def merge_bases(git, one, two):
    """Every commit where the branch at two left one, git's merge bases
    (`merge-base --all`, so criss-cross history, which has several, hides
    nothing behind the one git would pick), deepening a shallow clone step
    by step until git finds one; () where it never does, which the check
    fails on rather than guess. one and two are shas, so a deepening
    fetches the pull request's own commits."""
    for deepening in ((),) + DEEPENINGS:
        if deepening:
            shallow = git("rev-parse", "--is-shallow-repository")
            if shallow is None or shallow.strip() != b"true":
                return ()
            git("fetch", "--quiet", "--no-tags", *deepening, "origin", one, two)
        found = git("merge-base", "--all", one, two)
        if found is not None and found.split():
            return tuple(line.decode("ascii", "replace") for line in found.split())
    return ()


def checkout_of(git, pull, base):
    """What the checkout is: the commit whose change is the pull request's
    own, the commit GitHub's test merge was made on where HEAD is that
    merge, and whether HEAD is the pull request at all. HEAD is the test
    merge where it is a merge whose second parent is the head the event
    names, whatever its first: GitHub makes the merge on the base branch
    as it stands when the event fires, while the event's base is a commit
    it set when the pull request was opened and leaves while the base
    branch alone moves (symphony-template#315's pushes while behind main
    carried the base it was opened on), so the two agree only while the
    base has not moved past the pull request, which GitHub's rule that a
    branch be up to date no longer prevents; such a pull request is the
    test merge by its second parent alone, and its first parent is what
    the merge was made on, the merge's own change being read against it.
    Where the
    event names no head, HEAD is the test merge where it is a merge whose
    first parent is the base. The tip is the event's head wherever it
    names one, so a checkout of something else, the base as
    pull_request_target leaves it among them, still has the pull request's
    own change judged; and where the event names none, HEAD itself, as a
    sha, a head commit that is itself a merge among them. None for the tip
    where git cannot read HEAD."""
    here = git("rev-parse", "--verify", "--quiet", "HEAD^{commit}")
    parents = parents_of(git, "HEAD") if here is not None else None
    if parents is None:
        return None, None, False
    here = here.decode("ascii", "replace").strip()
    head = event_sha(pull.get("head"))
    test_merge = len(parents) == 2 and (parents[1] == head if head is not None else parents[0] == base)
    made_on = parents[0] if test_merge else None
    tip = parents[1] if test_merge else (head or here)
    return tip, made_on, test_merge or head is None or here == head


def change_of(git, pull):
    """The diffs a pull request's change is judged by, read with git from
    the checkout the check runs in, over every path, since a path that
    folds into the claim folder may be spelt in letters git's own icase
    pathspec does not fold (F31's second review), or why they could not
    be read: the
    branch's own change, from each commit where it left the base to its
    tip (git's three dots, from every merge base), so a file the base
    also holds, which the merge leaves alone, still counts where the
    branch added it; and, in GitHub's test merge, the merge's change
    against the commit it was made on, its first parent, as well, which
    is the base branch's tip when the event fired and not the event's
    base, a commit the base branch may have moved past, so that what the
    base branch did since is never read as the pull request's; the diffs
    from every merge base already hold whatever the merge brings, so that
    second reading is a second reading of the same change, never the only
    one. A checkout that is neither
    the test merge nor the event's head is no reading of the pull
    request's claim file, so it is a reason too, beside the change
    judged."""
    base = event_sha(pull.get("base"))
    if base is None:
        return None, "the event names no base commit"
    if not fetched(git, base):
        return None, git.said(f"the base commit {base[:12]} could not be fetched")
    tip, made_on, is_pull = checkout_of(git, pull, base)
    if tip is None:
        return None, git.said("git could not read the checked-out commit")
    if not fetched(git, tip):
        return None, git.said(f"the head commit {tip[:12]} could not be fetched")
    if made_on is not None and not fetched(git, made_on):
        return None, git.said(f"the commit the test merge was made on, {made_on[:12]}, could not be fetched")
    lefts = merge_bases(git, base, tip)
    if not lefts:
        return None, git.said(f"git found no commit where the branch left {base[:12]}")
    pairs = [(left, tip) for left in lefts] + ([(made_on, "HEAD")] if made_on is not None else [])
    diffs = [git("diff", "--no-renames", "--raw", "-z", older, newer) for older, newer in pairs]
    if any(d is None for d in diffs):
        return None, git.said(f"git could not compare {base[:12]} with the checkout")
    away = "" if is_pull else (f"the checkout is neither the head the event names, {tip[:12]}, nor GitHub's test "
                               f"merge of it, so the claim file looked for there is not the pull request's")
    return diffs, away


def changed_claims_in(root, pull):
    """The issues whose claim files the pull request adds or changes and
    the problems of what it does to the claim folder, read from the
    checkout (change_of), as (numbers, problems): ((), []) where root holds
    no .git, as a test's scratch folder does, where nothing is judged;
    and where it holds one, every problem, a reason the change could not
    be read among them, since a claim file the check cannot see may be the
    one it is there to stop (F31 and its reviews): git missing or refusing
    the checkout, a checkout whose top is another folder, a base or head
    that cannot be fetched, a branch whose merge base is never found, or a
    checkout that is not the pull request."""
    if not (Path(root) / ".git").exists():
        return (), []
    git = Git(root)
    top = git("rev-parse", "--show-toplevel")
    if top is None:
        return (), [UNREAD + git.said("git is missing or will not read the checkout")]
    if Path(os.fsdecode(top.strip())).resolve() != Path(root).resolve():
        return (), [UNREAD + f"the checkout's top is {os.fsdecode(top.strip())}, not the folder the check runs in"]
    diffs, why = change_of(git, pull)
    if diffs is None:
        return (), [UNREAD + why]
    numbers = sorted({n for diff in diffs for n in protocol_check.claim_files_changed(diff)})
    problems = list(dict.fromkeys(p for diff in diffs for p in protocol_check.claim_folder_problems(diff)))
    return tuple(numbers), problems + ([why] if why else [])


def run_check(argv, environ=None, root="."):
    """The protocol check: a run outside a pull request event passes, so
    the check never breaks a push build; the claim file is looked for in
    the checkout, root, and the claim files the pull request adds or
    changes are read from it with git (changed_claims_in)."""
    argparse.ArgumentParser(
        prog=f"{PROG} check",
        description="Hold the claude/<issue> pull request in GITHUB_EVENT_PATH to the claim protocol.",
    ).parse_args(argv)
    event, why_not = read_event(environment(environ))
    if event is None:
        print(why_not)
        return 0
    pull = event.get("pull_request") or {}
    changed, unread = changed_claims_in(root, pull) if pull else ((), [])
    problems = protocol_check.problems_for(
        pull, lambda number: (Path(root) / claims.claim_path(number)).is_file(), changed)
    problems += unread
    for problem in problems:
        print(problem)
    if problems:
        return 1
    print("protocol ok")
    return 0


# --- team ---------------------------------------------------------------------

def team_arguments(argv, environ=None):
    environ = environment(environ)
    parser = argparse.ArgumentParser(prog=f"{PROG} team", description="The clerk's steps (docs/agents/WORKERS.md).")
    parser.add_argument("--repo", default=environ.get("GITHUB_REPOSITORY"))
    sub = parser.add_subparsers(dest="command", required=True)
    seeding = sub.add_parser("seed")
    seeding.add_argument("--epic", type=int, required=True)
    seeding.add_argument("--plan", required=True, help="the coordinator's plan as a JSON file")
    seeding.add_argument("--team", default=environ.get("CLAUDE_TEAM", ""))
    seeding.add_argument("--registry", help="give the tasks to the workers of this registry folder, of any kind")
    routing = sub.add_parser("route")
    routing.add_argument("--event", required=True)
    routing.add_argument("--issue", type=int, required=True)
    routing.add_argument("--label", default=None)
    nexting = sub.add_parser("next")
    nexting.add_argument("--agent", required=True)
    nexting.add_argument("--team", default=environ.get("CLAUDE_TEAM", ""))
    nexting.add_argument("--registry", help="deal the shared work among this registry's workers in GitHub "
                                            "Actions, each by its area")
    nexting.add_argument("--label", default=None,
                         help="the worker's label, made before a take so the session can add it")
    dealing = sub.add_parser("deal", help="the next move of every worker in Actions of this registry, from one "
                                          "look at the repository, as a JSON object by worker on one line")
    dealing.add_argument("--registry", required=True, help="the registry whose workers in Actions are dealt to")
    dealing.add_argument("--workers", default=None,
                         help="a JSON list of the workers of this run (default: every worker in Actions): "
                              "theirs are the moves printed and the labels made before a take")
    schema = sub.add_parser("schema", help="the JSON schema of the coordinator's plan, on one line")
    schema.add_argument("--registry", help="the registry whose workers a task may go to")
    return parser.parse_args(argv)


def team_answer(client, args):
    """The one line the workflow reads, from the clerk's rules. With a
    registry, seed gives tasks to every registered worker, next deals the
    shared work among the workers in Actions, each by its area, from the
    asking worker's own look, and deal deals it once for the whole run,
    from one look, and answers for every worker of the run."""
    registered = workers.load(args.registry) if getattr(args, "registry", None) else None
    if args.command == "seed":
        names = [w.name for w in registered] if registered is not None else team.team_members(args.team)
        return team.seed(client, args.epic, load_json(args.plan), names)
    if args.command == "route":
        return team.route(client, args.event, args.issue, args.label)
    if args.command == "deal":
        return deal_answer(client, registered, args.workers)
    if registered is None:
        return team.next_move(client, args.agent, team.team_members(args.team), args.label)
    in_actions = workers.in_actions(registered)
    areas = {w.name: w.area for w in registered}
    return team.next_move(client, args.agent, [w.name for w in in_actions], args.label, areas)


def deal_answer(client, registered, run_json):
    """The deal's one line: a JSON object, worker to move, for the workers
    of the run, the registry's workers in Actions unless --workers names
    them; a name that is no worker here, or a value that is not JSON, is a
    ValueError naming the flag and the shape it takes, since the workflow
    took its list from the same registry; so is an empty --registry, which
    argparse lets through and which names no folder."""
    if registered is None:
        raise ValueError("--registry must name the registry folder")
    in_actions = workers.in_actions(registered)
    names = [w.name for w in in_actions]
    shape = f"--workers must be a JSON list of the registry's workers in Actions, {json.dumps(names)}"
    try:
        run = names if run_json is None else json.loads(run_json)
    except ValueError as problem:
        raise ValueError(f"{shape}; the value given is not JSON: {problem}") from None
    if not isinstance(run, list) or any(name not in names for name in run):
        raise ValueError(shape)
    moves = team.deal_moves(client, names, {w.name: w.area for w in registered},
                            {w.name: w.label for w in in_actions}, run)
    return json.dumps(moves, sort_keys=True, separators=(",", ":"))


def team_schema(args):
    """The plan's schema, which needs no GitHub: the registry's names, or
    none, for the owners."""
    names = [w.name for w in workers.load(args.registry)] if args.registry else []
    print(json.dumps(team.plan_schema(names), separators=(",", ":")))
    return 0


def run_team(argv, environ=None):
    environ = environment(environ)
    args = team_arguments(argv, environ)
    if args.command == "schema":
        try:
            return team_schema(args)
        except workers.WorkerError as problem:
            print(f"team: {problem}", file=sys.stderr)
            return 1
    token = environ.get("GITHUB_TOKEN") or environ.get("GH_TOKEN")
    if not token or not args.repo:
        print("GITHUB_TOKEN and --repo (or GITHUB_REPOSITORY) are needed", file=sys.stderr)
        return 1
    try:
        print(team_answer(github_adapter().GitHub(token, args.repo), args))
    except (team.ClientError, team.TeamError, workers.WorkerError, ValueError, KeyError) as problem:
        print(f"team: {problem}", file=sys.stderr)
        return 1
    return 0


# --- workers ------------------------------------------------------------------

def workers_arguments(argv):
    parser = argparse.ArgumentParser(
        prog=f"{PROG} workers",
        description="The registry of workers, one JSON file per worker (docs/agents/WORKERS.md); no network.")
    parser.add_argument("--dir", default=workers.REGISTRY, help=f"the registry's folder (default: {workers.REGISTRY})")
    sub = parser.add_subparsers(dest="command", required=True)
    sub.add_parser("check", help="hold every worker's file to the registry's rules; exit 1 on any problem")
    listing = sub.add_parser("list", help="the workers' names as a JSON list, as a workflow's matrix reads it")
    listing.add_argument("--runs", choices=workers.RUNS, help="only the workers that run this way")
    sub.add_parser("describe", help="one line per worker: its name, where it runs, its area, what it is for")
    getting = sub.add_parser("get", help="one setting of one worker")
    getting.add_argument("name")
    getting.add_argument("setting", choices=workers.KEYS + ("label",))
    forming = sub.add_parser("from-form", help="write the worker an \"Add a worker\" issue's answers describe")
    forming.add_argument("--body", help="a file holding the issue's body (default: stdin)")
    return parser.parse_args(argv)


def workers_check(args):
    found = workers.load(args.dir)
    listed = ", ".join(f"{w.name} ({w.runs})" for w in found) or "none"
    print(f"{len(found)} worker{'s' if len(found) != 1 else ''} in {args.dir}: {listed}")
    return 0


def workers_list(args):
    found = workers.load(args.dir)
    print(json.dumps([w.name for w in found if not args.runs or w.runs == args.runs]))
    return 0


def workers_describe(args):
    found = workers.load(args.dir)
    print("\n".join(workers.describe(w) for w in found) if found else "(no worker is registered)")
    return 0


def workers_get(args):
    worker = workers.named(workers.load(args.dir), args.name)
    if worker is None:
        raise CommandError(f"no worker named {args.name} in {args.dir}")
    value = worker.label if args.setting == "label" else worker.to_dict()[args.setting]
    print(value if isinstance(value, str) else json.dumps(value))
    return 0


def workers_from_form(args, stdin):
    """Write the worker the form's answers describe into the registry and
    print its name and how it runs, for the workflow that opens the pull
    request; a name already registered is refused, since a worker's file is
    changed by a pull request of its own, never by a second form."""
    body = Path(args.body).read_text(encoding="utf-8") if args.body else (stdin or sys.stdin).read()
    worker = workers.from_dict(workers.from_form(body))
    if workers.named(workers.load(args.dir), worker.name) is not None:
        raise CommandError(f"a worker named {worker.name} exists already; change {args.dir}/{worker.file_name} "
                           f"by pull request, or pick another name")
    folder = Path(args.dir)
    folder.mkdir(parents=True, exist_ok=True)
    (folder / worker.file_name).write_text(worker.to_json(), encoding="ascii")
    print(f"{worker.name} {worker.runs}")
    return 0


def run_workers(argv, stdin=None):
    args = workers_arguments(argv)
    commands = {"check": workers_check, "list": workers_list, "describe": workers_describe, "get": workers_get,
                "from-form": lambda given: workers_from_form(given, stdin)}
    try:
        return commands[args.command](args)
    except workers.WorkerError as problem:
        for line in str(problem).splitlines():
            print(f"workers: {line}", file=sys.stderr)
        return 1
    except (CommandError, OSError, UnicodeError) as problem:
        print(f"workers: {problem}", file=sys.stderr)
        return 1


# --- guard --------------------------------------------------------------------

GUARD_HELP = (f"usage: {PROG} guard < payload.json\n\n"
              "The PreToolUse hook on Bash: reads the tool call's JSON on stdin and exits 2, "
              "with the reason on stderr, to refuse a push the board rules forbid.")


def run_guard(argv, stdin=None):
    """The push guard as Claude Code's hook runs it. Its arguments are not
    parsed: exit 2 stops the tool call, so an argument error would refuse
    every Bash call of the session, and the guard fails open instead, as
    it does for a payload that is not JSON. Help is the one argument it
    answers, and it lets the call pass too."""
    if "-h" in argv or "--help" in argv:
        print(GUARD_HELP)
        return guard.PASS
    try:
        payload = json.load(sys.stdin if stdin is None else stdin)
    except ValueError:
        return guard.PASS
    code, reason = guard.run(payload)
    if reason:
        print(reason, file=sys.stderr)
    return code


locks = git_locks = None  # bound by run_locks


# --- locks --------------------------------------------------------------------

def locks_arguments(argv):
    parser = argparse.ArgumentParser(
        prog=f"{PROG} locks",
        description="Locks on files, symbols and named resources, kept in LOCKS.json on the branch locks "
                    "(sdlc_protocol.locks). verify obeys the switch SYMPHONY_LOCKS or .github/symphony-locks "
                    "(off, advise, enforce); the other commands do what they are asked.")
    parser.add_argument("--repo", default=".", help="the clone to work in")
    parser.add_argument("--remote", default="origin")
    parser.add_argument("--default-branch", default="main")
    sub = parser.add_subparsers(dest="command", required=True)
    sub.add_parser("status", help="the table, one entry a line")
    planning = sub.add_parser("plan", help="what a claim asks for to change KEYS: a write on each, reads on "
                                           "what each Python symbol's code depends on")
    planning.add_argument("keys", nargs="+", metavar="KEY")
    acquiring = sub.add_parser("acquire", help="take KEYS for --claim, all or none")
    acquiring.add_argument("keys", nargs="+", metavar="KEY")
    acquiring.add_argument("--claim", type=int, required=True)
    acquiring.add_argument("--worker", default="")
    acquiring.add_argument("--read", action="append", default=[], metavar="KEY", help="a read lock as well")
    acquiring.add_argument("--no-reads", action="store_true", help="take no read locks from the plan")
    acquiring.add_argument("--with-claim", action="store_true",
                           help="push the local branch claude/<claim> in the same atomic push: the claim itself")
    releasing = sub.add_parser("release", help="give up every lock of --claim")
    releasing.add_argument("--claim", type=int, required=True)
    sweeping = sub.add_parser("sweep", help="drop the locks of claims that have ended")
    sweeping.add_argument("--ended", type=int, action="append", default=[], metavar="N",
                          help="a claim the caller knows has ended (its issue closed, its release counted)")
    verifying = sub.add_parser("verify", help="hold the change from --base to --head to the locks: what it "
                                              "touches that another claim holds, and, for a claim, what it "
                                              "touches without holding it")
    verifying.add_argument("--base", required=True)
    verifying.add_argument("--head", default="HEAD")
    verifying.add_argument("--claim", type=int, default=0, help="the claim's issue, 0 for a pull request of none")
    overlapping = sub.add_parser("overlap", help="whether what the default branch --base changed since the "
                                                 "change's fork meets what --head changes or builds on, so the "
                                                 "change must bring main in; no lock table is read")
    overlapping.add_argument("--base", required=True)
    overlapping.add_argument("--head", default="HEAD")
    return parser.parse_args(argv)


def entry_line(entry):
    return f"{entry.key}\t{entry.mode}\t#{entry.claim}\t{entry.worker or 'unnamed'}"


def locks_status(args):
    _, entries = git_locks.fetch_table(args.repo, args.remote)
    print("\n".join(entry_line(entry) for entry in entries) or "no locks")
    return 0


def tracked_files(repo):
    listed = subprocess.run(["git", "ls-files", "-z"], cwd=repo, capture_output=True, text=True)
    return [name for name in listed.stdout.split("\0") if name]


def file_reader(repo):
    def read(path):
        try:
            return (Path(repo) / path).read_text(encoding="utf-8", errors="replace")
        except OSError:
            return None
    return read


def planned(args, claim=0, worker=""):
    return locks.plan(args.keys, file_reader(args.repo), tracked_files(args.repo), claim, worker)


def locks_plan(args):
    print("\n".join(entry_line(entry) for entry in planned(args)))
    return 0


def locks_acquire(args):
    wanted = planned(args, args.claim, args.worker)
    if args.no_reads:
        wanted = [entry for entry in wanted if entry.mode == locks.WRITE]
    wanted += [locks.Entry(key, args.claim, locks.READ, args.worker) for key in args.read]
    branch = claims.branch_of(args.claim) if args.with_claim else ""
    after = git_locks.acquire(args.repo, args.remote, wanted, args.default_branch, branch)
    print("\n".join(entry_line(entry) for entry in locks.held_by(after, args.claim)))
    return 0


def locks_release(args):
    git_locks.release(args.repo, args.remote, args.claim, args.default_branch)
    print(f"#{args.claim} holds no locks")
    return 0


def locks_sweep(args):
    after = git_locks.sweep(args.repo, args.remote, args.default_branch, args.ended)
    print(f"{len(after)} locks held")
    return 0


def switch_text(repo):
    try:
        return (Path(repo) / locks.SWITCH_FILE).read_text(encoding="utf-8")
    except OSError:
        return None


def locks_verify(args, environ):
    mode = locks.mode_of(environ, switch_text(args.repo))
    if mode == locks.OFF:
        print("locks are off")
        return 0
    try:
        problems = verify_problems(args)
    except git_locks.LocksFailed as problem:
        # At advise nothing fails, a clone too shallow to read included:
        # actions/checkout's default depth must not turn a run red.
        print(f"locks: {problem}")
        return 1 if mode == locks.ENFORCE else 0
    for line in problems:
        print(f"locks: {line}")
    if not problems:
        print("locks: the change touches nothing it may not")
    return 1 if problems and mode == locks.ENFORCE else 0


def locks_overlap(args, environ):
    mode = locks.mode_of(environ, switch_text(args.repo))
    if mode == locks.OFF:
        print("locks are off")
        return 0
    try:
        met = git_locks.overlap(args.repo, args.base, args.head)
    except git_locks.LocksFailed as problem:
        print(f"locks: {problem}")
        return 1 if mode == locks.ENFORCE else 0
    for ours, theirs in met:
        print(f"locks: {meeting_line(args.base, ours, theirs)}: bring {args.base} in")
    if not met:
        print(f"locks: nothing {args.base} changed since the base meets this change; it need not bring it in")
    return 1 if met and mode == locks.ENFORCE else 0


def meeting_line(base, ours, theirs):
    """Why a change must bring the default branch in: what each side did
    to the key, edited it or built on it, since a read on either side
    means that side changed something that depends on the key."""
    did = "edited" if theirs.mode == locks.WRITE else "built on"
    does = "edits" if ours.mode == locks.WRITE else "builds on"
    return f"{base} {did} {theirs.key} since this change's base, which this change {does}"


def verify_problems(args):
    _, entries = git_locks.fetch_table(args.repo, args.remote)
    numbers = {entry.claim for entry in entries} - {args.claim}
    alive = git_locks.live_claims(args.repo, args.remote, numbers, args.default_branch)
    held = [entry for entry in entries if entry.claim in alive or entry.claim == args.claim]
    touched = git_locks.touched(args.repo, args.base, args.head)
    problems = [f"{want.key} is held by #{other.claim} ({other.worker or 'unnamed'}, {other.mode})"
                for want, other in locks.conflicts(held, [locks.Entry(key, args.claim) for key in touched])]
    if args.claim:
        # A file the change adds, its own claim file among them, was no one
        # else's to edit; another claim that means to add it holds it, and
        # the conflicts above say so.
        new = set(git_locks.added(args.repo, args.base, args.head))
        problems += [f"{key} is changed but #{args.claim} holds no write lock on it"
                     for key in locks.uncovered(touched, held, args.claim) if locks.split_key(key)[0] not in new]
    return problems


def run_locks(argv, environ=None):
    # Imported here, as the GitHub adapter is, so the guard, which answers
    # before every Bash call, never pays for ast and the git side.
    global locks, git_locks
    from . import locks
    from .adapters import git_locks
    args = locks_arguments(argv)
    environ = os.environ if environ is None else environ
    commands = {"status": locks_status, "plan": locks_plan, "acquire": locks_acquire, "release": locks_release,
                "sweep": locks_sweep, "verify": lambda given: locks_verify(given, environ),
                "overlap": lambda given: locks_overlap(given, environ)}
    try:
        return commands[args.command](args)
    except locks.Held as problem:
        print(f"locks: refused: {problem}", file=sys.stderr)
        return 1
    except (git_locks.LocksFailed, git_locks.ClaimTaken) as problem:
        print(f"locks: {problem}", file=sys.stderr)
        return 1


# --- verify-vendored ----------------------------------------------------------

def run_verify_vendored(argv):
    """vendoring.py's check under the package's own command line, for a
    consumer that would rather name the module than the file."""
    parser = argparse.ArgumentParser(
        prog=f"{PROG} verify-vendored",
        description="Hold VENDOR_DIR/sdlc_protocol to the record VENDOR_DIR/VENDORED: exit 1 on any difference.")
    parser.add_argument("vendor_dir", metavar="VENDOR_DIR", help="the directory that holds sdlc_protocol/ and VENDORED")
    return vendoring.verify(parser.parse_args(argv).vendor_dir)


COMMANDS = {"board": run_board, "check": run_check, "team": run_team, "workers": run_workers, "guard": run_guard,
            "locks": run_locks, "verify-vendored": run_verify_vendored}


def main(argv=None):
    argv = sys.argv[1:] if argv is None else list(argv)
    if argv[:1] in (["-h"], ["--help"]):
        print(USAGE)
        return 0
    if not argv or argv[0] not in COMMANDS:
        print(USAGE, file=sys.stderr)
        return 2
    return COMMANDS[argv[0]](argv[1:])
