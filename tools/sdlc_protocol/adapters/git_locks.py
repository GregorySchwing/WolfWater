"""The lock table's git side: reading LOCKS.json from the branch `locks`,
writing a new table as one commit on the tip it read, and pushing it, with
the claim branch beside it when the locks are taken at the claim.

Every change is read, decided and pushed again from the start when the
push is refused because the table moved, so the decision is always made
on the tip the push replaces; git refuses a push that is not a fast-forward
of the remote's tip, which makes the push a compare-and-swap. Nothing here
forces a push or deletes a ref: the table only moves forward, and a lock
ends by an entry leaving it (sdlc_protocol.locks).

A claim is live while its branch `claude/<n>` is on the remote and its
head is not yet in the default branch. That is what git alone can tell;
the board's run, which reads the issues, also ends the locks of a claim
whose issue is closed or whose release it counted (live_claims's
`ended`).
"""

import subprocess

from .. import claims, locks
from ..gate import GateFailed
from .bare_git import GIT_TIMEOUT_SECONDS, git_environment, git_options, last_line, run_git

LOCKS_NAME = "symphony locks"
LOCKS_EMAIL = "locks@symphony.local"
# Every path this module hands git is a file's name, never a pattern: a
# file named a[1].py must not match a1.py.
OPTIONS = (*git_options(LOCKS_NAME, LOCKS_EMAIL), "--literal-pathspecs")
TRIES = 5

# Where the liveness reads keep the heads they fetch, apart from the
# clone's own remote-tracking refs, which a lease reads.
PRIVATE = "refs/symphony-locks"

# The refusals git prints for a ref that some other ref of an atomic push
# took down with it, as against the ref that was itself refused: the
# client's own words, and receive-pack's when the server refused another.
SWEPT_ALONG = ("atomic push failed", "atomic push failure", "atomic transaction failed")

# The refusals that mean the remote moved under the push, so a new read
# may succeed; any other refusal, a rule or a hook, will refuse again.
MOVED = ("non-fast-forward", "fetch first", "stale info", "atomic transaction failed")


class LocksFailed(Exception):
    """The table could not be read or moved; the message says why."""


class ClaimTaken(Exception):
    """The claim branch was refused: another worker holds the issue."""


def git(repo, *args):
    try:
        return run_git(list(args), cwd=repo, options=OPTIONS)
    except GateFailed as problem:
        raise LocksFailed(str(problem))


def tracking_ref(remote):
    return f"refs/remotes/{remote}/{locks.LOCKS_BRANCH}"


def fetch_table(repo, remote):
    """The table's tip on the remote and its entries; no tip and no
    entries while the branch `locks` is not there yet."""
    fetched = git(repo, "fetch", "--quiet", "--no-tags", remote,
                  f"+refs/heads/{locks.LOCKS_BRANCH}:{tracking_ref(remote)}")
    if fetched.returncode != 0:
        if "couldn't find remote ref" in fetched.stderr:
            return "", []
        raise LocksFailed(f"could not fetch {locks.LOCKS_BRANCH}: {last_line(fetched)}")
    tip = git(repo, "rev-parse", tracking_ref(remote)).stdout.strip()
    shown = git(repo, "show", f"{tip}:{locks.TABLE_FILE}")
    if shown.returncode != 0:
        raise LocksFailed(f"{locks.LOCKS_BRANCH} at {tip[:12]} holds no {locks.TABLE_FILE}")
    try:
        return tip, locks.read_table(shown.stdout)
    except (ValueError, KeyError, TypeError) as problem:
        raise LocksFailed(f"{locks.TABLE_FILE} at {tip[:12]} cannot be read: {problem}")


def remote_heads(repo, remote, pattern):
    listed = git(repo, "ls-remote", remote, pattern)
    if listed.returncode != 0:
        raise LocksFailed(f"could not list {pattern} on {remote}: {last_line(listed)}")
    heads = {}
    for line in listed.stdout.splitlines():
        sha, _, ref = line.partition("\t")
        heads[ref.removeprefix("refs/heads/")] = sha
    return heads


def live_claims(repo, remote, numbers, default_branch="main", ended=()):
    """Of the claims `numbers`, those still live: the branch on the
    remote, its head not in the default branch, and not among `ended`. A
    head whose ancestry cannot be read, as in a shallow clone, counts as
    live, since a lock dropped wrongly is worse than one kept a run
    longer."""
    wanted = {claims.branch_of(number): number for number in numbers if number not in set(ended)}
    if not wanted:
        return set()
    heads = remote_heads(repo, remote, "refs/heads/claude/*")
    present = {branch: heads[branch] for branch in wanted if branch in heads}
    if not present:
        return set()
    main = f"{PRIVATE}/{remote}/{default_branch}"
    refspecs = [f"+refs/heads/{default_branch}:{main}"]
    refspecs += [f"+refs/heads/{branch}:{PRIVATE}/{remote}/{branch}" for branch in present]
    # A fetch that fails leaves the ancestry below unreadable, and an
    # unreadable head counts as live, so its result needs no reading here.
    git(repo, "fetch", "--quiet", "--no-tags", remote, *refspecs)
    live = set()
    for branch, sha in present.items():
        landed = git(repo, "merge-base", "--is-ancestor", sha, main)
        if landed.returncode != 0:
            live.add(wanted[branch])
    return live


def commit_table(repo, parent, entries, message):
    """A commit whose tree is LOCKS.json alone, on `parent` or as the
    branch's first commit."""
    blob = git_input(repo, ["hash-object", "-w", "--stdin"], locks.table_text(entries))
    tree = git_input(repo, ["mktree"], f"100644 blob {blob}\t{locks.TABLE_FILE}\n")
    parents = ["-p", parent] if parent else []
    return git_input(repo, ["commit-tree", tree, *parents], message)


def git_input(repo, args, text):
    try:
        done = subprocess.run(["git", *OPTIONS, *args], cwd=repo, input=text, capture_output=True, text=True,
                              errors="replace", timeout=GIT_TIMEOUT_SECONDS, env=git_environment())
    except (subprocess.TimeoutExpired, FileNotFoundError) as problem:
        raise LocksFailed(f"git {args[0]} failed: {problem}")
    if done.returncode != 0:
        raise LocksFailed(f"git {args[0]} failed: {last_line(done)}")
    return done.stdout.strip()


def push(repo, remote, refspecs, options=()):
    """Push the refspecs all or none, never forced; each destination ref
    with what git says of refusing it, such as "[rejected] (fetch
    first)", or "" where it was taken. A ref the remote holds already at
    the commit pushed is "[up to date]": git sends nothing for it and
    checks no lease on it, so the caller decides what that means."""
    pushed = git(repo, "push", "--atomic", "--porcelain", *options, remote, *refspecs)
    outcome = {}
    for line in pushed.stdout.splitlines():
        parts = line.split("\t")
        if len(parts) >= 3 and parts[0] in ("*", " ", "+", "-", "!", "="):
            ref = parts[1].partition(":")[2]
            if parts[0] == "!":
                outcome[ref] = parts[2] or "refused"
            else:
                outcome[ref] = "[up to date]" if parts[0] == "=" else ""
    if pushed.returncode != 0 and not outcome:
        raise LocksFailed(f"could not push to {remote}: {last_line(pushed)}")
    return outcome


def change(repo, remote, decide, default_branch="main", claim_branch="", ended=(), tries=TRIES, owner=None):
    """Read the table, let `decide(entries, live)` give the new entries,
    and push them as one commit on the tip read, with `claim_branch`'s
    local head pushed as that branch in the same atomic push when given.
    The claim `owner` counts as live whatever its branch, since its
    entries and its branch arrive together. A refusal of the table alone
    starts over from a new read; a claim branch already on the remote is
    ClaimTaken. `decide` may raise locks.Held, which is passed on with
    nothing pushed."""
    for _ in range(tries):
        if claim_branch and remote_heads(repo, remote, f"refs/heads/{claim_branch}"):
            # Asked before the push, since git pushes the rest of an atomic
            # push beside a branch the remote holds at this very commit.
            raise ClaimTaken(f"{claim_branch} is on {remote} already: another worker holds the issue")
        tip, entries = fetch_table(repo, remote)
        numbers = {entry.claim for entry in entries}
        alive = live_claims(repo, remote, numbers, default_branch, ended) | ({owner} if owner else set())
        after = decide(entries, lambda entry: entry.claim in alive)
        if after == entries and not claim_branch:
            return entries
        # The table is pushed even unchanged beside a claim, so the claim
        # lands on the very table it was decided on.
        sha = commit_table(repo, tip, after, locks.change_message(entries, after))
        refspecs = [f"{sha}:refs/heads/{locks.LOCKS_BRANCH}"]
        options = []
        if claim_branch:
            # The claim is a new branch or nothing: a lease that expects no
            # branch at all, so commits built on another's claim are never
            # pushed onto it.
            options.append(f"--force-with-lease=refs/heads/{claim_branch}:")
            refspecs.append(f"refs/heads/{claim_branch}:refs/heads/{claim_branch}")
        outcome = push(repo, remote, refspecs, options)
        if claim_branch and outcome.get(f"refs/heads/{claim_branch}") == "[up to date]":
            # The remote has this very claim already, so the push made no
            # claim: a worker that checked out another's claim, or a
            # retry of a push whose answer was lost, which reads the same.
            raise ClaimTaken(f"{claim_branch} is on {remote} already at this commit: the push claimed nothing")
        refused = {ref: why for ref, why in outcome.items() if why}
        if not refused:
            return after
        if claim_branch and claim_held(repo, remote, claim_branch, refused):
            raise ClaimTaken(f"{claim_branch} is on {remote} already: another worker holds the issue")
        if not any(word in why for why in refused.values() for word in MOVED):
            raise LocksFailed("the remote refused the push: " + "; ".join(
                f"{ref.removeprefix('refs/heads/')} {why}" for ref, why in sorted(refused.items())))
    raise LocksFailed(f"{locks.LOCKS_BRANCH} kept moving: {tries} tries were refused")


def claim_held(repo, remote, claim_branch, refused):
    """Whether the claim branch was refused because it is there already.
    The wording alone cannot say so, since receive-pack marks every ref of
    an atomic push refused when it refuses one, so where it is not plainly
    the table alone the remote is asked."""
    why = refused.get(f"refs/heads/{claim_branch}", "")
    if not why or swept_along(why) and f"refs/heads/{locks.LOCKS_BRANCH}" in refused \
            and not swept_along(refused[f"refs/heads/{locks.LOCKS_BRANCH}"]):
        return False
    return claim_branch in remote_heads(repo, remote, f"refs/heads/{claim_branch}")


def swept_along(reason):
    return any(words in reason for words in SWEPT_ALONG)


def acquire(repo, remote, wanted, default_branch="main", claim_branch=""):
    """Take `wanted`, all for one claim, all or none. Without the claim
    branch in the same push the claim must be on the remote already, or
    the next change anyone makes would count its locks as ended."""
    owners = {entry.claim for entry in wanted}
    if len(owners) != 1 or min(owners) < 1:
        raise LocksFailed("locks are taken for one claim, an issue's number from 1")
    owner = owners.pop()
    branch = claims.branch_of(owner)
    if claim_branch and claim_branch != branch:
        raise LocksFailed(f"the locks of #{owner} go with {branch}, not {claim_branch}")
    if claim_branch and landed(repo, remote, f"refs/heads/{claim_branch}", default_branch):
        raise LocksFailed(f"{claim_branch} has no commit of its own beyond {default_branch}: commit the claim first")
    if not claim_branch and owner not in live_claims(repo, remote, {owner}, default_branch):
        raise LocksFailed(f"{branch} is not a live claim on {remote}: take the locks with the claim itself")
    return change(repo, remote, lambda entries, live: locks.acquire(entries, wanted, live),
                  default_branch, claim_branch, owner=owner)


def release(repo, remote, claim, default_branch="main"):
    return change(repo, remote, lambda entries, live: locks.release(locks.sweep(entries, live), claim),
                  default_branch)


def sweep(repo, remote, default_branch="main", ended=()):
    return change(repo, remote, lambda entries, live: locks.sweep(entries, live), default_branch, ended=ended)


def landed(repo, remote, ref, default_branch):
    """Whether the local `ref` is already in the remote's default branch,
    so a push of it claims nothing that would stay live."""
    main = f"{PRIVATE}/{remote}/{default_branch}"
    git(repo, "fetch", "--quiet", "--no-tags", remote, f"+refs/heads/{default_branch}:{main}")
    return git(repo, "merge-base", "--is-ancestor", ref, main).returncode == 0


def blob_text(repo, spec):
    """A blob's text with its line endings as git counts lines, or None
    when it is not there: read as bytes, since text mode would turn a lone
    carriage return into a line of its own."""
    try:
        done = subprocess.run(["git", *OPTIONS, "show", spec], cwd=repo, capture_output=True,
                              timeout=GIT_TIMEOUT_SECONDS, env=git_environment(), stdin=subprocess.DEVNULL)
    except (subprocess.TimeoutExpired, FileNotFoundError) as problem:
        raise LocksFailed(f"git show failed: {problem}")
    return done.stdout.decode("utf-8", errors="replace") if done.returncode == 0 else None


# How far a shallow clone is deepened, step by step, to find where a
# branch left the default branch; actions/checkout's default depth is 1.
DEEPENINGS = (50, 500, 5000)


def fork_points(repo, base, head):
    """Every commit where `head` left `base`, git's merge bases
    (`merge-base --all`, so criss-cross history hides nothing behind the
    one git would pick), deepening a shallow clone step by step until git
    finds one."""
    for deepening in (0,) + DEEPENINGS:
        if deepening:
            shallow = git(repo, "rev-parse", "--is-shallow-repository").stdout.strip() == "true"
            if not shallow and resolvable(repo, base) and resolvable(repo, head):
                break
            deepen(repo, deepening, base, head)
        found = git(repo, "merge-base", "--all", base, head)
        if found.returncode == 0 and found.stdout.split():
            return found.stdout.split()
    shallow = git(repo, "rev-parse", "--is-shallow-repository").stdout.strip() == "true"
    why = "the clone is shallow; check out with fetch-depth: 0" if shallow else last_line(found)
    raise LocksFailed(f"{base} and {head} share no history here: {why}")


def deepen(repo, depth, *revs):
    """Deepen each of `revs` by `depth`: a branch on a remote (origin/main)
    by its name, which also brings it where the clone has it not at all,
    since actions/checkout on a pull_request event fetches the test merge
    alone, one commit deep, and no branch; anything else by the commit it
    names, so the test merge's own parent is fetched as the protocol check
    fetches the pull request's own commits. Nothing else is fetched: a
    repository with hundreds of branches would otherwise be brought whole."""
    for rev in revs:
        remote, _, branch = rev.partition("/")
        if branch and git(repo, "remote", "get-url", remote).returncode == 0:
            git(repo, "fetch", "--quiet", "--no-tags", f"--deepen={depth}", remote,
                f"+refs/heads/{branch}:refs/remotes/{remote}/{branch}")
            continue
        sha = git(repo, "rev-parse", "--verify", "--quiet", f"{rev}^{{commit}}").stdout.strip()
        if sha:
            git(repo, "fetch", "--quiet", "--no-tags", f"--deepen={depth}", "origin", sha)


def fork_point(repo, base, head):
    return fork_points(repo, base, head)[0]


def resolvable(repo, rev):
    return bool(git(repo, "rev-parse", "--verify", "--quiet", f"{rev}^{{commit}}").stdout.strip())


def added(repo, base, head):
    """The files the change adds, which no claim could have been editing
    before it."""
    names = git(repo, "diff", "--no-color", "--no-textconv", "--no-renames", "--diff-filter=A", "--name-only", "-z",
                fork_point(repo, base, head), head)
    return sorted(filter(None, names.stdout.split("\0")))


def touched(repo, base, head):
    """The keys the change from the merge base of `base` and `head` to
    `head` touches, file by file; with several merge bases, from each."""
    keys = set()
    for fork in fork_points(repo, base, head):
        keys.update(touched_from(repo, fork, head))
    return sorted(keys)


def touched_from(repo, fork, head):
    names = git(repo, "diff", "--no-color", "--no-textconv", "--no-renames", "--name-only", "-z", fork, head)
    keys = set()
    for path in filter(None, names.stdout.split("\0")):
        base_text = blob_text(repo, f"{fork}:{path}")
        if git(repo, "cat-file", "-e", f"{head}:{path}").returncode != 0:
            base_text = None  # gone: the whole file
        diff = git(repo, "diff", "--no-color", "--no-textconv", "--no-renames", "--no-ext-diff", "-U0", fork, head, "--", path).stdout
        keys.update(locks.touched_keys(path, base_text, locks.hunks_of(diff)))
    return sorted(keys)


def tree_files(repo, rev):
    """The files `rev` holds, for the reads a plan takes from imports."""
    names = git(repo, "ls-tree", "-r", "--name-only", "-z", rev)
    if names.returncode != 0:
        raise LocksFailed(f"cannot list {rev}: {last_line(names)}")
    return [name for name in names.stdout.split("\0") if name]


def plan_at(repo, rev, keys):
    """The plan for `keys` as the files stand at `rev`, not in the work
    tree, since the default branch's side is read where it is."""
    return locks.plan(keys, lambda path: blob_text(repo, f"{rev}:{path}"), tree_files(repo, rev),
                      whole_files=True)


def own_change(repo, base, head):
    """The commit whose change is the branch's own: `head`, or, where
    `head` is a merge whose first parent is `base`'s tip, as GitHub's
    test merge (refs/pull/N/merge, actions/checkout's default on a
    pull_request event) is, its second parent, since judged at the test
    merge the branch has never left the base and nothing ever meets."""
    tip = git(repo, "rev-parse", "--verify", "--quiet", f"{base}^{{commit}}").stdout.strip()
    parents = git(repo, "rev-list", "--parents", "-n", "1", head).stdout.split()[1:]
    if len(parents) == 2 and parents[0] == tip:
        return parents[1], tip
    return git(repo, "rev-parse", "--verify", "--quiet", f"{head}^{{commit}}").stdout.strip(), tip


def overlap(repo, base, head):
    """What the default branch `base` changed since `head`'s fork from it
    that meets what `head` changes or builds on (locks.overlap), or [] when
    the default branch has not moved past the fork at all. The reads of a
    Python file written whole come from the whole file, since such an
    edit builds on whatever the file does."""
    # The base is read first, since the checkout may not have it at all
    # until fork_points has fetched it, and own_change tells the test
    # merge by the base's tip; a merge made against a tip the base has
    # since moved past is judged as itself, which gives what the base
    # did since the merge was made against the merge's own change.
    forks = fork_points(repo, base, head)
    own, tip = own_change(repo, base, head)
    if own != head:
        head, forks = own, fork_points(repo, base, own)
    if tip in forks:
        return []
    moved = plan_at(repo, base, touched(repo, head, base))
    change = plan_at(repo, head, touched(repo, base, head))
    return locks.overlap(moved, change)
