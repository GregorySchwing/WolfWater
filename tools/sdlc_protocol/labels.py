"""The labels the claim protocol reads and writes.

GitHub Issues is the queue, and labels are how an issue says where it
stands: todo in the queue, claimed by a branch, in-review with a pull
request waiting, blocker stopping everyone. One more label says whose an
issue is: worker:<name>, the worker that holds it, or that the
coordinator gave it before any claim. A worker is anything that claims
issues, a person's session, a routine on someone's Claude account or a run
in GitHub Actions, and its name is the one on its claim files. The two
older spellings, account:<name> and agent:<name>, are read as the same
label, so an issue labelled before workers had one name still says whose
it is. The board, the protocol check, the coordinator and the clerk all
read the same names, so they are spelled once, here.
"""

TODO = "todo"
CLAIMED = "claimed"
IN_REVIEW = "in-review"
BLOCKER = "blocker"
EPIC = "epic"
BLOCKED = "blocked"
# A goal is a request a person breaks into epics, each a sub-issue of the
# goal, which the coordinator then plans: the board leaves it out of its
# Queue as it leaves an epic out and draws it above its epics. Only a
# person puts it on, so it is not among the labels a repository makes up
# front (LABELS), and a repository whose label goal means something else
# renames that label first, since an issue carrying it is dealt to nobody.
GOAL = "goal"

# A worker's label: the name after the colon is the worker that holds the
# issue, or that the coordinator gave it. The older spellings, an
# account's label and a team member's, are read as the same label.
WORKER_PREFIX = "worker:"
ACCOUNT_PREFIX = "account:"
AGENT_PREFIX = "agent:"
HOLDER_PREFIXES = (WORKER_PREFIX, ACCOUNT_PREFIX, AGENT_PREFIX)
WORKER_COLOUR = "c5def5"
# The name the main account's label carried before workers had their own;
# kept for the readers of old labels and claim files.
MAIN_ACCOUNT = "main"

# The labels a repository needs for GitHub Issues to be the queue, and the
# two the coordinator's planning reads: name, colour (six hex digits, as the
# API takes them), description. A repository is set up by making every one
# of these; a worker's own label is made when the worker first needs it.
LABELS = (
    (TODO, "0e8a16", "In the queue, open for a claim"),
    (CLAIMED, "fbca04", "A claude/<issue> branch holds it"),
    (IN_REVIEW, "1d76db", "A pull request is open and waiting"),
    (BLOCKER, "b60205", "Stops everyone until it is settled"),
    (EPIC, "8b5cf6", "A request the coordinator plans into tasks (docs/agents/WORKERS.md)"),
    (BLOCKED, "e99695", "A worker waits on a peer or the coordinator (docs/agents/WORKERS.md)"),
)


def label_names(issue):
    """The names of an issue's labels, sorted, whether the API gave label
    objects or a fixture gave plain strings. Sorted, so a board rendered
    from the same issue twice shows the labels in the same order."""
    return sorted(label if isinstance(label, str) else label["name"]
                  for label in (issue or {}).get("labels") or [])


def name_after(label, prefix):
    """The name a label carries after the prefix, or "" when it does not
    start with the prefix: a bare prefix, worker: with nothing after the
    colon, names nobody."""
    return label[len(prefix):] if label.startswith(prefix) else ""


def names_with_prefix(labels, prefix):
    """What follows the prefix in every label that carries it with a name
    after it, in order (name_after)."""
    return [name_after(label, prefix) for label in labels if name_after(label, prefix)]


def worker_label(name):
    """The label of the worker with this name."""
    return WORKER_PREFIX + name


def holder_names(names):
    """The workers these label names name, in order, under any of the three
    spellings, each once."""
    found = []
    for name in names:
        for prefix in HOLDER_PREFIXES:
            if name_after(name, prefix) and name_after(name, prefix) not in found:
                found.append(name_after(name, prefix))
    return found


def is_holder_label(name):
    """Whether a label names a worker: a prefix of any spelling with a name
    after it. A bare prefix names nobody, as holder_names reads it, so a
    release leaves it on (F15's fourth review: until then a release took
    worker: off as a worker's label, which the board read as no worker's)."""
    return any(name_after(name, prefix) for prefix in HOLDER_PREFIXES)


def holder_labels(issue):
    """The labels on an issue that name a worker, under any spelling: what
    a release takes off with the label claimed."""
    return [name for name in label_names(issue) if is_holder_label(name)]


def worker_of(issue):
    """The worker an issue's label names, worker:<name> before the older
    spellings, or None when no label does."""
    names = label_names(issue)
    for prefix in HOLDER_PREFIXES:
        found = names_with_prefix(names, prefix)
        if found:
            return found[0]
    return None


def account_label(issue):
    """The account an issue's label names (account:main gives main), or
    None when no label does."""
    accounts = names_with_prefix(label_names(issue), ACCOUNT_PREFIX)
    return accounts[0] if accounts else None


def agent_label(issue):
    """The team member an issue's label names (agent:w1 gives w1), or None."""
    agents = names_with_prefix(label_names(issue), AGENT_PREFIX)
    return agents[0] if agents else None


def holders(issue):
    """Every worker an issue's labels name, under any spelling, as a set."""
    return set(holder_names(label_names(issue)))


def planned(names):
    """Whether these label names mark an issue the coordinator plans and
    nobody takes whole: an epic, planned into tasks, or a goal, planned
    into epics."""
    return EPIC in names or GOAL in names


def reserved(labels):
    """Whether these labels keep an issue off the open queue: an epic or a
    goal, which the coordinator plans (planned), or an issue a worker's
    label names, which that worker holds or was given and only that worker
    takes."""
    return planned(labels) or bool(holder_names(labels))


def blocked_by_counts(issue):
    """How many issues GitHub's dependencies say block this one, open and
    in all, from the dependency summary GitHub carries on every issue it
    lists: (0, 0) where the summary is missing or a count is not a whole
    number, since an issue whose summary GitHub does not give is read as
    blocked by nothing (spec/PROTOCOL.md, Pick)."""
    summary = (issue or {}).get("issue_dependencies_summary")
    if not isinstance(summary, dict):
        return (0, 0)
    return (whole_count(summary.get("blocked_by")), whole_count(summary.get("total_blocked_by")))


def whole_count(value):
    """A count as the summary gives it, or 0 where it is not a whole number
    of zero or more; a bool is not a count, though Python counts it."""
    return value if isinstance(value, int) and not isinstance(value, bool) and value >= 0 else 0


# The name team mode gave the same rule, kept for its callers.
reserved_for_team = reserved
