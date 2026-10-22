"""Hold a claude/<issue> pull request to the claim protocol.

The protocol says what a claim looks like; this module makes the checkable
parts of it a list of problems, for CI to fail on, so a pull request that
forgot them cannot merge however the session that made it was feeling:

- A claim branch belongs to the repository, never a fork.
- The body closes its own issue, so the merge closes it.
- A draft still carries docs/agents/claims/<issue>.md: the claim stays
  visible exactly as long as it holds.
- A pull request marked ready has deleted it: finishing lets go.
- A head branch that reads as a claim, claude/ and digits, is the
  issue's claim branch, its number as GitHub writes it; claude/03 claims
  nothing, while its author meant a claim on #3.
- No pull request adds or changes, against its base, a claim file other
  than its own issue's on its own claim branch: a claim file of another
  issue that reached the default branch would read as that issue's claim
  (Pick) and make the board and the push guard read a later claim's
  branch, once it brought the default branch in, as that file's worker's
  (F31). Deleting a claim file passes, so a stale one on the default
  branch can always be taken off.
- Nor does one make the claim folder, or a folder above it, anything but
  a folder, make a claim file a link, add a path below a folder in the
  claim folder, or add a path that is the claim folder's or a claim
  file's but for its letters' case: each would have a checkout, or git
  show, give claim files git does not hold (claim_folder_problems).

Any other head branch, dependabot's or a person's, is not a claim and
passes but for the claim files it adds or changes. The pull request is
the object GitHub's event carries; whether the claim file exists in the
checkout, and which claim files the change adds or changes, are the
caller's to say, so nothing here reads a file; a caller that names no
changed claim files (problems_for's two-argument form, as a consumer
written before F31 calls it) is judged on the head's own claim file alone.
"""

import re

from .claims import CLAIM_DIR, branch_of, claim_path, key_of, meant_number, number_in, number_of, reads_as_claim

# The forms GitHub links a pull request to an issue with; the protocol
# writes "Closes #N" and the rest are the same promise.
CLOSES = r"(?i)\b(?:close[sd]?|fix(?:e[sd])?|resolve[sd]?)\s+#%d\b"
PROTOCOL = "docs/agents/PROTOCOL.md"


def head_branch(pull):
    return ((pull or {}).get("head") or {}).get("ref") or ""


def claim_number(pull):
    """The issue a pull request claims, from its head branch; None when the
    branch is no claim."""
    return number_of(head_branch(pull))


def repository_of(side):
    """The full name of the repository a pull request's head or base is in."""
    return ((side or {}).get("repo") or {}).get("full_name")


def closes_its_issue(pull, number):
    return re.search(CLOSES % number, pull.get("body") or "") is not None


def claim_file_number(path):
    """The issue whose claim file a repository path is,
    docs/agents/claims/<issue>.md with the issue's number as claim_path
    writes it, or None for any other path, the folder's README.md among
    them."""
    folder, _, name = str(path).rpartition("/")
    if folder != CLAIM_DIR or not name.endswith(".md"):
        return None
    return number_in(name[:-len(".md")])


# The folders a claim file's path runs through, the claim folder last: a
# change that makes one of them anything but a folder, a symbolic link
# above all, makes the files the default branch shows there other than
# the ones git holds (F31's review).
CLAIM_FOLDERS = tuple("/".join(CLAIM_DIR.split("/")[:n]) for n in range(1, CLAIM_DIR.count("/") + 2))
# The modes of a file git holds as its own text: a claim file is one of
# them, never a symbolic link (120000) or a submodule (160000).
FILE_MODES = ("100644", "100755")


def raw_changes(raw):
    """Each entry of `git diff --raw --no-renames -z` output (bytes or
    text) as (status, old mode, new mode, path): a colon, the two modes,
    the two object names and the status, then the path, each ended by a
    NUL. An entry that does not hold together is left out."""
    text = raw.decode("utf-8", "surrogateescape") if isinstance(raw, bytes) else raw
    fields = text.split("\0")
    found = []
    for meta, path in zip(fields[0::2], fields[1::2]):
        parts = meta[1:].split() if meta.startswith(":") else []
        if len(parts) == 5 and path:
            found.append((parts[4], parts[0], parts[1], path))
    return found


def claim_files_changed(raw):
    """The issues whose claim files a diff adds or changes, in number
    order, from `git diff --raw --no-renames -z` output (raw_changes). A
    deleted file is left out, since taking a claim file off frees what it
    held."""
    return sorted({claim_file_number(path) for status, _, _, path in raw_changes(raw)
                   if not status.startswith("D") and claim_file_number(path) is not None})


def claim_folder_problems(raw):
    """The problems of a diff that changes what the claim files are rather
    than which (F31's reviews), each but a deletion, which passes as
    deleting a claim file does:
    - the claim folder, or a folder above it, made a file, a symbolic
      link or a submodule, so the default branch's checkout would show
      claim files git does not hold;
    - a claim file made a symbolic link or a submodule;
    - a path more than one level below the claim folder, since a folder
      there named as a claim file is (claims/9.md/<anything>) is read by
      `git show` as a listing that names a worker, and no claim file
      lies deeper than the folder itself;
    - a name in the claim folder that is a claim file's but for the case
      of its letters (9.MD), and a path that is a claim folder's, or lies
      in it, but for the case of its letters (Docs/Agents/claims), each
      of which a checkout on a file system that folds case reads as the
      claim file or the claim folder; case is folded as Python's
      str.upper and str.casefold together fold it, so both the dotless i
      and the long s fold as well."""
    problems = []
    for status, _, mode, path in raw_changes(raw):
        if status.startswith("D"):
            continue
        if path in CLAIM_FOLDERS:
            problems.append(f"{path} is made something other than a folder here (mode {mode}), so the claim files "
                            f"the default branch shows would not be the ones git holds ({PROTOCOL}, Rules)")
        elif claim_file_number(path) is not None and mode not in FILE_MODES:
            problems.append(f"{path} is made a link or a submodule here (mode {mode}), and a claim file is a file "
                            f"git holds ({PROTOCOL}, Rules)")
        elif path.startswith(CLAIM_DIR + "/") and "/" in path[len(CLAIM_DIR) + 1:]:
            problems.append(f"{path} lies below a folder in {CLAIM_DIR}, where only claim files and the README "
                            f"belong; a folder named as a claim file reads as one ({PROTOCOL}, Rules)")
        elif folds_into_claims(path):
            problems.append(f"{path} is a claim folder's or a claim file's path but for the case of its letters, "
                            f"which a checkout that folds case reads as one ({PROTOCOL}, Rules)")
    return problems


def folded(text):
    """A path as a file system that folds case compares it."""
    return text.upper().casefold()


def folds_into_claims(path):
    """Whether a path is a claim folder's, lies in the claim folder, or
    names a claim file there, once case is folded, and is not, as it is
    written."""
    own = path in CLAIM_FOLDERS or path.startswith(CLAIM_DIR + "/")
    folder = folded(CLAIM_DIR) + "/"
    inside = folded(path) in {folded(f) for f in CLAIM_FOLDERS} or folded(path).startswith(folder)
    if not own:
        return inside
    name = path[len(CLAIM_DIR) + 1:]
    return claim_file_number(path) is None and claim_file_number(CLAIM_DIR + "/" + folded(name)) is not None


def stray_claim_files(pull, changed_claims):
    """The problems of the claim files a pull request adds or changes that
    are not its own: every one, where its head is no claim branch here,
    and every one but its own issue's, where it is (F31)."""
    number = claim_number(pull)
    own = number if number is not None and repository_of(pull.get("head")) == repository_of(pull.get("base")) else None
    return [f"{claim_path(n)} is added or changed here, and a claim file is only its own claim branch's: once "
            f"merged it would read as a claim on #{n} on the default branch ({PROTOCOL}, Rules)"
            for n in sorted(set(changed_claims)) if n != own]


def problems_for(pull, claim_file_exists, changed_claims=()):
    """Every way this pull request breaks the protocol; [] when it keeps
    the rules. claim_file_exists says whether docs/agents/claims/<issue>.md
    is in the checkout: a bool, or a callable given the issue number, for
    a caller that would rather not parse the branch itself.
    changed_claims are the issues whose claim files the pull request adds
    or changes against its base (claim_files_changed), each of which but
    its own claim's is a problem (stray_claim_files)."""
    branch = head_branch(pull)
    number = claim_number(pull)
    strays = stray_claim_files(pull, changed_claims)
    if number is None:
        return ([misnamed(branch)] if reads_as_claim(branch) else []) + strays
    if callable(claim_file_exists):
        claim_file_exists = claim_file_exists(number)
    problems = []
    if repository_of(pull.get("head")) != repository_of(pull.get("base")):
        problems.append(f"{branch}: a claim branch lives in this repository, not in a fork "
                        f"({PROTOCOL}, Claim)")
    if not closes_its_issue(pull, number):
        problems.append(f"the body does not close #{number}, so the merge would leave the issue open; "
                        f'write "Closes #{number}" ({PROTOCOL}, Finish)')
    path = claim_path(number)
    if pull.get("draft"):
        if not claim_file_exists:
            problems.append(f"{path} is missing: a draft carries its claim file ({PROTOCOL}, Claim)")
    elif claim_file_exists:
        problems.append(f"{path} is still here: a pull request marked ready has deleted its claim file "
                        f"({PROTOCOL}, Finish)")
    return problems + strays


def misnamed(branch):
    """The problem with a head branch that reads as a claim and claims
    nothing (claims.reads_as_claim): the issue it seems to mean has its own
    claim branch, which is where the work goes."""
    meant = meant_number(key_of(branch))
    there = f"; the claim branch of #{meant} is {branch_of(meant)}, so push the work there" if meant is not None else ""
    return (f"{branch} reads as a claim branch and claims nothing, since a claim branch is claude/ and the "
            f"issue's number in one to nine ASCII digits with no leading zero{there} ({PROTOCOL}, Claim)")
