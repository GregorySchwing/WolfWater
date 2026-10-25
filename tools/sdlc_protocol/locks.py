"""Locks on files, on parts of files and on named resources, kept in git.

The claim locks an issue: the first push of `claude/<n>` wins. This module
locks what a claim works on, so two claims never edit one function, one
section or one Claude Doc at once, the rule a coordinator otherwise applies
by hand by naming one owner per resource. It is off unless the switch says
otherwise (mode_of), and nothing else in the package reads it yet.

A lock is an entry in one table, LOCKS.json on the branch `locks`. The
branch moves only by fast-forward pushes, so a push built on the tip a
worker read is accepted once and refused to anyone who read an older tip:
a compare-and-swap that git and GitHub decide, with no server and no
deletion. adapters/git_locks.py makes those commits and pushes; this
module holds the rules, as plain functions over what they are handed.

A key names what is locked:
  path          the whole file
  path#name     one symbol of it: a Python function or class, or
                Class.method; a Markdown section by its heading, a
                subsection as "Heading > Sub"
  path#<top>    everything in the file outside its named symbols
  scheme:id     a resource outside git, such as doc:a1cb9900, always
                whole: what follows a # in it is not read

Two keys meet when they name one resource and either is the whole of it,
or they name one symbol, or one symbol holds the other. Two entries
conflict when their keys meet, at least one is a write, and they belong to
different claims. A read lock is how a claim says what its edits depend
on, so a change to what another claim is building on waits.
"""

import ast
import json
import re
from dataclasses import asdict, dataclass

LOCKS_BRANCH = "locks"
TABLE_FILE = "LOCKS.json"
TABLE_VERSION = 1

WRITE = "write"
READ = "read"
MODES = (WRITE, READ)

TOP = "<top>"

OFF = "off"
ADVISE = "advise"
ENFORCE = "enforce"
SWITCH_ENV = "SYMPHONY_LOCKS"
SWITCH_FILE = ".github/symphony-locks"

# The kinds of file this module reads symbols of, and the separator by
# which one symbol holds another in each: Class.method in Python,
# "Heading > Sub" in Markdown, where a heading may itself hold a dot.
PYTHON = (".py",)
MARKDOWN = (".md", ".markdown")
SUFFIXES = PYTHON + MARKDOWN

SCHEME = re.compile(r"[a-z][a-z0-9-]*:[^/]")
HEADING = re.compile(r" {0,3}(#{1,6})(?:[ \t]+(.*?))?(?:[ \t]+#+)?[ \t]*$")
LONE_CR = re.compile(r"\r(?!\n)")
FENCE = re.compile(r" {0,3}(`{3,}|~{3,})")


class Held(Exception):
    """What a claim asked for is held by another: `conflicts` pairs each
    wanted entry with the entry that holds it."""

    def __init__(self, conflicts):
        self.conflicts = conflicts
        super().__init__("; ".join(f"{want.key} ({want.mode}) is held by #{held.claim} "
                                   f"({held.worker or 'unnamed'}, {held.mode})" for want, held in conflicts))


@dataclass(frozen=True, order=True)
class Entry:
    """One lock: the key, its mode, the claim's issue number and the
    worker. It belongs to the issue's claim and lives as long as that
    claim does (adapters/git_locks.live_claims), whatever its commits,
    since a worker may rebase its own claim branch."""
    key: str
    claim: int
    mode: str = WRITE
    worker: str = ""


@dataclass(frozen=True)
class Symbol:
    """A named region of a file, lines first to last, counted from 1."""
    name: str
    first: int
    last: int


def mode_of(environ=None, switch_text=None):
    """The switch: the environment's SYMPHONY_LOCKS first, else the text
    of .github/symphony-locks, else off. Anything unreadable is off, so a
    typo never turns enforcement on."""
    for value in ((environ or {}).get(SWITCH_ENV), switch_text):
        word = (value or "").strip().lower()
        if word in (OFF, ADVISE, ENFORCE):
            return word
    return OFF


def split_key(key):
    """A key's resource and its symbol, or None for the whole resource.
    The symbol follows the first # that ends a name this module reads
    symbols of, so a path that holds a # itself, notes#1.md, stays whole."""
    if is_named_resource(key):
        return key.partition("#")[0], None
    for index, mark in enumerate(key):
        if mark == "#" and key[:index].endswith(SUFFIXES):
            return key[:index], key[index + 1:] or None
    return key, None


def is_named_resource(resource):
    """A resource outside git, such as doc:a1cb9900, as against a path."""
    return bool(SCHEME.match(resource))


def holds(outer, inner, resource=""):
    """Whether symbol `outer` is `inner` or holds it, by the nesting of
    the resource's kind of file."""
    mark = " > " if resource.endswith(MARKDOWN) else "."
    return outer == inner or inner.startswith(outer + mark)


def keys_meet(first, second):
    """Whether two keys name overlapping parts of one resource."""
    resource, symbol = split_key(first)
    other_resource, other_symbol = split_key(second)
    if resource != other_resource:
        return False
    if symbol is None or other_symbol is None:
        return True
    return holds(symbol, other_symbol, resource) or holds(other_symbol, symbol, resource)


def covers(granted, touched):
    """Whether a lock on `granted` lets its holder change `touched`."""
    resource, symbol = split_key(granted)
    other_resource, other_symbol = split_key(touched)
    if resource != other_resource:
        return False
    if symbol is None:
        return True
    return other_symbol is not None and holds(symbol, other_symbol, resource)


def conflict(held, wanted):
    """Whether `held` stands in the way of `wanted`."""
    if held.claim == wanted.claim:
        return False
    if held.mode == READ and wanted.mode == READ:
        return False
    return keys_meet(held.key, wanted.key)


def conflicts(entries, wanted):
    """Each wanted entry paired with each held entry in its way."""
    return [(want, held) for want in wanted for held in entries if conflict(held, want)]


def overlap(moved, change):
    """Where what the default branch did since a change's base meets what
    the change does: each pair of an entry of `moved`, the plan of the
    default branch's own edits since the base, and an entry of `change`,
    the change's plan, whose keys meet with a write on either side. Two
    reads never meet, so a change that only builds on what main also
    builds on needs nothing from main; any other meeting is one the
    change's own check never saw, and is why the change brings main in."""
    return [(ours, theirs) for theirs in moved for ours in change
            if WRITE in (ours.mode, theirs.mode) and keys_meet(ours.key, theirs.key)]


def acquire(entries, wanted, live=None):
    """The table once `wanted` is granted, with the entries `live` says
    are dead dropped first; Held when any wanted entry is in the way of a
    live one. All or nothing: a refusal grants none of `wanted`."""
    kept = sweep(entries, live)
    standing = conflicts(kept, wanted)
    if standing:
        raise Held(standing)
    return merged(kept, wanted)


def merged(entries, wanted):
    """The entries with `wanted` added, a claim's write on a key taking
    the place of its read on the same key, and no entry twice."""
    table = {(entry.key, entry.claim): entry for entry in entries}
    for want in wanted:
        before = table.get((want.key, want.claim))
        if before is None or before.mode == READ:
            table[(want.key, want.claim)] = want
    return sorted(table.values())


def release(entries, claim):
    """The table without the claim's entries."""
    return sorted(entry for entry in entries if entry.claim != claim)


def sweep(entries, live=None):
    """The table without the entries whose claim has ended."""
    return sorted(entry for entry in entries if live is None or live(entry))


def held_by(entries, claim):
    return sorted(entry for entry in entries if entry.claim == claim)


def uncovered(touched, entries, claim):
    """The touched keys no write lock of the claim covers."""
    writes = [entry.key for entry in entries if entry.claim == claim and entry.mode == WRITE]
    return sorted({key for key in touched if not any(covers(granted, key) for granted in writes)})


def table_text(entries):
    """LOCKS.json's text: sorted, one entry a line, so a diff of the
    branch reads as what changed."""
    lines = ",\n".join("    " + json.dumps(asdict(entry), sort_keys=True) for entry in sorted(entries))
    body = f"[\n{lines}\n  ]" if entries else "[]"
    return f'{{\n  "version": {TABLE_VERSION},\n  "entries": {body}\n}}\n'


def read_table(text):
    """The entries of LOCKS.json's text; ValueError for a table this
    version cannot read, since a lock misread is worse than none."""
    if not (text or "").strip():
        return []
    data = json.loads(text)
    if not isinstance(data, dict) or data.get("version") != TABLE_VERSION:
        raise ValueError(f"{TABLE_FILE} is not version {TABLE_VERSION}")
    entries = []
    for item in data.get("entries", []):
        entry = Entry(key=str(item["key"]), claim=int(item["claim"]), mode=str(item.get("mode", WRITE)),
                      worker=str(item.get("worker", "")))
        if entry.mode not in MODES:
            raise ValueError(f"{TABLE_FILE}: {entry.key} has the mode {entry.mode!r}")
        entries.append(entry)
    return sorted(entries)


def change_message(before, after):
    """The commit message for one change of the table: who took or gave
    up what, one line each, under a summary line."""
    taken = sorted(set(after) - set(before))
    given = sorted(set(before) - set(after))
    lines = [f"#{entry.claim} ({entry.worker or 'unnamed'}) takes {entry.key} {entry.mode}" for entry in taken]
    lines += [f"#{entry.claim} gives up {entry.key} {entry.mode}" for entry in given]
    claims = sorted({entry.claim for entry in taken + given})
    summary = "Locks: " + (", ".join(f"#{claim}" for claim in claims) if claims else "no change")
    return summary + ("\n\n" + "\n".join(lines) if lines else "") + "\n"


def symbols_of(path, text):
    """The named regions of a file by its kind, or None where this module
    reads no symbols, so the file is locked whole."""
    if path.endswith(PYTHON):
        return python_symbols(text)
    if path.endswith(MARKDOWN):
        return markdown_symbols(text)
    return None


def python_symbols(text):
    """Top-level functions and classes, and the functions of each class,
    each from its first decorator to its last line; None when the file
    does not parse, so it is locked whole."""
    try:
        tree = ast.parse(git_lines(text))
    except (SyntaxError, ValueError):
        return None
    found = []
    for node in tree.body:
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)):
            found.append(symbol_of(node, node.name))
            if isinstance(node, ast.ClassDef):
                found += [symbol_of(inner, f"{node.name}.{inner.name}") for inner in node.body
                          if isinstance(inner, (ast.FunctionDef, ast.AsyncFunctionDef))]
    return found


def git_lines(text):
    """The text with each lone carriage return made a space, so ast counts
    lines as git does, at LF alone."""
    return LONE_CR.sub(" ", text)


def symbol_of(node, name):
    first = min([node.lineno] + [decorator.lineno for decorator in node.decorator_list])
    return Symbol(name, first, node.end_lineno)


def markdown_symbols(text):
    """Each section, from its heading to the line before the next heading
    of the same or a higher level, named by its headings from the top as
    "Heading > Sub". A heading inside a fenced block is no heading, and a
    fence closes only on the same mark at least as long; a setext heading,
    underlined, is read as text of the section above it. Lines are split
    as git splits them, at LF alone, so a line's number is git's."""
    lines = text.split("\n")
    if lines and lines[-1] == "":
        lines.pop()
    headings = []
    fence = ""
    for number, line in enumerate(lines, 1):
        marks = FENCE.match(line)
        if marks and not fence:
            fence = marks.group(1)
            continue
        if fence:
            if marks and marks.group(1)[0] == fence[0] and len(marks.group(1)) >= len(fence) \
                    and not line.strip().strip(fence[0]):
                fence = ""
            continue
        found = HEADING.match(line.rstrip("\r"))
        if found and found.group(2):
            headings.append((number, len(found.group(1)), found.group(2)))
    found = []
    path = []
    for index, (number, level, title) in enumerate(headings):
        path = [(lvl, name) for lvl, name in path if lvl < level] + [(level, title)]
        last = len(lines)
        for later_number, later_level, _ in headings[index + 1:]:
            if later_level <= level:
                last = later_number - 1
                break
        found.append(Symbol(" > ".join(name for _, name in path), number, last))
    return found


def innermost(symbols, line):
    """The narrowest symbol holding the line, else TOP."""
    holding = [symbol for symbol in symbols if symbol.first <= line <= symbol.last]
    return min(holding, key=lambda symbol: symbol.last - symbol.first).name if holding else TOP


def inserted_into(symbols, after):
    """The regions a line inserted after line `after` touches: the
    narrowest symbol it lands strictly inside; or, outside every symbol,
    TOP with the symbol that ends just before it and the one that starts
    just after it, since a three-way merge conflicts on edits that only
    abut, as an edit of a function's last line and a line added below it
    do."""
    holding = [symbol for symbol in symbols if symbol.first <= after and after + 1 <= symbol.last]
    if holding:
        return {min(holding, key=lambda symbol: symbol.last - symbol.first).name}
    names = {TOP}
    for edge in (after, after + 1):
        name = innermost(symbols, edge)
        if name != TOP:
            names.add(name)
    return names


def touched_keys(path, base_text, hunks):
    """The keys a change to one file touches, from its hunks on the base
    version, each (first line, count) as `git diff -U0` gives them; a
    count of 0 is an insertion after that line. A file new or gone, or
    one this module reads no symbols of, is touched whole."""
    if base_text is None or not hunks:
        return [path]  # new, gone, or changed with no lines to read: binary, or its mode alone
    symbols = symbols_of(path, base_text)
    if symbols is None:
        return [path]
    names = set()
    for first, count in hunks:
        if count == 0:
            names.update(inserted_into(symbols, first))
        else:
            names.update(innermost(symbols, line) for line in range(first, first + count))
    return sorted(f"{path}#{name}" for name in names)


HUNK = re.compile(r"^@@ -(\d+)(?:,(\d+))? \+\d+(?:,\d+)? @@")


def hunks_of(diff_text):
    """The base-side hunks of one file's `git diff -U0`."""
    found = []
    for line in diff_text.splitlines():
        match = HUNK.match(line)
        if match:
            found.append((int(match.group(1)), int(match.group(2) if match.group(2) is not None else 1)))
    return found


def python_reads(path, text, written, repository_files=()):
    """The keys a Python symbol's code depends on: the file's other
    symbols it names, a method's calls through self to its own class's
    methods, and the names it takes from the repository's own modules by
    a module-level `from ... import`. It sees direct references only; a
    dynamic call or a dependency through data is the required check's to
    catch."""
    symbols = python_symbols(text)
    if not symbols:
        return []
    tree = ast.parse(git_lines(text))
    node = find_node(tree, written)
    if node is None:
        return []
    own = {symbol.name for symbol in symbols}
    owner = written.split(".")[0] if "." in written else None
    imported = imported_names(tree, path, set(repository_files))
    reads = set()
    for inner in ast.walk(node):
        if isinstance(inner, ast.Name):
            if inner.id in own:
                reads.add(f"{path}#{inner.id}")
            elif inner.id in imported:
                reads.add(imported[inner.id])
        elif (isinstance(inner, ast.Attribute) and owner and isinstance(inner.value, ast.Name)
              and inner.value.id in ("self", "cls") and f"{owner}.{inner.attr}" in own):
            reads.add(f"{path}#{owner}.{inner.attr}")
    mine = f"{path}#{written}"
    return sorted(key for key in reads if not keys_meet(key, mine))


def python_file_reads(path, text, repository_files=()):
    """The keys a whole Python file depends on: what every symbol of it
    reads (python_reads), and what its module-level code names from the
    repository's own modules, for a file written whole or at `<top>`,
    where no one symbol says what the edit builds on. A file that does not
    parse reads nothing that can be seen."""
    try:
        tree = ast.parse(git_lines(text))
    except (SyntaxError, ValueError):
        return []
    reads = set()
    for symbol in python_symbols(text):
        reads.update(python_reads(path, text, symbol.name, repository_files))
    imported = imported_names(tree, path, set(repository_files))
    own = {symbol.name for symbol in python_symbols(text)}
    for node in tree.body:
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)):
            continue
        for inner in ast.walk(node):
            if not isinstance(inner, ast.Name):
                continue
            if inner.id in own:
                reads.add(f"{path}#{inner.id}")
            elif inner.id in imported:
                reads.add(imported[inner.id])
    return sorted(reads)


def find_node(tree, name):
    head, _, tail = name.partition(".")
    for node in tree.body:
        if getattr(node, "name", None) == head and isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef,
                                                                       ast.ClassDef)):
            if not tail:
                return node
            for inner in getattr(node, "body", []):
                if getattr(inner, "name", None) == tail:
                    return inner
    return None


def imported_names(tree, path, repository_files):
    """Each name a module-level `from X import Y` binds, mapped to the key
    of Y in the repository's file for X, where that file is in it."""
    package = path.rpartition("/")[0]
    names = {}
    for node in tree.body:
        if not isinstance(node, ast.ImportFrom) or (not node.module and not node.level):
            continue
        base = package
        for _ in range(max(node.level - 1, 0)):
            base = base.rpartition("/")[0]
        if not node.module:
            # from . import mod: each name is a module of the package,
            # read whole where the repository holds it.
            for alias in node.names:
                candidate = "/".join(part for part in (base, f"{alias.name}.py") if part)
                if candidate in repository_files:
                    names[alias.asname or alias.name] = candidate
            continue
        module = node.module.replace(".", "/")
        stem = "/".join(part for part in ((base if node.level else ""), module) if part)
        for candidate in (f"{stem}.py", f"{stem}/__init__.py"):
            if candidate in repository_files:
                for alias in node.names:
                    if alias.name != "*":
                        names[alias.asname or alias.name] = f"{candidate}#{alias.name}"
                break
    return names


def plan(targets, read_file, repository_files=(), claim=0, worker="", whole_files=False):
    """The entries a claim asks for to change `targets`: a write on each,
    and for each Python symbol a read on what its code depends on that the
    claim does not write itself. `read_file(path)` gives a file's text, or
    None when it is not there. With `whole_files`, a Python file written
    whole or at `<top>` reads what the whole file depends on too, which an
    overlap reading needs, since an edit there builds on whatever the file
    does, while a claim taking a whole file asks for no more than it
    names."""
    writes = sorted(set(targets))
    reads = set()
    for key in writes:
        resource, symbol = split_key(key)
        if not resource.endswith(PYTHON):
            continue
        text = read_file(resource)
        if text is None:
            continue
        if symbol and symbol != TOP:
            reads.update(python_reads(resource, text, symbol, repository_files))
        elif whole_files:
            reads.update(python_file_reads(resource, text, repository_files))
    reads = {key for key in reads if not any(covers(write, key) for write in writes)}
    return ([Entry(key, claim, WRITE, worker) for key in writes]
            + [Entry(key, claim, READ, worker) for key in sorted(reads)])
