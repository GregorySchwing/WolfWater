"""The workers a repository runs, kept as files: one small JSON file per
worker under .github/symphony/workers/, named for the worker.

A worker is anything that claims issues under the protocol: a person's
sessions, local or in the cloud, with whatever helpers they start; a
routine on someone's Claude account, which runs on its own schedule; or
runs in GitHub Actions on the repository's own Claude token. All of them
follow the same protocol, and the board, the coordinator and the clerk
treat them alike. Each has a name, unique in the repository, which is its
label (worker:<name>) and the name on its claim files' Worker line, so two
workers on one Claude account still tell their claims apart.

A worker's file says where it runs, what it is for, the area it prefers,
how its pull requests land and, for a worker in Actions, the model, the
turns one run may take and how many runs it may take in a row. A person
needs no file: their sessions take their GitHub login as their name. A
routine or a worker in Actions has one, made by hand or from the answers
to the repository's "Add a worker" issue form (from_form). The file holds
no credential, only a note of whose Claude account runs the worker. Nothing
here calls GitHub; load is the one function that reads files.
"""

import json
import re
from dataclasses import dataclass, field
from pathlib import Path

from . import labels

REGISTRY = ".github/symphony/workers"
RUNS = ("actions", "routine", "session")
LANDS = ("owner", "auto-merge")
# A name is lower case letters, digits and inner dashes, as a label, a
# file name and a branch-safe word all allow.
NAME = re.compile(r"^[a-z0-9](?:[a-z0-9-]{0,37}[a-z0-9])?$")
# A model as Claude Code's --model takes it: an alias or a full name, with
# no space a workflow line could split on.
MODEL = re.compile(r"^[A-Za-z0-9._:\[\]-]{0,80}$")
TURNS = range(1, 501)
ROUNDS = range(1, 51)
# A whole number of a worker's settings, as the form's answer or a file's
# JSON gives it, read only in ASCII digits, at most nine of them, far more
# than any setting takes: int() refuses a run of more than 4,300 digits
# with a ValueError, where every other bad value is named in one line, and
# it reads another script's digits as a number.
WHOLE_NUMBER = re.compile(r"[0-9]{1,9}")
ABOUT_LENGTH = 200
ACCOUNT_LENGTH = 120
PATTERN_LENGTH = 200


class WorkerError(Exception):
    """A worker's file, or the registry, that does not hold together; the
    message names every problem, one per line."""


@dataclass(frozen=True)
class Worker:
    """One worker as its file describes it, with the defaults filled in."""

    name: str
    runs: str
    about: str = ""
    area: tuple = field(default_factory=tuple)
    lands: str = "owner"
    model: str = ""
    max_turns: int = 80
    rounds: int = 6
    account: str = ""

    @property
    def label(self):
        return labels.worker_label(self.name)

    @property
    def file_name(self):
        return f"{self.name}.json"

    def to_dict(self):
        return {"name": self.name, "runs": self.runs, "about": self.about, "area": list(self.area),
                "lands": self.lands, "model": self.model, "max_turns": self.max_turns, "rounds": self.rounds,
                "account": self.account}

    def to_json(self):
        """The worker's file, every field written so a reader sees every
        setting there is, in the order to_dict gives."""
        return json.dumps(self.to_dict(), indent=2) + "\n"


KEYS = tuple(Worker("x", "session").to_dict())


# --- reading and checking -----------------------------------------------------

def is_plain_text(value, longest):
    return isinstance(value, str) and len(value) <= longest and value.isascii() and value.isprintable()


def pattern_problem(pattern):
    """What is wrong with one area pattern, or "": a path or glob relative
    to the repository's root, in ASCII, with no way out of it."""
    if not is_plain_text(pattern, PATTERN_LENGTH) or not pattern.strip():
        return f"an area entry must be a path of up to {PATTERN_LENGTH} printable ASCII characters"
    if pattern.startswith("/") or "\\" in pattern or ".." in pattern.split("/"):
        return f"the area entry {pattern!r} must be relative to the repository's root, with no .. in it"
    return ""


def field_problems(data):
    """What is wrong with each field of a worker's settings, one line each."""
    problems = []
    if not isinstance(data.get("name"), str) or not NAME.match(data["name"]):
        problems.append("name must be 1 to 39 lower case letters, digits or inner dashes, as nightly or alice-api")
    if data.get("runs") not in RUNS:
        problems.append("runs must be one of " + ", ".join(RUNS))
    if not is_plain_text(data.get("about", ""), ABOUT_LENGTH):
        problems.append(f"about must be one line of up to {ABOUT_LENGTH} printable ASCII characters")
    area = data.get("area", [])
    if not isinstance(area, list):
        problems.append("area must be a list of paths or globs")
    else:
        problems += [found for found in (pattern_problem(p) for p in area) if found]
    if data.get("lands", "owner") not in LANDS:
        problems.append("lands must be one of " + ", ".join(LANDS))
    if not isinstance(data.get("model", ""), str) or not MODEL.match(data.get("model", "")):
        problems.append("model must be a model's alias or name with no spaces, or empty for the default")
    for key, allowed in (("max_turns", TURNS), ("rounds", ROUNDS)):
        value = data.get(key, allowed.start)
        if not isinstance(value, int) or isinstance(value, bool) or value not in allowed:
            problems.append(f"{key} must be a whole number from {allowed.start} to {allowed.stop - 1}")
    if not is_plain_text(data.get("account", ""), ACCOUNT_LENGTH):
        problems.append(f"account must be a note of up to {ACCOUNT_LENGTH} printable ASCII characters, "
                        f"never a credential")
    return problems


def problems_of(data):
    """Every problem with a worker's settings as a dict, unknown keys first,
    since a misspelt key is the likeliest mistake and its value is lost."""
    if not isinstance(data, dict):
        return ["a worker's file must hold one JSON object"]
    unknown = sorted(set(data) - set(KEYS))
    found = [f"{key} is not a setting a worker has; the settings are " + ", ".join(KEYS) for key in unknown]
    return found + field_problems(data)


def from_dict(data):
    """The Worker the settings describe, or WorkerError naming every problem."""
    problems = problems_of(data)
    if problems:
        raise WorkerError("\n".join(problems))
    given = dict(data)
    given["area"] = tuple(given.get("area", ()))
    return Worker(**given)


def parse(text):
    """The Worker a file's text describes, or WorkerError."""
    try:
        data = json.loads(text, parse_int=json_integer)
    except ValueError as problem:
        raise WorkerError(f"not JSON: {problem}") from None
    return from_dict(data)


def json_integer(literal):
    """An integer of a worker's file as parse reads it: an int when
    WHOLE_NUMBER takes it, else the float it stands for, which the check
    names as no whole number as it names 1.5, never a run of digits
    int() refuses, which made the file read as no JSON."""
    return int(literal) if WHOLE_NUMBER.fullmatch(literal) else float(literal)


def load(directory=REGISTRY):
    """Every worker the folder holds, sorted by name: an empty list for a
    folder that is not there, and WorkerError naming every problem in every
    file, so the check says all of them at once."""
    folder = Path(directory)
    if not folder.is_dir():
        return []
    workers, problems = [], []
    for path in sorted(folder.glob("*.json")):
        try:
            worker = parse(path.read_text(encoding="utf-8"))
        except (OSError, UnicodeDecodeError, WorkerError) as problem:
            problems += [f"{path.name}: {line}" for line in str(problem).splitlines()]
            continue
        if path.name != worker.file_name:
            problems.append(f"{path.name}: the file must be named {worker.file_name}, after the worker")
            continue
        workers.append(worker)
    problems += [f"{path.name}: a worker's file ends in .json" for path in sorted(folder.iterdir())
                 if path.is_file() and path.suffix != ".json" and path.name != "README.md"]
    if problems:
        raise WorkerError("\n".join(problems))
    return workers


def describe(worker):
    """One line on a worker, for the coordinator planning work for it: its
    name, where it runs, its area and what it is for."""
    area = ", ".join(worker.area) or "anywhere"
    about = f"; {worker.about}" if worker.about else ""
    return f"{worker.name} ({worker.runs}): area {area}{about}"


def in_actions(workers):
    """The workers that run in GitHub Actions, in name order."""
    return [worker for worker in workers if worker.runs == "actions"]


def named(workers, name):
    """The worker with this name, or None."""
    return next((worker for worker in workers if worker.name == name), None)


# --- areas --------------------------------------------------------------------

def pattern_regex(pattern):
    """A compiled match for one area pattern: a directory ending in / takes
    everything under it, ** any run of directories, * and ? stay inside one
    directory, and a plain path names a file or the directory of that name."""
    if not any(mark in pattern for mark in "*?") and not pattern.endswith("/"):
        return re.compile(re.escape(pattern) + r"(?:/.*)?\Z")
    if pattern.endswith("/"):
        pattern += "**"
    out, i = [], 0
    while i < len(pattern):
        if pattern.startswith("**/", i):
            out.append("(?:.*/)?")
            i += 3
        elif pattern.startswith("**", i):
            out.append(".*")
            i += 2
        else:
            out.append({"*": "[^/]*", "?": "[^/]"}.get(pattern[i], re.escape(pattern[i])))
            i += 1
    return re.compile("".join(out) + r"\Z")


def in_area(area, path):
    """Whether a path falls in an area; every path falls in an empty one."""
    return not area or any(pattern_regex(pattern).match(path) for pattern in area)


def fits(area, files):
    """Whether work that touches these files fits a worker with this area:
    always for a worker with no area or work that names no files, and
    otherwise only when every file falls in the area, so a worker never
    takes work that reaches into another's."""
    return not area or not files or all(in_area(area, path) for path in files)


# --- the "Add a worker" form ----------------------------------------------------

# The form's questions, as GitHub writes each answer under a heading of the
# question's label, and the setting each answer fills. The repository's
# issue form (.github/ISSUE_TEMPLATE/worker.yml in the template) uses
# exactly these labels.
FORM = (
    ("Name", "name"),
    ("Where it runs", "runs"),
    ("What it is for", "about"),
    ("Area", "area"),
    ("How its pull requests land", "lands"),
    ("Model", "model"),
    ("Turns per run", "max_turns"),
    ("Rounds in a row", "rounds"),
    ("Claude account", "account"),
)
# What GitHub writes for a question left empty.
NO_RESPONSE = "_No response_"
HEADING = re.compile(r"^###\s+(.+?)\s*$")


def form_answers(body):
    """The form's answers by question label, each as the text under its
    heading, with an empty answer as ""."""
    answers, current = {}, None
    for line in (body or "").replace("\r\n", "\n").splitlines():
        heading = HEADING.match(line)
        if heading:
            current = heading.group(1)
            answers[current] = []
        elif current is not None:
            answers[current].append(line)
    return {label: "" if "\n".join(lines).strip() == NO_RESPONSE else "\n".join(lines).strip()
            for label, lines in answers.items()}


def first_word(answer):
    """A dropdown's value: the word before the option's colon, as in
    "routine: a routine on someone's Claude account"."""
    return answer.split(":", 1)[0].strip().lower()


def as_number(answer):
    """A number answer as an int when WHOLE_NUMBER takes it, or the text
    itself, which the check then names as no number."""
    text = answer.strip()
    return int(text) if WHOLE_NUMBER.fullmatch(text) else answer


def from_form(body):
    """The settings the form's answers give, as the dict a worker's file
    holds; only the answers given, so the defaults fill the rest. The
    check (from_dict) says what is wrong with them."""
    answers = form_answers(body)
    data = {}
    for label, key in FORM:
        answer = answers.get(label, "")
        if not answer:
            continue
        if key in ("runs", "lands"):
            data[key] = first_word(answer)
        elif key == "area":
            data[key] = [line.strip().lstrip("-* ").strip("`") for line in answer.splitlines() if line.strip()]
        elif key in ("max_turns", "rounds"):
            data[key] = as_number(answer)
        elif key == "name":
            data[key] = answer.strip().strip("`")
        else:
            data[key] = " ".join(answer.split())
    return data
