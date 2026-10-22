"""What Claude spends on this repository's work, recorded beside the work,
for the metrics workflow to add up (docs/setup/metrics.md, Economics).

Two kinds of run spend tokens here, and each leaves its usage where the
metrics workflow can read it with GITHUB_TOKEN alone:

  A session, a person's or a routine's, local or in the cloud. Claude
  Code runs `python3 tools/usage.py hook` as a PreToolUse hook on Bash
  (.claude/settings.json). When the call is a git commit, it reads the
  session's transcript, and its helpers' (the subagents' transcripts
  beside it), and leaves the session's usage so far in the git directory;
  the commit-msg hook, which the session-start hook installs where no
  other commit-msg hook is, then runs `python3 tools/usage.py trailer
  FILE`, which adds the usage since the session's last commit to the
  commit as one trailer line:

      Symphony-Usage: v1 session=6e4cb20b seconds=842 claude-opus-5-5=120,3400,0,52000,9100

  per model, the uncached input tokens, the cache writes kept five
  minutes and an hour, the cache reads and the output tokens; seconds is
  the time the session was at work between the two commits, gaps of more
  than five minutes left out. The line travels with the commit, so it
  needs no API, which a cloud session or a routine does not have, and the
  metrics workflow reads it from the pull request's commits. What a
  session spends after its last commit is not on any commit, and so not
  counted. The same two hooks add a second line, Symphony-Profile, the
  session's profile since its last commit by tools/profile.py, whose
  docstring gives it: counts of calls, reads, the context and seconds in
  tools, for the research programme in docs/research/beyond-cooperation/.
  It never stands in the usage line's way: the hook leaves the usage
  before it reads the profile, which takes longer, and a profile that
  cannot be read leaves that one commit without the line, the next
  commit's line counting from the last commit that carried one; a last
  commit whose profile is not known, made before the line existed, leaves
  the next commit without it. The counts written on a session's last
  commit are kept in the clone's common git directory, which its linked
  worktrees share, so a session that commits from two of them counts
  nothing twice.

  A run in GitHub Actions: a worker's, the coordinator's or @claude's.
  After claude-code-action the workflow runs `python3 tools/usage.py
  actions`, which reads the action's execution file, the Agent SDK's
  messages, and writes the run's record as JSON, which the workflow
  uploads as an artifact whose name starts with symphony-usage-; the
  metrics workflow lists and reads those artifacts. The workflow runs a
  copy of this file made before Claude starts, with --root naming the
  checkout, since the run may leave the checkout on a branch older than
  this file; profile.py, economics.py and the protocol library are copied
  beside it where the checkout has them, and the record then carries the
  run's profile too, the counts of a commit's Symphony-Profile line over
  the whole run. Inside Actions the hook and the trailer record nothing,
  so a run is never counted twice.

Nothing is recorded while the metrics switch is off (the word off in
.github/symphony-metrics, or SDLC_METRICS=off in the environment). No
word of the conversation is kept: token counts by model, seconds, turns,
the profile's counts, the first eight characters of the session's id,
and, for a run in Actions, which job ran and for which worker, issue and
pull request. A hook never stops the session or the commit: whatever it
cannot read, it leaves alone, and it always exits 0.
"""

import functools
import importlib.util
import json
import os
import re
import shlex
import stat
import subprocess
import sys
import time
from collections import deque
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import NamedTuple

TOOLS = Path(__file__).resolve().parent
ROOT = TOOLS.parent
sys.path.insert(0, str(TOOLS))

TRAILER_KEY = "Symphony-Usage"
TRAILER_VERSION = "v1"
# The trailer's line, its value read to the line's end and split on
# blanks: a lazy value followed by optional blanks before the end of the
# line would be tried at every blank of a long run of them, and anyone who
# can push writes a commit's message, which the metrics run reads.
TRAILER_LINE = re.compile(r"(?im)^[ \t]*symphony-usage:[ \t]*v1[ \t]+(.+)$")
# The five counts of a model's usage, in the order the trailer writes them.
FIELDS = ("input", "cache_write_5m", "cache_write_1h", "cache_read", "output")
# A model's name as a trailer or an Actions run's record may give it:
# letters, digits, dots, colons, dashes and a tag in brackets such as [1m],
# at most 81 characters, as the SDK names its models. A record anyone
# uploads may name anything, and reading a name for the price table
# (economics.normalized) takes time with the square of a run of brackets.
MODEL_KEY = re.compile(r"[A-Za-z0-9][A-Za-z0-9._:\[\]-]{0,80}")
# A pause longer than this between two lines of a transcript is time
# nobody was at work, and is left out of a session's seconds.
IDLE_GAP = 300
# A usage left by the hook belongs to the commit that follows it within
# this long; an older one is a commit that never happened.
PENDING_FRESH = 600
# How far past a << its delimiter is looked for: a delimiter is one word,
# and reading the rest of a long line for each << on it took time with the
# square of the line.
DELIMITER_REACH = 1024
# git's scissors line, below which a commit made with -v or commit.verbose
# carries its diff in the editor's template and where git cuts the message:
# the comment string, then 24 dashes, >8 and 24 dashes, and LF alone.
SCISSORS = re.compile(r"(?m)^(\S{1,16}) -{24} >8 -{24}\n")
# The comment characters git commit picks from for core.commentChar auto, in
# its order: the first that begins no line of the message.
AUTO_COMMENTS = "#;@!$%^&|:"
# The start of a line git counts as nothing in a message, as it counts blanks.
SIGNED_OFF = "Signed-off-by: "
# git commit's options as its option parser reads them, for the hook's
# reading of a command's words (commit_options): the short ones that take a
# value, the rest of their word or else the next word; the short flags; and
# the long ones, which may be shortened to any prefix that names one alone,
# with those of them that take a value.
SHORT_VALUED, SHORT_FLAGS = "mFcCt", "qvsaeinopz"
LONG_VALUED = ("file", "author", "date", "message", "reedit-message", "reuse-message", "fixup", "squash", "trailer",
               "template", "cleanup", "pathspec-from-file")
LONG_OPTIONS = LONG_VALUED + ("quiet", "verbose", "reset-author", "signoff", "edit", "status", "gpg-sign", "all",
                              "include", "interactive", "patch", "only", "no-verify", "dry-run", "short", "branch",
                              "ahead-behind", "porcelain", "long", "null", "amend", "no-post-rewrite",
                              "untracked-files", "pathspec-file-nul", "allow-empty", "allow-empty-message")
SWITCH_FILES = (".github/symphony-metrics", ".github/sdlc-metrics")
SWITCH_FILE = SWITCH_FILES[0]
SWITCH_VARIABLE = "SDLC_METRICS"
HOOK_MARKER = "symphony-usage"
SYNTHETIC = "<synthetic>"
COMMIT_MSG_HOOK = """#!/bin/sh
# symphony-usage: adds the Claude Code session's usage and profile since
# its last commit to this commit as the Symphony-Usage and Symphony-Profile
# trailers (tools/usage.py). It never stops a commit.
root=$(git rev-parse --show-toplevel 2>/dev/null) || exit 0
[ -f "$root/tools/usage.py" ] || exit 0
python3 "$root/tools/usage.py" trailer "$1" >/dev/null 2>&1 || true
exit 0
"""


# --- the switch ---------------------------------------------------------------

def first_word(text):
    words = (text or "").split()
    return words[0].lower() if words else ""


def sharing_off(root=ROOT, environ=None):
    """Whether the metrics switch is off, by the rule tools/metrics.py
    reads it with (its sharing): the variable SDLC_METRICS, in the
    environment, is off with any word but on; the file, the new name or
    the old, is on only with the word on first, says nothing when missing,
    and is off when it says anything else or cannot be read, an off in
    either file winning over an on in the other. A person who wants no
    usage of theirs on their commits sets SDLC_METRICS=off in their own
    environment. Not imported from metrics, since this runs before every
    git commit and must stay quick."""
    environ = os.environ if environ is None else environ
    said = first_word(environ.get(SWITCH_VARIABLE))
    if said and said != "on":
        return True
    for name in SWITCH_FILES:
        try:
            text = (Path(root) / name).read_text(encoding="utf-8")
        except FileNotFoundError:
            continue
        except (OSError, UnicodeDecodeError):
            return True
        if first_word(text) != "on":
            return True
    return False


def in_actions(environ=None):
    environ = os.environ if environ is None else environ
    return environ.get("GITHUB_ACTIONS") == "true"


# --- reading usage ------------------------------------------------------------

def zero():
    return [0] * len(FIELDS)


def counts_of(usage):
    """The five counts of one message's usage, as the API reports it: the
    cache writes split by how long they are kept where the API says so,
    and all at five minutes where it does not."""
    usage = usage or {}
    written = usage.get("cache_creation_input_tokens") or 0
    split = usage.get("cache_creation") or {}
    hour = split.get("ephemeral_1h_input_tokens") or 0
    five = split.get("ephemeral_5m_input_tokens")
    five = written - hour if five is None else five
    return [int(usage.get("input_tokens") or 0), int(five or 0), int(hour or 0),
            int(usage.get("cache_read_input_tokens") or 0), int(usage.get("output_tokens") or 0)]


def add(totals, model, counts):
    row = totals.setdefault(model, zero())
    for index, value in enumerate(counts):
        row[index] += value


def transcript_lines(path):
    """Every JSON object of a transcript, one a line; a line that is not
    one, a transcript still being written, is passed over, and so is a
    line with a byte that is not UTF-8, as where a read ends inside a
    character Claude Code is still writing, which would otherwise stop the
    whole read."""
    try:
        with open(path, "rb") as stream:
            for line in stream:
                try:
                    entry = json.loads(line)
                except ValueError:
                    continue
                if isinstance(entry, dict):
                    yield entry
    except OSError:
        return


def message_usage(entries):
    """Each assistant message's model and counts, by the message's id: a
    message written as several lines, one per block, carries the same usage
    on each, so it counts once."""
    found = {}
    for entry in entries:
        message = entry.get("message") if entry.get("type") == "assistant" else None
        if not isinstance(message, dict) or not message.get("usage"):
            continue
        model = message.get("model") or "unknown"
        if model == SYNTHETIC:
            continue
        key = message.get("id") or entry.get("requestId") or entry.get("uuid")
        found[key] = (model, counts_of(message.get("usage")))
    return found


def helper_transcripts(transcript):
    """The subagents' transcripts of a session, kept in a folder named for
    the session beside its own transcript."""
    folder = Path(transcript).with_suffix("") / "subagents"
    return sorted(folder.glob("*.jsonl")) if folder.is_dir() else []


def active_seconds(entries):
    """The seconds a session was at work: the time between consecutive
    lines, pauses longer than IDLE_GAP left out."""
    times = sorted(t for t in (moment(e.get("timestamp")) for e in entries) if t is not None)
    return int(sum(gap for gap in (b - a for a, b in zip(times, times[1:])) if gap <= IDLE_GAP))


def moment(text):
    try:
        return datetime.fromisoformat(str(text).replace("Z", "+00:00")).timestamp()
    except ValueError:
        return None


def session_usage(transcript):
    """A session's usage so far: by model, the five counts summed over its
    messages and its helpers', and its seconds at work."""
    main = list(transcript_lines(transcript))
    messages = message_usage(main)
    for helper in helper_transcripts(transcript):
        messages.update(message_usage(transcript_lines(helper)))
    totals = {}
    for model, counts in messages.values():
        add(totals, model, counts)
    return {"models": totals, "seconds": active_seconds(main)}


def since(now, before):
    """What a session spent since the usage last written on a commit."""
    models = {}
    for model, counts in now["models"].items():
        earlier = (before.get("models") or {}).get(model) or zero()
        delta = [max(0, a - b) for a, b in zip(counts, earlier)]
        if any(delta):
            models[model] = delta
    return {"models": models, "seconds": max(0, now["seconds"] - int(before.get("seconds") or 0))}


# --- the trailer ----------------------------------------------------------------

def session_token(session):
    """The eight characters of a session's id the trailer carries, enough to
    tell a repository's sessions apart and no more: the id without a
    leading session_, as a claude.ai one has, and only its letters, digits
    and dashes."""
    return re.sub(r"[^A-Za-z0-9-]", "", re.sub(r"^session_", "", str(session or "")))[:8]


def trailer(session, spent):
    """The trailer's value for what a session spent, models in name order."""
    parts = [TRAILER_VERSION, f"session={session_token(session)}", f"seconds={spent['seconds']}"]
    parts += [f"{model}={','.join(str(n) for n in counts)}" for model, counts in sorted(spent["models"].items())]
    return " ".join(parts)


def count_in(text):
    """A trailer's count as a whole number, or None when it spells none: ASCII
    digits below 2 ** 53, the bound the Actions records are read with
    (economics.actions_record), since anyone who can push writes a commit's
    message and int() refuses more than 4,300 digits with an error."""
    return int(text) if re.fullmatch(r"[0-9]{1,16}", text) and int(text) < 2 ** 53 else None


def parse_trailer(text):
    """Every Symphony-Usage trailer in a commit message, as (session,
    seconds, {model: counts}); a part that cannot be read is passed over."""
    found = []
    for match in TRAILER_LINE.finditer(text or ""):
        session, seconds, models = "", 0, {}
        for part in match.group(1).split():
            key, _, value = part.partition("=")
            if key == "session":
                session = value
            elif key == "seconds" and count_in(value) is not None:
                seconds = count_in(value)
            elif MODEL_KEY.fullmatch(key):
                counts = [count_in(n) for n in value.split(",")]
                if len(counts) == len(FIELDS) and None not in counts:
                    models[key] = counts
        found.append((session, seconds, models))
    return found


def git(*args, cwd=None):
    done = subprocess.run(["git", *args], cwd=cwd, capture_output=True, text=True, timeout=30)
    if done.returncode != 0:
        raise OSError(done.stderr.strip() or f"git {args[0]} failed")
    return done.stdout.strip()


def state_dir(cwd=None):
    """Where the hook leaves what it read for the commit that follows,
    inside the git directory, so no clone ever commits it."""
    folder = Path(git("rev-parse", "--absolute-git-dir", cwd=cwd)) / "symphony" / "usage"
    folder.mkdir(parents=True, exist_ok=True)
    return folder


def kept_dir(cwd=None):
    """Where the counts written on a session's last commit are kept: the
    clone's common git directory, the same folder as state_dir's in the
    clone's own worktree, and shared by its linked worktrees, each of
    which has a git directory of its own, so that a session committing
    from two of them counts each call once."""
    try:
        common = git("rev-parse", "--git-common-dir", cwd=cwd)
    except OSError:
        return state_dir(cwd)
    folder = (Path(cwd or os.getcwd()) / common).resolve() / "symphony" / "usage"
    folder.mkdir(parents=True, exist_ok=True)
    return folder


def write_json(path, data):
    temporary = Path(str(path) + ".part")
    temporary.write_text(json.dumps(data), encoding="utf-8")
    temporary.replace(path)


def read_json(path):
    try:
        data = json.loads(Path(path).read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return None
    return data if isinstance(data, dict) else None


# --- the hook and the commit-msg hook -------------------------------------------

def without_heredocs(command):
    """The command with its here-document bodies left out, since they are
    data the command reads, not commands the shell runs: a commit message
    written in one, with an apostrophe in it, would read as a quote left
    open and hide the git commit that follows it. (A body whose delimiter
    is not quoted can still run a $( ) of its own, left out with it.) A
    command with no << has none and is not read, since reading a long one a
    character at a time is most of what a here-document costs
    tools/profile.py's reading of a session."""
    text = str(command)
    if "<<" not in text:
        return text
    kept, at = [], 0
    for start, end in heredoc_bodies(text):
        kept.append(text[at:start])
        at = end
    kept.append(text[at:])
    return "".join(kept)


@dataclass(slots=True)
class Frame:
    """What the reading of a command holds open at one depth of $( ), <( )
    or >( ): the quote, how deep in arithmetic it is, how many ( it has
    open, and the here-documents opened at this depth, each a delimiter and
    whether its body's leading tabs go (<<-), whose bodies wait for a
    newline."""
    quote: str | None = None
    arithmetic: int = 0
    parens: int = 0
    waiting: deque = field(default_factory=deque)


def heredoc_bodies(text):
    """Where the bodies of a command's here-documents lie, each span from
    the newline after the operators to the end of the last delimiter. The
    command is read whole, as bash 5.2 reads it: a << in quotes, in a <<<
    here-string, in arithmetic, $(( )), $[ ] or (( )), or in a comment opens
    nothing, on a line a string carries on from an earlier one too; each
    $(, <( and >( reads quotes and operators anew, $( inside double quotes
    too, where Claude Code writes its commit messages, and the bodies of
    the here-documents opened in it follow the first newline inside it;
    those of one that closes on the line that opens them follow the very
    next newline, wherever it falls, before any other's. A backquoted
    command is read by itself, its here-documents ending with it."""
    from sdlc_protocol import guard
    frames, handed, i, word_start = [Frame()], deque(), 0, True
    while i < len(text):
        frame, ch = frames[-1], text[i]
        if handed and (newline := newline_at(text, i, frame.quote)) is not None:
            end = past_bodies(text, newline + 1, handed, len(frames) > 1)
            # A backslash before the newline joins the next line to this
            # one, the bodies between; with nothing after them it would be
            # left joining nothing, which no reading of the rest takes.
            yield (i if end == len(text) else newline), end
            i, word_start = end, True
            continue
        if ch == "\\" and frame.quote != "'":
            # An escaped character is part of a word; a backslash-newline is
            # the line continuation the shell takes out, word and all.
            word_start = word_start and text.startswith("\n", i + 1)
            i += 2
            continue
        if ch == "\n" and frame.waiting and not frame.quote:
            end = past_bodies(text, i + 1, frame.waiting, len(frames) > 1)
            yield i, end
            i, word_start = end, True
            continue
        close = backquote_end(text, i) if ch == "`" and frame.quote != "'" else None
        if close is not None and handed and "\n" in text[i:close]:
            # The bodies handed on wait for the newline inside it, so it is
            # read as it comes, up to that newline.
            step = 1
        elif close is not None:
            inner = text[i + 1:close]
            if "<<" in inner:
                yield from ((i + 1 + start, i + 1 + end) for start, end in heredoc_bodies(inner))
            step = min(close + 1, len(text)) - i
        elif frame.quote:
            step = inside_quotes(text, i, frames)
        elif ch == "#" and word_start:
            step = line_end(text, i) - i
        else:
            step = outside_quotes(text, i, frames, handed)
        word_start = text[i + step - 1] in guard.WORD_BREAKS
        i += step


def newline_at(text, i, quote):
    """The newline the reading meets at i: text[i] itself, or the one a
    backslash outside single quotes takes out of the line; None where the
    reading meets none there."""
    if text[i] == "\n":
        return i
    if text[i] == "\\" and quote != "'" and text.startswith("\n", i + 1):
        return i + 1
    return None


def backquote_end(text, i):
    """Where the backquoted command that opens at i ends: at the first
    backquote after it that no backslash escapes, or at the text's end."""
    at = i + 1
    while at < len(text):
        if text[at] == "`":
            return at
        at += 2 if text[at] == "\\" else 1
    return len(text)


def inside_quotes(text, i, frames):
    """One step of the reading inside quotes, and its length: the quote's
    own character closes it, and inside double quotes $( opens a frame,
    while $(( is arithmetic."""
    frame = frames[-1]
    if frame.quote == '"' and text.startswith("$((", i):
        return 3
    if frame.quote == '"' and text.startswith("$(", i):
        frames.append(Frame())
        return 2
    if text[i] == frame.quote:
        frame.quote = None
    return 1


def outside_quotes(text, i, frames, handed):
    """One step of the reading outside quotes, and its length: a quote
    opens; $(, <( and >( open a frame, and the ) that closes one hands the
    here-documents it left waiting to be read at the next newline;
    arithmetic and ( nest, and the old $[ ] is passed over whole; and a <<
    that is neither a <<< nor in arithmetic opens a here-document, whose
    delimiter is looked for within reach on its line."""
    frame, ch = frames[-1], text[i]
    if ch in "'\"":
        frame.quote = ch
    elif text.startswith("$((", i):
        frame.arithmetic += 1
        return 3
    elif text.startswith("$(", i) or (text.startswith(("<(", ">("), i) and not frame.arithmetic):
        frames.append(Frame())
        return 2
    elif text.startswith("$[", i):
        return bracket_end(text, i) - i
    elif text.startswith("((", i):
        frame.arithmetic += 1
        return 2
    elif frame.arithmetic and text.startswith("))", i):
        frame.arithmetic -= 1
        return 2
    elif ch == "(":
        frame.parens += 1
    elif ch == ")" and frame.parens:
        frame.parens -= 1
    elif ch == ")" and len(frames) > 1:
        frames.pop()
        handed.extend(frame.waiting)
    elif text.startswith("<<<", i):
        return 3
    elif text.startswith("<<", i) and not frame.arithmetic:
        heredoc = delimiter_of(text[i + 2:i + 2 + DELIMITER_REACH].split("\n", 1)[0])
        if heredoc:
            frame.waiting.append(heredoc)
        return 2
    return 1


def bracket_end(text, i):
    """Where the old arithmetic $[ ] that opens at i ends: past the ] that
    closes it, on a later line too, or at the text's end, where none does."""
    depth, stop = 0, len(text)
    for at in range(i + 1, stop):
        if text[at] == "[":
            depth += 1
        elif text[at] == "]":
            depth -= 1
            if not depth:
                return at + 1
    return stop


def line_end(text, i):
    """Where the line holding text[i] ends: its newline, or the text's end."""
    end = text.find("\n", i)
    return len(text) if end == -1 else end


def past_bodies(text, at, waiting, in_substitution):
    """Where the command goes on after the bodies of the here-documents
    waiting at a newline, the first body starting at at, each taken from
    waiting as it is read: the bodies follow in order, each up to its
    delimiter's line, the last ending at its delimiter, and one the command
    ends inside runs to its end; inside $( ) a body may also end at its
    delimiter where a ) follows on the line, as bash 5.2 takes it, and the
    command goes on there, those not yet read left waiting."""
    end = at
    while waiting:
        word, dashed = waiting[0]
        while True:
            if at >= len(text):
                waiting.clear()
                return len(text)
            stop = line_end(text, at)
            closing = delimiter_length(text[at:stop], word, dashed, in_substitution)
            if closing is not None:
                break
            at = stop + 1
        waiting.popleft()
        end = at + closing
        if end < stop:
            return end
        at = stop + 1
    return end


def delimiter_length(line, word, dashed, in_substitution):
    """How much of a line closes a here-document's body, or None: the whole
    line where it is the delimiter alone, after its leading tabs for <<-,
    and inside $( ) the delimiter that begins the line where a ) follows it
    anywhere on the line, which bash 5.2 keeps for older scripts."""
    body = line.lstrip("\t") if dashed else line
    if body == word:
        return len(line)
    if in_substitution and body.startswith(word) and ")" in body[len(word):]:
        return len(line) - len(body) + len(word)
    return None


def delimiter_of(rest):
    """The word after a <<, as the shell reads it, and whether a - joined to
    the << (<<-) takes the body's leading tabs away; None where no word
    follows. Its quotes or a backslash only make the body literal, and a
    carriage return, as a command written with CRLF line ends has, is part
    of the word, as it is to bash."""
    dashed = rest.startswith("-")
    lexer = shlex.shlex(rest[1:] if dashed else rest, posix=True, punctuation_chars=True)
    lexer.whitespace_split = True
    lexer.whitespace = " \t"
    lexer.commenters = ""
    try:
        word = next(iter(lexer), "")
    except ValueError:
        return None
    if not word or all(ch in "();<>|&" for ch in word):
        return None
    return word, dashed


def is_commit(command):
    """Whether a shell command runs git commit, in any simple command of it,
    its here-documents' bodies left out as data."""
    return commit_arguments(command) is not None


def commit_arguments(command):
    """The words after commit in the first simple command of a shell
    command that runs git commit, or None where none does. git may stand
    after assignments, a reserved word of the shell's, as in if git commit;
    then, or { git commit; }, and a wrapper with its options, as in time -p
    or nice -n 5, which the guard reads past as it does for a push."""
    from sdlc_protocol import guard
    for argv in guard.simple_commands(without_heredocs(command or "")):
        start = guard.command_start(argv)
        if start is None or os.path.basename(argv[start]) != "git":
            continue
        rest = argv[start + 1:]
        while rest and rest[0].startswith("-"):
            rest = rest[2:] if rest[0] in ("-c", "-C") else rest[1:]
        if rest and rest[0] == "commit":
            return rest[1:]
    return None


class Words(NamedTuple):
    """What a git commit's words say of its message (commit_options): the
    --cleanup they give, None where they give none; whether git opens the
    editor; the template -t gives, "" for --no-template and None where they
    give none; and whether --fixup is given, in any of its forms."""

    cleanup: str | None
    editor: bool
    template: str | None
    fixup: bool


def commit_options(arguments):
    """What a git commit's words say of its message, as Words. git opens
    the editor unless -m, -F, -C or a plain --fixup gives the message or
    --no-edit says not to, and does for -e however the message is given.
    The words are read as git's option parser reads them: a long option by
    any prefix that names it alone, short options run together, a value
    the rest of its word or else the next word, and nothing after -- an
    option."""
    cleanup, given, edit, template, fixup = None, False, None, None, False
    words = iter(arguments)
    for word in words:
        if word == "--":
            break
        if word.startswith("--"):
            name, joined, value = word[2:].partition("=")
            option, negated = long_option(name) or (None, False)
            if option in LONG_VALUED and not joined and not negated:
                value = next(words, "")
            if option == "cleanup":
                cleanup = value
            elif option == "edit":
                edit = not negated
            elif option == "template":
                template = "" if negated else value
            elif option in ("message", "file", "reuse-message") and not negated:
                given = True
            elif option == "fixup" and not negated:
                fixup = True
                given = given or not re.match(r"(amend|reword):", value)
        elif word.startswith("-") and len(word) > 1:
            for index, flag in enumerate(word[1:], 2):
                if flag == "e":
                    edit = True
                elif flag in SHORT_VALUED:
                    value = word[index:] if index < len(word) else next(words, "")
                    given = given or flag in "mFC"
                    template = value if flag == "t" else template
                    break
                elif flag not in SHORT_FLAGS:
                    break
    return Words(cleanup, not given if edit is None else edit, template, fixup)


def long_option(name):
    """The long option of git commit that --name gives, and whether it is
    given as --no-name, as git's option parser reads the options the hook
    reads: the option of that name, else the one option the name begins;
    None where it names none or more than one."""
    if name in LONG_OPTIONS:
        return name, False
    found = [(option, False) for option in LONG_OPTIONS if option.startswith(name)]
    if name.startswith("no-") and len(name) > 3:
        found += [(option, True) for option in LONG_OPTIONS if option.startswith(name[3:])]
    return found[0] if len(found) == 1 else None


def run_hook(payload, environ=None):
    """The PreToolUse hook: before a git commit, leave the session's usage
    so far for the commit-msg hook, then its profile so far beside it.
    Silent, and never in the way: the usage is left first, so that a hook
    stopped while it reads the profile, at Claude Code's timeout for a
    hook, still leaves the usage for the commit."""
    if not isinstance(payload, dict) or in_actions(environ):
        return
    command = (payload.get("tool_input") or {}).get("command") if isinstance(payload.get("tool_input"), dict) else ""
    transcript, session = payload.get("transcript_path"), payload.get("session_id")
    arguments = commit_arguments(command) if transcript and session else None
    if arguments is None:
        return
    cwd = payload.get("cwd") or None
    root = git("rev-parse", "--show-toplevel", cwd=cwd)
    if sharing_off(root, environ):
        return
    pending = state_dir(cwd) / "pending.json"
    words = commit_options(arguments)
    left = {"session": str(session), "at": time.time(), "cleanup": words.cleanup, "editor": words.editor,
            "template": given_template(words.template, cwd), "fixup": words.fixup, **session_usage(transcript)}
    write_json(pending, left)
    profiled = profile_so_far(transcript)
    # Another session's hook in the same clone may have left its own usage
    # while this one read the profile; that one's commit comes next.
    again = read_json(pending) or {}
    if profiled is not None and (again.get("session"), again.get("at")) == (left["session"], left["at"]):
        write_json(pending, {**left, "profile": profiled})


def given_template(template, cwd=None):
    """The template a command's -t names, as a path from where the command
    runs, since git reads it from there and the commit-msg hook runs at the
    top of the tree; "" for --no-template, None where the words name none."""
    if not template:
        return template
    return str(Path(cwd or os.getcwd()) / os.path.expanduser(template))


def run_trailer(message_file, cwd=None, environ=None, now=None):
    """The commit-msg hook's step: add to the commit the usage and the
    profile the hook left since the session's last commit, and remember
    them as written. A commit with nothing spent since the last carries no
    line. It remembers nothing when its hook read nothing at all, a
    transcript not there for that one hook, so the next line counts from
    the last; and it remembers the counts read when they fell below those
    remembered, a transcript replaced under the session, so the lines after
    it count on from there rather than wait for the session to spend past
    the old counts."""
    if in_actions(environ):
        return
    environ = os.environ if environ is None else environ
    folder = state_dir(cwd)
    pending = read_json(folder / "pending.json")
    now = time.time() if now is None else now
    if not pending or now - float(pending.get("at") or 0) > PENDING_FRESH:
        return
    editor = editor_opened(pending, environ)
    template = pending.get("template")
    template = commit_template(cwd) if template is None else str(template)
    if says_nothing(message_file, cwd, editor, cleanup_mode(pending, cwd), template, pending.get("fixup") is True):
        return
    (folder / "pending.json").unlink(missing_ok=True)
    session = re.sub(r"[^A-Za-z0-9-]", "", str(pending.get("session")))
    if not session:
        return
    written = kept_dir(cwd) / f"{session}.json"
    before = read_json(written) or {}
    spent = since(pending, before)
    if not spent["models"] and not pending.get("models"):
        return
    profiled, profile_kept = profile_line(session, pending.get("profile"), before)
    state = {"models": pending.get("models") or {}, "seconds": pending.get("seconds") or 0}
    if profile_kept is not None:
        state["profile"] = profile_kept
    write_json(written, state)
    if not spent["models"]:
        return
    lines = [f"{TRAILER_KEY}: {trailer(session, spent)}"]
    if profiled:
        lines.append(profiled)
    add_trailers(message_file, lines, cwd, editor)


def editor_opened(pending, environ):
    """Whether git opened the editor for the commit whose commit-msg hook
    this is. git tells its hooks by GIT_EDITOR, which it sets to : where it
    opens none; where : is the person's own editor, which git leaves as it
    is, the words the hook read say instead (commit_options)."""
    return environ.get("GIT_EDITOR") != ":" or pending.get("editor") is True


def cleanup_mode(pending, cwd=None):
    """How git cleans up the commit's message: as the command's --cleanup
    says, else commit.cleanup in git's config, a git -c's included, which
    git hands its hooks; None for git's default."""
    if pending.get("cleanup"):
        return str(pending["cleanup"])
    try:
        return git("config", "--get", "commit.cleanup", cwd=cwd) or None
    except OSError:
        return None


def says_nothing(message_file, cwd=None, editor=True, cleanup=None, template="", fixup=False):
    """Whether git, once it has cleaned the commit's message up, refuses it
    for saying nothing, as it refuses a bare git commit's where the editor
    is true, as in Claude Code's shell, or would commit nothing but its own
    comments: the lines added to it would give it a message, the usage line
    its subject, so nothing is added and the usage waits for the commit
    that follows. git cuts the message at its scissors line, below which -v
    and commit.verbose put the diff, and cutting where git does not only
    holds the usage back; it takes the comments out as cleanup says
    (cleanup_mode), by default where it opened the editor (editor_opened),
    and counts Signed-off-by lines as nothing; and it refuses a message
    that is the commit template and nothing more (untouched), and an
    amend! commit of --fixup with nothing past its subject. The text goes
    to git as UTF-8 whatever the locale, and with its line ends as they
    are, since git cuts only at a scissors line that LF alone ends."""
    try:
        text = Path(message_file).read_bytes().decode("utf-8", errors="replace")
    except OSError:
        return False
    cut = SCISSORS.search(text)
    if cut:
        text = text[:cut.start()]
    if cleanup == "verbatim":
        return not text or (fixup and amend_says_nothing(text))
    strips = cleanup == "strip" or (cleanup in (None, "default") and editor)
    if strips or editor:
        alone = comments_alone(text, cut, cwd)
        if alone is None or alone:
            return bool(alone)
    elif nothing_but(text, (SIGNED_OFF,)):
        return True
    if not template and not fixup:
        return False
    cleaned = stripped(text, cut, cwd, comments=strips)
    return cleaned is not None and (untouched(cleaned, template, strips, cut, cwd)
                                    or (fixup and amend_says_nothing(cleaned)))


def comments_alone(text, cut, cwd=None):
    """Whether the text is nothing but blanks, comments and Signed-off-by
    lines, the comment string the scissors line's where there is one, as
    git picks it for commentChar auto, else the config's; None where git
    cannot say. For auto with no scissors line, any string git can pick
    that takes it all will do: wherever git's own comments are in the file,
    only git's pick takes them, and with none there the worst is a usage
    line held back a commit."""
    if not cut and (comment_setting(cwd) or "").lower() == "auto":
        return any(nothing_but(text, (mark, SIGNED_OFF)) for mark in AUTO_COMMENTS)
    left = stripped(text, cut, cwd)
    return None if left is None else nothing_but(left, (SIGNED_OFF,))


def untouched(cleaned, template, comments, cut, cwd=None):
    """Whether a cleaned-up message is the commit template and nothing
    more, as git's template_untouched reads it: the template, read from its
    file and cleaned up as the message was, begins the message, and the rest
    is blanks and Signed-off-by lines. A template git cannot read, or an
    empty one, is none."""
    try:
        raw = Path(template).read_bytes().decode("utf-8", errors="replace") if template else ""
    except OSError:
        return False
    model = stripped(raw, cut, cwd, comments) if raw else None
    return model is not None and cleaned.startswith(model) and nothing_but(cleaned[len(model):], (SIGNED_OFF,))


def amend_says_nothing(cleaned):
    """Whether an amend! commit's cleaned-up message has nothing past its
    subject but blanks and Signed-off-by lines, which git refuses for
    --fixup=amend: and reword: as an empty body. Kept verbatim, only no body
    at all is empty to git, but the lines would be the body all the same,
    which an autosquash makes the amended commit's message, so they wait."""
    if not cleaned.startswith("amend! "):
        return False
    return nothing_but(cleaned[subject_length(cleaned):], (SIGNED_OFF,))


def subject_length(text):
    """How long a message's subject is, as git's commit_subject_length reads
    it: its lines up to the first that is blank or only whitespace."""
    at = 0
    while at < len(text):
        end = text.find("\n", at)
        end = len(text) if end == -1 else end + 1
        if not text[at:end].strip():
            break
        at = end
    return at


def commit_template(cwd=None):
    """commit.template in git's config, a git -c's included, as a path from
    where the hook runs, ~ taken as git takes it; "" where none is set."""
    try:
        path = git("config", "--path", "--get", "commit.template", cwd=cwd)
    except OSError:
        return ""
    return str(Path(cwd or ".") / path) if path else ""


def nothing_but(text, starts):
    """Whether every line of the text is blank or begins with one of starts."""
    return all(not line.strip() or line.startswith(starts) for line in text.splitlines())


def comment_setting(cwd=None):
    """core.commentChar or core.commentString as git reads them, two names
    for one setting of which the last given counts, and auto in any case;
    None where neither is set or git cannot say."""
    try:
        listed = git("config", "--get-regexp", r"^core\.comment(char|string)$", cwd=cwd)
    except OSError:
        return None
    last = listed.splitlines()[-1] if listed else ""
    return last.partition(" ")[2] or None


def stripped(text, cut, cwd=None, comments=True):
    """The text cleaned up by git stripspace, its comment lines taken out
    where comments is true, the comment string the scissors line's where it
    was cut at one; None where git cannot."""
    config = ["-c", f"core.commentChar={cut[1]}"] if cut else []
    options = ["--strip-comments"] if comments else []
    try:
        done = subprocess.run(["git", *config, "stripspace", *options], input=text, cwd=cwd,
                              capture_output=True, encoding="utf-8", errors="replace", timeout=30)
    except (OSError, subprocess.SubprocessError):
        return None
    return done.stdout if done.returncode == 0 else None


def add_trailers(message_file, lines, cwd=None, editor=True):
    """Add the trailer lines to the commit's message by git's own rules for
    a trailer block, or at its end where git cannot. A message whose last
    line has no newline, as git merge leaves a merge's, is given one first:
    git interpret-trailers would join the lines to that line, with no blank
    line between, the usage line part of the merge's subject. Where git
    opened the editor and a scissors line cuts the message, git is told the
    comment string that line begins with, which git commit may have picked
    for commentChar auto and git interpret-trailers would not know: the
    lines would go below the scissors line, where git cuts them away."""
    text = end_last_line(message_file)
    cut = SCISSORS.search(text) if editor else None
    config = ["-c", f"core.commentChar={cut[1]}"] if cut else []
    try:
        git(*config, "interpret-trailers", "--in-place", *[arg for line in lines for arg in ("--trailer", line)],
            str(message_file), cwd=cwd)
    except OSError:
        with open(message_file, "a", encoding="utf-8") as stream:
            stream.write("\n" + "".join(f"{line}\n" for line in lines))


def end_last_line(path):
    """The text of the file at path, its last line ended with a newline
    where it had none; empty where the file cannot be read."""
    try:
        data = Path(path).read_bytes()
        if data and not data.endswith(b"\n"):
            with open(path, "ab") as stream:
                stream.write(b"\n")
            data += b"\n"
    except OSError:
        return ""
    return data.decode("utf-8", errors="replace")


# --- the profile beside the usage -------------------------------------------------

@functools.lru_cache(maxsize=None)
def profiler():
    """tools/profile.py, which reads a session's profile for its line,
    loaded once by its path beside this file, since the standard library
    has a module of the same name, which a test run under cProfile has
    already imported; loaded only when a profile is read, so that this
    file alone still records a run in Actions, from the copy each workflow
    makes of it."""
    spec = importlib.util.spec_from_file_location("symphony_profile", TOOLS / "profile.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def profile_so_far(transcript):
    """The session's profile counts so far, or None when they cannot be read."""
    try:
        return profiler().session_counts(transcript)
    except Exception:
        return None


def run_profile(execution_file):
    """A run's profile, the same counts as a commit's line over the whole
    execution file, or None when they cannot be read, as where the workflow
    copied this file without profile.py; the file has no clock, so its
    seconds in tools are nothing."""
    try:
        profile = profiler()
        return profile.since(profile.session_counts(execution_file), {})
    except Exception:
        return None


def profile_line(session, now, before):
    """The Symphony-Profile line for the profile since the session's last
    commit that carried one, and the counts to keep for the next commit.
    No profile read now gives no line and keeps the counts kept before, so
    the next line counts this commit's calls too. A last commit recorded
    without a profile, by a hook older than the line, gives no line, since
    the counts since it are not known and this line would claim them all,
    and starts the count from now; so do counts kept before that cannot be
    read."""
    earlier = before.get("profile") if before else {}
    earlier = earlier if isinstance(earlier, dict) else None
    if not isinstance(now, dict):
        return None, earlier
    try:
        profile = profiler()
        kept = profile.kept(now)
        if earlier is None or not profile.known(earlier):
            return None, kept
        return f"{profile.TRAILER_KEY}: {profile.trailer(session, profile.since(now, earlier))}", kept
    except Exception:
        return None, earlier


def hooks_dir(cwd=None):
    """Where git looks for this clone's hooks: core.hooksPath when set,
    else the git directory's hooks."""
    try:
        configured = git("config", "core.hooksPath", cwd=cwd)
    except OSError:
        configured = ""
    if configured:
        top = git("rev-parse", "--show-toplevel", cwd=cwd)
        return (Path(top) / Path(configured).expanduser()).resolve()
    return Path(git("rev-parse", "--absolute-git-dir", cwd=cwd)) / "hooks"


def install(cwd=None, environ=None):
    """Put the commit-msg hook in place, unless the switch is off, the run is
    in GitHub Actions, whose runs record their usage as an artifact, or the
    clone has a commit-msg hook of its own, which is left as it is. What
    was done, in a line."""
    if in_actions(environ):
        return "usage: in GitHub Actions the run's record is uploaded instead; no commit-msg hook"
    root = git("rev-parse", "--show-toplevel", cwd=cwd)
    if sharing_off(root, environ):
        return "usage: the metrics switch is off; no commit-msg hook"
    path = hooks_dir(cwd) / "commit-msg"
    if path.exists() and HOOK_MARKER not in path.read_text(encoding="utf-8", errors="replace"):
        return f"usage: {path} is the clone's own commit-msg hook; left as it is, and no usage is recorded"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(COMMIT_MSG_HOOK, encoding="ascii")
    path.chmod(path.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)
    return f"usage: commit-msg hook in place at {path}"


# --- a run in Actions -------------------------------------------------------------

def result_of(messages):
    """The SDK's result message, the last of the run's messages, or {}."""
    for message in reversed(messages if isinstance(messages, list) else []):
        if isinstance(message, dict) and message.get("type") == "result":
            return message
    return {}


def run_models(messages, result):
    """A run's counts by model: the result's modelUsage where the SDK gives
    it, with its cache writes at five minutes, else the assistant messages'
    own usage, summed once per message."""
    by_model = result.get("modelUsage")
    if isinstance(by_model, dict) and by_model:
        totals = {}
        for model, used in by_model.items():
            if isinstance(used, dict):
                add(totals, model, [int(used.get("inputTokens") or 0), int(used.get("cacheCreationInputTokens") or 0), 0,
                                    int(used.get("cacheReadInputTokens") or 0), int(used.get("outputTokens") or 0)])
        return totals
    totals = {}
    for model, counts in message_usage(m for m in messages if isinstance(m, dict)).values():
        add(totals, model, counts)
    return totals


def actions_record(messages, job, workflow, worker="", issue=None, pull=None, at=None):
    """The record of one run in Actions, from the execution file's messages."""
    result = result_of(messages)
    return {
        "v": 1, "source": "actions", "workflow": workflow, "job": job, "worker": worker or "",
        "issue": issue, "pull": pull,
        "at": at or datetime.now(timezone.utc).replace(microsecond=0).isoformat().replace("+00:00", "Z"),
        "seconds": round((result.get("duration_ms") or 0) / 1000, 1),
        "api_seconds": round((result.get("duration_api_ms") or 0) / 1000, 1),
        "turns": int(result.get("num_turns") or 0),
        "reported_usd": float(result.get("total_cost_usd") or 0),
        "error": bool(result.get("is_error")) or not result,
        "models": run_models(messages if isinstance(messages, list) else [], result),
    }


def number_or_none(text):
    return int(text) if str(text or "").isdecimal() else None


# --- the command line -------------------------------------------------------------

def main(argv=None, stdin=None, environ=None):
    argv = sys.argv[1:] if argv is None else list(argv)
    command = argv[0] if argv else ""
    try:
        if command == "hook":
            try:
                payload = json.load(stdin or sys.stdin)
            except ValueError:
                return 0
            run_hook(payload, environ)
        elif command == "trailer" and len(argv) == 2:
            run_trailer(argv[1], environ=environ)
        elif command == "install":
            print(install(environ=environ))
        elif command == "actions":
            return run_actions(argv[1:], environ)
        else:
            print("usage: python3 tools/usage.py hook|trailer FILE|install|actions --execution-file F --out F "
                  "--workflow W --job J [--worker W] [--issue N] [--pull N] [--root DIR]", file=sys.stderr)
            return 2
    # A hook never stands in a session's way, whatever went wrong.
    except Exception as problem:
        if command == "install":
            print(f"usage: {problem}")
    return 0


def run_actions(argv, environ=None):
    """Write the record of the run whose execution file is named, unless
    the switch is off; never a failure of the workflow."""
    import argparse
    parser = argparse.ArgumentParser(prog="python3 tools/usage.py actions")
    parser.add_argument("--execution-file", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--workflow", required=True)
    parser.add_argument("--job", required=True)
    parser.add_argument("--worker", default="")
    parser.add_argument("--issue", default="")
    parser.add_argument("--pull", default="")
    parser.add_argument("--root", default=str(ROOT),
                        help="the checkout whose metrics switch counts, when this file runs from a copy")
    args = parser.parse_args(argv)
    if sharing_off(args.root, environ):
        print("usage: the metrics switch is off; no record")
        return 0
    try:
        messages = json.loads(Path(args.execution_file).read_text(encoding="utf-8", errors="replace"))
    except (OSError, ValueError):
        print("usage: no execution file to read; no record")
        return 0
    record = actions_record(messages, args.job, args.workflow, args.worker,
                            number_or_none(args.issue), number_or_none(args.pull))
    profiled = run_profile(args.execution_file)
    if profiled is not None:
        record["profile"] = profiled
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(record, indent=2) + "\n", encoding="ascii")
    tokens = sum(sum(counts) for counts in record["models"].values())
    print(f"usage: {args.workflow}/{args.job}: {tokens} tokens in {record['seconds']} s, recorded")
    return 0


if __name__ == "__main__":
    sys.exit(main())
