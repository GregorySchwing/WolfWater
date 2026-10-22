"""What a Claude Code session spent, read from the transcript it leaves:
its context, call by call; its tool calls, by tool; what it read, and
read again; and where its time went. The first instrument of the research
programme in docs/research/beyond-cooperation/README.md, which measures a
way of working before changing it.

  python3 tools/profile.py TRANSCRIPT [--json FILE] [--paths N] [--tools N] [--no-helpers]

prints a report and, with --json, writes the record docs/research/
beyond-cooperation/README.md describes, aggregates only: counts, bytes,
seconds and tokens, with model names and the names of the harness's own
tools, every MCP tool folded into one row and any other tool into another,
and never a path, a command, a prompt or any text of the conversation, so
a record can be committed beside the question it answers. --paths N
prints the N files read most often and --tools N the N tools called most
often by their own names, which stay on the screen.

The transcript is the one Claude Code writes under ~/.claude/projects/,
one JSON object a line, the same file tools/usage.py reads for the usage
trailer, and its helpers' transcripts (the subagents') are read beside it
as usage.py reads them: their calls join the per-tool rows and their
tokens the totals, with their part given in the record's helpers block,
while the context, the time and the prompts are the session's own. The
Agent SDK's execution file of a run in Actions, one JSON array of the
same messages, is read too, though its messages carry no timestamps, so
its time measures are empty and it has no session token; its helpers'
messages are in the same array, each marked with the tool call that
started its helper (parent_tool_use_id), and are read as the helpers'
transcripts are. Each API call is
counted once by its message id: a message written as several lines, one
per block, carries its usage on each, the input and cache counts the same
on every line and the output count as far as the message had got when
the line was written, so the largest is taken; a third-party counter that
sums every line over-counts.

The time measures rest on the transcript's own clock, line by line: the
active seconds are the gaps between consecutive lines with pauses of more
than five minutes left out, as usage.py counts them; the seconds in tools
are the union of each tool call's interval, from the line that asked to
the line that answered, so two tools running at once count once; the
seconds outside tools are the active seconds not covered by that union,
which is the model's time, its thinking and its answer, and the harness's
own time between lines, which the transcript cannot tell apart; the
harness's own count of the thinking is given beside them. Per tool, the
seconds are summed over its calls, so parallel calls each count. The
tokens each tool brought into the context are an estimate, and the key
says so: the context's growth from one call to the next, less the call's
own output, shared among the call's answered tools by the size of their
answers; a reminder the harness added with the answers is counted with
them. A repeated call is one the same as an earlier call in the same
transcript, the same tool with the same input: MAST's step repetition
(arXiv:2503.13657) read by code, and an upper bound on waste, since a
git status after a change is a repeat and no waste.

The same counts ride on a session's commits, as a second line beside the
Symphony-Usage trailer, written by the same hooks of tools/usage.py: the
profile since the session's last commit,

    Symphony-Profile: v1 session=6e4cb20b calls=42 failed=3 repeated=1 reads=12 rereads=4 read_bytes=58210 reread_bytes=9120 peak=181234 compactions=0 drops=0 tool_seconds=95

each count the rise of the record's own since that commit, the calls and
the reads with the helpers' in them, the compactions, the drops and the
seconds in tools the session's own, and peak the largest context of the
session's own calls made since; the session is the usage line's, so the
two lines pair up. It carries counts alone, no tool's name, since an MCP
server names its tools, so it stands on a public commit as the usage line
does, and the metrics switch that stops the one stops the other.
tools/metrics.py sums it per pull request (pull_profile), with the
profile the record of each Actions run carries, per work setting and in
the coordination block's two sides (group, coordination), into
PROFILE.json, which stays on the repository's own branch metrics.

Standard library only, like every tool here; the shell commands are split
by the vendored protocol library's reader, the one the push guard uses.
"""

import argparse
import functools
import json
import posixpath
import re
import sys
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

TOOLS = Path(__file__).resolve().parent
sys.path.insert(0, str(TOOLS))

import economics  # noqa: E402
import usage  # noqa: E402
from sdlc_protocol import guard  # noqa: E402

SCHEMA = "symphony-profile/1"
# Programs that print a file: a shell command of one of them naming a file
# is a read of it, as a call of the Read tool is; one naming none is a
# filter on another command's output.
READERS = ("cat", "head", "tail", "sed", "less", "more", "nl")
# Those that print part of a file, as a Read with an offset or a limit does.
PARTIAL_READERS = ("head", "tail", "sed")
SEARCHERS = ("grep", "rg", "ag", "ack", "egrep", "fgrep")
LISTERS = ("ls", "find", "tree", "wc", "stat", "du", "file")
GIT = ("git", "gh")
# Builtins that read and run nothing worth counting.
QUIET = ("cd", "echo", "printf", "true", "false", "exit", "export", "set", "test", "[", "pwd", "read", "shift",
         "local", "declare", "unset", "break", "continue", "return", "wait", ":")
# Words the shell reads before a command and that name none: its reserved
# words, the guard's and the ones that close a construct, and the wrappers
# that run the command after them as it is (timeout and nice take a value
# first and are not read); a loop or case header names no command at all.
KEYWORDS = guard.RESERVED | {"}", "done", "fi", "esac"}
WRAPPERS = ("time", "sudo", "nohup", "exec", "command", "builtin", "env", "xargs")
HEADERS = ("for", "select", "case")
SHELL_USES = ("read", "filter", "search", "list", "git", "run")
# The readers' options whose value is the next word, so that it is not
# taken for a file: head -n 20 FILE, sed -e SCRIPT FILE, nl -w 3 FILE; the
# short ones by their letter, the long ones by name (a value joined with =
# is in the same word).
SHORT_VALUED = {"head": "nc", "tail": "nc", "sed": "efl", "nl": "bdfhilnsvw"}
LONG_VALUED = {"head": ("--lines", "--bytes"), "tail": ("--lines", "--bytes"),
               "sed": ("--expression", "--file"),
               "nl": ("--body-numbering", "--footer-numbering", "--header-numbering", "--line-increment",
                      "--number-format", "--number-separator", "--number-width", "--section-delimiter",
                      "--starting-line-number")}
# The sed options that give the script, after which every plain word is a file.
SCRIPT_OPTIONS = ("--expression", "--file")
# A fall of the context between two consecutive calls by more than this
# share is a drop: a compaction, or the harness trimming what it sends.
DROP = 0.5
GLOB = "*?[$~`"
# A transcript is named for its session id; any other name is a run's file.
SESSION_STEM = re.compile(r"[0-9a-f]{8}(-[0-9a-f]+)*")
BY_TOOL_KEYS = ("calls", "failed", "repeated", "input_bytes", "result_bytes", "seconds_in_tool",
                "context_tokens_estimated")
# The record names the harness's own tools, these, and folds every other:
# each MCP tool into MCP_TOOLS, since an MCP server names its tools and a
# server's name can say what a person uses, and any other name, a tool an
# Agent SDK program defines for itself, into OTHER_TOOLS. --tools shows the
# names on the screen alone.
HARNESS_TOOLS = frozenset({
    "Agent", "Task", "Bash", "BashOutput", "KillBash", "KillShell", "Read", "Write", "Edit", "MultiEdit", "Glob",
    "Grep", "LS", "NotebookEdit", "NotebookRead", "WebFetch", "WebSearch", "TodoRead", "TodoWrite", "ExitPlanMode",
    "EnterPlanMode", "AskUserQuestion", "SlashCommand", "Skill", "ToolSearch", "TaskCreate", "TaskGet", "TaskList",
    "TaskUpdate", "TaskStop", "TaskOutput", "SendMessage", "ListAgents", "ReadNotifications", "SubagentHandback",
    "Monitor", "Workflow", "ScheduleWakeup", "CronCreate", "CronDelete", "CronList", "EnterWorktree", "ExitWorktree",
    "ListMcpResourcesTool", "ReadMcpResourceTool", "ReadMcpResourceDirTool", "PushNotification", "Artifact",
    "ArtifactComments", "ArtifactData", "ReportFindings", "ListSkills", "SearchSkills", "SuggestSkills",
    "ListPlugins", "SearchPlugins", "SuggestPluginInstall", "ListConnectors", "SuggestConnectors",
    "SearchMcpRegistry", "DesignSync",
})
MCP_TOOLS = "mcp"
OTHER_TOOLS = "other"
SECONDS_KEYS = ("span", "active", "in_tools", "outside_tools", "thinking")
HELPERS_KEYS = ("transcripts", "calls", "seconds_in_tool", "output_tokens")
READS_KEYS = ("calls", "files_named", "files", "rereads", "partial", "bytes", "reread_bytes")
# The line a commit carries, read to its end and split on blanks, as
# usage.py reads its own, since anyone who can push writes a commit's
# message: a lazy value before optional blanks would be tried at every
# blank of a long run of them.
TRAILER_KEY = "Symphony-Profile"
TRAILER_VERSION = "v1"
TRAILER_LINE = re.compile(r"(?im)^[ \t]*symphony-profile:[ \t]*v1[ \t]+(.+)$")
TRAILER_FIELDS = ("calls", "failed", "repeated", "reads", "rereads", "read_bytes", "reread_bytes", "peak",
                  "compactions", "drops", "tool_seconds")
SUMS = tuple(field for field in TRAILER_FIELDS if field != "peak")
# What PROFILE.json (tools/metrics.py) gives for each work setting and each
# side of the comparison, and the ratios, over 1 where working in parallel
# did better, as the coordination block's are (tools/economics.py).
GROUP_KEYS = ("pull_requests", "merged", "profiled", "profiled_merged", *SUMS, "peak_median")
RATIO_KEYS = ("calls_per_merge", "reads_per_merge", "rereads_per_merge", "peak")
# The session as the line writes it, at most usage.session_token's eight characters.
TRAILER_SESSION = re.compile(r"[A-Za-z0-9-]{1,8}")
# A blank line, which ends a paragraph of a commit's message.
PARAGRAPH_BREAK = re.compile(r"\n[ \t\r]*\n")


# --- reading the transcript -----------------------------------------------------

def entries_of(path):
    """A transcript's entries: one JSON array, as the Agent SDK writes a
    run's execution file in Actions, which claude-code-action writes over
    many lines, or one JSON object a line, as Claude Code writes a
    session's. The array is read whole first, since a line of it can be a
    whole object by itself, an empty one, and read as a line it would be
    taken for the file's one entry."""
    if opens_an_array(path):
        try:
            data = json.loads(Path(path).read_text(encoding="utf-8", errors="replace"))
        except (OSError, ValueError):
            data = None
        if isinstance(data, list):
            return [e for e in data if isinstance(e, dict)]
    return list(usage.transcript_lines(path))


def opens_an_array(path):
    try:
        with open(path, encoding="utf-8", errors="replace") as stream:
            return stream.read(4096).lstrip()[:1] == "["
    except OSError:
        return False


def threads_of(entries):
    """A transcript's entries parted by whose they are: the session's own,
    and each helper's, in the order each helper first speaks, where the
    entries name the tool call that started their helper, as an execution
    file's do (parent_tool_use_id); a session's transcript keeps its
    helpers in files of their own and names none."""
    own, helpers = [], {}
    for entry in entries:
        parent = entry.get("parent_tool_use_id")
        if parent:
            helpers.setdefault(str(parent), []).append(entry)
        else:
            own.append(entry)
    return own, list(helpers.values())


def blocks(message):
    """The content blocks of a message; a string content has none."""
    content = message.get("content") if isinstance(message, dict) else None
    return [b for b in content if isinstance(b, dict)] if isinstance(content, list) else []


def text_bytes(content):
    """The size of what a tool returned as the model saw it: a string, the
    text of a text block, and any other block by the size of its JSON (a
    tool reference, or an image by its data, which is not what an image
    costs in tokens)."""
    if isinstance(content, str):
        return len(content.encode("utf-8"))
    if not isinstance(content, list):
        return 0
    total = 0
    for block in content:
        if isinstance(block, dict) and block.get("type") == "text":
            total += len(str(block.get("text") or "").encode("utf-8"))
        else:
            total += len(json.dumps(block, sort_keys=True).encode("utf-8"))
    return total


def api_calls(entries):
    """Each API call of a transcript once, in order of its first line, with
    the tool calls it made: a message's lines share its id, and its
    tool_use blocks are told apart by theirs. Its usage is the largest its
    lines carry, since a helper's lines carry the output count as it grew."""
    calls, by_id = [], {}
    for entry in entries:
        message = entry.get("message") if entry.get("type") == "assistant" else None
        if not isinstance(message, dict) or message.get("model") == usage.SYNTHETIC:
            continue
        key = message.get("id") or entry.get("requestId") or entry.get("uuid")
        at = usage.moment(entry.get("timestamp"))
        call = by_id.get(key)
        counts = usage.counts_of(message.get("usage")) if message.get("usage") else None
        thinking = entry.get("thinkingDurationMs") or 0
        if call is None:
            call = {"model": message.get("model") or "unknown", "at": at, "until": at, "counts": counts,
                    "tools": [], "thinking_ms": thinking}
            by_id[key] = call
            calls.append(call)
        else:
            if at is not None and (call["until"] is None or at > call["until"]):
                call["until"] = at
            if counts:
                call["counts"] = [max(a, b) for a, b in zip(call["counts"], counts)] if call["counts"] else counts
            call["thinking_ms"] = max(call["thinking_ms"], thinking)
        seen = {t["id"] for t in call["tools"]}
        for block in blocks(message):
            if block.get("type") == "tool_use" and block.get("id") and block["id"] not in seen:
                call["tools"].append({"id": block["id"], "name": str(block.get("name") or "unknown"),
                                      "input": block.get("input") if isinstance(block.get("input"), dict) else {},
                                      "at": at, "cwd": entry.get("cwd") or "", "result": None})
    return calls


def tool_results(entries):
    """Each tool's answer by the id of the call it answers: when it came,
    how big it was and whether it was an error."""
    found = {}
    for entry in entries:
        if entry.get("type") != "user":
            continue
        for block in blocks(entry.get("message")):
            if block.get("type") == "tool_result" and block.get("tool_use_id"):
                found[block["tool_use_id"]] = {"at": usage.moment(entry.get("timestamp")),
                                               "bytes": text_bytes(block.get("content")),
                                               "error": bool(block.get("is_error"))}
    return found


def context_of(counts):
    """The tokens an API call sent: uncached input, cache writes and cache reads."""
    return counts[0] + counts[1] + counts[2] + counts[3]


def is_prompt(entry):
    """A person's own message: a user line that is no tool result, no
    attachment the harness added and no summary a compaction wrote."""
    if entry.get("type") != "user" or entry.get("isMeta") or entry.get("isCompactSummary"):
        return False
    message = entry.get("message")
    if isinstance(message, dict) and isinstance(message.get("content"), str):
        return True
    return bool(blocks(message)) and all(b.get("type") == "text" for b in blocks(message))


def compactions(entries):
    """How many times the context was compacted: the boundary lines the
    harness writes, or, where it writes none, the summaries it starts the
    continued conversation with."""
    boundaries = sum(1 for e in entries if e.get("type") == "system" and e.get("subtype") == "compact_boundary")
    summaries = sum(1 for e in entries if e.get("type") == "user" and e.get("isCompactSummary"))
    return max(boundaries, summaries)


# --- what was read --------------------------------------------------------------

def resolve(path, base):
    """A path as the transcript names it, made absolute against the line's
    working directory where it is a plain path; a pattern or a variable is
    kept as written, since the file it means cannot be told."""
    if not path or path == "-" or any(ch in path for ch in GLOB):
        return None
    if posixpath.isabs(path) or (len(path) > 1 and path[1] == ":"):
        return posixpath.normpath(path)
    return posixpath.normpath(posixpath.join(base or "", path)) if base else path


def program_of(argv):
    """The program a simple command runs, after any VAR=value words, the
    shell's keywords and the wrappers that run it as it is; a loop or case
    header runs none."""
    words = list(argv)
    while words and (guard.ASSIGNMENT.match(words[0]) or words[0] in KEYWORDS or words[0] in WRAPPERS):
        words.pop(0)
    if not words or words[0] in HEADERS:
        return "", []
    return guard.program(words[0]), words[1:]


# The here-document reader is usage.py's, whose hook leaves the bodies out
# too before it looks for a git commit.
without_heredocs = usage.without_heredocs


@functools.lru_cache(maxsize=256)
def commands_of(command):
    """The simple commands of a shell command, its here-documents left
    out: read once for the three readers of each call (what it is for,
    what it reads, and whether it reads part of a file), since splitting a
    long command is most of what profiling a session costs."""
    return tuple(tuple(argv) for argv in guard.simple_commands(without_heredocs(command)))


def is_write(program, arguments):
    """sed editing in place writes the file and reads it for nobody."""
    return program == "sed" and any(
        word == "--in-place" or word.startswith("--in-place=") or (word.startswith("-") and not word.startswith("--")
                                                                     and "i" in word)
        for word in arguments)


def files_of(program, arguments):
    """The files a reading command names: its words that are no option
    (less and more take +options too), no option's value and not the
    standard input (-), and for sed, whose script is its first plain word
    unless an option gave it, the words after the script."""
    files, skip, script_seen = [], False, program != "sed"
    for word in arguments:
        if skip:
            skip = False
        elif word == "-" or (word.startswith("+") and program in ("less", "more")):
            continue
        elif word.startswith("--"):
            name, joined, _ = word.partition("=")
            if name in LONG_VALUED.get(program, ()):
                skip = not joined
                script_seen = script_seen or name in SCRIPT_OPTIONS
        elif word.startswith("-") and len(word) > 1:
            if word[-1] in SHORT_VALUED.get(program, "") and not word[1:-1].isdigit():
                skip = True
                script_seen = script_seen or word[-1] in "ef"
        elif not script_seen:
            script_seen = True
        else:
            files.append(word)
    return files


def shell_reads(command, base):
    """The files a shell command prints, as resolve() names them, following
    a cd inside the command."""
    paths = []
    for argv in commands_of(command):
        program, arguments = program_of(argv)
        if program == "cd" and len(arguments) == 1:
            base = resolve(arguments[0], base) or base
        elif program in READERS and not is_write(program, arguments):
            paths.extend(p for p in (resolve(f, base) for f in files_of(program, arguments)) if p)
    return paths


def shell_uses(command):
    """What each simple command of a shell command is for."""
    uses = Counter()
    for argv in commands_of(command):
        program, arguments = program_of(argv)
        if not program or program in QUIET:
            continue
        if program in READERS and not is_write(program, arguments):
            uses["read" if files_of(program, arguments) else "filter"] += 1
        elif program in SEARCHERS:
            uses["search"] += 1
        elif program in LISTERS:
            uses["list"] += 1
        elif program in GIT:
            uses["git"] += 1
        else:
            uses["run"] += 1
    return uses


def reads_of(tool):
    """The files one tool call reads, and whether it reads part of them."""
    name, given = tool["name"], tool["input"]
    if name == "Read":
        path = resolve(str(given.get("file_path") or ""), tool["cwd"])
        return ([path] if path else []), bool(given.get("offset") or given.get("limit"))
    if name == "Bash":
        command = str(given.get("command") or "")
        partial = any(program_of(argv)[0] in PARTIAL_READERS for argv in commands_of(command))
        return shell_reads(command, tool["cwd"]), partial
    return [], False


# --- time -------------------------------------------------------------------------

def tool_seconds(tool):
    """From the line that asked to the line that answered; nothing when
    either is missing or the answer came first."""
    result = tool["result"]
    if result is None or result["at"] is None or tool["at"] is None or result["at"] < tool["at"]:
        return 0
    return result["at"] - tool["at"]


def union(intervals):
    """The intervals merged, so time two tools shared counts once."""
    merged = []
    for start, end in sorted(intervals):
        if merged and start <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], end))
        else:
            merged.append((start, end))
    return merged


def active_outside(times, intervals):
    """The active seconds, each gap between consecutive lines of at most
    IDLE_GAP, less the part of each a tool was running. The lines and the
    merged intervals are walked once, side by side, as both are in order:
    a long session has thousands of each, and the hook reads it before
    every commit."""
    times, merged = sorted(times), union(intervals)
    total, first = 0.0, 0
    for a, b in zip(times, times[1:]):
        if b - a > usage.IDLE_GAP:
            continue
        while first < len(merged) and merged[first][1] <= a:
            first += 1
        covered, later = 0.0, first
        while later < len(merged) and merged[later][0] < b:
            covered += max(0.0, min(b, merged[later][1]) - max(a, merged[later][0]))
            later += 1
        total += (b - a) - covered
    return total


# --- the profile -------------------------------------------------------------------

def empty_tool():
    return dict.fromkeys(BY_TOOL_KEYS, 0)


def tool_kind(name):
    """The name a tool's row has in the record: its own for the harness's
    tools, else one name for every MCP tool and one for every other."""
    if name in HARNESS_TOOLS:
        return name
    return MCP_TOOLS if name.startswith("mcp__") else OTHER_TOOLS


class Profile:
    """The record as it is built, one transcript at a time."""

    def __init__(self):
        self.models = {}
        self.by_tool = {}
        self.calls = Counter()
        self.shell = Counter()
        self.reads = Counter()
        self.files = set()
        self.per_path = Counter()
        self.context = {}
        self.sizes = []
        self.seconds = dict.fromkeys(SECONDS_KEYS, 0)
        self.prompts = 0
        self.attachments = 0
        self.helpers = dict.fromkeys(HELPERS_KEYS, 0)
        self.lines = 0
        self.started = None
        self.ended = None

    def add_transcript(self, entries, main):
        calls, results = api_calls(entries), tool_results(entries)
        for call in calls:
            if call["counts"]:
                usage.add(self.models, call["model"], call["counts"])
            for tool in call["tools"]:
                tool["result"] = results.get(tool["id"])
        self.add_tools(calls)
        if not main:
            self.add_helper(calls)
        if main:
            self.lines = len(entries)
            self.add_context(calls)
            self.context["compactions"] = compactions(entries)
            self.add_time(entries, calls)
            self.prompts = sum(1 for e in entries if is_prompt(e))
            self.attachments = sum(1 for e in entries if e.get("type") == "attachment")

    def add_tools(self, calls):
        seen, made = set(), set()
        for call in calls:
            for tool in call["tools"]:
                row = self.by_tool.setdefault(tool["name"], empty_tool())
                row["calls"] += 1
                given = json.dumps(tool["input"], sort_keys=True)
                row["input_bytes"] += len(given.encode("utf-8"))
                self.calls["total"] += 1
                if (tool["name"], given) in made:
                    row["repeated"] += 1
                    self.calls["repeated"] += 1
                made.add((tool["name"], given))
                result = tool["result"]
                if result is None:
                    self.calls["unanswered"] += 1
                else:
                    row["result_bytes"] += result["bytes"]
                    row["failed"] += int(result["error"])
                    self.calls["failed"] += int(result["error"])
                    row["seconds_in_tool"] += tool_seconds(tool)
                if tool["name"] == "Bash":
                    self.shell.update(shell_uses(str(tool["input"].get("command") or "")))
                self.add_reads(tool, seen)

    def add_helper(self, calls):
        """What a helper's transcript added to the rows and the totals, so
        that a reader can take it back out."""
        self.helpers["transcripts"] += 1
        self.helpers["calls"] += sum(len(c["tools"]) for c in calls)
        self.helpers["output_tokens"] += sum(c["counts"][4] for c in calls if c["counts"])
        self.helpers["seconds_in_tool"] += sum(tool_seconds(t) for c in calls for t in c["tools"])

    def add_reads(self, tool, seen):
        paths, partial = reads_of(tool)
        if not paths:
            return
        new = [p for p in paths if p not in seen]
        self.reads["calls"] += 1
        self.reads["files_named"] += len(paths)
        self.reads["rereads"] += len(paths) - len(new)
        self.reads["partial"] += int(partial)
        size = tool["result"]["bytes"] if tool["result"] else 0
        self.reads["bytes"] += size
        if not new:
            self.reads["reread_bytes"] += size
        seen.update(paths)
        self.files.update(paths)
        self.per_path.update(paths)

    def add_context(self, calls):
        """The context call by call, and its growth shared out to the tools
        whose answers it carried; a call whose usage says it sent nothing
        is an error line, not a call."""
        measured = [c for c in calls if c["counts"] and context_of(c["counts"]) > 0]
        sizes = [context_of(c["counts"]) for c in measured]
        self.sizes = sizes
        cached = sum(c["counts"][3] for c in measured)
        self.context = {"calls": len(sizes), "first": sizes[0] if sizes else 0, "peak": max(sizes, default=0),
                        "last": sizes[-1] if sizes else 0, "mean": round(sum(sizes) / len(sizes)) if sizes else 0,
                        "cached_share": round(cached / sum(sizes), 3) if sum(sizes) else 0,
                        "compactions": 0, "drops": 0, "growth": 0}
        for call, size, following in zip(measured, sizes, sizes[1:]):
            if following < size * (1 - DROP):
                self.context["drops"] += 1
                continue
            growth = max(0, following - size - call["counts"][4])
            self.context["growth"] += growth
            self.share_growth(call, growth)

    def share_growth(self, call, growth):
        answered = [t for t in call["tools"] if t["result"] is not None]
        if not answered:
            return
        weights = [t["result"]["bytes"] for t in answered]
        if not sum(weights):
            weights = [1] * len(answered)
        for tool, weight in zip(answered, weights):
            self.by_tool[tool["name"]]["context_tokens_estimated"] += round(growth * weight / sum(weights))

    def add_time(self, entries, calls):
        times = [t for t in (usage.moment(e.get("timestamp")) for e in entries) if t is not None]
        if not times:
            return
        self.started, self.ended = min(times), max(times)
        intervals = union([(t["at"], t["result"]["at"]) for c in calls for t in c["tools"]
                           if t["result"] and t["at"] is not None and t["result"]["at"] is not None
                           and t["result"]["at"] >= t["at"]])
        self.seconds["span"] = self.ended - self.started
        self.seconds["active"] = usage.active_seconds(entries)
        self.seconds["in_tools"] = sum(end - start for start, end in intervals)
        self.seconds["outside_tools"] = active_outside(times, intervals)
        self.seconds["thinking"] = sum(c["thinking_ms"] for c in calls) / 1000

    def record(self, session):
        usd, unpriced = economics.priced(self.models)
        tokens = economics.by_kind(self.models)
        tokens["total"] = sum(tokens.values())
        folded = {}
        for name, row in self.by_tool.items():
            into = folded.setdefault(tool_kind(name), empty_tool())
            for key, value in row.items():
                into[key] += value
        by_tool = {name: {k: round(v) for k, v in row.items()} for name, row in sorted(folded.items())}
        return {
            "schema": SCHEMA, "session": session, "helpers": {k: round(v) for k, v in self.helpers.items()},
            "started_at": stamp(self.started), "ended_at": stamp(self.ended),
            "prompts": self.prompts, "attachments": self.attachments,
            "seconds": {k: round(self.seconds[k]) for k in SECONDS_KEYS},
            "models": {m: list(c) for m, c in sorted(self.models.items())},
            "tokens": tokens, "estimated_usd": round(usd, 4), "tokens_unpriced": unpriced,
            "prices_as_of": economics.PRICES_AS_OF,
            "context": dict(self.context),
            "calls": {"total": self.calls["total"], "failed": self.calls["failed"],
                      "unanswered": self.calls["unanswered"], "repeated": self.calls["repeated"], "by_tool": by_tool},
            "shell": {k: self.shell[k] for k in SHELL_USES},
            "reads": {k: (len(self.files) if k == "files" else self.reads[k]) for k in READS_KEYS},
        }


def stamp(moment):
    if moment is None:
        return None
    return datetime.fromtimestamp(moment, timezone.utc).replace(microsecond=0).isoformat().replace("+00:00", "Z")


def profile(transcript, helpers=True):
    """The record of one session: its transcript and, by default, its
    helpers', those its transcript holds and those in files beside it."""
    built = Profile()
    own, held = threads_of(entries_of(transcript))
    built.add_transcript(own, main=True)
    if helpers:
        for entries in held:
            built.add_transcript(entries, main=False)
        for helper in usage.helper_transcripts(transcript):
            built.add_transcript(entries_of(helper), main=False)
    stem = Path(transcript).stem
    return built.record(usage.session_token(stem) if SESSION_STEM.fullmatch(stem) else None), built


# --- the trailer --------------------------------------------------------------------

def session_counts(transcript):
    """A session's counts so far, as its commits' trailers add them up: the
    sums, from the record with its helpers', and the context of each of the
    session's own calls in order, from which the peak since a commit is
    taken."""
    record, built = profile(transcript)
    calls, reads, context = record["calls"], record["reads"], record["context"]
    sums = {"calls": calls["total"], "failed": calls["failed"], "repeated": calls["repeated"],
            "reads": reads["files_named"], "rereads": reads["rereads"], "read_bytes": reads["bytes"],
            "reread_bytes": reads["reread_bytes"], "compactions": context.get("compactions", 0),
            "drops": context.get("drops", 0), "tool_seconds": record["seconds"]["in_tools"]}
    return {"sums": sums, "sizes": list(built.sizes)}


def kept(now):
    """What is kept of a session's counts once they are on a commit: the
    sums, and how many of its calls the peak has looked at."""
    return {"sums": dict(now["sums"]), "api_calls": len(now["sizes"])}


def known(before):
    """Whether counts kept at a commit can be counted from: as kept() keeps
    them, every sum and the count of calls a whole number, none below
    nothing; or empty, as kept after a session's first commit whose
    profile was not read, from which the next line counts everything.
    Anything else, a file damaged or written by another version, could
    make a line claim calls an earlier line carried."""
    if before == {}:
        return True
    sums, calls = before.get("sums"), before.get("api_calls")
    return (isinstance(sums, dict) and all(is_count(sums.get(key)) for key in SUMS) and is_count(calls))


def is_count(value):
    return isinstance(value, int) and not isinstance(value, bool) and value >= 0


def since(now, before):
    """The profile since the counts kept() kept at the session's last
    commit, or since its start when before is empty: each sum's rise, and
    the largest context of the calls made since."""
    earlier = before.get("sums") or {}
    spent = {key: max(0, int(now["sums"].get(key) or 0) - int(earlier.get(key) or 0)) for key in SUMS}
    spent["peak"] = max(now["sizes"][int(before.get("api_calls") or 0):], default=0)
    return spent


def trailer(session, spent):
    """The trailer's value for the profile since a commit, its session as
    the usage line gives it."""
    counts = [f"{key}={int(spent[key])}" for key in TRAILER_FIELDS]
    return " ".join([TRAILER_VERSION, f"session={usage.session_token(session)}", *counts])


def parse_trailer(text):
    """Every Symphony-Profile trailer in a commit message, as (session,
    {field: count}): read in the message's last paragraph alone, where git
    and the hook put a commit's trailers, so that a line quoted in the
    body counts for nothing. So it is read from the commits the hook
    wrote, a pull request's own, and never from a squash commit, which
    moves each commit's trailers into its body. A part that cannot be read
    is passed over: a count past usage.count_in's bound, as the usage
    line's are, a session of more than the token's eight characters or of
    others, and a key said again, of which the first counts."""
    found = []
    for match in TRAILER_LINE.finditer(trailer_block(text)):
        session, counts, said = "", {}, set()
        for part in match.group(1).split():
            key, _, value = part.partition("=")
            if key in said:
                continue
            said.add(key)
            if key == "session":
                session = value if TRAILER_SESSION.fullmatch(value) else ""
            elif key in TRAILER_FIELDS and usage.count_in(value) is not None:
                counts[key] = usage.count_in(value)
        found.append((session, counts))
    return found


def trailer_block(text):
    """A commit message's last paragraph."""
    return PARAGRAPH_BREAK.split((text or "").strip())[-1]


# --- the sums over pull requests ----------------------------------------------------

def measured(counts):
    """Whether a line's or a run's counts are a profile: any call made or
    any context held. A run that made no API call records every count 0,
    and a line or a run of zeros counted as profiled would pull the median
    peak down and the per-merge ratios with it, where the economics count
    a pull request as measured only with tokens above 0."""
    return bool(counts) and bool(counts.get("calls") or counts.get("peak"))


def commit_lines(commit):
    """The Symphony-Profile lines of a commit as GitHub lists it, those
    that are a profile (measured)."""
    message = ((commit or {}).get("commit") or {}).get("message") or ""
    return [(session, counts) for session, counts in parse_trailer(message) if measured(counts)]


def pull_profile(work, records):
    """What one pull request's profile adds up to (an economics.Work): the
    lines on its commits and the profiles of the Actions runs that name it
    or its claim's issue, as economics.spent_on takes their usage; each sum
    added up and peak the largest; None where none of them has a profile."""
    found = [counts for commit in work.commits for _, counts in commit_lines(commit)]
    found += [record["profile"] for record in records if record.get("profile") and economics.names(record, work)]
    if not found:
        return None
    total = {key: sum(counts.get(key, 0) for counts in found) for key in SUMS}
    total["peak"] = max(counts.get("peak", 0) for counts in found)
    return total


def run_sums(records):
    """The sums of the profiles of some Actions runs."""
    found = [record["profile"] for record in records if record.get("profile")]
    return {key: sum(counts.get(key, 0) for counts in found) for key in SUMS}


def group(works, profiles, extra=None):
    """One work setting, or one side of the comparison, as GROUP_KEYS: its
    sums are its profiled pull requests' (those with any profile) and what
    it spends beside them (extra), so they are divided by its profiled
    merges and never by merges whose profile is unknown; peak_median is the
    median of its profiled pull requests' peaks, rounded half up as every
    median of the field metrics is (economics.median_seconds)."""
    profiled = [w for w in works if profiles.get(w.number) is not None]
    row = {"pull_requests": len(works), "merged": sum(1 for w in works if w.merged is not None),
           "profiled": len(profiled), "profiled_merged": sum(1 for w in profiled if w.merged is not None)}
    row.update({key: sum(profiles[w.number][key] for w in profiled) + (extra or {}).get(key, 0) for key in SUMS})
    row["peak_median"] = economics.median_seconds([profiles[w.number]["peak"] for w in profiled])
    return row


def ratios(alone, parallel):
    """Each comparison as a ratio over 1 where working in parallel did
    better: fewer calls, reads and re-reads per profiled merge, and a
    smaller median peak; None where either side has nothing to compare."""
    def per_merge(side, key):
        return economics.per(side[key], side["profiled_merged"])

    per_merges = [economics.inverse(per_merge(alone, key), per_merge(parallel, key))
                  for key in ("calls", "reads", "rereads")]
    return dict(zip(RATIO_KEYS, (*per_merges, economics.inverse(alone["peak_median"], parallel["peak_median"]))))


def coordination(works, profiles, now, others=(), window_records=()):
    """The coordination block's siblings: the window's pull requests split
    into the same two sides as economics.coordination splits them, alone
    and parallel (economics.sides), each as group gives it, with the
    profiles of the coordinator's runs that name none of them on the
    parallel side (economics.overhead), as the price of coordinating, and
    the ratios."""
    alone, parallel = economics.sides(works, now, others)
    alone_side = group(alone, profiles)
    parallel_side = group(parallel, profiles, run_sums(economics.overhead(works, window_records)))
    return {"alone": alone_side, "parallel": parallel_side, "ratios": ratios(alone_side, parallel_side)}


# --- the report ---------------------------------------------------------------------

def share(part, whole):
    return f"{100 * part / whole:.0f}%" if whole else "n/a"


def report(record):
    """The record as a person reads it."""
    seconds, context, calls, reads, tokens = (record[k] for k in ("seconds", "context", "calls", "reads", "tokens"))
    helpers = record["helpers"]
    who = f"session {record['session']}" if record["session"] else "a run"
    when = f"{record['started_at']} to {record['ended_at']}" if record["started_at"] else "no clock"
    lines = [f"Profile of {who} ({when}, {helpers['transcripts']} helper transcript"
             f"{'s' if helpers['transcripts'] != 1 else ''}, {record['prompts']} prompts, "
             f"{record['attachments']} attachments)",
             f"Time: span {seconds['span']:,} s, active {seconds['active']:,} s; in tools {seconds['in_tools']:,} s "
             f"({share(seconds['in_tools'], seconds['active'])} of active), outside tools "
             f"{seconds['outside_tools']:,} s ({share(seconds['outside_tools'], seconds['active'])}), "
             f"the model and the harness between lines; the harness counted {seconds['thinking']:,} s of thinking",
             f"Tokens: input {tokens['input']:,}, cache writes {tokens['cache_write']:,}, cache reads "
             f"{tokens['cache_read']:,}, output {tokens['output']:,}; estimated USD {record['estimated_usd']:.2f} "
             f"at list prices of {record['prices_as_of']} ({', '.join(record['models']) or 'no model'})",
             f"Context: {context.get('calls', 0):,} API calls, first {context.get('first', 0):,} tokens, peak "
             f"{context.get('peak', 0):,}, last {context.get('last', 0):,}, mean {context.get('mean', 0):,}; "
             f"{share(context.get('cached_share', 0), 1)} served from cache; {context.get('compactions', 0)} "
             f"compactions, {context.get('drops', 0)} drops; grew by {context.get('growth', 0):,} tokens, "
             f"shared out to the tools below (an estimate)",
             f"Tool calls: {calls['total']:,}, {calls['failed']:,} failed, {calls['unanswered']:,} unanswered, "
             f"{calls['repeated']:,} repeated (the same tool with the same input again)"]
    lines.append(f"  {'tool':34} {'calls':>6} {'failed':>6} {'in KB':>8} {'out KB':>8} {'seconds':>8} {'context':>10}")
    for name, row in calls["by_tool"].items():
        lines.append(f"  {name[:34]:34} {row['calls']:>6,} {row['failed']:>6,} {row['input_bytes'] / 1024:>8.1f} "
                     f"{row['result_bytes'] / 1024:>8.1f} {row['seconds_in_tool']:>8,} "
                     f"{row['context_tokens_estimated']:>10,}")
    if helpers["transcripts"]:
        lines.append(f"Helpers: {helpers['transcripts']} transcripts added {helpers['calls']:,} calls, "
                     f"{helpers['seconds_in_tool']:,} s and {helpers['output_tokens']:,} output tokens to the rows "
                     f"and totals above; the context and the time are the session's own")
    shell = record["shell"]
    lines.append("Shell commands: " + ", ".join(f"{k} {shell[k]:,}" for k in SHELL_USES))
    lines.append(f"Reads: {reads['calls']:,} calls naming {reads['files_named']:,} files ({reads['files']:,} "
                 f"distinct), {reads['rereads']:,} re-reads ({share(reads['rereads'], reads['files_named'])}), "
                 f"{reads['partial']:,} partial; {reads['bytes'] / 1024:,.1f} KB returned, "
                 f"{reads['reread_bytes'] / 1024:,.1f} KB of it by calls that only re-read")
    return "\n".join(lines) + "\n"


def most_read(built, count):
    rows = sorted(built.per_path.items(), key=lambda item: (-item[1], item[0]))[:count]
    return "".join(f"{n:>5}  {path}\n" for path, n in rows)


def most_called(built, count):
    rows = sorted(built.by_tool.items(), key=lambda item: (-item[1]["calls"], item[0]))[:count]
    return "".join(f"{row['calls']:>5}  {name}\n" for name, row in rows)


def main(argv=None):
    parser = argparse.ArgumentParser(prog="python3 tools/profile.py", description=__doc__.split("\n\n")[0])
    parser.add_argument("transcript", help="a session's transcript (.jsonl), or a run's execution file")
    parser.add_argument("--json", help="write the record here, aggregates only")
    parser.add_argument("--paths", type=int, default=0, help="also print the N files read most often")
    parser.add_argument("--tools", type=int, default=0, help="also print the N tools called most often, by name")
    parser.add_argument("--no-helpers", action="store_true", help="leave the helpers' transcripts out")
    args = parser.parse_args(argv)
    if not Path(args.transcript).is_file():
        print(f"profile: no transcript at {args.transcript}", file=sys.stderr)
        return 2
    record, built = profile(args.transcript, helpers=not args.no_helpers)
    if not built.lines:
        print(f"profile: no entries read from {args.transcript}", file=sys.stderr)
        return 2
    sys.stdout.write(report(record))
    if args.paths:
        sys.stdout.write(f"Files read most often (not in the record):\n{most_read(built, args.paths)}")
    if args.tools:
        sys.stdout.write(f"Tools called most often, by name (not in the record):\n{most_called(built, args.tools)}")
    if args.json:
        out = Path(args.json)
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text(json.dumps(record, indent=2, sort_keys=True) + "\n", encoding="ascii")
        print(f"record written to {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
