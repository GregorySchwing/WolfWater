"""How much of each Claude account's limits is used, and which account a
run should take: the pool of slots the token wizard stores
(scripts/tokens.py, docs/setup/claude-tokens.md).

A subscription has a five-hour window and a seven-day window, and some
plans a weekly limit for one model, such as Fable. GitHub never gives a
secret back, so nothing can read a slot's usage by looking at its token;
three kinds of reading are kept instead, one file per slot and kind on the
branch `pool`:

  run    a Claude run in Actions (the workflows record it once step 2 of
         docs/research/claude-token-pool.md wires them). Claude Code reads the
         anthropic-ratelimit-unified-* headers of each response and, in
         its stream-json output, emits a rate_limit_event whose
         unifiedWindows carry the five-hour and seven-day windows'
         utilization and reset; `observe` takes the last one from the
         run's execution file, at no cost.
  probe  `probe` asks Claude for one word with a slot's token and reads
         the same event: a fresh reading for a small part of the window,
         run by .github/workflows/pool.yml when a person asks for one.
  cswap  claude-swap on the owner's machine reads each of its accounts'
         usage page, which a token from claude setup-token cannot (it is
         inference-only), and that page has the per-model weekly limits;
         `python3 scripts/tokens.py usage` maps cswap's accounts to slots
         as the wizard remembered them and sends the readings here.

No reading holds a token, an e-mail or a word of a conversation: a slot's
number, when, the windows' fractions and resets, and per-model limits by
the model's name. `pick` is the auto mode: of the slots in the pool, the
one with the most room in its five-hour window, then its seven-day
window, among those no window holds at or past the threshold; when every
slot is held, the one whose hold ends first, which the caller is told.
A window whose reset has passed counts as empty, since it is.

  python3 tools/pool.py observe --execution-file F --slot N --out DIR [--source probe]
  python3 tools/pool.py probe --slot N --out FILE
  python3 tools/pool.py publish FILE...        onto the branch pool (Actions)
  python3 tools/pool.py show DIR               the slots as the panel shows them
  python3 tools/pool.py pick DIR --slots 1,2,3 [--threshold 0.9]
  python3 tools/pool.py choose --worker NAME --pool 1,2,3   the slot a worker's run takes
  python3 tools/pool.py settings --slot main --settings JSON   a slot's model, effort and fallback

Which account a worker in Actions runs on is the file
.github/symphony/accounts.json, {"default": "main", "workers": {"NAME":
"auto"}}: main is CLAUDE_CODE_OAUTH_TOKEN, a number from 1 to 32 is that
slot's SYMPHONY_CLAUDE_<slot>, auto is pick over the pool, the slots
of SYMPHONY_CLAUDE_ACCOUNTS, and a list of slots, ["2", "5"], is pick
over those of them the pool holds, main only where the list names it, so
the Tokens tab can give each worker or agent its own tokens. With no file, or no line for the worker and
no default, the worker runs on main, as every worker did before.

How each slot runs is the variable SYMPHONY_CLAUDE_SETTINGS, which the
control panel's Tokens tab writes when a token's tile is set: a JSON
object by slot, {"1": {"model": "claude-opus-5-5", "effort": "high",
"fallback": "claude-sonnet-5-5", "resume": true}}. The model, when set,
is the one the slot's runs take in place of the worker's; the effort is
Claude Code's --effort; the fallback is the lesser model a run takes when
the model's own weekly limit is spent while the slot's other windows have
room, and Claude Code's --fallback-model when the model is busy. resume
false keeps a slot out of auto, once a reading after resumedAt shows it at
its limit, until a person presses Resume, which sets resumedAt anew;
paused true keeps it out of auto whatever its readings say. A setting that
does not hold together is left out with a word on stderr, so one bad
value never stops a worker; a slot with no settings runs as before.
"""

import argparse
import json
import re
import subprocess
import sys
import tempfile
from datetime import datetime, timedelta, timezone
from pathlib import Path

VERSION = 1
BRANCH = "pool"
SOURCES = ("run", "probe", "cswap")
# An organization's pool holds every member's tokens, so it has room for 32;
# the workflows name each slot's secret by its fixed name (workers.yml).
SLOTS = tuple(str(slot) for slot in range(1, 33)) + ("main",)
WINDOWS = ("five_hour", "seven_day")
THRESHOLD = 0.9
UNKNOWN_HOLD = timedelta(hours=5)
FILE_NAME = re.compile(r"^slot-(main|[1-9]|[12][0-9]|3[0-2])\.(run|probe|cswap)\.json$")
NAME_SHAPE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9 ._-]{0,39}$")
PROBE = ("claude", "-p", "Reply with the one word ok.", "--model", "haiku", "--max-turns", "1",
         "--output-format", "stream-json", "--verbose")
PROBE_TIMEOUT = 120
PUSH_ATTEMPTS = 5
ASSIGNMENTS = ".github/symphony/accounts.json"
CHOICES = SLOTS + ("auto",)
EFFORTS = ("low", "medium", "high", "xhigh", "max")
# The models a token's tile offers. A setting names one of them or is
# left out: a name of the right shape that no model answers to would fail
# every run on the token, and a setting must never stop a run. The control
# panel's tile list (assets/js/symphony-panel.js, TILE_MODELS) is the same.
KNOWN_MODELS = ("claude-fable-5-1", "claude-opus-5-5", "claude-sonnet-5-5", "claude-haiku-5-5")
# The worker's own model, as the registry takes it (tools/sdlc_protocol/workers.py).
WORKER_MODEL = re.compile(r"[A-Za-z0-9._:\[\]-]{0,80}")
# A token's own id, which stays when its slot changes, and the GitHub login
# of the person it belongs to, as the Tokens tab writes them: the id labels
# the issues the token worked (token:<id>) for the leaderboard.
TOKEN_ID = re.compile(r"[0-9a-f]{8}")
LOGIN = re.compile(r"[A-Za-z0-9](?:[A-Za-z0-9]|-(?=[A-Za-z0-9])){0,38}")
# The per-model weekly limits by the family a model's name holds.
FAMILIES = ("fable", "opus", "sonnet", "haiku")


class PoolError(Exception):
    """What could not be read or written; never holds a token."""


# Times.

def iso(moment):
    return moment.astimezone(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def parse_time(value):
    """An ISO time or unix seconds, as Claude Code and cswap write them, or None."""
    if isinstance(value, bool):
        return None
    if isinstance(value, (int, float)):
        try:
            return datetime.fromtimestamp(value, tz=timezone.utc)
        except (OverflowError, OSError, ValueError):
            return None
    if isinstance(value, str) and len(value) <= 40:
        try:
            moment = datetime.fromisoformat(value.replace("Z", "+00:00"))
        except ValueError:
            return None
        return moment if moment.tzinfo else moment.replace(tzinfo=timezone.utc)
    return None


def fraction(value):
    """A window's utilization as a fraction, 0 and up, or None."""
    if isinstance(value, bool) or not isinstance(value, (int, float)) or value != value:
        return None
    return max(0.0, min(float(value), 10.0))


def window(utilization, resets):
    used = fraction(utilization)
    if used is None:
        return None
    found = {"utilization": round(used, 4)}
    moment = parse_time(resets)
    if moment is not None:
        found["resetsAt"] = iso(moment)
    return found


# Readings.

def reading(slot, source, windows, scoped=(), status=None, now=None):
    record = {"version": VERSION, "slot": slot, "source": source, "at": iso(now or datetime.now(timezone.utc)),
              "windows": {name: found for name, found in windows.items() if found}}
    if scoped:
        record["scoped"] = list(scoped)
    if status in ("allowed", "allowed_warning", "rejected"):
        record["status"] = status
    return record


def messages_in(text):
    """The messages of an execution file: one JSON array, or a JSON object
    a line, as stream-json writes them."""
    text = text.strip()
    if text.startswith("["):
        try:
            found = json.loads(text)
        except json.JSONDecodeError:
            found = None
        if isinstance(found, list):
            return [item for item in found if isinstance(item, dict)]
    found = []
    for line in text.splitlines():
        try:
            item = json.loads(line)
        except json.JSONDecodeError:
            continue
        if isinstance(item, dict):
            found.append(item)
    return found


# The per-model weekly windows a rate_limit_event can name, by the model.
MODEL_WINDOWS = {"seven_day_opus": "Opus", "seven_day_sonnet": "Sonnet"}


def windows_from_event(info):
    """The windows a rate_limit_event's rate_limit_info carries, and its
    per-model ones: every window under unifiedWindows, and the top-level
    one where it names one. A refusal names the window that refused and
    may give no utilization: that window is full until its reset."""
    windows, scoped = {}, []
    unified = info.get("unifiedWindows")
    if isinstance(unified, dict):
        for name in WINDOWS:
            each = unified.get(name)
            if isinstance(each, dict):
                windows[name] = window(each.get("utilization"), each.get("resetsAt"))
    kind = info.get("rateLimitType")
    used = info.get("utilization")
    if info.get("status") == "rejected":
        used = 1.0 if fraction(used) is None else max(fraction(used), 1.0)
    if kind in WINDOWS and (info.get("status") == "rejected" or not windows.get(kind)):
        windows[kind] = window(used, info.get("resetsAt"))
    elif kind in MODEL_WINDOWS:
        found = window(used, info.get("resetsAt"))
        if found:
            found["name"] = MODEL_WINDOWS[kind]
            scoped.append(found)
    return windows, scoped


def observed(messages, slot, source, now=None):
    """The reading of the last rate_limit_event that names a window, or None."""
    for message in reversed(messages):
        info = message.get("rate_limit_info") if message.get("type") == "rate_limit_event" else None
        if not isinstance(info, dict):
            continue
        windows, scoped = windows_from_event(info)
        if any(windows.values()) or scoped:
            return reading(slot, source, windows, scoped, status=info.get("status"), now=now)
    return None


def from_cswap(row, slot, now=None):
    """A reading from one account row of `cswap list --json`: its usage, or
    the last good one when the newest fetch failed, with when it was read."""
    usage = row.get("usage") if isinstance(row.get("usage"), dict) else row.get("lastGoodUsage")
    if not isinstance(usage, dict):
        return None
    windows = {}
    for name, key in (("five_hour", "fiveHour"), ("seven_day", "sevenDay")):
        each = usage.get(key)
        if isinstance(each, dict) and fraction(each.get("pct")) is not None:
            windows[name] = window(each["pct"] / 100, each.get("resetsAt"))
    scoped = []
    for each in usage.get("scoped") or []:
        if isinstance(each, dict) and isinstance(each.get("name"), str) and NAME_SHAPE.fullmatch(each["name"]) \
                and fraction(each.get("pct")) is not None:
            found = window(each["pct"] / 100, each.get("resetsAt"))
            found["name"] = each["name"]
            scoped.append(found)
    if not windows and not scoped:
        return None
    when = parse_time(row.get("usageFetchedAt") or row.get("lastGoodFetchedAt")) or now
    return reading(slot, "cswap", windows, scoped, now=when)


def file_name(record):
    return f"slot-{record['slot']}.{record['source']}.json"


def write(path, record):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    Path(path).write_text(json.dumps(record, indent=2, sort_keys=True) + "\n", encoding="ascii")


# The slots as they stand.

def load(folder):
    """Every reading in a checkout of the branch pool, by slot."""
    found = {}
    for path in sorted(Path(folder).glob("slot-*.json")):
        if not FILE_NAME.match(path.name):
            continue
        try:
            record = json.loads(path.read_text(encoding="utf-8"))
        except (OSError, ValueError):
            continue
        if well_formed(record, FILE_NAME.match(path.name)):
            found.setdefault(record["slot"], []).append(record)
    return found


def well_formed(record, name):
    """A reading as this file writes it, of the slot and kind its file's
    name says; anyone who may write the branch could put anything there,
    and one stray file must not stop the panel or the picker."""
    if not isinstance(record, dict) or record.get("version") != VERSION or not parse_time(record.get("at")):
        return False
    if (record.get("slot"), record.get("source")) != name.groups():
        return False
    windows, scoped = record.get("windows", {}), record.get("scoped", [])
    if not isinstance(windows, dict) or not isinstance(scoped, list):
        return False
    return all(isinstance(each, dict) for each in list(windows.values()) + scoped)


def current(records, now):
    """A slot's windows from its newest reading of each, a window past its
    reset counting as empty, and its per-model limits from the newest
    reading that has them."""
    windows, scoped = {}, []
    for record in sorted(records, key=lambda each: parse_time(each["at"])):
        for name, each in (record.get("windows") or {}).items():
            if name in WINDOWS and isinstance(each, dict) and fraction(each.get("utilization")) is not None:
                windows[name] = dict(each, at=record["at"], source=record["source"])
        if record.get("scoped"):
            scoped = [dict(each, at=record["at"]) for each in record["scoped"] if isinstance(each, dict)]
    for each in list(windows.values()) + scoped:
        resets = parse_time(each.get("resetsAt"))
        if resets is not None and resets <= now:
            each["utilization"] = 0.0
            each["reset"] = True
    return windows, scoped


def held_until(windows, scoped, threshold, now):
    """When the slot is free again: the latest reset of a window at or past
    the threshold, or None when nothing holds it. A full window with no
    reset known holds it for the length of a five-hour window from its
    reading, so a slot is tried again rather than held for good on a
    reading nobody may renew."""
    until = None
    for each in list(windows.values()) + scoped:
        used = fraction(each.get("utilization"))
        if used is not None and used >= threshold:
            resets = parse_time(each.get("resetsAt"))
            if resets is None:
                read = parse_time(each.get("at")) or now
                resets = read + UNKNOWN_HOLD
            if resets > now and (until is None or resets > until):
                until = resets
    return until


def family(model):
    """The per-model limit's name a model's runs count against, or None
    where the name says no family, as the subscription's default does."""
    for each in FAMILIES:
        if each in (model or ""):
            return each.capitalize()
    return None


def scoped_for(scoped, model):
    """The per-model limits that hold a run of the model: its family's
    alone, or every one where the family is not known, as before."""
    name = family(model)
    return [each for each in scoped if name is None or each.get("name") == name]


def tripped(records, since, threshold):
    """Whether a reading taken at or after since saw a window at its limit,
    for a slot that does not resume by itself."""
    since = parse_time(since)
    for record in records:
        if since is not None and parse_time(record["at"]) < since:
            continue
        seen = list((record.get("windows") or {}).values()) + list(record.get("scoped") or [])
        if any(isinstance(each, dict) and (fraction(each.get("utilization")) or 0.0) >= threshold for each in seen):
            return True
    return False


def slot_state(records, setting, model, threshold, now):
    """(until, model, why) for one slot: until is when a held slot is free
    again (datetime.max for one a person must resume), model the one its run
    takes, and why a word on the fallback or the hold."""
    setting = setting or {}
    model = setting.get("model") or model
    if setting.get("paused"):
        return datetime.max.replace(tzinfo=timezone.utc), model, "paused"
    if setting.get("resume") is False and tripped(records, setting.get("resumedAt"), threshold):
        return datetime.max.replace(tzinfo=timezone.utc), model, "at its limit, waiting for Resume"
    windows, scoped = current(records, now)
    general = held_until(windows, [], threshold, now)
    if general:
        return general, model, "at its limit"
    own = held_until({}, scoped_for(scoped, model), threshold, now)
    if not own:
        return None, model, ""
    fallback = setting.get("fallback")
    if fallback and fallback != model and not held_until({}, scoped_for(scoped, fallback), threshold, now):
        return None, fallback, f"{model or 'the model'} is at its weekly limit, so {fallback}"
    return own, model, "the model is at its weekly limit"


def pick_run(found, slots, settings=None, model="", threshold=THRESHOLD, now=None, spread=0):
    """(slot, until, model, why): the slot with the most room, the model
    its run takes and why it is not the usual one; for every slot held,
    the one whose hold ends first, with until set. Workers of one run all
    read the pool before any of them writes it, so spread, each worker's
    place in the run, deals them the free slots in turn, the one with the
    most room to the first."""
    now = now or datetime.now(timezone.utc)
    settings = settings or {}
    if not slots:
        raise PoolError("the pool names no slot")
    free, held = [], []
    for order, slot in enumerate(slots):
        records = found.get(slot, [])
        until, taken, why = slot_state(records, settings.get(slot), model, threshold, now)
        windows, _ = current(records, now)
        room = tuple(fraction((windows.get(name) or {}).get("utilization")) or 0.0 for name in WINDOWS)
        # A slot no reading names comes after every slot one does: its
        # secret may be unset or revoked, which no reading would show.
        (held if until else free).append((until, (not windows, room), order, slot, taken, why))
    if free:
        free.sort(key=lambda each: (each[1], each[2]))
        _, _, _, slot, taken, why = free[spread % len(free)]
        return slot, None, taken, why
    until, _, _, slot, taken, why = min(held, key=lambda each: (each[0], each[2]))
    return slot, until, taken, why


def pick(found, slots, threshold=THRESHOLD, now=None, spread=0):
    """(slot, None) for the slot with the most room, or (slot, until) for
    the one whose hold ends first when every slot is held."""
    slot, until, _, _ = pick_run(found, slots, threshold=threshold, now=now, spread=spread)
    return slot, until


def parse_slots(text):
    """A list of slots as SYMPHONY_CLAUDE_ACCOUNTS holds it, in its order, each once."""
    slots = []
    for word in re.split(r"[,\s]+", text or ""):
        if word and word not in SLOTS:
            raise PoolError(f"a slot is a number from 1 to 32 or main, not {word[:20]!r}")
        if word and word not in slots:
            slots.append(word)
    return slots


def parse_settings(text, complain=None):
    """SYMPHONY_CLAUDE_SETTINGS as {slot: setting}, each value checked and
    one that is not what it should be left out and named to complain."""
    complain = complain or (lambda words: print(f"pool: {words}", file=sys.stderr))
    if not (text or "").strip():
        return {}
    try:
        data = json.loads(text)
    except json.JSONDecodeError:
        complain("SYMPHONY_CLAUDE_SETTINGS is not JSON; every slot runs as it would with none")
        return {}
    if not isinstance(data, dict):
        complain("SYMPHONY_CLAUDE_SETTINGS is not an object by slot; every slot runs as it would with none")
        return {}
    settings = {}
    for slot, each in data.items():
        if slot not in SLOTS or not isinstance(each, dict):
            complain(f"SYMPHONY_CLAUDE_SETTINGS: {str(slot)[:20]!r} is no slot's settings; left out")
            continue
        settings[slot] = setting_of(slot, each, complain)
    return settings


def setting_of(slot, each, complain):
    found = {}
    for key in ("model", "fallback"):
        value = each.get(key)
        if value in (None, ""):
            continue
        if value in KNOWN_MODELS:
            found[key] = value
        else:
            complain(f"slot {slot}: the {key} is one of {', '.join(KNOWN_MODELS)}; left out")
    effort = each.get("effort")
    if effort not in (None, ""):
        if effort in EFFORTS:
            found["effort"] = effort
        else:
            complain(f"slot {slot}: the effort is one of {', '.join(EFFORTS)}; left out")
    for key, shape in (("id", TOKEN_ID), ("owner", LOGIN)):
        value = each.get(key)
        if value in (None, ""):
            continue
        if isinstance(value, str) and shape.fullmatch(value):
            found[key] = value
        else:
            complain(f"slot {slot}: the {key} is not of its shape; left out")
    for key in ("resume", "paused"):
        if isinstance(each.get(key), bool):
            found[key] = each[key]
    if each.get("resumedAt") is not None:
        if parse_time(each["resumedAt"]):
            found["resumedAt"] = each["resumedAt"]
        else:
            complain(f"slot {slot}: resumedAt is not a time; left out")
    return found


def summary(found, threshold=THRESHOLD, now=None):
    """The slots as the control panel shows them, as JSON."""
    now = now or datetime.now(timezone.utc)
    rows = []
    for slot in SLOTS:
        if slot not in found:
            continue
        windows, scoped = current(found[slot], now)
        until = held_until(windows, scoped, threshold, now)
        rows.append({"slot": slot, "windows": windows, "scoped": scoped,
                     "heldUntil": iso(until) if until else None,
                     "readAt": max(record["at"] for record in found[slot])})
    return {"version": VERSION, "at": iso(now), "threshold": threshold, "slots": rows}


# Which account a worker takes.

def assignments(path):
    """The file's assignments, checked: {worker: choice} and the default."""
    try:
        text = Path(path).read_text(encoding="utf-8")
    except FileNotFoundError:
        return {}, "main"
    try:
        data = json.loads(text)
    except json.JSONDecodeError as error:
        raise PoolError(f"{path} is not JSON ({error.msg}, line {error.lineno})") from None
    if not isinstance(data, dict) or set(data) - {"default", "workers"}:
        raise PoolError(f"{path} holds an object with default and workers, and nothing else")
    default = data.get("default", "main")
    workers = data.get("workers", {})
    if not is_choice(default):
        raise PoolError(f"{path}: the default is main, auto, a slot from 1 to 32 or a list of slots")
    if not isinstance(workers, dict) or not all(isinstance(name, str) and is_choice(choice) for name, choice in workers.items()):
        raise PoolError(f"{path}: each worker takes main, auto, a slot from 1 to 32 or a list of slots")
    return workers, default


def is_choice(choice):
    """main, auto, one slot, or a list of slots, main among them perhaps,
    each named once: the tokens a worker or an agent may take, which the
    Tokens tab gives it per repository."""
    if isinstance(choice, list):
        return 0 < len(choice) == len(set(choice)) and all(isinstance(slot, str) and slot in SLOTS for slot in choice)
    return isinstance(choice, str) and choice in CHOICES


def branch_readings(folder, cwd=None):
    """The readings on origin's branch pool, read into the folder given,
    which stays empty where the branch is not there yet."""
    folder = Path(folder)
    fetched = subprocess.run(["git", "fetch", "--quiet", "--no-tags", "origin",
                              f"+refs/heads/{BRANCH}:refs/remotes/origin/{BRANCH}"], cwd=cwd, capture_output=True)
    if fetched.returncode != 0:
        return folder
    for name in git("ls-tree", "--name-only", f"refs/remotes/origin/{BRANCH}", cwd=cwd).splitlines():
        if FILE_NAME.match(name):
            shown = subprocess.run(["git", "show", f"refs/remotes/origin/{BRANCH}:{name}"], cwd=cwd, capture_output=True)
            if shown.returncode == 0:
                (folder / name).write_bytes(shown.stdout)
    return folder


def choose_run(worker, pool_text, path=ASSIGNMENTS, folder=None, now=None, spread=0, settings=None, model=""):
    """What a worker's run takes, as a dict: the slot, how it was chosen,
    for a held auto pick when the first slot is free (until), and the
    model, effort and fallback the slot's settings give, the model the
    worker's own where the slot names none."""
    settings = settings or {}
    workers, default = assignments(path)
    choice = workers.get(worker, default)
    how = "assigned" if worker in workers else "default"
    if isinstance(choice, list):
        # A list is auto over its own slots, those the pool holds, and
        # falls back to main only where it names main.
        pooled = parse_slots(pool_text)
        slots = [slot for slot in choice if slot != "main" and slot in pooled]
        if slots:
            choice = "auto"
        elif "main" in choice:
            choice, how = "main", how + ", with none of its slots in SYMPHONY_CLAUDE_ACCOUNTS"
        else:
            raise PoolError(f"{worker} is given slots {', '.join(choice)}, and SYMPHONY_CLAUDE_ACCOUNTS holds none of them")
    elif choice == "auto":
        slots = [slot for slot in parse_slots(pool_text) if slot != "main"]
        if not slots:
            choice, how = "main", "auto, with no slot in SYMPHONY_CLAUDE_ACCOUNTS"
    if choice != "auto":
        return run_of(choice, how, None, (settings.get(choice) or {}).get("model") or model, settings)
    if folder is None:
        with tempfile.TemporaryDirectory(prefix="pool-") as scratch:
            found = load(branch_readings(scratch))
    else:
        found = load(folder)
    slot, until, taken, why = pick_run(found, slots, settings, model, now=now, spread=spread)
    return run_of(slot, "auto" + (f", {why}" if why and not until else ""), until, taken, settings)


def run_of(slot, how, until, model, settings):
    setting = settings.get(slot) or {}
    model = model or setting.get("model", "")
    fallback = setting.get("fallback", "")
    return {"slot": slot, "how": how, "until": until, "model": model or "",
            "effort": setting.get("effort", ""), "fallback": fallback if fallback != model else "",
            "id": setting.get("id", "")}


def choose(worker, pool_text, path=ASSIGNMENTS, folder=None, now=None, spread=0):
    """(slot, how, until): the slot a worker's run takes, how it was
    chosen, and, for auto with every slot held, when the first is free."""
    found = choose_run(worker, pool_text, path, folder, now, spread)
    return found["slot"], found["how"], found["until"]


# Probing and publishing, in Actions.

def probe(slot, run=subprocess.run):
    """One word from Claude on the token in CLAUDE_CODE_OAUTH_TOKEN, for the
    rate-limit event its answer brings; the token stays in the environment
    claude reads it from."""
    try:
        done = run(list(PROBE), capture_output=True, text=True, timeout=PROBE_TIMEOUT, stdin=subprocess.DEVNULL)
    except (OSError, subprocess.TimeoutExpired) as error:
        raise PoolError(f"claude could not be run for slot {slot} ({type(error).__name__})") from None
    record = observed(messages_in(done.stdout or ""), slot, "probe")
    if record is None:
        raise PoolError(f"claude answered for slot {slot} with no rate-limit reading (exit {done.returncode})")
    return record


def git(*args, cwd=None, stdin=None):
    done = subprocess.run(["git", *args], cwd=cwd, input=stdin, capture_output=True, text=True)
    if done.returncode != 0:
        raise PoolError(f"git {args[0]}: {(done.stderr or '').strip().splitlines()[-1:] or 'failed'}")
    return done.stdout.strip()


def newest_by_name(paths):
    """One path per reading's name, the newest by its own time: jobs of
    one run that took the same slot each leave a reading of that name,
    and the last to finish says most."""
    chosen = {}
    for path in paths:
        name = Path(path).name
        if not FILE_NAME.match(name):
            raise PoolError(f"{name} is not a reading's name")
        try:
            at = parse_time(json.loads(Path(path).read_text(encoding="utf-8")).get("at"))
        except (OSError, ValueError, AttributeError):
            at = None
        at = at or datetime.min.replace(tzinfo=timezone.utc)
        if name not in chosen or at > chosen[name][0]:
            chosen[name] = (at, path)
    return {name: path for name, (_, path) in chosen.items()}


def publish(paths, cwd=None):
    """The readings onto the branch pool, beside the ones there: one commit
    with no parent, as the branches board and metrics are kept, pushed
    with a lease on the commit read, and tried again when another run
    pushed first."""
    entries = {}
    for name, path in newest_by_name(paths).items():
        entries[name] = git("hash-object", "-w", str(path), cwd=cwd)
    for _ in range(PUSH_ATTEMPTS):
        previous = ""
        if subprocess.run(["git", "fetch", "--quiet", "--no-tags", "origin", f"+refs/heads/{BRANCH}:refs/remotes/origin/{BRANCH}"],
                          cwd=cwd, capture_output=True).returncode == 0:
            previous = git("rev-parse", f"refs/remotes/origin/{BRANCH}", cwd=cwd)
        tree = {}
        if previous:
            for line in git("ls-tree", previous, cwd=cwd).splitlines():
                meta, name = line.split("\t", 1)
                if FILE_NAME.match(name):
                    tree[name] = meta.split()[2]
        tree.update(entries)
        listing = "".join(f"100644 blob {blob}\t{name}\n" for name, blob in sorted(tree.items()))
        made = git("mktree", cwd=cwd, stdin=listing)
        commit = git("commit-tree", made, "-m", "Usage of the pool's slots", cwd=cwd)
        pushed = subprocess.run(["git", "push", "--quiet", f"--force-with-lease=refs/heads/{BRANCH}:{previous}",
                                 "origin", f"{commit}:refs/heads/{BRANCH}"], cwd=cwd, capture_output=True)
        if pushed.returncode == 0:
            return commit
        refused = (pushed.stderr or b"").decode("utf-8", "replace").strip().splitlines()
        if not any("stale info" in line or "fetch first" in line or "non-fast-forward" in line for line in refused):
            raise PoolError(f"the push to {BRANCH} was refused: {refused[-1] if refused else 'no reason given'}")
    raise PoolError(f"the branch {BRANCH} moved under every one of {PUSH_ATTEMPTS} pushes")


# The command line.

def out_path(out, record):
    """--out names a file when it ends in .json, else the folder the
    reading goes into under its own name, made if need be."""
    out = Path(out)
    return out if out.suffix == ".json" else out / file_name(record)


def build_parser():
    parser = argparse.ArgumentParser(prog="python3 tools/pool.py", description="The usage of the pool's Claude slots.")
    commands = parser.add_subparsers(dest="command", required=True)
    observe = commands.add_parser("observe", help="the last rate-limit reading of a run's execution file")
    observe.add_argument("--execution-file", required=True)
    observe.add_argument("--slot", required=True, choices=SLOTS)
    observe.add_argument("--source", default="run", choices=("run", "probe"), help="probe for a run that pool.yml started to read the slot")
    observe.add_argument("--out", required=True)
    probing = commands.add_parser("probe", help="one word from Claude, for a fresh reading")
    probing.add_argument("--slot", required=True, choices=SLOTS)
    probing.add_argument("--out", required=True)
    publishing = commands.add_parser("publish", help="readings onto the branch pool")
    publishing.add_argument("paths", nargs="+")
    showing = commands.add_parser("show", help="the slots as they stand, as JSON")
    showing.add_argument("folder")
    showing.add_argument("--threshold", type=float, default=THRESHOLD)
    picking = commands.add_parser("pick", help="the slot a run should take")
    picking.add_argument("folder")
    picking.add_argument("--slots", required=True)
    picking.add_argument("--threshold", type=float, default=THRESHOLD)
    choosing = commands.add_parser("choose", help="the slot a worker's run takes, as step outputs")
    choosing.add_argument("--worker", required=True)
    choosing.add_argument("--pool", default="")
    choosing.add_argument("--assignments", default=ASSIGNMENTS)
    choosing.add_argument("--spread", type=int, default=0, help="the worker's place in its run, from 0")
    choosing.add_argument("--settings", default="", help="SYMPHONY_CLAUDE_SETTINGS, the slots' settings as JSON")
    choosing.add_argument("--model", default="", help="the worker's own model, taken where the slot names none")
    showing_settings = commands.add_parser("settings", help="one slot's model, effort, fallback and id, as step outputs")
    showing_settings.add_argument("--slot", required=True, choices=SLOTS)
    showing_settings.add_argument("--settings", default="", help="SYMPHONY_CLAUDE_SETTINGS, the slots' settings as JSON")
    return parser


def main(argv=None):
    args = build_parser().parse_args(argv)
    try:
        if args.command == "observe":
            try:
                text = Path(args.execution_file).read_text(encoding="utf-8", errors="replace")
            except OSError:
                print(f"pool: no execution file at {args.execution_file}; nothing observed")
                return 0
            record = observed(messages_in(text), args.slot, args.source)
            if record is None:
                print("pool: the run carried no rate-limit reading; nothing observed")
                return 0
            write(out_path(args.out, record), record)
        elif args.command == "probe":
            record = probe(args.slot)
            write(out_path(args.out, record), record)
        elif args.command == "publish":
            print(f"pool: written as {publish(args.paths)}")
        elif args.command == "show":
            print(json.dumps(summary(load(args.folder), args.threshold), indent=2, sort_keys=True))
        elif args.command == "pick":
            slot, until = pick(load(args.folder), parse_slots(args.slots), args.threshold)
            print(slot)
            if until:
                print(f"pool: every slot is held; slot {slot} is free first, at "
                      f"{iso(until)}", file=sys.stderr)
        elif args.command == "choose":
            model = args.model if WORKER_MODEL.fullmatch(args.model or "") else ""
            run = choose_run(args.worker, args.pool, args.assignments, spread=max(args.spread, 0),
                             settings=parse_settings(args.settings), model=model)
            print(f"slot={run['slot']}")
            for key in ("model", "effort", "fallback", "id"):
                if run[key]:
                    print(f"{key}={run[key]}")
            until = run["until"]
            if until:
                print("wait=true")
                print(f"until={'when a person presses Resume on its tile' if until.year == 9999 else 'at ' + iso(until)}")
            slot = run["slot"]
            print(f"pool: {args.worker} runs on {'slot ' + slot if slot != 'main' else 'the main token'} ({run['how']})"
                  + (f" with {run['model']}" if run["model"] else "")
                  + (", though every slot of the pool is held" if until else ""), file=sys.stderr)
        elif args.command == "settings":
            run = run_of(args.slot, "", None, "", parse_settings(args.settings))
            for key in ("model", "effort", "fallback", "id"):
                if run[key]:
                    print(f"{key}={run[key]}")
    except PoolError as error:
        print(f"pool: {error}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
