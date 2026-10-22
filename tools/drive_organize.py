"""Sort a whole Google Drive or OneDrive into projects, each an epic with
its tasks and the directions it could take, proposed for a person to accept.

The drive sync (tools/drive_sync.py) keeps a few picked files in step with
the repository. This reads far more, a whole drive if the person shares
its root, without bringing it in: .github/workflows/drive-organize.yml runs
it in three steps, each started by hand, since only someone who may write
here can start a workflow, and nothing of it runs on a schedule.

  catalog   list every folder and file the items in drive.json's
            "organize" (else its "items") hold, with their size, date and
            link, and digest a sample of them (tools/file_digest.py), to
            the branch drive-catalog: catalog.json, TREE.md, digests.jsonl
  organize  Claude reads that and returns projects as JSON (schema below);
            open makes each a proposed epic with its tasks under it, and
            one index issue, "Your drive, organized"
  accept    a person's yes to one project: it is labelled tracking and its
            folders join drive.json's items, by a pull request, so the sync
            brings them in; ready, which the sync's own run calls, then
            gives its tasks todo and their workers' labels once they are here

  python3 tools/drive_organize.py catalog --config drive.json --out DIR
  python3 tools/drive_organize.py schema --registry .github/symphony/workers
  python3 tools/drive_organize.py open --plan PLAN.json --catalog DIR/catalog.json --config drive.json
  python3 tools/drive_organize.py accept --epic N --config drive.json --write NEW.json
  python3 tools/drive_organize.py ready --config drive.json --root .

drive.json, as the control panel's Folders tab writes it, with "organize"
optional:

  {"provider": "microsoft", "folder": "drive", "items": [...],
   "organize": {"items": [{"id": "root", "name": "OneDrive"}]}}

A proposed project is a plain issue labelled epic, proposed, drive and
project:<slug>, opened by the workflow's own token, so it starts no
workflow; its tasks carry proposed, no todo and no "Epic: #N" line, so no
worker takes one and the coordinator routes none until the project's files
have arrived after accept. The organizer never proposes a slug it proposed before, open
or closed, so closing a proposal is how a person says no. Standard library
only; every name the drive gives is data, cleaned before it is a path, and
a catalog names a person's files, so the workflow runs only in a private
repository and logs counts alone.
"""

import argparse
import glob
import json
import os
import re
import subprocess
import sys
import tempfile
import time
import urllib.parse
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import drive_sync  # noqa: E402
import file_digest  # noqa: E402
from sdlc_protocol import team, workers as registry  # noqa: E402

MAX_NODES = 20000
MAX_DEPTH = 16
DIGEST_FILES = 400
DIGEST_BYTES = 120 << 20
DIGEST_EACH = 5 << 20
TREE_DEPTH = 4
TREE_LINES = 1500
NAMES_A_FOLDER = 10
# The folders Microsoft 365 keeps for its own apps at a OneDrive's top: moving
# one breaks the app that writes there, and what is in them (chat files,
# recordings) is no project's, so they are listed but never read into or moved.
MANAGED = {"apps", "attachments", "meetings", "microsoft copilot chat files", "microsoft teams chat files",
           "microsoft teams data", "recordings", "notebooks", "desktop", "documents", "pictures"}
MAX_PROJECTS = 30
MAX_TASKS = 12
READABLE = {"document", "spreadsheet", "table", "slides", "notebook", "text", "code", "data", "pdf", "image"}
DOMAINS = ["data analysis", "accounting", "finance", "business", "marketing", "sales", "software",
           "research", "writing", "teaching", "design and art", "media", "operations", "personal", "other"]
SLUG = re.compile(r"[a-z0-9][a-z0-9-]{0,39}")
OWNER = re.compile(r"[A-Za-z0-9][A-Za-z0-9_.-]{0,38}")
REPO_NAME = re.compile(r"[A-Za-z0-9][A-Za-z0-9_.-]{0,99}")
TEMPLATE = "OpenSpineConsortium/symphony-template"
BOT = "github-actions[bot]"
LABEL_DRIVE, LABEL_PROPOSED, LABEL_TRACKING = "drive", "proposed", "tracking"
LABEL_EPIC, LABEL_TODO = "epic", "todo"
# A folder larger than this is not tracked whole: it would fill the
# repository and outlast the sync's run.
MAX_TRACK_BYTES = 2 << 30
MAX_TRACK_FILES = 2000
INDEX_MARKER = "<!-- symphony-drive-index -->"
# The organizer's own lines are read only where it writes them, the head
# of a body or its end, and every text Claude gives is cleaned of comment
# marks (plain), so nothing a drive's file said can pose as one.
HEAD = re.compile(r"\A<!-- symphony-drive-project: ([a-z0-9-]{1,40}) -->\n<!-- symphony-drive-items: (\[.*?\]) -->\n")
OWNER_LINE = re.compile(r"\n<!-- symphony-drive-owner: ([A-Za-z0-9_.-]{1,39}) -->\n\Z")
PROJECT_LINE = re.compile(r"\AProject: #(\d+)$", re.MULTILINE)


class OrganizeError(Exception):
    """A step that cannot go on, said in one line with no drive's words in it."""


# The catalog

def node(provider, meta, path, drive=""):
    """One folder or file, the same whatever the provider."""
    if provider == "google":
        mime = meta.get("mimeType", "")
        return {"id": meta["id"], "path": path, "drive": "", "folder": mime == drive_sync.GOOGLE_FOLDER,
                "name": meta.get("name", ""), "mime": mime, "size": as_int(meta.get("size")),
                "modified": meta.get("modifiedTime", ""), "web": meta.get("webViewLink", "")}
    drive = drive or meta.get("parentReference", {}).get("driveId", "")
    return {"id": meta["id"], "path": path, "drive": drive, "folder": "folder" in meta,
            "name": meta.get("name", ""), "mime": meta.get("file", {}).get("mimeType", ""),
            "size": as_int(meta.get("size")), "modified": meta.get("lastModifiedDateTime", ""),
            "web": meta.get("webUrl", "")}


def as_int(value):
    try:
        return max(0, int(value))
    except (TypeError, ValueError):
        return 0


def children(provider_name, provider, folder):
    if provider_name == "google":
        return provider.children(folder["id"])
    return provider.children(folder["drive"], folder["id"])


def top(provider_name, provider, item):
    """The item a walk starts from. A OneDrive's top keeps the id root,
    which Graph takes for it, so the moves can name it as their home's
    parent and tell the folders Microsoft 365's apps keep there."""
    if provider_name == "google":
        meta = provider.meta(item["id"])
    else:
        meta = provider.meta(item.get("drive", ""), item["id"])
    name = item.get("name") or meta.get("name", "") or item["id"]
    found = {**node(provider_name, meta, drive_sync.segment(name), item.get("drive", "")), "parent": ""}
    if provider_name == "microsoft" and (item["id"] == "root" or "root" in meta):
        found["id"] = "root"
    return found


def kept_by_apps(provider_name, entry):
    return (provider_name == "microsoft" and entry["folder"] and entry.get("parent") == "root"
            and entry["name"].strip().lower() in MANAGED)


def walk(provider_name, provider, items):
    """Every folder and file the items hold, each once, depth first and
    by name, with what the drive would not list named, and whether the
    walk stopped at MAX_NODES."""
    found, seen, problems, stopped = [], set(), [], False
    for item in items:
        try:
            stack = [(top(provider_name, provider, item), 0)]
        except drive_sync.DriveError as problem:
            problems.append(f"{item.get('name') or item.get('id', '?')}: {problem}")
            continue
        while stack:
            current, depth = stack.pop()
            if (current["drive"], current["id"]) in seen or skipped(current):
                continue
            if len(found) >= MAX_NODES:
                stopped = True
                break
            seen.add((current["drive"], current["id"]))
            found.append(current)
            if not current["folder"] or kept_by_apps(provider_name, current):
                continue
            if depth >= MAX_DEPTH:
                problems.append(f"{current['path']}: deeper than {MAX_DEPTH} levels, not listed")
                continue
            try:
                kids = [{**node(provider_name, meta, f"{current['path']}/{drive_sync.segment(meta.get('name', ''))}",
                                current["drive"]), "parent": current["id"]}
                        for meta in children(provider_name, provider, current)]
            except drive_sync.DriveError as problem:
                problems.append(f"{current['path']}: {problem}")
                continue
            stack += [(kid, depth + 1) for kid in sorted(kids, key=lambda k: k["path"], reverse=True)]
    return found, problems, stopped


def skipped(entry):
    name = entry["name"]
    return name.startswith("~$") or drive_sync.SIBLING in name or entry["mime"] == "application/vnd.google-apps.shortcut"


def catalog_of(provider_name, found, problems, stopped):
    """The catalog: the folders with what they hold, rolled up to every
    folder above, and the files with their kind."""
    folders = {n["path"]: {"path": n["path"], "id": n["id"], "drive": n["drive"], "web": n["web"],
                           "files": 0, "bytes": 0, "kinds": {}, "latest": ""}
               for n in found if n["folder"]}
    files = []
    for entry in found:
        if entry["folder"]:
            continue
        kind = file_digest.kind_of(entry["name"], entry["mime"])
        files.append({**{k: entry[k] for k in ("path", "id", "drive", "name", "mime", "size", "modified", "web")},
                      "kind": kind, "native": entry["mime"] in drive_sync.GOOGLE_EXPORTS})
        for folder in ancestors(entry["path"], folders):
            folder["files"] += 1
            folder["bytes"] += entry["size"]
            folder["kinds"][kind] = folder["kinds"].get(kind, 0) + 1
            folder["latest"] = max(folder["latest"], entry["modified"])
    totals = {}
    for entry in files:
        totals[entry["kind"]] = totals.get(entry["kind"], 0) + 1
    return {"provider": provider_name, "made": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
            "complete": not stopped and not problems, "problems": problems, "totals": totals,
            "folders": sorted(folders.values(), key=lambda f: f["path"]), "files": files,
            "parents": {n["id"]: n["parent"] for n in found if n.get("parent")},
            "managed": sorted(n["path"] for n in found if kept_by_apps(provider_name, n))}


def ancestors(path, folders):
    parts = path.split("/")
    return [folders[p] for p in ("/".join(parts[:i]) for i in range(1, len(parts))) if p in folders]


def tree(catalog):
    """TREE.md: the folders to TREE_DEPTH with what each holds and the
    names of a few of its files, the first thing Claude reads."""
    direct = {}
    for entry in catalog["files"]:
        direct.setdefault(entry["path"].rpartition("/")[0], []).append(entry)
    lines = [f"# The drive, {sum(catalog['totals'].values())} files",
             "kinds: " + ", ".join(f"{n} {k}" for k, n in sorted(catalog["totals"].items(), key=lambda x: -x[1])),
             ""]
    loose = [e for e in catalog["files"] if "/" not in e["path"]]
    lines += [f"- {e['path']} ({e['kind']}, {file_digest.human(e['size'])})" for e in loose[:NAMES_A_FOLDER]]
    for folder in catalog["folders"]:
        depth = folder["path"].count("/")
        if depth > TREE_DEPTH:
            continue
        if len(lines) >= TREE_LINES:
            lines.append("... the tree goes on; catalog.json lists every folder and file")
            break
        kinds = ", ".join(f"{n} {k}" for k, n in sorted(folder["kinds"].items(), key=lambda x: -x[1])[:4])
        indent = "  " * depth
        lines.append(f"{indent}- {folder['path'].rpartition('/')[2]}/  {folder['files']} files, "
                     f"{file_digest.human(folder['bytes'])}{'; ' + kinds if kinds else ''}"
                     f"{'; latest ' + folder['latest'][:10] if folder['latest'] else ''}")
        here = sorted(direct.get(folder["path"], []), key=lambda e: e["modified"], reverse=True)
        lines += [f"{indent}  - {e['name']} ({e['kind']})" for e in here[:NAMES_A_FOLDER]]
        if len(here) > NAMES_A_FOLDER:
            lines.append(f"{indent}  - and {len(here) - NAMES_A_FOLDER} more")
    return "\n".join(lines) + "\n"


def choose(files):
    """Which files to digest: the newest few of every folder in turn, so a
    folder of ten thousand photographs takes no more of the budget than a
    folder of three spreadsheets, until DIGEST_FILES or DIGEST_BYTES."""
    groups = {}
    for entry in files:
        if entry["kind"] in READABLE and (entry["native"] or 0 < entry["size"] <= DIGEST_EACH):
            groups.setdefault(entry["path"].rpartition("/")[0], []).append(entry)
    queues = [sorted(group, key=lambda e: e["modified"], reverse=True) for _, group in sorted(groups.items())]
    chosen, spent = [], 0
    while queues and len(chosen) < DIGEST_FILES:
        for queue in list(queues):
            if not queue:
                queues.remove(queue)
                continue
            entry = queue.pop(0)
            if spent + entry["size"] > DIGEST_BYTES:
                continue
            chosen.append(entry)
            spent += entry["size"]
            if len(chosen) >= DIGEST_FILES:
                break
    return chosen


def download_entry(entry):
    """What drive_sync's providers take to download a file."""
    return {"id": entry["id"], "drive": entry["drive"], "native": entry["native"], "mime": entry["mime"]}


def digests(provider, chosen):
    """Each chosen file's digest; the file itself is gone once read."""
    with tempfile.TemporaryDirectory() as scratch:
        for entry in chosen:
            name = entry["name"]
            if entry["native"]:
                name += drive_sync.GOOGLE_EXPORTS[entry["mime"]][1]
            target = Path(scratch) / ("file" + Path(drive_sync.segment(name)).suffix.lower()[:12])
            try:
                data = provider.download(download_entry(entry))
                if len(data) > DIGEST_EACH:
                    raise drive_sync.DriveError("larger than the digest takes")
                target.write_bytes(data)
                found = file_digest.digest(target, entry["kind"], name)
            except drive_sync.DriveError:
                found = {"kind": entry["kind"], "size": entry["size"], "text": "could not be downloaded"}
            finally:
                target.unlink(missing_ok=True)
            yield {"path": entry["path"], "kind": entry["kind"], "text": found["text"]}


def write_catalog(provider_name, provider, items, out):
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    found, problems, stopped = walk(provider_name, provider, items)
    if not found:
        raise OrganizeError("the drive listed none of the items: " + "; ".join(problems[:3]))
    catalog = catalog_of(provider_name, found, problems, stopped)
    (out / "catalog.json").write_text(json.dumps(catalog, indent=1) + "\n", encoding="utf-8")
    (out / "TREE.md").write_text(tree(catalog), encoding="utf-8")
    count = 0
    with open(out / "digests.jsonl", "w", encoding="utf-8") as stream:
        for line in digests(provider, choose(catalog["files"])):
            stream.write(json.dumps(line) + "\n")
            count += 1
    return {"folders": len(catalog["folders"]), "files": len(catalog["files"]), "digested": count,
            "unlisted": len(problems), "stopped": stopped}


# What Claude returns

def schema(owners):
    """The JSON schema of the organizer's answer, as Claude Code's
    --json-schema takes it: each project's tasks are the coordinator's
    tasks (team.plan_schema), so the clerk's own code opens them."""
    task = team.plan_schema(owners)["properties"]["tasks"]["items"]
    text = {"type": "string"}
    texts = {"type": "array", "items": text}
    project = {"type": "object",
               "properties": {"slug": {"type": "string", "pattern": f"^{SLUG.pattern}$"},
                              "title": text, "domain": {"type": "string", "enum": DOMAINS},
                              "summary": text, "folders": texts, "files": texts,
                              "directions": {**texts, "minItems": 1, "maxItems": 8},
                              "plan": text,
                              "repository": {"type": "object",
                                             "properties": {"name": {"type": "string",
                                                                     "pattern": f"^{REPO_NAME.pattern}$"},
                                                            "team": {"type": "boolean"}, "why": text},
                                             "required": ["name", "team", "why"]},
                              "tasks": {"type": "array", "items": task, "minItems": 1, "maxItems": MAX_TASKS}},
               "required": ["slug", "title", "domain", "summary", "folders", "files", "directions", "plan",
                            "tasks"]}
    return {"type": "object",
            "properties": {"summary": text,
                           "projects": {"type": "array", "items": project, "minItems": 1,
                                        "maxItems": MAX_PROJECTS},
                           "unsorted": texts, "directions": texts},
            "required": ["summary", "projects", "unsorted", "directions"]}


# GitHub, through gh

def gh(args, payload=None):
    """One gh call; its standard error is left out, since it can quote an
    issue's text, and the log of the run is kept."""
    result = subprocess.run(["gh", *args], input=None if payload is None else json.dumps(payload),
                            capture_output=True, text=True)
    if result.returncode:
        raise OrganizeError(f"gh {' '.join(a for a in args[:2] if not a.startswith('-'))} failed")
    return result.stdout


class GitHub:
    """The few calls the organizer makes, on the repository GITHUB_REPOSITORY names."""

    def __init__(self, repo, run=gh):
        self.repo, self.run = repo, run

    def api(self, method, path, payload=None):
        args = ["api", "-X", method, f"repos/{self.repo}/{path}"] + (["--input", "-"] if payload is not None else [])
        out = self.run(args, payload)
        return json.loads(out) if out.strip() else {}

    def issues(self):
        out = self.run(["api", "--paginate", f"repos/{self.repo}/issues?state=all&labels={LABEL_DRIVE}&per_page=100",
                        "--jq", ".[] | {number, id, state, body, user: .user.login, "
                                "labels: [.labels[].name], pr: (.pull_request != null)}"])
        return [json.loads(line) for line in out.splitlines() if line.strip()]

    def issue(self, number):
        found = self.api("GET", f"issues/{number}")
        return {"number": found.get("number"), "id": found.get("id"), "state": found.get("state"),
                "body": found.get("body") or "", "user": (found.get("user") or {}).get("login", ""),
                "labels": [label.get("name") for label in found.get("labels", [])]}

    def create(self, title, body, labels):
        return self.api("POST", "issues", {"title": title[:250], "body": body, "labels": labels})

    def set_body(self, number, body):
        self.api("PATCH", f"issues/{number}", {"body": body})

    def add_labels(self, number, labels):
        self.api("POST", f"issues/{number}/labels", {"labels": labels})

    def remove_label(self, number, name):
        try:
            self.api("DELETE", f"issues/{number}/labels/{urllib.parse.quote(name, safe='')}")
        except OrganizeError:
            pass  # a label already gone counts as taken off

    def label(self, name, color, description):
        self.run(["label", "create", name, "--force", "--color", color, "--description", description,
                  "--repo", self.repo])

    def sub_issue(self, parent, child_id):
        try:
            self.api("POST", f"issues/{parent}/sub_issues", {"sub_issue_id": child_id})
        except OrganizeError:
            pass  # the link is for people; the epic's task list is what the clerk reads


# Opening the projects

def ours(issue):
    return issue.get("user") == BOT and not issue.get("pr")


def head_of(issue):
    """The project's slug and items an issue of the organizer's own begins
    with, or None."""
    mark = HEAD.match(issue.get("body") or "") if ours(issue) else None
    if not mark:
        return None
    try:
        items = json.loads(mark.group(2))
    except ValueError:
        items = []
    return mark.group(1), items if isinstance(items, list) else []


def known_projects(issues):
    """The slugs the organizer has proposed before, open or closed."""
    found = {}
    for issue in issues:
        head = head_of(issue)
        if head:
            found.setdefault(head[0], issue)
    return found


def plain(text, limit=4000):
    """Text Claude wrote, with no comment marks left in it to pose as the
    organizer's own lines."""
    return str(text).replace("<!--", "&lt;!--").replace("-->", "--&gt;")[:limit]


def hidden(value):
    """JSON for an HTML comment: no > in it, so nothing ends the comment."""
    return json.dumps(value, ensure_ascii=True, separators=(",", ":")).replace(">", "\\u003e")


def tick(text):
    return "`" + str(text).replace("`", "'").replace("\n", " ")[:200] + "`"


def link(web):
    return f" ([open]({web}))" if isinstance(web, str) and web.startswith("https://") and ")" not in web else ""


def items_for(project, lookup):
    """The drive's folders and files a project names, as drive.json items,
    each named by the catalog, never by Claude, and the paths left out: one
    the catalog lacks, a top item, which is the whole of what was shared,
    and a folder too large to track."""
    items, left_out, seen = [], [], set()
    for path in list(project.get("folders", [])) + list(project.get("files", [])):
        entry = lookup.get(str(path).strip().strip("/"))
        if not entry:
            continue
        if "/" not in entry["path"] or entry.get("bytes", entry.get("size", 0)) > MAX_TRACK_BYTES \
                or entry.get("files", 0) > MAX_TRACK_FILES:
            left_out.append(entry["path"])
        elif (entry["drive"], entry["id"]) not in seen:
            seen.add((entry["drive"], entry["id"]))
            items.append({"id": entry["id"], "name": entry["path"].rpartition("/")[2], "drive": entry["drive"],
                          "path": entry["path"], "web": entry.get("web", "")})
    return items, left_out


def project_body(project, items, left_out, folder):
    places = "\n".join(f"- {tick(item['path'])}{link(item['web'])}" for item in items) or "- (none the catalog knows)"
    if left_out:
        places += ("\n\nNot tracked whole, being all that was shared or too large for the sync: "
                   + ", ".join(tick(path) for path in left_out[:20]))
    directions = "\n".join(f"- {plain(line, 600)}" for line in project["directions"][:8])
    saved = [{k: item[k] for k in ("id", "name", "drive")} for item in items]
    return (f"<!-- symphony-drive-project: {project['slug']} -->\n<!-- symphony-drive-items: {hidden(saved)} -->\n"
            f"Domain: {plain(project['domain'], 40)}\n\n{plain(project['summary']).strip()}\n\n"
            f"## Where it lives in the drive\n\n{places}\n\n"
            f"## Directions you could take\n\n{directions}\n\n"
            f"{repository_section(project.get('repository'))}"
            "## What happens next\n\n"
            "Proposed by the drive organizer: no worker takes its tasks yet. To go ahead, run the Drive "
            "organizer workflow with the step accept and this issue's number. The folders above then join "
            f"the drive sync, which brings them into `{folder}/`, and once they are there the tasks get todo "
            "and go to the workers. To say no, close this issue; the organizer never proposes it again.")


def repository_section(repository):
    """A Symphony repository of its own, where the project would sooner live
    in one that people join than in this one: the link makes it from the
    template, and a team's is opened to volunteers by crowd mode, each on
    their own Claude account. The owner makes it; nothing here does."""
    if not isinstance(repository, dict) or not REPO_NAME.fullmatch(str(repository.get("name", ""))):
        return ""
    owner, _, name = TEMPLATE.partition("/")
    url = f"https://github.com/new?template_owner={owner}&template_name={name}&name={repository['name']}"
    team = (" For people who join with their own Claude accounts, its owner then runs "
            "`python3 scripts/crowd.py setup` in it (docs/setup/crowdsourcing.md).") if repository.get("team") else ""
    return (f"## A repository of its own\n\n{plain(repository.get('why', ''), 600).strip()}\n\n"
            f"Make it as {tick(repository['name'])} with [this link]({url}), choosing its owner there, then move "
            f"this issue's tasks to it.{team}\n\n")


def clean_task(task):
    return {"title": plain(task.get("title", ""), 250), "owner": task.get("owner", ""),
            "files": [plain(path, 300) for path in task.get("files", [])[:50]], "body": plain(task.get("body", ""))}


def task_issue_body(epic, task, owner):
    """A task's body as the clerk writes it, but under "Project: #N": the
    "Epic: #N" line, which the coordinator's routing reads, is written only
    once the task starts (start_tasks), so closing a proposal's task starts
    no integration of a project nobody accepted."""
    body = team.task_body(epic, task).replace(f"Epic: #{epic}", f"Project: #{epic}", 1)
    note = ("\nThe files above are paths in the drive. The drive sync brings their folders into the repository "
            "once the project is accepted. `python3 tools/file_digest.py FILE` reads a spreadsheet, a deck or "
            "a notebook as text.\n")
    return body + note + (f"\n<!-- symphony-drive-owner: {owner} -->\n" if owner else "")


def ensure_labels(github, slugs):
    github.label(LABEL_DRIVE, "1d76db", "Opened by the drive organizer")
    github.label(LABEL_PROPOSED, "fbca04", "Proposed by the drive organizer; accept to start it")
    github.label(LABEL_TRACKING, "bfd4f2", "Accepted; its tasks start once the drive sync brings its files in")
    github.label(LABEL_EPIC, "3e4b9e", "A request the coordinator plans into tasks")
    for slug in slugs:
        github.label(f"project:{slug}", "c5def5", "A project the drive organizer found")


def open_project(github, project, catalog_lookup, folder, owners):
    """The epic first, with no head, so a run that fails halfway leaves
    nothing the organizer reads as a project, then its tasks, then the
    head, the plan and the label epic."""
    slug = project["slug"]
    items, left_out = items_for(project, catalog_lookup)
    epic = github.create(plain(project["title"], 250), "Being opened by the drive organizer.",
                         [LABEL_DRIVE, LABEL_PROPOSED])
    tasks = [clean_task(task) for task in project["tasks"][:MAX_TASKS]]
    made = []
    for task, owner in zip(tasks, team.owners_for(tasks, list(owners))):
        owner = owner if owner and OWNER.fullmatch(owner) else None
        issue = github.create(task["title"], task_issue_body(epic["number"], task, owner),
                              [LABEL_DRIVE, LABEL_PROPOSED, f"project:{slug}"])
        github.sub_issue(epic["number"], issue["id"])
        made.append((issue["number"], task["title"], owner))
    body = project_body(project, items, left_out, folder)
    github.set_body(epic["number"], body + team.plan_section(plain(project.get("plan", "")), made))
    github.add_labels(epic["number"], [LABEL_EPIC, f"project:{slug}"])
    return epic["number"]


def index_body(plan, numbers, folder):
    rows = "\n".join(f"- #{number} {plain(title, 250)} ({plain(domain, 40)}){' (proposed before)' if old else ''}"
                     for number, title, domain, old in numbers)
    unsorted = "\n".join(f"- {tick(plain(path, 300))}" for path in plan.get("unsorted", [])[:100]) or "- nothing"
    directions = "\n".join(f"- {plain(line, 600)}" for line in plan.get("directions", [])[:20]) or "- none"
    return (f"{INDEX_MARKER}\n{plain(plan.get('summary', '')).strip()}\n\n## Projects\n\n{rows}\n\n"
            "Each is an epic with its tasks under it, proposed: accept one with the Drive organizer "
            f"workflow's step accept and its number, and its folders join the sync into `{folder}/`.\n\n"
            f"## What fits no project\n\n{unsorted}\n\n## Directions across the drive\n\n{directions}\n")


def whole(project):
    """Whether a project has every field the schema requires, of its type,
    before anything of it is opened: a plan a session wrote is held to the
    schema here, as the workflow's Claude run is held by --json-schema."""
    texts = ("title", "summary", "plan")
    lists = ("folders", "files", "directions", "tasks")
    return (all(isinstance(project.get(key), str) for key in texts)
            and all(isinstance(project.get(key), list) for key in lists)
            and project.get("domain") in DOMAINS and bool(project["tasks"]) and bool(project["directions"])
            and all(isinstance(task, dict) for task in project["tasks"]))


def open_plan(github, plan, catalog, folder, owners=()):
    """Every new project as a proposed epic, then the index issue, opened
    or rewritten; returns the epics' numbers, new ones first."""
    lookup = {entry["path"]: entry for entry in catalog["folders"] + catalog["files"]}
    issues = github.issues()
    known = known_projects(issues)
    projects, slugs = [], set()
    for project in plan.get("projects", [])[:MAX_PROJECTS]:
        slug = str(project.get("slug", ""))
        if SLUG.fullmatch(slug) and slug not in slugs and whole(project):
            slugs.add(slug)
            projects.append(project)
    ensure_labels(github, sorted(slugs - set(known)))
    numbers, opened = [], []
    for project in projects:
        old = known.get(project["slug"])
        number = old["number"] if old else open_project(github, project, lookup, folder, owners)
        if not old:
            opened.append(number)
        numbers.append((number, project["title"], project["domain"], bool(old)))
    body = index_body(plan, numbers, folder)
    index = next((i for i in issues if ours(i) and (i.get("body") or "").startswith(INDEX_MARKER)), None)
    if index:
        github.set_body(index["number"], body)
    else:
        index = github.create("Your drive, organized", body, [LABEL_DRIVE])
    return {"index": index["number"], "opened": opened, "kept": [n for n, *_, old in numbers if old]}


# Accepting a project, and starting its tasks

def accept(github, number, config):
    """A person's yes to a proposed project: proposed comes off, tracking
    goes on, and the drive.json it returns tracks the project's folders.
    Its tasks wait for start_tasks, once their files are here. Running it
    again on a project it accepted is safe, so a failed pull request can
    be made again; one the organizer did not open is refused."""
    epic = github.issue(number)
    head = head_of(epic)
    if not head or LABEL_EPIC not in epic["labels"]:
        raise OrganizeError(f"#{number} is no project the drive organizer proposed")
    if epic["state"] != "open":
        raise OrganizeError(f"#{number} is closed; reopen it to accept it")
    if LABEL_PROPOSED in epic["labels"] or LABEL_TRACKING not in epic["labels"]:
        github.label(LABEL_TRACKING, "bfd4f2", "Accepted; its tasks start once the drive sync brings its files in")
        github.add_labels(number, [LABEL_TRACKING])
        github.remove_label(number, LABEL_PROPOSED)
    return tracked(config, head[1])


def start_tasks(github, epic):
    """Give an accepted project's open tasks to their workers: the
    "Epic: #N" line the coordinator reads, todo and the worker's label."""
    found = team.TASKS_LINE.findall(epic["body"])
    started = []
    for task_number in (int(n) for n in (found[-1].split() if found else [])):
        task = github.issue(task_number)
        if not ours(task) or task["state"] != "open" or LABEL_PROPOSED not in task["labels"]:
            continue
        body = PROJECT_LINE.sub(f"Epic: #{epic['number']}", task["body"], count=1)
        github.set_body(task_number, body)
        owner = OWNER_LINE.search(task["body"])
        labels = [LABEL_TODO]
        if owner:
            name = team.worker_label(owner.group(1))
            github.label(name, team.WORKER_COLOUR, "A worker's own label")
            labels.append(name)
        github.add_labels(task_number, labels)
        github.remove_label(task_number, LABEL_PROPOSED)
        started.append(task_number)
    github.remove_label(epic["number"], LABEL_TRACKING)
    return started


def arrived(items, config, base):
    """Whether every item is in drive.json and the sync has brought it in,
    a Google file under its exported name."""
    have = {(item.get("drive", ""), item.get("id")) for item in config.get("items", [])}
    for item in items:
        key = (str(item.get("drive", "")), item.get("id"))
        name = drive_sync.segment(str(item.get("name", "")))
        if key not in have or not ((base / name).exists() or any(base.glob(glob.escape(name) + ".*"))):
            return False
    return True


def ready(github, config, root):
    """Start every accepted project whose files the sync has brought in;
    the drive sync's run calls it, and so does accept."""
    base = Path(root) / config.get("folder", "papers")
    started = {}
    for issue in github.issues():
        head = head_of(issue)
        if head and issue["state"] == "open" and LABEL_TRACKING in issue["labels"] \
                and LABEL_EPIC in issue["labels"] and arrived(head[1], config, base):
            started[issue["number"]] = start_tasks(github, issue)
    return started


def tracked(config, items):
    """drive.json with the project's items added, and how many were new."""
    config = {**config, "items": list(config.get("items", []))}
    have = {(item.get("drive", ""), item.get("id")) for item in config["items"]}
    added = 0
    for item in items:
        if not (isinstance(item, dict) and isinstance(item.get("id"), str) and item["id"]):
            continue
        key = (str(item.get("drive", "")), item["id"])
        if key not in have:
            have.add(key)
            config["items"].append({"id": item["id"], "name": str(item.get("name", ""))[:120],
                                    **({"drive": key[0]} if key[0] else {})})
            added += 1
    return config, added


# The command line

def organize_items(config):
    items = (config.get("organize") or {}).get("items") or config.get("items") or []
    if not items:
        raise OrganizeError("drive.json names no items to organize")
    return items


def provider_for(config):
    kind = drive_sync.PROVIDERS.get(config.get("provider"))
    if not kind:
        raise OrganizeError("drive.json names no provider this knows")
    token = drive_sync.access_token(os.environ["SYMPHONY_DRIVE_TOKEN_URL"], config["provider"],
                                    os.environ["SYMPHONY_DRIVE_REFRESH"])
    return kind(token)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = parser.add_subparsers(dest="command", required=True)
    making = sub.add_parser("catalog")
    making.add_argument("--config", required=True)
    making.add_argument("--out", required=True)
    shape = sub.add_parser("schema")
    shape.add_argument("--registry", default=".github/symphony/workers")
    opening = sub.add_parser("open")
    opening.add_argument("--plan", required=True)
    opening.add_argument("--catalog", required=True)
    opening.add_argument("--config", required=True)
    opening.add_argument("--registry", default=".github/symphony/workers")
    taking = sub.add_parser("accept")
    taking.add_argument("--epic", required=True, type=int)
    taking.add_argument("--config", required=True)
    taking.add_argument("--write", required=True)
    starting = sub.add_parser("ready")
    starting.add_argument("--config", required=True)
    starting.add_argument("--root", default=".")
    args = parser.parse_args(argv)
    try:
        return run_command(args)
    except (OrganizeError, drive_sync.DriveError, registry.WorkerError, KeyError, ValueError) as problem:
        print(f"drive_organize: {problem}", file=sys.stderr)
        return 1


def registered(directory):
    """The registry's workers, any of whom a task may go to."""
    return [worker.name for worker in registry.load(directory)]


def run_command(args):
    if args.command == "schema":
        print(json.dumps(schema(registered(args.registry)), separators=(",", ":")))
        return 0
    config = json.loads(Path(args.config).read_text(encoding="utf-8"))
    if args.command == "catalog":
        counts = write_catalog(config["provider"], provider_for(config), organize_items(config), args.out)
        print(json.dumps(counts))
        return 0
    github = GitHub(os.environ["GITHUB_REPOSITORY"])
    folder = config.get("folder", "papers")
    if args.command == "open":
        plan = json.loads(Path(args.plan).read_text(encoding="utf-8"))
        catalog = json.loads(Path(args.catalog).read_text(encoding="utf-8"))
        result = open_plan(github, plan, catalog, folder, registered(args.registry))
        print(json.dumps({"index": result["index"], "opened": len(result["opened"]), "kept": len(result["kept"])}))
        return 0
    if args.command == "ready":
        started = ready(github, config, args.root)
        print(json.dumps({"projects": len(started), "tasks": sum(len(t) for t in started.values())}))
        return 0
    new_config, added = accept(github, args.epic, config)
    Path(args.write).write_text(json.dumps(new_config, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"added": added}))
    return 0


if __name__ == "__main__":
    sys.exit(main())
