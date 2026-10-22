"""Keep files a person picked in Google Drive or OneDrive in step with a
folder of this repository, both ways, from GitHub Actions.

The person connects a drive once, on the control panel's Folders page
(docs/research/drive-sync.md): one press signs in to Google or Microsoft,
a picker chooses the files and folders to track, and the page writes
.github/symphony/drive.json, the list below, and saves the drive's refresh
token as the repository's secret SYMPHONY_DRIVE_REFRESH. No client secret
is kept here: the run trades the refresh token for an access token at the
sign-in service, which holds the client secrets and keeps nothing, the way
the panel's own sign-in does.

  {"provider": "google" | "microsoft", "folder": "papers",
   "items": [{"id": "...", "name": "...", "drive": "<OneDrive drive id>"}]}

Each run of .github/workflows/drive-sync.yml lists what the items hold now
(a folder's files, its subfolders' too), compares it with the run before,
kept as state.json on the branch drive-state, and with the files on the
default branch, and for each file:

  changed in the drive only       pull: it rides the pull request drive/sync
  changed on the default branch   push: the new copy goes up in place (a
                                  native Google Doc, Sheet or Slides, which
                                  takes no upload, gets a copy beside it)
  changed on both, or seen first  the default branch's copy goes up beside
                                  the drive's as "<name> (Symphony
                                  revisions)", and the drive's comes down

so nothing either side wrote is lost. A native Google file comes down
exported (a Doc as .docx with its comments, a Sheet as .xlsx, Slides as
.pptx). A file the drive no longer has stays in the repository. A pulled
file is pending until drive/sync merges, so a later run neither pulls it
twice nor pushes the older copy up, and a copy it once staged that merges
after the drive moved on is pulled over, never pushed. Files are known by
the drive's id, so a rename or two files of one name mix nothing up.

  python3 tools/drive_sync.py folder --config .github/symphony/drive.json \\
      --state STATE.json --root .
  python3 tools/drive_sync.py run --config .github/symphony/drive.json \\
      --state STATE.json --root . --out STAGING --summary SUMMARY.json

The first names the repository's folder the run may write, or refuses
one that is no plain path inside the repository or that already holds
the repository's own files, such as scripts or tools, since whoever may
edit the drive's folder would then rewrite them and the sync's pull
request merges by itself.

Standard library only; every value read from the drive is data, never a
path outside the folder: a name is cleaned to one path segment.
"""

import argparse
import hashlib
import json
import os
import re
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

GOOGLE = "https://www.googleapis.com"
GRAPH = "https://graph.microsoft.com/v1.0"
DOCX = "application/vnd.openxmlformats-officedocument.wordprocessingml.document"
GOOGLE_EXPORTS = {
    "application/vnd.google-apps.document": (DOCX, ".docx"),
    "application/vnd.google-apps.spreadsheet":
        ("application/vnd.openxmlformats-officedocument.spreadsheetml.sheet", ".xlsx"),
    "application/vnd.google-apps.presentation":
        ("application/vnd.openxmlformats-officedocument.presentationml.presentation", ".pptx"),
}
GOOGLE_FOLDER = "application/vnd.google-apps.folder"
SIBLING = " (Symphony revisions)"
MAX_DEPTH = 6
# Graph's simple upload takes 4 MB; a larger file needs an upload session.
SIMPLE_UPLOAD = 4 * 1024 * 1024
FOLDER = re.compile(r"[A-Za-z0-9_][A-Za-z0-9_. -]*(/[A-Za-z0-9_][A-Za-z0-9_. -]*)*")
REASON = re.compile(r"[A-Za-z0-9_.-]{1,64}")
# Graph throttles a long walk with 429 (and a read with 503), as Google does
# with 429 (its 403 rateLimitExceeded is not asked again),
# naming the wait in Retry-After; a call is asked again that many times, so a
# throttled folder is listed late instead of left out of the catalog.
THROTTLE_TRIES = 4
THROTTLE_WAIT_MOST = 120
# Cloudflare, in front of the sign-in service, refuses urllib's own
# User-Agent (Python-urllib/3.x) with a plain-text 403, "error code: 1010",
# before the service sees the call, so every call names this program.
USER_AGENT = "symphony-drive/1 (+https://openspineconsortium.com/symphony/)"
EDGE_CODE = re.compile(rb"error code: ?([0-9]{3,5})")


class DriveError(Exception):
    """A call the drive refused or a file this run will not move, for the summary."""


# HTTP, with no redirect followed while a token rides the request

class NoRedirect(urllib.request.HTTPRedirectHandler):
    def redirect_request(self, *args, **kwargs):
        return None


OPENER = urllib.request.build_opener(NoRedirect)


def call(method, url, token=None, data=None, headers=None, opener=None, hops=3, sleep=time.sleep):
    """One request; a redirect is followed without the token, since urllib
    would send the Authorization header to whatever host it names, and only
    for a GET, so a write is never turned into a read that looks like it
    worked. A throttled answer is asked again after the wait it names (see
    throttled). Every failure the network or the drive gives is a DriveError."""
    for attempt in range(THROTTLE_TRIES):
        request = urllib.request.Request(url, data=data, method=method, headers={"User-Agent": USER_AGENT, **(headers or {})})
        if token:
            request.add_header("Authorization", f"Bearer {token}")
        try:
            with (opener or OPENER).open(request, timeout=60) as response:
                return response.status, dict(response.headers), response.read()
        except urllib.error.HTTPError as error:
            location = error.headers.get("Location") if error.headers else None
            if error.code in (301, 302, 303, 307, 308) and location and method == "GET" and hops > 0:
                return call("GET", urllib.parse.urljoin(url, location), None, None, None, opener, hops - 1, sleep)
            wait = throttled(method, error)
            if wait is not None and attempt + 1 < THROTTLE_TRIES:
                sleep(wait)
                continue
            raise DriveError(f"{method} {url.split('?')[0]} answered {error.code}: {reason_in(error.read())}") from None
        except (urllib.error.URLError, OSError, ValueError) as error:
            raise DriveError(f"{method} {url.split('?')[0]} failed: {error}") from None


def throttled(method, error):
    """The seconds to wait before asking again, or None when the answer is
    not one to ask again. A 429 did nothing, so any call is asked again; a
    503 is asked again only for a read, since a write it answered may have
    been done. The wait is Retry-After's, bounded, else a few seconds."""
    if error.code != 429 and not (error.code == 503 and method == "GET"):
        return None
    named = str(error.headers.get("Retry-After", "") if error.headers else "").strip()
    seconds = int(named) if re.fullmatch(r"[0-9]{1,6}", named) else 5
    return min(max(seconds, 1), THROTTLE_WAIT_MOST)


def reason_in(body):
    """The short error code of an answer, and nothing else of it: the log
    of a public repository is public, and a drive's message can name a
    file, a person or an address."""
    try:
        answer = json.loads(body) if isinstance(body, (bytes, str)) else body
    except ValueError:
        answer = None
    if answer is None and isinstance(body, bytes):
        edge = EDGE_CODE.search(body[:512])
        if edge:
            return f"cloudflare_{edge.group(1).decode('ascii')}"
    found = answer.get("error") if isinstance(answer, dict) else None
    if isinstance(found, dict):
        found = found.get("status") or found.get("code")
    return found if isinstance(found, str) and REASON.fullmatch(found) else "no reason given"


def call_json(method, url, token=None, payload=None, opener=None):
    data = json.dumps(payload).encode("utf-8") if payload is not None else None
    headers = {"Content-Type": "application/json"} if payload is not None else {}
    _, _, body = call(method, url, token, data, headers, opener)
    try:
        return json.loads(body or b"{}")
    except ValueError:
        raise DriveError(f"{method} {url.split('?')[0]} answered with no JSON") from None


def access_token(token_url, provider, refresh_token, opener=None):
    """Trade the refresh token for an access token at the sign-in service."""
    answer = call_json("POST", token_url, None, {"provider": provider, "refresh_token": refresh_token}, opener)
    if not answer.get("access_token"):
        raise DriveError(f"the sign-in service gave no access token: {reason_in(answer)}")
    return answer["access_token"]


# Names

def segment(name):
    """A drive's name as one safe path segment."""
    cleaned = re.sub(r'[\x00-\x1f<>:"/\\|?*]', "_", name).strip().strip(".")
    return cleaned[:120] or "_"


def blob(data):
    return hashlib.sha256(data).hexdigest()


# Google Drive

class Google:
    def __init__(self, token, opener=None):
        self.token, self.opener = token, opener

    def get(self, url):
        return call_json("GET", url, self.token, opener=self.opener)

    def meta(self, item_id):
        fields = "id,name,mimeType,version,modifiedTime,parents,size,webViewLink"
        return self.get(f"{GOOGLE}/drive/v3/files/{q(item_id)}?fields={fields}&supportsAllDrives=true")

    def children(self, folder_id):
        query = q(f"'{folder_id}' in parents and trashed=false")
        return self.pages(f"{GOOGLE}/drive/v3/files?q={query}&pageSize=200&supportsAllDrives=true"
                          "&includeItemsFromAllDrives=true&corpora=allDrives&orderBy=name"
                          "&fields=nextPageToken,files(id,name,mimeType,version,modifiedTime,parents,size,"
                          "webViewLink)", "files")

    def pages(self, url, key):
        token, found = "", []
        while True:
            page = self.get(url + (f"&pageToken={q(token)}" if token else ""))
            found += page.get(key, [])
            token = page.get("nextPageToken")
            if not token:
                return found

    def comments_mark(self, file_id):
        """A short hash of the file's Drive comments, since a comment on a
        Google Doc changes no version of the file itself."""
        comments = self.pages(f"{GOOGLE}/drive/v3/files/{q(file_id)}/comments?pageSize=100"
                              "&fields=nextPageToken,comments(id,modifiedTime,resolved,replies(id,modifiedTime))",
                              "comments")
        return blob(json.dumps(comments, sort_keys=True).encode("utf-8"))[:12]

    def entry(self, meta, rel):
        native = meta["mimeType"] in GOOGLE_EXPORTS
        version = f"{meta.get('version', '')}+{self.comments_mark(meta['id'])}"
        return {"id": meta["id"], "rel": rel, "version": version, "native": native,
                "mime": meta["mimeType"], "parent": (meta.get("parents") or [""])[0]}

    def walk(self, meta, prefix, depth=0):
        if meta["mimeType"] == GOOGLE_FOLDER:
            if depth >= MAX_DEPTH:
                return
            base = f"{prefix}{segment(meta['name'])}/"
            for child in self.children(meta["id"]):
                yield from self.walk(child, base, depth + 1)
        elif meta["mimeType"] in GOOGLE_EXPORTS or not meta["mimeType"].startswith("application/vnd.google-apps."):
            name = segment(meta["name"])
            suffix = GOOGLE_EXPORTS.get(meta["mimeType"], (None, ""))[1]
            if suffix and not name.lower().endswith(suffix):
                name += suffix
            if SIBLING not in name and not name.startswith("~$"):
                yield self.entry(meta, prefix + name)

    def listing(self, items):
        for item in items:
            yield from self.walk(self.meta(item["id"]), "")

    def download(self, entry):
        if entry["native"]:
            mime = GOOGLE_EXPORTS[entry["mime"]][0]
            url = f"{GOOGLE}/drive/v3/files/{q(entry['id'])}/export?mimeType={q(mime)}"
        else:
            url = f"{GOOGLE}/drive/v3/files/{q(entry['id'])}?alt=media&supportsAllDrives=true"
        return call("GET", url, self.token, opener=self.opener)[2]

    def upload(self, entry, data):
        url = f"{GOOGLE}/upload/drive/v3/files/{q(entry['id'])}?uploadType=media&supportsAllDrives=true"
        call("PATCH", url, self.token, data, {"Content-Type": entry["mime"]}, self.opener)
        return self.entry(self.meta(entry["id"]), entry["rel"])["version"]

    def upload_beside(self, entry, name, data):
        """Write the copy beside the drive's, replacing an earlier one of the
        same name, so revisions do not pile up."""
        mime = mime_of(name)
        parent = entry["parent"]
        query = q(f"name = '{name.replace(chr(92), chr(92) * 2).replace(chr(39), chr(92) + chr(39))}'"
                  f" and '{parent}' in parents and trashed=false")
        found = self.get(f"{GOOGLE}/drive/v3/files?q={query}&supportsAllDrives=true&includeItemsFromAllDrives=true"
                         "&corpora=allDrives&fields=files(id)").get("files", []) if parent else []
        if found:
            target = found[0]["id"]
        else:
            target = call_json("POST", f"{GOOGLE}/drive/v3/files?supportsAllDrives=true", self.token,
                               {"name": name, "mimeType": mime, "parents": [parent] if parent else []},
                               self.opener)["id"]
        url = f"{GOOGLE}/upload/drive/v3/files/{q(target)}?uploadType=media&supportsAllDrives=true"
        call("PATCH", url, self.token, data, {"Content-Type": mime}, self.opener)


def mime_of(name):
    suffix = Path(name).suffix.lower()
    for mime, ext in GOOGLE_EXPORTS.values():
        if ext == suffix:
            return mime
    return "application/octet-stream"


# OneDrive, through Microsoft Graph

class Microsoft:
    def __init__(self, token, opener=None):
        self.token, self.opener = token, opener

    def item_url(self, drive, item_id):
        return f"{GRAPH}/drives/{q(drive)}/items/{q(item_id)}" if drive else f"{GRAPH}/me/drive/items/{q(item_id)}"

    def meta(self, drive, item_id):
        return call_json("GET", self.item_url(drive, item_id) + "?$select=id,name,cTag,eTag,file,folder,root,"
                         "parentReference,size,lastModifiedDateTime,webUrl", self.token, opener=self.opener)

    def children(self, drive, item_id):
        url = (self.item_url(drive, item_id) + "/children?$top=200&$select=id,name,cTag,eTag,file,folder,"
               "parentReference,size,lastModifiedDateTime,webUrl")
        while url:
            page = call_json("GET", url, self.token, opener=self.opener)
            yield from page.get("value", [])
            url = page.get("@odata.nextLink")

    def walk(self, drive, meta, prefix, depth=0):
        drive = drive or meta.get("parentReference", {}).get("driveId", "")
        if "folder" in meta:
            if depth >= MAX_DEPTH:
                return
            base = f"{prefix}{segment(meta['name'])}/"
            for child in self.children(drive, meta["id"]):
                yield from self.walk(drive, child, base, depth + 1)
        elif "file" in meta and SIBLING not in meta["name"] and not meta["name"].startswith("~$"):
            yield {"id": meta["id"], "rel": prefix + segment(meta["name"]), "version": meta.get("cTag", ""),
                   "native": False, "drive": drive, "mime": meta.get("file", {}).get("mimeType", ""),
                   "parent": meta.get("parentReference", {}).get("id", ""), "etag": meta.get("eTag", "")}

    def listing(self, items):
        for item in items:
            yield from self.walk(item.get("drive", ""), self.meta(item.get("drive", ""), item["id"]), "")

    def download(self, entry):
        return call("GET", self.item_url(entry["drive"], entry["id"]) + "/content", self.token,
                    opener=self.opener)[2]

    def put(self, url, data, headers=None):
        if len(data) > SIMPLE_UPLOAD:
            raise DriveError("a file over 4 MB is not uploaded yet; it stays in the repository")
        _, _, body = call("PUT", url, self.token, data, {"Content-Type": "application/octet-stream",
                                                         **(headers or {})}, self.opener)
        return json.loads(body or b"{}")

    def upload(self, entry, data):
        answer = self.put(self.item_url(entry["drive"], entry["id"]) + "/content", data,
                          {"If-Match": entry["etag"]} if entry.get("etag") else None)
        return answer.get("cTag", "")

    def upload_beside(self, entry, name, data):
        url = self.item_url(entry["drive"], entry["parent"]) + f":/{q(name)}:/content"
        self.put(url, data)


def q(text):
    return urllib.parse.quote(text, safe="")


PROVIDERS = {"google": Google, "microsoft": Microsoft}


# What to do with each file

def plan(remote, main, records):
    """Decide, for each file the drive or the last run knows, what this run
    does: first (no record yet), pull, push, conflict (ours beside, then
    theirs down), forget, wait (pulled, not merged yet) or nothing. remote
    and records are keyed by the drive's id, since a path can change or be
    shared; main maps a path to the default branch's blob."""
    actions = {}
    for file_id in sorted(set(remote) | set(records)):
        entry, record = remote.get(file_id), records.get(file_id)
        here = main.get(record["rel"]) if record else main.get(entry["rel"])
        actions[file_id] = decide(entry, here, record)
    return actions


def decide(entry, main_blob, record):
    if entry is None:
        return "forget" if record else "nothing"
    if record is None:
        return "first"
    drive_moved = entry["version"] != record["version"]
    renamed = entry["rel"] != record["rel"]
    if record.get("pending"):
        if main_blob == record.get("base"):
            return "pull" if drive_moved or renamed else "wait"
        return "conflict"
    # A copy this sync once staged and the drive has since moved past is no
    # revision of ours: it is a pull request that merged late.
    main_moved = (main_blob is not None and main_blob != record["blob"]
                  and main_blob not in record.get("staged", []))
    stale = main_blob is not None and main_blob != record["blob"] and not main_moved
    if main_moved and (drive_moved or renamed):
        return "conflict"
    if drive_moved or renamed or stale:
        return "pull"
    return "push" if main_moved else "nothing"


def settle(records, main):
    """A pulled file whose copy is now on the default branch has merged."""
    for record in records.values():
        if record.get("pending") and main.get(record["rel"]) == record["blob"]:
            record["pending"] = False
            record.pop("base", None)


def sibling_name(rel):
    path = Path(rel)
    return f"{path.stem}{SIBLING}{path.suffix}"


# A run

def main_blobs(root, folder):
    base = Path(root) / folder
    if not base.is_dir():
        return {}
    return {p.relative_to(base).as_posix(): blob(p.read_bytes()) for p in base.rglob("*")
            if p.is_file() and not p.is_symlink() and not p.name.startswith("~$")}


def listing(provider, items):
    """The drive's files by id, each at one path, and the picked items the
    drive would not list: an id two picked items both reach is listed once,
    and two files of one name are told apart by a hash of the id, in the
    same way whatever order the drive lists them. One item deleted or no
    longer shared is named and the others still sync."""
    found, problems = {}, []
    for item in items:
        try:
            for entry in provider.listing([item]):
                found.setdefault(entry["id"], entry)
        except DriveError as problem:
            problems.append(f"{item.get('name') or item.get('id', '?')}: {problem}")
    taken = {}
    for file_id in sorted(found):
        rel = found[file_id]["rel"]
        if rel in taken:
            path = Path(rel)
            rel = f"{path.with_suffix('').as_posix()}-{blob(file_id.encode())[:8]}{path.suffix}"
        taken[rel] = file_id
        found[file_id] = {**found[file_id], "rel": rel}
    return found, problems


def held_back(action, entry, record, records):
    """What a run with part of the drive unlisted leaves for the next: a
    file it did not see may only sit in the item it could not list, and a
    path may have moved, or a new file taken an unseen file's path, only
    because a file of the same name went unseen."""
    if action == "forget":
        return True
    if record is None:
        return any(r.get("rel") == entry["rel"] for r in records.values())
    return entry is not None and entry["rel"] != record["rel"]


def run(config, state, root, out, provider):
    folder = config.get("folder", "papers")
    summary = {"pulled": [], "pushed": [], "beside": [], "waiting": [], "problems": []}
    remote, unlisted = listing(provider, config["items"])
    summary["problems"] += unlisted
    if unlisted and len(unlisted) == len(config["items"]):
        raise DriveError("the drive listed none of the picked items: " + "; ".join(unlisted[:3]))
    main = main_blobs(root, folder)
    if "files" in state and any("rel" not in r for r in state["files"].values()):
        state["files"] = {}
    # Records are paths inside one folder; another folder starts afresh.
    if state.get("folder", folder) != folder:
        state["files"] = {}
    state["folder"] = folder
    records = state.setdefault("files", {})
    settle(records, main)
    for file_id, action in plan(remote, main, records).items():
        if unlisted and held_back(action, remote.get(file_id), records.get(file_id), records):
            continue
        try:
            act(action, file_id, remote.get(file_id), main, records, Path(root) / folder,
                Path(out) / folder, provider, summary)
        except (DriveError, OSError) as problem:
            rel = (remote.get(file_id) or records.get(file_id) or {}).get("rel", file_id)
            summary["problems"].append(f"{rel}: {problem}")
    return summary


def act(action, file_id, entry, main, records, source, staging, provider, summary):
    record = records.get(file_id)
    if action == "forget":
        records.pop(file_id, None)
    elif action == "wait":
        stage(staging / entry["rel"], provider.download(entry))
        summary["waiting"].append(entry["rel"])
    elif action == "push":
        push(entry, record, main, source, provider, summary)
    elif action in ("first", "pull", "conflict"):
        ours = record["rel"] if record else entry["rel"]
        data = provider.download(entry)
        digest = blob(data)
        if action != "pull" and main.get(ours) not in (None, digest):
            put_beside(entry, ours, source, provider, summary)
        if digest == main.get(entry["rel"]):
            records[file_id] = keep(entry, digest, record)
            return
        stage(staging / entry["rel"], data)
        summary["pulled"].append(entry["rel"])
        records[file_id] = {**keep(entry, digest, record), "pending": True, "base": main.get(entry["rel"])}


def keep(entry, digest, record):
    staged = ((record or {}).get("staged", []) + [digest])[-5:]
    return {"rel": entry["rel"], "version": entry["version"], "blob": digest, "pending": False,
            "staged": staged}


def push(entry, record, main, source, provider, summary):
    data = (source / record["rel"]).read_bytes()
    if entry["native"]:
        provider.upload_beside(entry, sibling_name(entry["rel"]), data)
        summary["beside"].append(entry["rel"])
        version = entry["version"]
    else:
        version = provider.upload(entry, data) or entry["version"]
        summary["pushed"].append(entry["rel"])
    record.update({"version": version, "blob": main[record["rel"]], "pending": False})


def put_beside(entry, rel, source, provider, summary):
    """Our copy goes up beside the drive's before the drive's comes down;
    where it cannot, the drive's still comes down, since ours stays in the
    repository's history."""
    try:
        provider.upload_beside(entry, sibling_name(entry["rel"]), (source / rel).read_bytes())
        summary["beside"].append(entry["rel"])
    except DriveError as problem:
        summary["problems"].append(f"{rel}: our copy could not go up beside the drive's ({problem}); "
                                   "it stays in the repository's history")


def stage(path, data):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(data)


def checked_folder(config, state, root):
    """The folder the run may write, or a DriveError saying why not."""
    folder = config.get("folder", "papers")
    if not isinstance(folder, str) or not FOLDER.fullmatch(folder) or ".." in folder:
        raise DriveError(f"drive.json names a folder this run will not use: {folder!r}")
    top = Path(root) / folder.split("/")[0]
    if (top.exists() or top.is_symlink()) and not tracks(state, root, folder):
        raise DriveError(f"drive.json names {folder!r}, but {top.name!r} already holds this repository's "
                         "own files; pick a folder the repository does not have yet, or, to take back "
                         "a folder the sync kept before, delete the branch drive-state")
    return folder


def tracks(state, root, folder):
    """Whether the last run kept this folder. A state from before runs
    wrote their folder counts only where every file under it is one the
    sync brought in, and there is at least one."""
    if "folder" in state:
        return state["folder"] == folder
    files = state.get("files")
    records = {r.get("rel") for r in files.values() if isinstance(r, dict)} if isinstance(files, dict) else set()
    here = set(main_blobs(root, folder))
    return bool(here) and here <= records


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = parser.add_subparsers(dest="command", required=True)
    go = sub.add_parser("run")
    go.add_argument("--config", required=True)
    go.add_argument("--state", required=True)
    go.add_argument("--root", default=".")
    go.add_argument("--out", required=True)
    go.add_argument("--summary", required=True)
    gate = sub.add_parser("folder")
    gate.add_argument("--config", required=True)
    gate.add_argument("--state", required=True)
    gate.add_argument("--root", default=".")
    args = parser.parse_args(argv)
    config = json.loads(Path(args.config).read_text(encoding="utf-8"))
    state_path = Path(args.state)
    state = json.loads(state_path.read_text(encoding="utf-8")) if state_path.is_file() else {}
    if args.command == "folder":
        try:
            print(f"folder={checked_folder(config, state if isinstance(state, dict) else {}, args.root)}")
        except DriveError as problem:
            print(f"drive_sync: {problem}", file=sys.stderr)
            return 1
        return 0
    summary = {"pulled": [], "pushed": [], "beside": [], "waiting": [], "problems": []}
    code = 0
    try:
        kind = PROVIDERS[config["provider"]]
        token = access_token(os.environ["SYMPHONY_DRIVE_TOKEN_URL"], config["provider"],
                             os.environ["SYMPHONY_DRIVE_REFRESH"])
        summary = run(config, state, args.root, args.out, kind(token))
    except (DriveError, KeyError) as problem:
        print(f"drive_sync: {problem}", file=sys.stderr)
        code = 1
    finally:
        # What went up is recorded even when a later call failed, so the
        # next run never sends it up again beside itself.
        state_path.write_text(json.dumps(state, indent=1, sort_keys=True) + "\n", encoding="utf-8")
        Path(args.summary).write_text(json.dumps(summary, indent=1) + "\n", encoding="utf-8")
    print(json.dumps({k: len(v) for k, v in summary.items()}))
    return code


if __name__ == "__main__":
    sys.exit(main())
