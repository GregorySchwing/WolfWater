"""Word documents under review: read their comments, queue the revisions as
tasks, and write a worker's revisions back as tracked changes.

A paper kept as a .docx file in the repository is reviewed in Word, by
comments, and revised by a worker like any other task (docs/setup/
word-revisions.md). Claude cannot read a .docx file by itself, since it is a
zip of XML, so this tool turns each one into a review sheet beside it,
<name>.review.md: every paragraph numbered, with the comments anchored to it
underneath, which a worker reads and the pull request shows as a diff. The
worker answers with a file of revisions, <name>.revisions.json, a list of

  {"paragraph": 12, "old": "the paragraph as it is", "text": "as revised",
   "comment": "3", "why": "one line for the author"}

and this tool writes them into the .docx as tracked changes, word by word,
each run keeping its formatting, so the author accepts or rejects each one
in Word, and a "why" becomes a comment of Claude's on the paragraph. "old"
must match the paragraph exactly, so a revision meant for one paragraph never
lands on another after the text moved.

  python3 tools/docx_review.py extract FILE.docx          the JSON it reads
  python3 tools/docx_review.py sheet FILE.docx [...]      write the review sheets
  python3 tools/docx_review.py tasks FOLDER [--gh]        one task per document
  python3 tools/docx_review.py watch FILE.docx [...]      queue new comments (papers.yml)
  python3 tools/docx_review.py apply FILE.docx REVISIONS [--out OUT.docx]
  python3 tools/docx_review.py panel FILE.docx --lens "Name=REVISIONS" [...] [--out OUT.docx]
  python3 tools/docx_review.py sample OUT.docx            a small paper to try

A paragraph holding anything but plain text runs (a field, a picture, a
hyperlink, a tracked change already there) is refused by apply, which names
it, since rewriting it could lose what it holds; the author revises it by
hand. Comments marked done in Word are left out of the sheet and the tasks,
and so are Claude's own, which answer comments and are no new work.

watch is what .github/workflows/papers.yml runs when a document changes on
the default branch: a document's task is one issue for good, found by a
marker in its body, which also lists the comments it has seen by a key that
survives Word renumbering them; new ones are told on the issue, which is
reopened or queued again, and the workflow then wakes the workers.
"""

import argparse
import copy
import datetime
import difflib
import hashlib
import io
import json
import re
import subprocess
import sys
import zipfile
from pathlib import Path
from xml.etree import ElementTree as ET
from xml.sax.saxutils import escape

W_URI = "http://schemas.openxmlformats.org/wordprocessingml/2006/main"
W = "{%s}" % W_URI
MC = "{http://schemas.openxmlformats.org/markup-compatibility/2006}"
W14 = "{http://schemas.microsoft.com/office/word/2010/wordml}"
W15 = "{http://schemas.microsoft.com/office/word/2012/wordml}"
XML_SPACE = "{http://www.w3.org/XML/1998/namespace}space"

DOCUMENT = "word/document.xml"
COMMENTS = "word/comments.xml"
COMMENTS_EXTENDED = "word/commentsExtended.xml"

TRACKED = {W + "ins": "ins", W + "moveTo": "ins", W + "del": "del", W + "moveFrom": "del"}

# What a paragraph and its runs may hold for apply to rewrite it safely.
PLAIN_PARAGRAPH = {W + "pPr", W + "r", W + "commentRangeStart", W + "commentRangeEnd",
                   W + "bookmarkStart", W + "bookmarkEnd", W + "proofErr"}
PLAIN_RUN = {W + "rPr", W + "t", W + "tab", W + "commentReference",
             W + "lastRenderedPageBreak"}
MARKERS = {W + "commentRangeStart", W + "commentRangeEnd", W + "bookmarkStart",
           W + "bookmarkEnd"}


class Refused(Exception):
    """A revision this tool will not write, with the reason for the person."""


# Reading

def read_parts(path):
    with zipfile.ZipFile(path) as archive:
        return {info.filename: archive.read(info.filename) for info in archive.infolist()}


def parse(data):
    """Parse a part, registering its prefixes so they survive a rewrite."""
    for _, (prefix, uri) in ET.iterparse(io.BytesIO(data), events=["start-ns"]):
        try:
            ET.register_namespace(prefix, uri)
        except ValueError:
            pass
    return ET.fromstring(data)


def serialize(root, original):
    """Write a part back with the original declaration and root start tag,
    since Word's mc:Ignorable names prefixes ElementTree would drop, adding
    every declaration ElementTree hoisted to the root from deeper down (a
    picture's a: and pic:), and refuse a result that does not parse."""
    text = ET.tostring(root, encoding="unicode")
    start_tag = re.compile(r"<[A-Za-z_][^\s>/]*[^>]*?/?>")
    declaration = re.match(r"<\?xml[^>]*\?>\s*", original.decode("utf-8"))
    rest = original.decode("utf-8")[declaration.end() if declaration else 0:]
    original_tag = start_tag.search(rest).group(0)
    written_tag = start_tag.search(text).group(0)
    declared = set(re.findall(r"xmlns(?::([\w.-]+))?=", original_tag))
    missing = [d for d in re.findall(r'\sxmlns(?::[\w.-]+)?="[^"]*"', written_tag)
               if (re.match(r"\sxmlns(?::([\w.-]+))?=", d).group(1) or "") not in declared]
    merged = re.sub(r"\s*/?>$", "", original_tag) + "".join(missing) + ">"
    text = text.replace(written_tag, merged, 1)
    head = declaration.group(0) if declaration else ""
    data = (head + text).encode("utf-8")
    ET.fromstring(data)
    return data


def paragraphs(document_root):
    """The body's paragraphs in order, tables' included, leaving out the
    copy Word keeps in mc:Fallback, which it never shows."""
    body = document_root.find(W + "body")
    return [p for p, _ in shown(body) if p.tag == W + "p"] if body is not None else []


def shown(element):
    """Every element Word shows, in document order, with its paragraph."""
    for child in element:
        if child.tag == MC + "Fallback":
            continue
        yield child, None
        yield from shown(child)


def paragraph_text(paragraph, view="accepted"):
    """The paragraph's text with tracked changes accepted, or rejected."""
    parts = []
    for element, inside in walk(paragraph):
        if view == "accepted" and inside == "del" or view == "rejected" and inside == "ins":
            continue
        if element.tag in (W + "t", W + "delText"):
            parts.append(element.text or "")
        elif element.tag == W + "tab":
            parts.append("\t")
        elif element.tag in (W + "br", W + "cr"):
            parts.append("\n")
    return "".join(parts)


def walk(element, inside=None):
    """Every descendant in document order, with the tracked change it is in."""
    for child in element:
        if child.tag in (W + "p", W + "pPr", W + "rPr", MC + "Fallback"):
            continue
        state = TRACKED.get(child.tag, inside)
        yield child, state
        yield from walk(child, state)


def paragraph_style(paragraph):
    style = paragraph.find(W + "pPr/" + W + "pStyle")
    return style.get(W + "val") if style is not None else ""


def comment_records(parts):
    if COMMENTS not in parts:
        return []
    done, parents = comment_threads(parts)
    records = []
    for comment in parse(parts[COMMENTS]).iter(W + "comment"):
        texts = [paragraph_text(p) for p in comment.iter(W + "p")]
        ids = [p.get(W14 + "paraId") for p in comment.iter(W + "p")]
        last = ids[-1] if ids else None
        records.append({
            "id": comment.get(W + "id"),
            "author": comment.get(W + "author", ""),
            "date": comment.get(W + "date", ""),
            "text": "\n".join(texts).strip(),
            "done": done.get(last, False),
            "para_id": last,
            "parent_para_id": parents.get(last),
        })
    by_para = {r["para_id"]: r["id"] for r in records if r["para_id"]}
    for record in records:
        record["reply_to"] = by_para.get(record.pop("parent_para_id"))
        record.pop("para_id")
    return records


def comment_threads(parts):
    """Which comments Word marks done, and which reply to which."""
    if COMMENTS_EXTENDED not in parts:
        return {}, {}
    done, parents = {}, {}
    for entry in parse(parts[COMMENTS_EXTENDED]).iter(W15 + "commentEx"):
        para = entry.get(W15 + "paraId")
        done[para] = entry.get(W15 + "done") == "1"
        if entry.get(W15 + "paraIdParent"):
            parents[para] = entry.get(W15 + "paraIdParent")
    return done, parents


def anchors(document_root):
    """For each comment id, the paragraphs it spans and the text it marks,
    with a range Word starts outside any paragraph, as before a table."""
    found, open_ids = {}, set()
    body = document_root.find(W + "body")
    index_of = {id(p): i for i, p in enumerate(paragraphs(document_root))}
    for element, index, inside in positions(body, index_of):
        cid = element.get(W + "id")
        if element.tag == W + "commentRangeStart":
            open_ids.add(cid)
            found.setdefault(cid, {"paragraphs": [], "anchor": ""})
        elif element.tag == W + "commentRangeEnd":
            open_ids.discard(cid)
        elif element.tag == W + "commentReference":
            found.setdefault(cid, {"paragraphs": [], "anchor": ""})
            open_ids.add(cid)
        elif inside != "del" and element.tag == W + "t":
            for open_id in open_ids:
                found[open_id]["anchor"] += element.text or ""
        elif inside != "del" and element.tag == W + "tab":
            for open_id in open_ids:
                found[open_id]["anchor"] += "\t"
        for open_id in open_ids:
            if index is not None and index not in found[open_id]["paragraphs"]:
                found[open_id]["paragraphs"].append(index)
        if element.tag == W + "commentReference":
            open_ids.discard(cid)
    return found


def positions(element, index_of, index=None, inside=None):
    """Every element Word shows, with its paragraph's number and the tracked
    change it is in, leaving out paragraph and run properties."""
    for child in element:
        if child.tag in (W + "pPr", W + "rPr", MC + "Fallback"):
            continue
        here = index_of.get(id(child), index)
        state = TRACKED.get(child.tag, inside)
        yield child, here, state
        yield from positions(child, index_of, here, state)


def extract(path):
    parts = read_parts(path)
    root = parse(parts[DOCUMENT])
    marks = anchors(root)
    comments = comment_records(parts)
    for comment in comments:
        comment.update(marks.get(comment["id"], {"paragraphs": [], "anchor": ""}))
    return {
        "file": str(path),
        "paragraphs": [{"index": i, "style": paragraph_style(p), "text": paragraph_text(p)}
                       for i, p in enumerate(paragraphs(root))],
        "comments": comments,
    }


def open_comments(review):
    """The comments still to act on: not done, and not in a done thread."""
    done = {c["id"] for c in review["comments"] if c["done"]}
    return [c for c in review["comments"] if c["id"] not in done and c["reply_to"] not in done]


# The review sheet

def sheet(review):
    name = Path(review["file"]).name
    lines = [f"# Review sheet: {name}", "",
             f"Written by `python3 tools/docx_review.py sheet` from {name}; do not edit it by hand.",
             "Revise by writing " + Path(name).stem + ".revisions.json beside it (docs/setup/word-revisions.md).", ""]
    by_paragraph = {}
    for comment in open_comments(review):
        last = comment["paragraphs"][-1] if comment["paragraphs"] else None
        by_paragraph.setdefault(last, []).append(comment)
    for paragraph in review["paragraphs"]:
        if not paragraph["text"].strip() and paragraph["index"] not in by_paragraph:
            continue
        style = f" ({paragraph['style']})" if paragraph["style"] else ""
        lines.append(f"[{paragraph['index']}]{style} {paragraph['text']}")
        for comment in by_paragraph.get(paragraph["index"], []):
            lines.extend(comment_lines(comment))
        lines.append("")
    for comment in by_paragraph.get(None, []):
        lines.extend(comment_lines(comment) + [""])
    return "\n".join(lines).rstrip() + "\n"


def comment_lines(comment):
    reply = f", replying to comment {comment['reply_to']}" if comment["reply_to"] else ""
    span = ", ".join(str(i) for i in comment["paragraphs"]) or "none"
    lines = [f"> Comment {comment['id']} by {comment['author']}{reply}, on paragraphs {span}:"]
    if comment["anchor"].strip():
        lines.append(f'> marks "{comment["anchor"].strip()}"')
    lines.extend("> " + line for line in comment["text"].splitlines() or [""])
    return lines


def write_sheets(paths):
    written = []
    for path in paths:
        review = extract(path)
        target = Path(path).with_suffix(".review.md")
        target.write_text(sheet(review), encoding="utf-8", newline="\n")
        written.append(str(target))
    return written


# Tasks

# Comments a worker wrote are answers, never new work, or each revision
# would queue the next.
REVISER = "Claude"
DOCUMENT_MARK = "<!-- symphony-document: {} -->"
COMMENTS_MARK = re.compile(r"<!-- symphony-comments: ([0-9a-f ]*) -->")


def comment_key(comment):
    """A comment's identity across saves: Word renumbers w:id each time it
    saves, while the author, the date and the words stay."""
    text = "\x1f".join((comment["author"], comment["date"], comment["text"]))
    return hashlib.sha1(text.encode("utf-8")).hexdigest()[:12]


def awaiting(path):
    """The open comments of a document that a person wrote, by key."""
    return {comment_key(c): c for c in open_comments(extract(path)) if c["author"] != REVISER}


def task_title(path):
    return f"Task: revise {Path(path).as_posix()}"


def task_body(path, comments):
    path = Path(path)
    stem = path.with_suffix("").as_posix()
    keys = " ".join(sorted(comments))
    return (
        "### Goal\n\n"
        f"{path.name} answers its open review comments, as tracked changes the author can accept "
        "or reject in Word:\n\n" + comment_list(comments.values()) + "\n\n"
        "### Acceptance\n\n"
        f"Read {stem}.review.md (run `python3 tools/docx_review.py sheet {path.as_posix()}` if it is "
        f"missing or older than the document). Answer the comments this issue lists, the body's and "
        f"those in later comments of this issue not answered yet. Write {stem}.revisions.json, one "
        "entry per paragraph changed, with the paragraph's exact old text, the new text, the comment it "
        f"answers and a one-line why; run `python3 tools/docx_review.py apply {path.as_posix()} "
        f"{stem}.revisions.json` and then `python3 tools/docx_review.py sheet {path.as_posix()}`; commit "
        "all three. A comment that asks for something the text cannot give (new data, a figure, the "
        "author's judgement) gets an entry with no change, only its why. Where the document changed on "
        "the default branch meanwhile, take that copy and apply the revisions again. "
        "`bash scripts/check.sh` passes.\n\n"
        "### Files likely touched\n\n"
        f"{path.as_posix()}\n{stem}.review.md\n{stem}.revisions.json\n\n"
        + DOCUMENT_MARK.format(path.as_posix()) + "\n"
        + f"<!-- symphony-comments: {keys} -->\n"
    )


def comment_list(comments):
    return "\n".join(f"- {c['author']}: {one_line(c['text'])}"
                      + (f' (on "{one_line(c["anchor"])}")' if c["anchor"].strip() else "")
                      for c in comments)


def one_line(text):
    flat = " ".join(text.split())
    return flat if len(flat) <= 120 else flat[:117] + "..."


def documents(folder):
    return [p for p in sorted(Path(folder).rglob("*.docx")) if not p.name.startswith("~$")]


def tasks(folder):
    found = []
    for path in documents(folder):
        comments = awaiting(path)
        if comments:
            found.append({"title": task_title(path), "body": task_body(path, comments), "file": str(path)})
    return found


def gh(*args):
    result = subprocess.run(["gh", *args], check=True, capture_output=True, text=True)
    return result.stdout


TRUSTED = {"OWNER", "MEMBER", "COLLABORATOR"}
WORKFLOW_BOT = "github-actions[bot]"
HOLDING = {"claimed", "in-review"}


def task_issues(run=gh):
    """Every issue that may be a document's task: opened by this repository's
    workflows or by someone who may write, so nobody else's issue, however
    it is marked, can be taken for one and steer a worker."""
    lines = run("api", "--paginate", "repos/{owner}/{repo}/issues?state=all&per_page=100", "--jq",
                ".[] | select(.pull_request | not) | {number, state, body, author: .user.login, "
                "association: .author_association, labels: [.labels[].name]}")
    issues = [json.loads(line) for line in lines.splitlines() if line.strip()]
    trusted = [i for i in issues if i["author"] == WORKFLOW_BOT or i["association"] in TRUSTED]
    return sorted(trusted, key=lambda i: i["number"])


def watch(paths, run=gh):
    """Queue the new comments of each document: open its task, or tell its
    task of the comments it has not listed yet and queue it again. Returns
    the documents queued, so the caller can wake the workers; a document
    that cannot be read is named on standard error and the rest go on."""
    issues = task_issues(run)
    queued = []
    for path in paths:
        try:
            if Path(path).is_file() and queue_document(path, issues, run):
                queued.append(str(path))
        except (Refused, KeyError, zipfile.BadZipFile, ET.ParseError, subprocess.CalledProcessError) as problem:
            print(f"docx_review: {path}: {problem}", file=sys.stderr)
    return queued


def queue_document(path, issues, run):
    comments = awaiting(path)
    mark = DOCUMENT_MARK.format(Path(path).as_posix())
    issue = next((i for i in issues if mark in (i.get("body") or "").splitlines()), None)
    if issue is None:
        if not comments:
            return False
        run("issue", "create", "--title", task_title(path), "--body", task_body(path, comments), "--label", "todo")
        return True
    listed = COMMENTS_MARK.search(issue["body"])
    known = set(listed.group(1).split()) if listed else set()
    new = {k: c for k, c in comments.items() if k not in known}
    if not new:
        return False
    keys = f"<!-- symphony-comments: {' '.join(sorted(known | set(new)))} -->"
    body = COMMENTS_MARK.sub(keys, issue["body"]) if listed else issue["body"].rstrip("\n") + "\n" + keys + "\n"
    number = str(issue["number"])
    run("issue", "edit", number, "--body", body)
    issue["body"] = body
    run("issue", "comment", number, "--body", "New review comments on this document:\n\n"
        + comment_list(new.values()))
    labels = set(issue.get("labels", []))
    if issue["state"].upper() == "CLOSED":
        # A finished task still carries the labels of its last claim; a
        # reopened one is queued afresh.
        run("issue", "reopen", number)
        stale = sorted(l for l in labels if l in HOLDING or l.startswith(("worker:", "account:", "agent:")))
        for label in stale:
            run("issue", "edit", number, "--remove-label", label)
        labels -= set(stale)
        issue["state"] = "OPEN"
    if not labels & (HOLDING | {"todo"}):
        run("issue", "edit", number, "--add-label", "todo")
    return True


# Applying revisions as tracked changes

def apply(path, revisions, out=None, author="Claude", now=None):
    parts = read_parts(path)
    original = parts[DOCUMENT]
    root = parse(original)
    comments_root = parse(parts[COMMENTS]) if COMMENTS in parts else None
    stamp = (now or datetime.datetime.now(datetime.timezone.utc)).strftime("%Y-%m-%dT%H:%M:%SZ")
    ids = IdSource(root, comments_root, *(parse(data) for name, data in parts.items()
                                          if name.startswith("word/") and name.endswith(".xml")
                                          and name not in (DOCUMENT, COMMENTS)))
    all_paragraphs = paragraphs(root)
    numbers = [revision["paragraph"] for revision in revisions]
    if len(numbers) != len(set(numbers)):
        raise Refused("two revisions name one paragraph; write one entry per paragraph")
    for revision in revisions:
        index = revision["paragraph"]
        if not 0 <= index < len(all_paragraphs):
            raise Refused(f"paragraph {index} does not exist")
        paragraph = all_paragraphs[index]
        current = paragraph_text(paragraph)
        if current != revision["old"]:
            raise Refused(f"paragraph {index} no longer reads as the revision's old text")
        bad = [ch for ch in revision["text"] if ord(ch) < 32 and ch != "\t"]
        if bad:
            raise Refused(f"paragraph {index}: the new text holds a line break or control character; "
                          "a paragraph's text is one line, tabs allowed")
        if revision.get("why") and comments_root is None:
            raise Refused(f"paragraph {index}: the document has no comments part to answer in")
        if revision["text"] != current:
            rewrite(paragraph, revision["text"], author, stamp, ids, index)
        if revision.get("why") and comments_root is not None:
            note = revision["why"]
            if revision.get("comment"):
                note = f"On comment {revision['comment']}: {note}"
            add_comment(paragraph, comments_root, note, author, stamp, ids)
    parts[DOCUMENT] = serialize(root, original)
    if comments_root is not None:
        parts[COMMENTS] = serialize(comments_root, parts[COMMENTS])
    write_parts(path, parts, out or path)


AUTHOR = re.compile(r"[^\x00-\x1f<>&\"]{1,60}")


def panel(path, lenses, out=None, now=None):
    """Several reviewers' revisions, each under its own author, one after
    another: each lens is a name and its revisions, written as apply writes
    them, so the author sees in Word whose every change and comment is. A
    revision an earlier lens has already changed the paragraph for, or that
    no longer reads as its old text, or that a paragraph's content refuses,
    becomes that lens's comment on the paragraph, "Suggests: ...", so no
    reviewer's view is lost and none overwrites another's. Returns what each
    lens did."""
    target = out or path
    report, source = [], path
    for author, revisions in lenses:
        if not AUTHOR.fullmatch(author) or author != author.strip():
            raise Refused(f"the author name {author[:20]!r} is not one line of up to 60 plain characters")
        kept, revised, suggested, skipped = sort_lens(source, revisions)
        apply(source, kept, target, author, now)
        source = target
        report.append({"author": author, "revised": revised, "suggested": suggested,
                       "commented": sum(1 for k in kept if k.get("why")) - suggested,
                       "skipped": skipped})
    return report


def sort_lens(path, revisions):
    """One lens's revisions as apply takes them: those that still fit, and
    the rest as comments, with what could be neither said."""
    parts = read_parts(path)
    all_paragraphs = paragraphs(parse(parts[DOCUMENT]))
    has_comments = COMMENTS in parts
    kept, seen, suggested, skipped, revised = [], set(), 0, [], 0
    for revision in revisions:
        index, text = revision.get("paragraph"), str(revision.get("text", ""))
        if not isinstance(index, int) or not 0 <= index < len(all_paragraphs) or index in seen:
            skipped.append(f"paragraph {index}: not in the document, or named twice")
            continue
        if any(ord(ch) < 32 and ch != "\t" for ch in text):
            skipped.append(f"paragraph {index}: the new text holds a line break")
            continue
        paragraph = all_paragraphs[index]
        current = paragraph_text(paragraph)
        why = str(revision.get("why", ""))
        try:
            old_atoms(paragraph, index)
            fits = current == revision.get("old")
        except Refused:
            fits = False
        comment_only = text in (current, revision.get("old"))
        if fits and (text != current or why and has_comments):
            kept.append({**revision, "text": text, "why": why if has_comments else ""})
            revised += text != current
        elif not has_comments:
            skipped.append(f"paragraph {index}: changed since the review sheet or not plain text, and the "
                           "document has no comments part to suggest it in")
            continue
        elif comment_only and why:
            # A note with no change of its own: a comment holds it, on a
            # paragraph an earlier lens has changed too.
            kept.append({"paragraph": index, "old": current, "text": current, "why": why,
                         "comment": revision.get("comment", "")})
        elif comment_only:
            skipped.append(f"paragraph {index}: nothing to change or say")
            continue
        else:
            note = f"Suggests: {text}" + (f" Why: {why}" if why else "")
            kept.append({"paragraph": index, "old": current, "text": current, "why": note,
                         "comment": revision.get("comment", "")})
            suggested += 1
        seen.add(index)
    return kept, revised, suggested, skipped


class IdSource:
    """Fresh w:id values above every one the document and comments use."""

    def __init__(self, *roots):
        used = [int(e.get(W + "id")) for r in roots if r is not None for e in r.iter()
                if (e.get(W + "id") or "").isdigit() and len(e.get(W + "id")) < 10]
        self.next = max(used, default=0) + 1

    def take(self):
        self.next += 1
        return str(self.next - 1)


def rewrite(paragraph, new_text, author, stamp, ids, index):
    atoms, markers, runs = old_atoms(paragraph, index)
    old_text = "".join(ch for ch, _ in atoms)
    pieces = diff_pieces(old_text, new_text)
    for child in list(paragraph):
        if child.tag != W + "pPr":
            paragraph.remove(child)
    emitted = []  # (kind, run index or marker element, char)
    pending = list(markers)

    def markers_through(position):
        while pending and pending[0][0] < position:
            emitted.append(("marker", pending.pop(0)[1], None))

    for kind, old_start, old_end, inserted in pieces:
        if kind == "ins":
            # A marker where the insertion stands, such as a comment's end
            # after the words it replaces, comes after it, so the comment
            # still covers the replacement.
            markers_through(old_start)
            source = atoms[old_start - 1][1] if old_start > 0 else (atoms[0][1] if atoms else None)
            emitted.extend(("ins", source, ch) for ch in inserted)
            continue
        for offset in range(old_start, old_end):
            markers_through(offset + 1)
            emitted.append((kind, atoms[offset][1], atoms[offset][0]))
    markers_through(len(atoms) + 1)
    for element in group(emitted, runs, author, stamp, ids):
        paragraph.append(element)


def old_atoms(paragraph, index):
    """The paragraph's characters, each with the run it came from, and the
    markers (comment ranges, bookmarks, comment references) by position."""
    atoms, markers, runs = [], [], []
    for child in paragraph:
        if child.tag not in PLAIN_PARAGRAPH:
            raise Refused(f"paragraph {index} holds {local(child.tag)}; revise it by hand")
        if child.tag in MARKERS:
            markers.append((len(atoms), child))
        elif child.tag == W + "r":
            kinds = {c.tag for c in child}
            if kinds - PLAIN_RUN:
                raise Refused(f"paragraph {index} holds {local(sorted(kinds - PLAIN_RUN)[0])}; revise it by hand")
            if child.find(W + "rPr/" + W + "rPrChange") is not None:
                raise Refused(f"paragraph {index} holds a tracked format change; accept it first")
            if W + "commentReference" in kinds:
                markers.append((len(atoms), child))
                continue
            runs.append(child.find(W + "rPr"))
            for element in child:
                text = element.text or "" if element.tag == W + "t" else "\t" if element.tag == W + "tab" else ""
                atoms.extend((ch, len(runs) - 1) for ch in text)
    return atoms, markers, runs


def local(tag):
    return tag.split("}")[-1]


def diff_pieces(old, new):
    """Word-level differences as (kind, old start, old end, inserted text)."""
    old_tokens, new_tokens = tokens(old), tokens(new)
    old_at, new_at = offsets(old_tokens), offsets(new_tokens)
    matcher = difflib.SequenceMatcher(a=old_tokens, b=new_tokens, autojunk=False)
    pieces = []
    for tag, i1, i2, j1, j2 in joined(matcher.get_opcodes(), old_tokens):
        if tag == "equal":
            pieces.append(("eq", old_at[i1], old_at[i2], ""))
            continue
        if i2 > i1:
            pieces.append(("del", old_at[i1], old_at[i2], ""))
        if j2 > j1:
            pieces.append(("ins", old_at[i2], old_at[i2], new[new_at[j1]:new_at[j2]]))
    return pieces


def joined(opcodes, old_tokens):
    """Changes kept apart only by spaces become one, so a rewritten phrase
    reads as one deletion and one insertion, not word by word."""
    def spaces(op):
        return op[0] == "equal" and all(not t.strip() for t in old_tokens[op[1]:op[2]])

    out = []
    for index, op in enumerate(opcodes):
        between = 0 < index < len(opcodes) - 1 and opcodes[index + 1][0] != "equal"
        changes = op[0] != "equal" or (spaces(op) and between)
        if changes and out and out[-1][0] == "replace":
            out[-1] = ("replace", out[-1][1], op[2], out[-1][3], op[4])
        elif changes:
            out.append(("replace", *op[1:]))
        else:
            out.append(op)
    return out


def tokens(text):
    return re.findall(r"\s+|\w+|[^\w\s]", text)


def offsets(token_list):
    at = [0]
    for token in token_list:
        at.append(at[-1] + len(token))
    return at


def group(emitted, runs, author, stamp, ids):
    """Consecutive characters of one kind and one run become one run; runs
    of a deletion or an insertion go inside one w:del or w:ins."""
    out, wrapper, current = [], None, None
    for kind, source, ch in emitted:
        if kind == "marker":
            out.append(source)
            wrapper = current = None
            continue
        if current is None or current["kind"] != kind or current["source"] != source:
            run = make_run(runs[source] if source is not None and runs else None)
            current = {"kind": kind, "source": source, "run": run, "text": ""}
            if kind == "eq":
                out.append(run)
                wrapper = None
            else:
                if wrapper is None or wrapper.tag != W + kind:
                    wrapper = ET.Element(W + kind, {W + "id": ids.take(), W + "author": author,
                                                    W + "date": stamp})
                    out.append(wrapper)
                wrapper.append(run)
        add_char(current, ch)
        if kind == "eq":
            wrapper = None
    return out


def make_run(rpr):
    run = ET.Element(W + "r")
    if rpr is not None:
        run.append(copy.deepcopy(rpr))
    return run


def add_char(current, ch):
    run = current["run"]
    text_tag = W + ("delText" if current["kind"] == "del" else "t")
    if ch == "\t":
        run.append(ET.Element(W + "tab"))
        return
    last = run[-1] if len(run) else None
    if last is None or last.tag != text_tag:
        last = ET.SubElement(run, text_tag, {XML_SPACE: "preserve"})
        last.text = ""
    last.text += ch


def add_comment(paragraph, comments_root, text, author, stamp, ids):
    cid = ids.take()
    comment = ET.SubElement(comments_root, W + "comment", {
        W + "id": cid, W + "author": author, W + "date": stamp, W + "initials": author[:2]})
    ET.SubElement(ET.SubElement(ET.SubElement(comment, W + "p"), W + "r"), W + "t").text = text
    start = ET.Element(W + "commentRangeStart", {W + "id": cid})
    first = 1 if len(paragraph) and paragraph[0].tag == W + "pPr" else 0
    paragraph.insert(first, start)
    paragraph.append(ET.Element(W + "commentRangeEnd", {W + "id": cid}))
    reference = ET.SubElement(paragraph, W + "r")
    ET.SubElement(reference, W + "commentReference", {W + "id": cid})


def write_parts(source, parts, target):
    with zipfile.ZipFile(source) as archive:
        order = [info for info in archive.infolist()]
    buffer = io.BytesIO()
    with zipfile.ZipFile(buffer, "w", zipfile.ZIP_DEFLATED) as archive:
        for info in order:
            archive.writestr(info, parts[info.filename])
    Path(target).write_bytes(buffer.getvalue())


# A sample paper

SAMPLE_PARAGRAPHS = [
    ("Title", [("A small paper for trying the review loop", None)]),
    (None, [("We measured lumbar levels in ", None), ("forty", "b"),
            (" radiographs and found the detector agrees with readers most of the time.", None)]),
    (None, [("The method is described in the next section, and it is very novel.", None)]),
    (None, [("Results are shown in Table 1.", None)]),
]
SAMPLE_COMMENTS = [
    (1, "agrees with readers most of the time", "Reviewer 1", "Give the agreement as a number, kappa 0.81."),
    (2, "it is very novel", "Reviewer 2", "Avoid claiming novelty; say what is new."),
    (3, "Table 1", "Reviewer 2", "Table 1 is missing from this draft."),
]


def sample(out):
    """Write a small .docx with a bold run and three comments, for trying the loop."""
    body, comments = [], []
    for index, (style, pieces) in enumerate(SAMPLE_PARAGRAPHS):
        xml = "<w:p>"
        if style:
            xml += f'<w:pPr><w:pStyle w:val="{style}"/></w:pPr>'
        for text, fmt in pieces:
            rpr = "<w:rPr><w:b/></w:rPr>" if fmt == "b" else ""
            for cid, (at, anchor, _, _) in enumerate(SAMPLE_COMMENTS):
                if at == index and anchor in text:
                    before, text = text.split(anchor, 1)
                    xml += run_xml(before, rpr) + f'<w:commentRangeStart w:id="{cid}"/>'
                    xml += run_xml(anchor, rpr) + f'<w:commentRangeEnd w:id="{cid}"/>'
                    xml += f'<w:r><w:commentReference w:id="{cid}"/></w:r>'
            xml += run_xml(text, rpr)
        body.append(xml + "</w:p>")
    for cid, (_, _, author, text) in enumerate(SAMPLE_COMMENTS):
        comments.append(f'<w:comment w:id="{cid}" w:author="{author}" w:date="2026-10-18T00:00:00Z" '
                        f'w:initials="R"><w:p><w:r><w:t>{escape(text)}</w:t></w:r></w:p></w:comment>')
    ns = f'xmlns:w="{W_URI}"'
    rel = "http://schemas.openxmlformats.org/officeDocument/2006/relationships"
    files = {
        "[Content_Types].xml": (
            '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
            '<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types">'
            '<Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/>'
            '<Default Extension="xml" ContentType="application/xml"/>'
            '<Override PartName="/word/document.xml" ContentType="application/vnd.openxmlformats-'
            'officedocument.wordprocessingml.document.main+xml"/>'
            '<Override PartName="/word/comments.xml" ContentType="application/vnd.openxmlformats-'
            'officedocument.wordprocessingml.comments+xml"/></Types>'),
        "_rels/.rels": (
            '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
            '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">'
            f'<Relationship Id="rId1" Type="{rel}/officeDocument" Target="word/document.xml"/>'
            '</Relationships>'),
        "word/_rels/document.xml.rels": (
            '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
            '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">'
            f'<Relationship Id="rId1" Type="{rel}/comments" Target="comments.xml"/></Relationships>'),
        DOCUMENT: (f'<?xml version="1.0" encoding="UTF-8" standalone="yes"?><w:document {ns}>'
                   f'<w:body>{"".join(body)}</w:body></w:document>'),
        COMMENTS: (f'<?xml version="1.0" encoding="UTF-8" standalone="yes"?><w:comments {ns}>'
                   f'{"".join(comments)}</w:comments>'),
    }
    with zipfile.ZipFile(out, "w", zipfile.ZIP_DEFLATED) as archive:
        for name, text in files.items():
            archive.writestr(name, text)


def run_xml(text, rpr):
    return f'<w:r>{rpr}<w:t xml:space="preserve">{escape(text)}</w:t></w:r>' if text else ""


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = parser.add_subparsers(dest="command", required=True)
    sub.add_parser("extract").add_argument("file")
    sub.add_parser("sheet").add_argument("files", nargs="+")
    queue = sub.add_parser("tasks")
    queue.add_argument("folder")
    queue.add_argument("--gh", action="store_true", help="open or update the issues with gh")
    watching = sub.add_parser("watch")
    watching.add_argument("files", nargs="*")
    write = sub.add_parser("apply")
    write.add_argument("file")
    write.add_argument("revisions")
    write.add_argument("--out")
    write.add_argument("--author", default="Claude")
    lensing = sub.add_parser("panel")
    lensing.add_argument("file")
    lensing.add_argument("--lens", action="append", required=True, help="Author=revisions.json, once per reviewer")
    lensing.add_argument("--out")
    sub.add_parser("sample").add_argument("out")
    args = parser.parse_args(argv)
    try:
        if args.command == "extract":
            print(json.dumps(extract(args.file), indent=2, ensure_ascii=False))
        elif args.command == "sheet":
            print("\n".join(write_sheets(args.files)))
        elif args.command == "tasks":
            if args.gh:
                print("\n".join(watch(documents(args.folder))))
            else:
                print(json.dumps(tasks(args.folder), indent=2))
        elif args.command == "watch":
            print("\n".join(watch(args.files)))
        elif args.command == "apply":
            revisions = json.loads(Path(args.revisions).read_text(encoding="utf-8"))
            apply(args.file, revisions, args.out, args.author)
            print(args.out or args.file)
        elif args.command == "panel":
            lenses = []
            for lens in args.lens:
                author, sep, file = lens.partition("=")
                if not sep or not author:
                    raise Refused(f"--lens takes Author=revisions.json, got {lens[:40]!r}")
                lenses.append((author, json.loads(Path(file).read_text(encoding="utf-8"))))
            print(json.dumps(panel(args.file, lenses, args.out), indent=1))
        elif args.command == "sample":
            sample(args.out)
            print(args.out)
    except (Refused, KeyError, zipfile.BadZipFile) as problem:
        print(f"docx_review: {problem}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
