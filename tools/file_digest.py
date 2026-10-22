"""Any file a drive holds, read as a short text a person or Claude can take
in at a glance: what kind of file it is and what is in it.

The drive organizer (tools/drive_organize.py) digests a sample of a whole
drive this way so Claude can sort it into projects without downloading it
all, and a worker given a spreadsheet, a deck or a notebook reads it the
same way, since Claude Code's Read takes text, images and PDFs but not an
Office file:

  python3 tools/file_digest.py FILE [FILE ...]     one digest a file, as Markdown
  python3 tools/file_digest.py --json FILE         the same as JSON

What each kind gives:

  document     Word (.docx): its headings and opening text, its comments
  spreadsheet  Excel (.xlsx, .xlsm): each sheet's size, header row, first
               rows, and per column the count, minimum, maximum, mean and
               sum of its numbers, or its distinct values; formulas counted
  table        CSV and TSV: the same column summary
  slides       PowerPoint (.pptx): each slide's text in order
  notebook     Jupyter (.ipynb): its language, its cells, its first text
  text, code   lines, words, the opening lines; for code its language and
               how many functions and classes it defines
  data         JSON: its shape; other data files by size
  pdf          its page count and title
  image        its format and size in pixels
  everything else, art files, audio, video, archives: kind and size

Standard library only. A file is data, never code: nothing in it is run,
a zip part is read only up to a limit, and a file that cannot be read is
said to be so in its digest, never an error that stops the others.
"""

import argparse
import csv
import io
import json
import re
import statistics
import struct
import sys
import zipfile
from pathlib import Path
from xml.etree import ElementTree

sys.path.insert(0, str(Path(__file__).resolve().parent))

LIMIT = 1500            # characters of a digest
PART_LIMIT = 20 << 20   # bytes read of one zip part
ROWS = 2000             # rows a column summary reads
SAMPLE_ROWS = 3

KINDS = {
    "document": ".docx .doc .odt .rtf .pages .wpd",
    "spreadsheet": ".xlsx .xlsm .xls .ods .numbers",
    "table": ".csv .tsv",
    "slides": ".pptx .ppt .odp .key",
    "notebook": ".ipynb",
    "pdf": ".pdf",
    "text": ".txt .md .markdown .tex .bib .rst .org .log",
    "code": (".py .js .mjs .ts .tsx .jsx .java .c .h .cpp .hpp .cc .cs .go .rs .rb .php .r .m .jl "
             ".sh .ps1 .bat .sql .swift .kt .scala .lua .pl .html .htm .css .scss .vue .dart .f90 "
             ".yml .yaml .toml .ini .cfg"),
    "data": ".json .xml .parquet .feather .sav .dta .sas7bdat .mat .h5 .hdf5 .npy .npz .sqlite .db .rds",
    "image": ".png .jpg .jpeg .gif .webp .svg .heic .tif .tiff .bmp .ico .raw .cr2 .nef",
    "art": (".psd .ai .sketch .fig .xd .indd .afdesign .afphoto .kra .procreate .clip .blend .obj "
            ".fbx .stl .3mf .dwg .dxf .skp"),
    "audio": ".mp3 .wav .m4a .flac .aac .ogg .aiff .mid .midi",
    "video": ".mp4 .mov .avi .mkv .webm .m4v .wmv",
    "archive": ".zip .tar .gz .tgz .7z .rar .bz2 .xz",
    "email": ".eml .msg .mbox",
}
BY_SUFFIX = {suffix: kind for kind, suffixes in KINDS.items() for suffix in suffixes.split()}
NATIVE = {
    "application/vnd.google-apps.document": "document",
    "application/vnd.google-apps.spreadsheet": "spreadsheet",
    "application/vnd.google-apps.presentation": "slides",
}
LANGUAGES = {
    ".py": "Python", ".js": "JavaScript", ".mjs": "JavaScript", ".ts": "TypeScript", ".tsx": "TypeScript",
    ".jsx": "JavaScript", ".java": "Java", ".c": "C", ".h": "C", ".cpp": "C++", ".hpp": "C++",
    ".cc": "C++", ".cs": "C#", ".go": "Go", ".rs": "Rust", ".rb": "Ruby", ".php": "PHP", ".r": "R",
    ".m": "MATLAB", ".jl": "Julia", ".sh": "shell", ".ps1": "PowerShell", ".bat": "batch", ".sql": "SQL",
    ".swift": "Swift", ".kt": "Kotlin", ".scala": "Scala", ".lua": "Lua", ".pl": "Perl", ".html": "HTML",
    ".htm": "HTML", ".css": "CSS", ".scss": "SCSS", ".vue": "Vue", ".dart": "Dart", ".f90": "Fortran",
    ".yml": "YAML", ".yaml": "YAML", ".toml": "TOML", ".ini": "INI", ".cfg": "INI",
}
# Spaces and tabs only, never \s, which crosses lines and makes a file of
# blank lines take quadratic time.
DEFINITION = re.compile(r"^[ \t]*(?:export[ \t]+)?(?:async[ \t]+)?(?:def|class|function|fn|func|sub|struct|"
                        r"interface|public[ \t]+(?:static[ \t]+)?[\w<>\[\]]+[ \t]+\w+[ \t]*\()\b", re.MULTILINE)
SCAN = 200 * 1024       # characters of code the definition count reads
SHEET = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"
RELS = "{http://schemas.openxmlformats.org/package/2006/relationships}"
OFFICE_REL = "{http://schemas.openxmlformats.org/officeDocument/2006/relationships}"
DRAWING = "{http://schemas.openxmlformats.org/drawingml/2006/main}"
PRESENTATION = "{http://schemas.openxmlformats.org/presentationml/2006/main}"


def kind_of(name, mime=""):
    """The kind of a file from its name, or a Google file's type."""
    if mime in NATIVE:
        return NATIVE[mime]
    return BY_SUFFIX.get(Path(name).suffix.lower(), "other")


def digest(path, kind=None, name=None):
    """A file's digest: {"kind", "size", "text"}, the text at most LIMIT."""
    path = Path(path)
    name = name or path.name
    kind = kind or kind_of(name)
    size = path.stat().st_size
    reader = READERS.get(kind)
    try:
        text = reader(path, name) if reader else ""
    except Exception as problem:  # a file is data: any fault in it is said, never raised
        text = f"could not be read ({type(problem).__name__})"
    head = f"{kind}, {human(size)}"
    return {"kind": kind, "size": size, "text": clip(f"{head}\n{text}".strip())}


def clip(text, limit=LIMIT):
    return text if len(text) <= limit else text[:limit - 3].rstrip() + "..."


def human(size):
    for unit in ("bytes", "KB", "MB", "GB"):
        if size < 1024 or unit == "GB":
            return f"{size:.0f} {unit}" if unit == "bytes" else f"{size:.1f} {unit}"
        size /= 1024
    return f"{size} bytes"


def part(archive, name):
    with archive.open(name) as stream:
        return stream.read(PART_LIMIT)


def xml(archive, name):
    return ElementTree.fromstring(part(archive, name))


# Word

def read_document(path, _name):
    import docx_review
    with zipfile.ZipFile(path) as archive:
        # docx_review reads every part whole, so a part that inflates past
        # the limit is refused before it is read.
        if any(info.file_size > PART_LIMIT for info in archive.infolist()):
            raise ValueError("a part inflates past the limit")
    review = docx_review.extract(path)
    lines, words = [], 0
    for paragraph in review["paragraphs"]:
        text = paragraph["text"].strip()
        if not text:
            continue
        words += len(text.split())
        style = paragraph["style"] or ""
        lines.append(("# " if style.lower().startswith(("heading", "title")) else "") + text)
    open_comments = len(docx_review.open_comments(review))
    head = f"{len(review['paragraphs'])} paragraphs, {words} words, {open_comments} open comments"
    return head + "\n" + "\n".join(lines[:40])


# Excel

def read_spreadsheet(path, _name):
    with zipfile.ZipFile(path) as archive:
        names = set(archive.namelist())
        strings = shared_strings(archive) if "xl/sharedStrings.xml" in names else []
        out = []
        for title, target in sheets(archive):
            if target not in names:
                continue
            rows, formulas, extent = sheet_rows(archive, target, strings)
            out.append(f"## sheet {title}: {extent or 'empty'}, {formulas} formulas")
            out.append(table_summary(rows))
        return "\n".join(out) or "no sheets"


def shared_strings(archive):
    root = xml(archive, "xl/sharedStrings.xml")
    return ["".join(t.text or "" for t in item.iter(f"{SHEET}t")) for item in root.iter(f"{SHEET}si")]


def sheets(archive):
    """Each sheet's title and the zip part that holds it, in the workbook's order."""
    book = xml(archive, "xl/workbook.xml")
    rels = xml(archive, "xl/_rels/workbook.xml.rels")
    targets = {rel.get("Id"): rel.get("Target", "") for rel in rels.iter(f"{RELS}Relationship")}
    for sheet in book.iter(f"{SHEET}sheet"):
        target = targets.get(sheet.get(f"{OFFICE_REL}id"), "")
        target = target.lstrip("/") if target.startswith("/") else "xl/" + target
        yield sheet.get("name", "?"), target


def sheet_rows(archive, target, strings):
    """The first ROWS rows of a sheet as lists of values, the count of its
    formulas and the extent it declares."""
    root = xml(archive, target)
    dimension = root.find(f"{SHEET}dimension")
    rows, formulas = [], len(list(root.iter(f"{SHEET}f")))
    for row in root.iter(f"{SHEET}row"):
        if len(rows) >= ROWS:
            break
        values = {}
        for cell in row.iter(f"{SHEET}c"):
            column = column_index(cell.get("r", ""), len(values))
            values[column] = cell_value(cell, strings)
        if values:
            width = max(values) + 1
            rows.append([values.get(i, "") for i in range(min(width, 200))])
    return rows, formulas, dimension.get("ref") if dimension is not None else ""


def column_index(ref, fallback):
    letters = re.match(r"[A-Z]+", ref)
    if not letters:
        return fallback
    index = 0
    for letter in letters.group():
        index = index * 26 + ord(letter) - 64
    return index - 1


def cell_value(cell, strings):
    kind = cell.get("t")
    value = cell.find(f"{SHEET}v")
    if kind == "inlineStr":
        return "".join(t.text or "" for t in cell.iter(f"{SHEET}t"))
    if value is None or value.text is None:
        return ""
    if kind == "s":
        index = int(value.text) if value.text.isdigit() and len(value.text) < 10 else -1
        return strings[index] if 0 <= index < len(strings) else ""
    return value.text


# CSV and TSV

def read_table(path, name):
    raw = path.read_bytes()[:PART_LIMIT].decode("utf-8-sig", errors="replace")
    delimiter = "\t" if name.lower().endswith(".tsv") else sniff(raw)
    rows = []
    for row in csv.reader(io.StringIO(raw), delimiter=delimiter):
        if len(rows) >= ROWS:
            break
        rows.append(row)
    total = raw.count("\n")
    return f"about {total} rows\n" + table_summary(rows)


def sniff(raw):
    try:
        return csv.Sniffer().sniff(raw[:4096], delimiters=",;\t|").delimiter
    except csv.Error:
        return ","


def table_summary(rows):
    """The header, a few rows and each column's summary of a table read as
    rows, the first non-empty row taken for the header."""
    rows = [row for row in rows if any(str(v).strip() for v in row)]
    if not rows:
        return "no rows"
    header, body = rows[0], rows[1:]
    lines = [f"{len(body)} rows read, {len(header)} columns",
             "header: " + " | ".join(str(v)[:30] for v in header[:20])]
    lines += ["row: " + " | ".join(str(v)[:30] for v in row[:20]) for row in body[:SAMPLE_ROWS]]
    for index, title in enumerate(header[:30]):
        column = [row[index] for row in body if index < len(row) and str(row[index]).strip()]
        lines.append(f"- {str(title)[:40] or f'column {index + 1}'}: {column_summary(column)}")
    return "\n".join(lines)


def column_summary(values):
    numbers = [n for n in (number(v) for v in values) if n is not None]
    if values and len(numbers) >= 0.8 * len(values):
        return (f"{len(numbers)} numbers, min {short(min(numbers))}, max {short(max(numbers))}, "
                f"mean {short(statistics.fmean(numbers))}, sum {short(sum(numbers))}")
    distinct = sorted({str(v)[:24] for v in values})
    shown = ", ".join(distinct[:6]) + (", ..." if len(distinct) > 6 else "")
    return f"{len(values)} values, {len(distinct)} distinct: {shown}" if values else "empty"


def number(value):
    text = str(value).strip().replace(",", "").replace("$", "").rstrip("%")
    if not text or len(text) > 40:
        return None
    try:
        found = float(text)
    except ValueError:
        return None
    return found if found == found and abs(found) != float("inf") else None


def short(value):
    return f"{value:.4g}" if abs(value) < 1e6 else f"{value:.3e}"


# PowerPoint

def read_slides(path, _name):
    with zipfile.ZipFile(path) as archive:
        names = set(archive.namelist())
        order = slide_order(archive) if "ppt/presentation.xml" in names else []
        order = [p for p in order if p in names] or sorted(
            (n for n in names if re.fullmatch(r"ppt/slides/slide\d+\.xml", n)),
            key=lambda n: int(re.search(r"\d+", n).group()))
        out = [f"{len(order)} slides"]
        for number_, name in enumerate(order, 1):
            words = [t.text.strip() for t in xml(archive, name).iter(f"{DRAWING}t") if t.text and t.text.strip()]
            out.append(f"{number_}. " + " / ".join(words)[:200])
        return "\n".join(out)


def slide_order(archive):
    deck = xml(archive, "ppt/presentation.xml")
    rels = xml(archive, "ppt/_rels/presentation.xml.rels")
    targets = {rel.get("Id"): rel.get("Target", "") for rel in rels.iter(f"{RELS}Relationship")}
    return ["ppt/" + targets.get(s.get(f"{OFFICE_REL}id"), "").lstrip("/")
            for s in deck.iter(f"{PRESENTATION}sldId")]


# Notebooks, text, code and data

def read_notebook(path, _name):
    book = json.loads(path.read_bytes()[:PART_LIMIT])
    cells = book.get("cells", []) if isinstance(book, dict) else []
    language = (book.get("metadata", {}).get("kernelspec", {}) or {}).get("language", "?")
    count = {}
    for cell in cells:
        count[cell.get("cell_type", "?")] = count.get(cell.get("cell_type", "?"), 0) + 1
    first = next(("".join(c.get("source", [])) for c in cells if c.get("cell_type") == "markdown"), "")
    return f"{language}, " + ", ".join(f"{n} {k} cells" for k, n in sorted(count.items())) + "\n" + first[:600]


def read_text(path, name):
    text = path.read_bytes()[:PART_LIMIT].decode("utf-8", errors="replace")
    lines = text.splitlines()
    head = f"{len(lines)} lines, {len(text.split())} words"
    suffix = Path(name).suffix.lower()
    if suffix in LANGUAGES:
        head = f"{LANGUAGES[suffix]}, {len(lines)} lines, {len(DEFINITION.findall(text[:SCAN]))} functions or classes"
    return head + "\n" + "\n".join(line[:160] for line in lines[:30])


def read_data(path, name):
    if Path(name).suffix.lower() != ".json":
        return ""
    found = json.loads(path.read_bytes()[:PART_LIMIT])
    if isinstance(found, dict):
        return "an object with keys: " + ", ".join(list(found)[:30])
    if isinstance(found, list):
        first = found[0] if found else None
        keys = ", ".join(list(first)[:20]) if isinstance(first, dict) else type(first).__name__
        return f"a list of {len(found)}, the first: {keys}"
    return type(found).__name__


def read_pdf(path, _name):
    raw = path.read_bytes()[:PART_LIMIT]
    pages = len(re.findall(rb"/Type\s*/Page(?![a-zA-Z])", raw))
    title = re.search(rb"/Title\s*\(([^)]{1,200})\)", raw)
    named = title.group(1).decode("latin-1", errors="replace") if title else ""
    return f"{pages or 'unknown'} pages" + (f", title {named}" if named else "")


# Images

def read_image(path, name):
    head = path.read_bytes()[:256 * 1024]
    found = image_size(head, Path(name).suffix.lower())
    return f"{found[0]}, {found[1]} x {found[2]} pixels" if found else ""


def image_size(head, suffix):
    if head[:8] == b"\x89PNG\r\n\x1a\n" and len(head) >= 24:
        return ("PNG", *struct.unpack(">II", head[16:24]))
    if head[:6] in (b"GIF87a", b"GIF89a") and len(head) >= 10:
        return ("GIF", *struct.unpack("<HH", head[6:10]))
    if head[:4] == b"RIFF" and head[8:12] == b"WEBP":
        return webp_size(head)
    if head[:2] == b"\xff\xd8":
        return jpeg_size(head)
    if suffix == ".svg":
        text = head.decode("utf-8", errors="replace")
        width, height = (re.search(rf'\b{a}="([\d.]+)', text) for a in ("width", "height"))
        if width and height:
            return ("SVG", width.group(1), height.group(1))
        box = re.search(r'viewBox="[\d.\-]+[ ,]+[\d.\-]+[ ,]+([\d.]+)[ ,]+([\d.]+)', text)
        return ("SVG", box.group(1), box.group(2)) if box else None
    return None


def jpeg_size(head):
    at = 2
    while at + 9 < len(head):
        if head[at] != 0xFF:
            at += 1
            continue
        marker = head[at + 1]
        length = struct.unpack(">H", head[at + 2:at + 4])[0]
        if marker in (0xC0, 0xC1, 0xC2, 0xC3) and at + 9 <= len(head):
            height, width = struct.unpack(">HH", head[at + 5:at + 9])
            return ("JPEG", width, height)
        at += 2 + max(length, 2)
    return None


def webp_size(head):
    chunk = head[12:16]
    if chunk == b"VP8X" and len(head) >= 30:
        return ("WebP", 1 + int.from_bytes(head[24:27], "little"), 1 + int.from_bytes(head[27:30], "little"))
    if chunk == b"VP8 " and len(head) >= 30:
        width, height = struct.unpack("<HH", head[26:30])
        return ("WebP", width & 0x3FFF, height & 0x3FFF)
    if chunk == b"VP8L" and len(head) >= 25:
        bits = int.from_bytes(head[21:25], "little")
        return ("WebP", (bits & 0x3FFF) + 1, ((bits >> 14) & 0x3FFF) + 1)
    return None


READERS = {
    "document": lambda p, n: read_document(p, n) if n.lower().endswith(".docx") else "",
    "spreadsheet": lambda p, n: read_spreadsheet(p, n) if n.lower().endswith((".xlsx", ".xlsm")) else "",
    "table": read_table,
    "slides": lambda p, n: read_slides(p, n) if n.lower().endswith(".pptx") else "",
    "notebook": read_notebook,
    "text": read_text,
    "code": read_text,
    "data": read_data,
    "pdf": read_pdf,
    "image": read_image,
}


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--json", action="store_true")
    parser.add_argument("files", nargs="+")
    args = parser.parse_args(argv)
    found = []
    for name in args.files:
        if not Path(name).is_file():
            print(f"file_digest: no file {name}", file=sys.stderr)
            return 1
        found.append({"file": name, **digest(name)})
    if args.json:
        print(json.dumps(found, indent=1))
    else:
        print("\n\n".join(f"# {item['file']}\n{item['text']}" for item in found))
    return 0


if __name__ == "__main__":
    sys.exit(main())
