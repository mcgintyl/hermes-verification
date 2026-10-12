"""Refresh paper7/file_manifest.csv and paper7/MANIFEST_SHA256.txt.

file_manifest.csv (relative_path, bytes, sha256, description) lists every file in paper7/ except
itself, MANIFEST_SHA256.txt, results/, __pycache__/ and operating-system or editor files (.DS_Store,
Thumbs.db, desktop.ini, *.pyc, *.pyo, *.swp, *.swo; reproduce.py skips the same ones). Descriptions
live in file_manifest.csv: this script keeps them, recomputes bytes and SHA-256, drops rows for files
that no longer exist, and stops with a list of any new file that has no row yet (add a row
"path,,,description", then rerun).
MANIFEST_SHA256.txt ("sha256  path", sorted by path) is then written for the same files plus
file_manifest.csv. reproduce.py checks both (checks H.1 to H.3).

Usage:  python paper7/make_manifest.py [--check]     (--check: report differences, write nothing)

Line endings: the hashes are for the files as committed, with LF line endings (.gitattributes marks
paper7/** -text). A refresh refuses to run if any text file contains CRLF, so a manifest is never written
from a checkout whose text files were converted.
"""
import argparse
import csv
import hashlib
import sys
from pathlib import Path

PKG = Path(__file__).resolve().parent
MANIFEST = PKG / "file_manifest.csv"
SHA_LIST = PKG / "MANIFEST_SHA256.txt"
SKIP_DIRS = {"results", "__pycache__", ".git"}
TEXT_SUFFIXES = {".md", ".py", ".csv", ".txt", ".json"}
JUNK_NAMES = {".DS_Store", "Thumbs.db", "desktop.ini"}        # the same lists as reproduce.py
JUNK_SUFFIXES = {".pyc", ".pyo", ".swp", ".swo"}


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def package_files():
    out = []
    for p in sorted(PKG.rglob("*")):
        parts = p.relative_to(PKG).parts
        if (not p.is_file() or p in (MANIFEST, SHA_LIST) or set(parts[:-1]) & SKIP_DIRS
                or p.name in JUNK_NAMES or p.suffix in JUNK_SUFFIXES):
            continue
        out.append(p.relative_to(PKG).as_posix())
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--check", action="store_true")
    args = ap.parse_args(argv)
    with open(MANIFEST, encoding="utf-8", newline="") as f:
        rows = list(csv.DictReader(f))
    desc = {r["relative_path"]: r["description"] for r in rows}
    files = package_files()
    missing = [f for f in files if not desc.get(f)]
    if missing:
        sys.exit("no description in file_manifest.csv for:\n  " + "\n  ".join(missing))
    new_rows, changed = [], []
    for r in rows:
        rel = r["relative_path"]
        if rel not in files:
            continue
        p = PKG / rel
        row = {"relative_path": rel, "bytes": str(p.stat().st_size), "sha256": sha256(p), "description": desc[rel]}
        if (r["bytes"], r["sha256"]) != (row["bytes"], row["sha256"]):
            changed.append(rel)
        new_rows.append(row)
    dropped = [r["relative_path"] for r in rows if r["relative_path"] not in files]
    crlf = [rel for rel in files if ((PKG / rel).suffix in TEXT_SUFFIXES or rel.endswith(".gitignore"))
            and b"\r\n" in (PKG / rel).read_bytes()]
    short = lambda xs: "none" if not xs else "%d (%s%s)" % (len(xs), ", ".join(xs[:8]), ", ..." if len(xs) > 8 else "")
    print("%d files; changed: %s; dropped: %s" % (len(new_rows), short(changed), short(dropped)))
    if args.check:
        listed = {}
        if SHA_LIST.is_file():
            for line in SHA_LIST.read_text(encoding="utf-8").splitlines():
                if line.strip():
                    h, rel = line.split("  ", 1)
                    listed[rel] = h
        want = {r["relative_path"]: r["sha256"] for r in new_rows}
        want["file_manifest.csv"] = sha256(MANIFEST)
        stale = sorted(k for k in set(want) | set(listed) if want.get(k) != listed.get(k))
        print("MANIFEST_SHA256.txt: %s" % ("up to date" if not stale else "differs for " + short(stale)))
        sys.exit(1 if (changed or dropped or stale) else 0)
    if crlf:
        sys.exit("refusing to refresh the manifests: these files contain CRLF (a converted checkout?):\n  "
                 + "\n  ".join(crlf))
    with open(MANIFEST, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["relative_path", "bytes", "sha256", "description"], lineterminator="\n")
        w.writeheader()
        w.writerows(new_rows)
    lines = {r["relative_path"]: r["sha256"] for r in new_rows}
    lines["file_manifest.csv"] = sha256(MANIFEST)
    with open(SHA_LIST, "w", encoding="utf-8", newline="\n") as f:
        f.writelines("%s  %s\n" % (lines[k], k) for k in sorted(lines))
    print("wrote", MANIFEST.name, "and", SHA_LIST.name)


if __name__ == "__main__":
    main()
