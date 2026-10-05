"""Refresh file_manifest.csv (relative_path, bytes, sha256, description) for every file in paper10/.

Descriptions live in file_manifest.csv itself: this script keeps them, recomputes bytes and
SHA-256, drops rows for files that no longer exist, and stops with a list of any new file that
has no row yet (add a row with its description, then rerun). Generated folders
(reproduce_outputs/, regenerated/, data/rebuilt/, __pycache__/) and file_manifest.csv are skipped.

Usage:  python make_manifest.py [--check]     (--check: report differences, write nothing)

Line endings: the hashes are for the files as committed, with LF line endings. In a checkout whose
text files were converted to CRLF (Git core.autocrlf=true without a `paper10/** -text` attribute),
--check reports files whose LF-normalized bytes match as CRLF-only differences and still exits 0,
and a refresh refuses to run, so a manifest is never written from converted files.
"""
import argparse
import csv
import hashlib
import sys
from pathlib import Path

PKG = Path(__file__).resolve().parent
MANIFEST = PKG / 'file_manifest.csv'
SKIP_DIRS = {'reproduce_outputs', 'regenerated', '__pycache__', '.git'}
TEXT_SUFFIXES = {'.md', '.py', '.csv', '.txt', '.json', '.diff', '.gitignore'}


def sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for chunk in iter(lambda: f.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()


def sha256_lf(path):
    """SHA-256 of the file with CRLF converted to LF (text files only), or None."""
    if Path(path).suffix not in TEXT_SUFFIXES and Path(path).name != '.gitignore':
        return None
    raw = Path(path).read_bytes()
    return hashlib.sha256(raw.replace(b'\r\n', b'\n')).hexdigest() if b'\r\n' in raw else None


def package_files():
    out = []
    for p in sorted(PKG.rglob('*')):
        rel = p.relative_to(PKG).as_posix()
        if not p.is_file() or p == MANIFEST or set(p.relative_to(PKG).parts[:-1]) & SKIP_DIRS or rel.startswith('data/rebuilt/'):
            continue
        out.append(rel)
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--check', action='store_true')
    args = ap.parse_args(argv)
    with open(MANIFEST, encoding='utf-8', newline='') as f:
        rows = list(csv.DictReader(f))
    desc = {r['relative_path']: r['description'] for r in rows}
    files = package_files()
    missing = [f for f in files if not desc.get(f)]
    if missing:
        sys.exit('no description in file_manifest.csv for:\n  ' + '\n  '.join(missing))
    order = [r['relative_path'] for r in rows if r['relative_path'] in files]
    new_rows, changed, crlf_only, crlf_any = [], [], [], []
    for rel in order:
        p = PKG / rel
        row = {'relative_path': rel, 'bytes': str(p.stat().st_size), 'sha256': sha256(p), 'description': desc[rel]}
        old = next(r for r in rows if r['relative_path'] == rel)
        lf = sha256_lf(p)
        if lf is not None:
            crlf_any.append(rel)
        if (old['bytes'], old['sha256']) != (row['bytes'], row['sha256']):
            (crlf_only if lf == old['sha256'] else changed).append(rel)
        new_rows.append(row)
    dropped = [r['relative_path'] for r in rows if r['relative_path'] not in files]
    print(f'{len(new_rows)} files; changed: {changed or "none"}; dropped: {dropped or "none"}')
    if crlf_only:
        print(f'CRLF-only differences (LF-normalized bytes match the manifest): {len(crlf_only)} files; '
              'mark paper10/** -text in .gitattributes or check out with core.autocrlf=false')
    if args.check:
        sys.exit(1 if (changed or dropped) else 0)
    if crlf_any:
        sys.exit('refusing to refresh the manifest: these files have CRLF line endings (a converted checkout?):\n  '
                 + '\n  '.join(crlf_any))
    with open(MANIFEST, 'w', encoding='utf-8', newline='') as f:
        w = csv.DictWriter(f, fieldnames=['relative_path', 'bytes', 'sha256', 'description'], lineterminator='\n')
        w.writeheader()
        w.writerows(new_rows)
    print('wrote', MANIFEST)


if __name__ == '__main__':
    main()
