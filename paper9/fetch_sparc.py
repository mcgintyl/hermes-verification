#!/usr/bin/env python3
"""Fetch the SPARC rotation curves this package needs, and verify them.

    python fetch_sparc.py                 download the official archive and extract
    python fetch_sparc.py --from-dir DIR  copy from a SPARC folder you already have
    python fetch_sparc.py --check         only verify what is already in data/sparc

The SPARC data are not redistributed with this package. They are published by
Lelli, McGaugh & Schombert (2016), AJ 152, 157, at http://astroweb.case.edu/SPARC/,
whose terms ask that users cite that paper. This script downloads the public archive
Rotmod_LTG.zip, takes only the 133 files the paper scores, and checks every one of
them against data/sparc_manifest.csv, so you can be sure you have the same bytes the
published result used.
"""
import argparse
import hashlib
import shutil
import sys
import urllib.request
import zipfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))
from hermes_clusters import galaxies  # noqa: E402

URL = "https://astroweb.case.edu/SPARC/Rotmod_LTG.zip"
ARCHIVE_SHA256 = "0a80cc90714828cc28b7dd57923576714d209f2490328c087c4a4ad607faf588"
ARCHIVE_BYTES = 110737
DEST = ROOT / "data" / "sparc"


def verify(manifest, quiet=False):
    """Check every file in data/sparc against the manifest. Returns the list of problems."""
    problems = []
    for galaxy, (name, sha, size) in sorted(manifest.items()):
        path = DEST / name
        if not path.is_file():
            problems.append((galaxy, name, "missing"))
            continue
        blob = path.read_bytes()
        if len(blob) != size:
            problems.append((galaxy, name, "wrong size: %d, expected %d" % (len(blob), size)))
        elif hashlib.sha256(blob).hexdigest() != sha:
            problems.append((galaxy, name, "checksum mismatch"))
    if not quiet:
        print("%d of %d files present and verified" % (len(manifest) - len(problems), len(manifest)))
        for p in problems[:10]:
            print("  problem: %s (%s) %s" % p)
    return problems


def from_archive(blob, manifest):
    with zipfile.ZipFile(blob) as z:
        by_name = {n.split("/")[-1]: n for n in z.namelist() if n.endswith("_rotmod.dat")}
        for _, (name, _, _) in sorted(manifest.items()):
            if name not in by_name:
                raise SystemExit("the archive does not contain %s" % name)
            (DEST / name).write_bytes(z.read(by_name[name]))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--from-dir", metavar="DIR", help="copy the files from an existing SPARC folder")
    ap.add_argument("--check", action="store_true", help="verify what is already in data/sparc and stop")
    args = ap.parse_args()
    manifest = galaxies.load_manifest()
    DEST.mkdir(parents=True, exist_ok=True)

    if args.check:
        return 1 if verify(manifest) else 0

    if args.from_dir:
        source = Path(args.from_dir)
        for galaxy, (name, _, _) in sorted(manifest.items()):
            found = galaxies.rotmod_path(source, galaxy, manifest)
            if found is None:
                print("  not found in %s: %s (%s)" % (source, galaxy, name))
                continue
            shutil.copyfile(found, DEST / name)
        print("copied from %s" % source)
    else:
        print("downloading %s (%d bytes)" % (URL, ARCHIVE_BYTES))
        with urllib.request.urlopen(URL, timeout=120) as response:
            blob = response.read()
        digest = hashlib.sha256(blob).hexdigest()
        print("  %d bytes, SHA-256 %s" % (len(blob), digest))
        if digest != ARCHIVE_SHA256:
            print("  WARNING: this is not the archive this package was built against.")
            print("  Expected SHA-256 %s. The per-file checks below are what matter;" % ARCHIVE_SHA256)
            print("  if they pass, the data are the same even though the archive was repackaged.")
        archive = DEST.parent / "Rotmod_LTG.zip"
        archive.write_bytes(blob)
        from_archive(archive, manifest)
        archive.unlink()

    problems = verify(manifest)
    if problems:
        print("\n%d file(s) did not verify. Do not use this data for reproduction." % len(problems))
        return 1
    print("\nAll %d rotation curves verified. Now run: python reproduce_galaxies.py" % len(manifest))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
