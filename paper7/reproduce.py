#!/usr/bin/env python3
"""Rebuild the Paper 7 M33 numbers from public inputs and check them.

    python paper7/reproduce.py

Paper 7: "Field Audit: Board-Complete Data and the Limits of External Gravity-Model
Testing" (McGinty 2026). Zenodo, all versions: 10.5281/zenodo.20550091. The checks
follow the version 4 text: Tables 2, 2b and 3, Appendix A, and every M33 number in the
prose. Nothing is downloaded; every input is in this repository:

  paper7/00_final_board/m33_corbelli2014_table1_import.csv   Corbelli et al. (2014) Table 1
  paper7/00_final_board/m33_corbelli_reconstruction_kernel.py the frozen thick-disk kernel
  paper7/00_final_board/m33_corbelli_reconstructed_board.csv  the board (58 radii)
  paper7/00_final_board/m33_age_bracketed_amber_lock.csv      t50 = 6.0 +/- 1.2 and 7.1 +/- 1.3 Gyr
  paper7/01b_hermes_result_canonical_gate/                   the canonical-gate predictions of v3
  paper1/hermes_gate_phi.py                                  the log-grid gate (canonical)
  paper9/hermes_clusters/gate.py                             the chain-rule gate of Paper 1's scores
  paper8/verify_eta.py                                       the Paper 8 operator (Section 6 only)
  paper1/Hermes_ConfigG_PerGalaxy_133_Export.csv and
  paper8/expected/paper8_sparc_133_model_comparison.csv      the SPARC medians quoted in Section 4
  paper7/expected/                                           the tables v4 prints, at full precision

Prints one line per check with the place the number appears in the paper, writes
results/ (the rebuilt board, the regenerated tables and checks.csv), compares the tables
with expected/, and exits 1 on any FAIL.

Conventions (each one changes a checked value if altered):
  board     g_bar = (Vdisk^2 + Vbulge^2 + sign(Vgas) Vgas^2) / R; g98 = 98th percentile of g_bar
  gate      log-grid shear, paper1/hermes_gate_phi.py (the chain-rule gate is checked separately)
  Hermes    beta = pi exp(-t50 g98 / 46654) - 1/sqrt(2 pi); V = sqrt(max(g_bar (1 + beta phi), 0) R)
  MOND      simple interpolating function, a0 = 1.2e-10 m/s^2 with the exact parsec
            (3702.813 (km/s)^2/kpc)
  chi2/N    mean over the N points of (Vobs - V)^2 / (errV^2 + 386); N is not reduced
  residual  V_model - Vobs
  Table 3   Vdisk multiplied by the scale; g98, beta and phi recomputed at each scale; the
            source-native disk has M*/L(3.6 um) = 4.9e9 / 10^(0.4 (24.60 + 3.24 - 3.2)) = 0.6826;
            sweep 0.5 to 3.0 in steps of 0.0025, then 0.0001 within 0.005 of each age's minimum;
            each printed row is evaluated at its printed scale (the minimum row at 1.368, the
            3-dp rounding of the fine-grid minimum 1.3677)
Printed numbers are checked to half a unit of their last printed digit (ties allowed).

Needs Python 3.10 or later (the frozen kernel's annotations use the X | Y syntax), NumPy, SciPy
and pandas.

Not rebuilt: the Section 5 numbers (Repetto case, SPARC-overlap THINGS and dwarf runs), whose
figure-extracted inputs are not in this repository; version 4 quotes them from the original runs.
"""
import argparse, ast, csv, hashlib, json, os, re, sys

sys.dont_write_bytecode = True        # keep __pycache__ out of the hashed package folders
import numpy as np                    # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
BOARD_DIR = os.path.join(HERE, "00_final_board")
CANON_DIR = os.path.join(HERE, "01b_hermes_result_canonical_gate")
EXP = os.path.join(HERE, "expected")
OUT = os.path.join(HERE, "results")
sys.path[0:0] = [BOARD_DIR, os.path.join(ROOT, "paper1"), os.path.join(ROOT, "paper9"), os.path.join(ROOT, "paper8")]
import m33_corbelli_reconstruction_kernel as kernel       # noqa: E402  (frozen kernel)
from hermes_gate_phi import hermes_phi                    # noqa: E402  (log-grid)
from hermes_clusters.gate import phi_chain                # noqa: E402  (chain-rule)
import verify_eta as ve                                   # noqa: E402  (Paper 8 operator)

PI = np.pi
F = 1.0 / np.sqrt(2.0 * PI)
SIG2 = 386.0                                   # (km/s)^2, McGinty (2026a) Appendix A
K_DIV = 46654.0                                # c/(2 pi) in Gyr (km/s)^2/kpc
KPC_M = 3.0856775814913673e19                  # exact parsec (IAU 2015 B2) x 1000
A0 = 1.2e-10 * KPC_M / 1.0e6                   # 3702.813 (km/s)^2/kpc
AGES = (6.0, 7.1)
# Source-native M*/L at 3.6 um (03_baryonic_convention_sensitivity report): M* = 4.9e9 Msun
# (Corbelli 2014), total [3.6] = 3.2 mag, distance modulus 24.60, solar 3.6 um magnitude 3.24.
ML0 = 4.9e9 / 10 ** (0.4 * (24.60 + 3.24 - 3.2))
UNIT_ML_SCALE = 1.210                          # SPARC-like unit-M/L analogue, sqrt(1/0.683)

IMPORT_CSV = os.path.join(BOARD_DIR, "m33_corbelli2014_table1_import.csv")
BOARD_CSV = os.path.join(BOARD_DIR, "m33_corbelli_reconstructed_board.csv")
KERNEL_PY = os.path.join(BOARD_DIR, "m33_corbelli_reconstruction_kernel.py")
AGE_CSV = os.path.join(BOARD_DIR, "m33_age_bracketed_amber_lock.csv")
CANON_PP = os.path.join(CANON_DIR, "m33_hermes_per_point_predictions_canonical.csv")
CANON_PY = os.path.join(CANON_DIR, "compute_m33_hermes_canonical.py")
CANON_MASKS = os.path.join(CANON_DIR, "m33_hermes_mask_summary_canonical.csv")
SUPERSEDED = os.path.join(HERE, "01_hermes_result_locked", "m33_hermes_age_bracket_results.csv")
BOARD_COPY = os.path.join(HERE, "01_hermes_result_locked", "m33_corbelli_reconstructed_board_input.csv")
MOND_CFG = os.path.join(HERE, "02_mond_comparator", "m33_mond_comparator_run_config.json")
SWEEP_03 = os.path.join(HERE, "03_baryonic_convention_sensitivity", "m33_vdisk_scaling_series_both_ages.csv")
PAPER_TEXT = os.path.join(HERE, "paper7_main.md")
P1_EXPORT = os.path.join(ROOT, "paper1", "Hermes_ConfigG_PerGalaxy_133_Export.csv")
P8_SPARC = os.path.join(ROOT, "paper8", "expected", "paper8_sparc_133_model_comparison.csv")

# SHA-256 of the inputs other papers and scripts read (LF bytes, as committed)
SHA = {
    IMPORT_CSV: "6dba81de462849123497a6e38207259b90d492b17b808e290af8392321986cfb",
    BOARD_CSV: "ecd29c4dec518008ea80952c6411bdb956d306312de78ed58d8fcb73a9cc089a",
    KERNEL_PY: "75c9ca0b0c6674283db3429144a50be688b975be8d335353ce32ed7ebf8740c1",
    AGE_CSV: "f408e30b863e3896a454e649520373ee1799bf037a8f18a8d751b56e8e8e50b0",
    CANON_PP: "42b2b94e190f8f361162442d40c7a13652dc04e120be6262542a489f98c4eef3",
}

# ------------------------------------------------------------------ the printed tables (v4)
MASKS = [  # short name, Table 2 label, condition
    ("full", "Full board", lambda R: R > 0),
    ("R<=8", "Main disk, R <= 8 kpc", lambda R: R <= 8),
    ("R>8", "Warp zone, R > 8 kpc", lambda R: R > 8),
    ("R<=10", "Main disk, R <= 10 kpc", lambda R: R <= 10),
    ("R>10", "Outer tail, R > 10 kpc", lambda R: R > 10),
    ("R>=9.1", "Barker break outer, R >= 9.1 kpc", lambda R: R >= 9.1),
    ("R>=11.6", "Barker S2 far outer, R >= 11.6 kpc", lambda R: R >= 11.6),
    ("R>=15", "Putman tidal-radius exposure, R >= 15 kpc", lambda R: R >= 15),
]
MODELS = [("H6", "Hermes primary t50=6.0"), ("H7", "Hermes conservative t50=7.1"), ("MOND", "MOND fixed simple")]
# Table 2: mask, model, N, R min, R max, chi2/N, chi2 fraction (%), mean, median |res|, RMS residual
TABLE2 = [
    ("full", "H6", 58, 0.24, 22.72, 3.22, 100.00, -28.02, 32.00, 38.99),
    ("full", "H7", 58, 0.24, 22.72, 3.39, 100.00, -29.49, 33.35, 39.91),
    ("full", "MOND", 58, 0.24, 22.72, 0.12, 100.00, -4.12, 3.69, 7.67),
    ("R<=8", "H6", 27, 0.24, 7.57, 0.78, 11.24, -4.18, 10.55, 17.46),
    ("R<=8", "H7", 27, 0.24, 7.57, 0.78, 10.75, -5.71, 10.37, 17.49),
    ("R<=8", "MOND", 27, 0.24, 7.57, 0.02, 9.31, 0.15, 2.58, 3.08),
    ("R>8", "H6", 31, 8.06, 22.72, 5.35, 88.76, -48.78, 49.61, 50.79),
    ("R>8", "H7", 31, 8.06, 22.72, 5.65, 89.25, -50.20, 51.03, 52.10),
    ("R>8", "MOND", 31, 8.06, 22.72, 0.20, 90.69, -7.84, 7.81, 10.09),
    ("R<=10", "H6", 31, 0.24, 9.53, 0.93, 15.41, -7.21, 14.61, 19.09),
    ("R<=10", "H7", 31, 0.24, 9.53, 0.96, 15.14, -8.73, 13.03, 19.39),
    ("R<=10", "MOND", 31, 0.24, 9.53, 0.02, 10.66, 0.47, 2.58, 3.07),
    ("R>10", "H6", 27, 10.02, 22.72, 5.86, 84.59, -51.92, 57.00, 53.37),
    ("R>10", "H7", 27, 10.02, 22.72, 6.17, 84.86, -53.33, 58.26, 54.68),
    ("R>10", "MOND", 27, 10.02, 22.72, 0.23, 89.34, -9.39, 8.55, 10.75),
    ("R>=9.1", "H6", 28, 9.53, 22.72, 5.71, 85.50, -50.99, 54.71, 52.63),
    ("R>=9.1", "H7", 28, 9.53, 22.72, 6.02, 85.84, -52.41, 56.04, 53.95),
    ("R>=9.1", "MOND", 28, 9.53, 22.72, 0.22, 89.34, -9.04, 8.25, 10.55),
    ("R>=11.6", "H6", 23, 11.97, 22.72, 6.39, 78.61, -55.21, 57.67, 56.15),
    ("R>=11.6", "H7", 23, 11.97, 22.72, 6.70, 78.51, -56.54, 59.05, 57.43),
    ("R>=11.6", "MOND", 23, 11.97, 22.72, 0.25, 84.93, -10.11, 8.90, 11.42),
    ("R>=15", "H6", 16, 15.39, 22.72, 7.16, 61.31, -60.83, 59.48, 61.13),
    ("R>=15", "H7", 16, 15.39, 22.72, 7.46, 60.80, -62.05, 60.78, 62.35),
    ("R>=15", "MOND", 16, 15.39, 22.72, 0.32, 74.83, -11.96, 13.22, 13.05),
]
# Table 2b: sigma_int^2; Hermes (6.0 Gyr) full, R <= 10, R > 10, share of chi2 at R > 10 (%); MOND full, R > 10
TABLE2B = [
    (0, "118.48", "87.49", "154.06", "60.5", "2.88", "4.14"),
    (100, "9.52", "3.44", "16.51", "80.7", "0.33", "0.61"),
    (386, "3.22", "0.93", "5.86", "84.6", "0.12", "0.23"),
    (900, "1.53", "0.40", "2.82", "85.9", "0.06", "0.11"),
]
TABLE2B_COLS = ["Hermes full", "Hermes R <= 10", "Hermes R > 10", "Hermes share R > 10 (%)", "MOND full", "MOND R > 10"]
# Table 3: Vdisk scale, M/L factor, global M*/L (3.6 um), t50; Hermes chi2/N, outer chi2/N, outer mean;
# MOND chi2/N, outer chi2/N, outer mean; note.  "Outer" is R > 10 kpc.
TABLE3 = [
    ("1.000", "1.00", "0.683", 6.0, "3.22", "5.86", "-51.92", "0.12", "0.23", "-9.39", "Baseline/source-native"),
    ("1.000", "1.00", "0.683", 7.1, "3.39", "6.17", "-53.33", "0.12", "0.23", "-9.39", "Conservative age bracket"),
    ("1.368", "1.87", "1.278", 6.0, "2.41", "4.00", "-42.57", "0.53", "0.13", "5.13",
     "Lowest Hermes chi2nu on a continuous scaling (not adopted)"),
    ("1.500", "2.25", "1.536", 6.0, "2.53", "3.56", "-40.05", "0.99", "0.32", "10.36", "Discrete sweep"),
    ("1.500", "2.25", "1.536", 7.1, "2.52", "4.06", "-43.01", "0.99", "0.32", "10.36", "Discrete sweep"),
    ("2.000", "4.00", "2.731", 6.0, "3.04", "2.38", "-32.39", "4.27", "2.14", "29.95", "Discrete sweep"),
    ("2.000", "4.00", "2.731", 7.1, "2.70", "2.89", "-36.01", "4.27", "2.14", "29.95", "Discrete sweep"),
    ("2.500", "6.25", "4.267", 6.0, "3.61", "1.43", "-23.60", "9.92", "5.57", "49.01", "Discrete sweep"),
    ("2.500", "6.25", "4.267", 7.1, "2.97", "1.75", "-26.69", "9.92", "5.57", "49.01", "Discrete sweep"),
    ("3.000", "9.00", "6.144", 6.0, "4.68", "0.73", "-11.95", "17.93", "10.49", "67.60", "Discrete sweep"),
    ("3.000", "9.00", "6.144", 7.1, "4.03", "0.82", "-13.86", "17.93", "10.49", "67.60", "Discrete sweep"),
]
SWEEP_COARSE = (5000, 30000, 25)      # Vdisk scale 0.5000 to 3.0000 in steps of 0.0025 (units of 1e-4)
SWEEP_FINE_HALF = 50                  # then steps of 0.0001 within +/- 0.0050 of each age's coarse minimum

# Operating-system and editor files that are never part of the package (the root .gitignore keeps
# them out of git); make_manifest.py skips the same ones.
JUNK_NAMES = {".DS_Store", "Thumbs.db", "desktop.ini"}
JUNK_SUFFIXES = (".pyc", ".pyo", ".swp", ".swo")

RESULTS, PROSE, INFO = [], [], []


# ------------------------------------------------------------------ helpers
def check(cid, where, what, got, want, tol):
    ok = bool(abs(float(got) - float(want)) <= tol * (1 + 1e-6))
    RESULTS.append((cid, where, what, float(got), float(want), tol, ok))
    print("  [%s] %-7s %-24s %-58s got %-12.7g want %-10.7g tol %g"
          % ("PASS" if ok else "FAIL", cid, where, what, float(got), float(want), tol))


def check_true(cid, where, what, cond, detail=""):
    RESULTS.append((cid, where, what, float(bool(cond)), 1.0, 0, bool(cond)))
    print("  [%s] %-7s %-24s %s %s" % ("PASS" if cond else "FAIL", cid, where, what, detail))


def printed(cid, where, what, got, shown, tol):
    """A number the paper prints: checked to `tol` and listed in prose_numbers.csv."""
    check(cid, where, what, got, float(shown), tol)
    PROSE.append((cid, where, what, float(got), shown, tol))


def info(cid, where, what, got):
    """A value printed for information only: no pass/fail, not counted as a check."""
    INFO.append((cid, where, what, float(got), "", "", "info"))
    print("  [INFO] %-7s %-24s %-58s got %-12.7g" % (cid, where, what, float(got)))


def emit_check(kind, cid, where, what, got, target, tol=0, detail=""):
    """Route a check from prose_checks() to printed(), check() or check_true()."""
    if kind == "printed":
        printed(cid, where, what, got, target, tol)
    elif kind == "check":
        check(cid, where, what, got, target, tol)
    else:
        check_true(cid, where, what, got, detail)


def sha256(path):
    with open(path, "rb") as fh:
        return hashlib.sha256(fh.read()).hexdigest()


def read_csv(path):
    with open(path, encoding="utf-8-sig", newline="") as fh:
        return list(csv.DictReader(fh))


def num(x, fmt="%.10f"):
    return fmt % (float(x) + 0.0)     # + 0.0 turns -0.0 into 0.0


def write_csv(name, header, rows):
    with open(os.path.join(OUT, name), "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, lineterminator="\n")
        w.writerow(header)
        w.writerows(rows)
    return name


def cells_match(a, b, tol=1e-9):
    """Compare two CSV row lists cell by cell: numbers to tol (relative or absolute), text exactly."""
    if len(a) != len(b):
        return False, "row count %d vs %d" % (len(a), len(b))
    worst = 0.0
    for ra, rb in zip(a, b):
        if list(ra) != list(rb):
            return False, "columns differ"
        for k in ra:
            try:
                x, y = float(ra[k]), float(rb[k])
                d = abs(x - y) / max(1.0, abs(y))
                worst = max(worst, d)
                if d > tol:
                    return False, "%s: %r vs %r" % (k, ra[k], rb[k])
            except ValueError:
                if ra[k] != rb[k]:
                    return False, "%s: %r vs %r" % (k, ra[k], rb[k])
    return True, "max rel diff %.1e" % worst


# ------------------------------------------------------------------ models
def gbar(B, scale=1.0, vgas=None, vdisk=None):
    Vg = B["Vg"] if vgas is None else vgas
    Vd = (B["Vd"] if vdisk is None else vdisk) * scale
    return (Vd ** 2 + B["Vb"] ** 2 + np.sign(Vg) * Vg ** 2) / B["R"]


def beta(t50, g98):
    return PI * np.exp(-t50 * g98 / K_DIV) - F


def v_hermes(R, gb, t50, gate=hermes_phi, phi=None):
    phi = gate(R, gb) if phi is None else phi
    g = gb * (1.0 + beta(t50, np.percentile(gb, 98)) * phi)
    return np.sqrt(np.maximum(g, 0.0) * R)


def v_mond(R, gb):
    x = np.maximum(np.abs(gb) / A0, 1e-30)
    return np.sqrt(np.maximum(gb * (1.0 + np.sqrt(1.0 + 4.0 / x)) / 2.0, 0.0) * R)


def stats(B, V, m, sig2=SIG2):
    res = V - B["Vo"]
    q = res ** 2 / (B["e"] ** 2 + sig2)
    return dict(N=int(m.sum()), Rmin=float(B["R"][m].min()), Rmax=float(B["R"][m].max()),
                chi2nu=float(q[m].mean()), frac=float(100.0 * q[m].sum() / q.sum()),
                mean=float(res[m].mean()), med=float(np.median(np.abs(res[m]))),
                rms=float(np.sqrt(np.mean(res[m] ** 2))))


def load_board(path=BOARD_CSV):
    rows = read_csv(path)
    col = lambda n: np.array([float(r[n]) for r in rows])
    return dict(rows=rows, R=col("R_kpc"), Vo=col("Vobs_kms"), e=col("errV_kms"), Vg=col("Vgas_kms"),
                Vd=col("Vdisk_kms"), Vb=col("Vbulge_kms"), Vbar=col("Vbar_kms"), col=col)


def sweep_row(B, n):
    """Everything Table 3 needs at Vdisk scale n / 10000."""
    s = n / 10000.0
    R, gb = B["R"], gbar(B, s)
    phi, g98 = hermes_phi(R, gb), float(np.percentile(gb, 98))
    out = dict(n=n, scale=s, g98=g98)
    for t in AGES:
        V = v_hermes(R, gb, t, phi=phi)
        out[t] = (stats(B, V, R > 0)["chi2nu"], stats(B, V, R <= 10)["chi2nu"],
                  stats(B, V, R > 10)["chi2nu"], stats(B, V, R > 10)["mean"])
    VM = v_mond(R, gb)
    out["mond"] = (stats(B, VM, R > 0)["chi2nu"], stats(B, VM, R > 10)["chi2nu"], stats(B, VM, R > 10)["mean"])
    return out


# ------------------------------------------------------------------ manifests
def manifest_files():
    """Every file in paper7/ except MANIFEST_SHA256.txt, results/, __pycache__/ and OS/editor junk."""
    out = []
    for dp, dn, fns in os.walk(HERE):
        rel_dir = os.path.relpath(dp, HERE).replace(os.sep, "/")
        dn[:] = [d for d in dn if d not in ("__pycache__", "results", ".git")]
        for fn in fns:
            rel = fn if rel_dir == "." else rel_dir + "/" + fn
            if rel != "MANIFEST_SHA256.txt" and fn not in JUNK_NAMES and not fn.endswith(JUNK_SUFFIXES):
                out.append(rel)
    return sorted(out)


def read_manifest():
    out, path = {}, os.path.join(HERE, "MANIFEST_SHA256.txt")
    if not os.path.isfile(path):
        return out
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            if line.strip() and not line.startswith("#"):
                h, rel = line.rstrip("\n").split("  ", 1)
                out[rel] = h
    return out


# ------------------------------------------------------------------ paper text
def appendix_a_rows(path):
    """The rows of the Appendix A table in a paper text, as lists of cell strings."""
    rows, inside = [], False
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            if re.match(r"^#+\s*Appendix A", line):
                inside = True
                continue
            if inside and re.match(r"^#+\s", line):
                break
            if inside and line.startswith("|"):
                cells = [c.strip() for c in line.strip().strip("|").split("|")]
                if re.match(r"^-?\d+(\.\d+)?$", cells[0]):
                    rows.append(cells)
    return rows


def decimals(cell):
    return len(cell.split(".")[1]) if "." in cell else 0


def is_v4_text(path):
    """True when the shipped paper text is version 4 (its version note opens "Version 4")."""
    with open(path, encoding="utf-8") as fh:
        return re.search(r"\bVersion 4\b", fh.read()) is not None


def module_constant(path, name):
    """The value of a top-level `NAME = <literal>` assignment in a script, read without running it."""
    with open(path, encoding="utf-8") as fh:
        tree = ast.parse(fh.read())
    for node in tree.body:
        if isinstance(node, ast.Assign) and any(isinstance(t, ast.Name) and t.id == name for t in node.targets):
            return ast.literal_eval(node.value)
    raise KeyError(name)


# ------------------------------------------------------------------ prose (Abstract, Sections 4 and 6)
def prose_checks(B, emit):
    """The M33 numbers of the Abstract and Sections 4 and 6 on board B (checks P.1 to P.40).

    emit(kind, id, where, what, value, target, tol, detail): kind "printed" (target is the printed
    string), "check" (a support value) or "true" (value is a condition). Section 5 runs this on the
    board; Section 7 runs it again on the 3-dp Appendix A table."""
    R, Vo, e = B["R"], B["Vo"], B["e"]
    gb = gbar(B)
    phi = hermes_phi(R, gb)
    g98 = float(np.percentile(gb, 98))
    V = {"H6": v_hermes(R, gb, 6.0, phi=phi), "H7": v_hermes(R, gb, 7.1, phi=phi), "MOND": v_mond(R, gb)}
    S = {k: {s: stats(B, V[k], f(R)) for s, _, f in MASKS} for k in V}
    emit("printed", "P.1", "Abstract", "Hermes chi2/N (approx. 3.2)", S["H6"]["full"]["chi2nu"], "3.2", 0.05)
    emit("printed", "P.2", "Abstract", "MOND chi2/N (approx. 0.1)", S["MOND"]["full"]["chi2nu"], "0.1", 0.05)
    emit("printed", "P.3", "Sec 4 para 3", "MOND chi2/N", S["MOND"]["full"]["chi2nu"], "0.12", 0.005)
    emit("printed", "P.4", "Sec 4 para 3", "MOND median |residual| (km/s)", S["MOND"]["full"]["med"], "3.7", 0.05)
    q0 = lambda k, f: stats(B, V[k], f(R), 0.0)
    emit("printed", "P.5", "Sec 4 para 3", "MOND chi2/N without the floor", q0("MOND", lambda r: r > 0)["chi2nu"], "2.88", 0.005)
    emit("printed", "P.6", "Sec 4 para 3", "Hermes chi2/N without the floor", q0("H6", lambda r: r > 0)["chi2nu"], "118", 0.5)
    emit("printed", "P.7", "Sec 4 para 4", "Hermes chi2/N, primary age", S["H6"]["full"]["chi2nu"], "3.22", 0.005)
    emit("printed", "P.8", "Sec 4 para 4", "Hermes chi2/N, conservative age", S["H7"]["full"]["chi2nu"], "3.39", 0.005)
    emit("printed", "P.9", "Sec 4 para 4", "Hermes median |residual| (approx. 32 km/s)", S["H6"]["full"]["med"], "32", 0.5)
    emit("printed", "P.10", "Sec 4 para 6", "Hermes R <= 10 chi2/N", S["H6"]["R<=10"]["chi2nu"], "0.93", 0.005)
    emit("printed", "P.11", "Sec 4 para 6", "Hermes R <= 10 mean residual (km/s)", S["H6"]["R<=10"]["mean"], "-7.2", 0.05)
    emit("printed", "P.12", "Sec 4 para 6", "Hermes R > 10 chi2/N", S["H6"]["R>10"]["chi2nu"], "5.86", 0.005)
    emit("printed", "P.13", "Sec 4 para 6", "Hermes R > 10 mean residual (km/s)", S["H6"]["R>10"]["mean"], "-51.9", 0.05)
    emit("printed", "P.14", "Sec 4 para 6", "share of Hermes chi2 at R > 10 (approx. 85%)", S["H6"]["R>10"]["frac"], "85", 0.5)
    res = V["H6"] - Vo
    hi = (R <= 10) & (res > 10.0)
    lo = (R <= 10) & (res < -10.0)
    emit("true", "P.15", "Sec 4 para 6", "main-disk points > 10 km/s high are the %d at 0.46-1.45 kpc" % hi.sum(),
         np.array_equal(hi, (R >= 0.46) & (R <= 1.45)), None)
    emit("printed", "P.16", "Sec 4 para 6", "overshoot at 0.46-1.45 kpc: smallest (km/s)", res[hi].min(), "11", 0.5)
    emit("printed", "P.17", "Sec 4 para 6", "overshoot at 0.46-1.45 kpc: largest (km/s)", res[hi].max(), "17", 0.5)
    emit("printed", "P.18", "Sec 4 para 6", "errV at 0.46-1.45 kpc: smallest (km/s)", e[hi].min(), "2.5", 0.05)
    emit("printed", "P.19", "Sec 4 para 6", "errV at 0.46-1.45 kpc: largest (km/s)", e[hi].max(), "5.0", 0.05)
    emit("true", "P.20", "Sec 4 para 6", "main-disk points > 10 km/s low are the %d at 4.64-9.53 kpc" % lo.sum(),
         np.array_equal(lo, (R >= 4.64) & (R <= 9.53)), None)
    emit("printed", "P.21", "Sec 4 para 6", "deficit zone: first radius (kpc)", R[lo].min(), "4.6", 0.05)
    emit("printed", "P.22", "Sec 4 para 6", "deficit zone: last radius (kpc)", R[lo].max(), "9.5", 0.05)
    emit("printed", "P.23", "Sec 4 para 6", "deficit at 4.64-9.53 kpc: smallest (km/s)", -res[lo].max(), "19", 0.5)
    emit("printed", "P.24", "Sec 4 para 6", "deficit at 4.64-9.53 kpc: largest (km/s)", -res[lo].min(), "35", 0.5)
    emit("printed", "P.25", "Sec 4 para 6", "errV at 4.64-9.53 kpc: smallest (km/s)", e[lo].min(), "1.4", 0.05)
    emit("printed", "P.26", "Sec 4 para 6", "errV at 4.64-9.53 kpc: largest (km/s)", e[lo].max(), "3.2", 0.05)
    emit("printed", "P.27", "Sec 4 para 6", "share of Hermes chi2 at R > 10 without the floor (%)", q0("H6", lambda r: r > 10)["frac"], "61", 0.5)
    for i, (t, lab) in enumerate(((4.8, "t50 6.0 - 1.2"), (7.2, "t50 6.0 + 1.2"))):
        Vt = v_hermes(R, gb, t, phi=phi)
        emit("printed", "P.%d" % (28 + 2 * i), "Sec 4 para 6", "%s Gyr: full chi2/N" % lab, stats(B, Vt, R > 0)["chi2nu"], ("3.05", "3.40")[i], 0.005)
        emit("printed", "P.%d" % (29 + 2 * i), "Sec 4 para 6", "%s Gyr: R > 10 chi2/N" % lab, stats(B, Vt, R > 10)["chi2nu"], ("5.51", "6.20")[i], 0.005)
    gu = gbar(B, UNIT_ML_SCALE)
    Vu = {t: v_hermes(R, gu, t) for t in AGES}
    emit("printed", "P.32", "Sec 4 para 7", "unit-M/L analogue (Vdisk x 1.210): Hermes chi2/N", stats(B, Vu[6.0], R > 0)["chi2nu"], "2.55", 0.005)
    emit("printed", "P.33", "Sec 4 para 7", "unit-M/L analogue: MOND chi2/N", stats(B, v_mond(R, gu), R > 0)["chi2nu"], "0.19", 0.005)
    emit("check", "P.34", "Sec 4 para 7 (support)", "unit-M/L analogue, t50 7.1: Hermes chi2/N (2.67)", stats(B, Vu[7.1], R > 0)["chi2nu"], 2.67, 0.005)
    ge = gbar(B, np.sqrt(1.0 / ML0))
    emit("true", "P.35", "Sec 4 para 7 (support)", "same 2 dp at the exact scale sqrt(1/0.682647) = %.4f" % np.sqrt(1.0 / ML0),
         round(stats(B, v_hermes(R, ge, 6.0), R > 0)["chi2nu"], 2) == 2.55 and round(stats(B, v_mond(R, ge), R > 0)["chi2nu"], 2) == 0.19, None)
    Vcap = np.sqrt(gb * (1.0 + PI - F) * R)            # beta at its ceiling (psi = 0) and phi = 1
    emit("check", "P.36", "Sec 4 para 10 (support)", "outer points under-predicted even at the cap (of 27)", int((Vcap[R > 10] < Vo[R > 10]).sum()), 27, 0)
    # Section 6: Paper 8's operator (paper8/verify_eta.py) with a younger effective age outside 10 kpc
    psi_sys = ve.compute_psi(6.0, g98)
    psi_outer = lambda to: np.where(R <= 10, psi_sys, ve.compute_psi(to, g98))
    for i, (to, shown) in enumerate(((1.0, "1.30"), (2.0, "1.53"))):
        Ve = np.sqrt(np.maximum(ve.compute_g_model_eta(gb, phi, psi_outer(to), psi_sys), 0.0) * R)
        emit("printed", "P.%d" % (37 + i), "Sec 6", "Paper 8 refinement, outer age %.0f Gyr: R > 10 chi2/N" % to, stats(B, Ve, R > 10)["chi2nu"], shown, 0.005)
    # "neither alone closes the residual": each piece alone leaves the R > 10 chi2/N above 1
    Va = np.sqrt(np.maximum(ve.compute_g_model_eta(gb, phi, psi_sys, psi_sys), 0.0) * R)
    emit("check", "P.39", "Sec 6 (support)", "eta alone (one age, 6.0 Gyr): R > 10 chi2/N above 1 (Paper 8: 2.60)",
         stats(B, Va, R > 10)["chi2nu"], 2.60, 0.005)
    Vb = np.sqrt(np.maximum(ve.compute_g_model_baseline(gb, phi, psi_outer(0.0)), 0.0) * R)
    emit("check", "P.40", "Sec 6 (support)", "component age alone (outer 0 Gyr, no eta): R > 10 chi2/N above 1 (Paper 8: 4.12)",
         stats(B, Vb, R > 10)["chi2nu"], 4.12, 0.005)


# ------------------------------------------------------------------ main
def main():
    argparse.ArgumentParser(description=__doc__.split("\n")[0]).parse_args()
    os.makedirs(OUT, exist_ok=True)

    print("0. Shipped files")
    man = read_manifest()
    files = manifest_files()
    bad = [rel for rel, h in man.items() if not os.path.isfile(os.path.join(HERE, rel)) or sha256(os.path.join(HERE, rel)) != h]
    if not man:
        bad = ["MANIFEST_SHA256.txt missing or empty"]
    check("H.1", "MANIFEST_SHA256", "listed files whose SHA-256 differs or is missing (of %d)" % len(man), len(bad), 0, 0)
    if bad:
        print("          " + ", ".join(bad[:5]) + (" ..." if len(bad) > 5 else ""))
    extra = sorted(set(files) - set(man))
    check_true("H.2", "MANIFEST_SHA256", "the manifest lists every file in paper7/ (results/, __pycache__/, OS/editor junk excluded)", not extra and set(man) <= set(files),
               "unlisted: %d (%s)" % (len(extra), ", ".join(extra[:5])) if extra else "")
    fm = read_csv(os.path.join(HERE, "file_manifest.csv"))
    fm_ok = ({r["relative_path"] for r in fm} == set(man) - {"file_manifest.csv"}
             and all(man.get(r["relative_path"]) == r["sha256"]
                     and os.path.isfile(os.path.join(HERE, r["relative_path"]))
                     and int(r["bytes"]) == os.path.getsize(os.path.join(HERE, r["relative_path"])) for r in fm))
    check_true("H.3", "file_manifest.csv", "same files, bytes and SHA-256 as MANIFEST_SHA256.txt (%d rows)" % len(fm), fm_ok)
    for cid, path, what in (("H.4", IMPORT_CSV, "Corbelli (2014) Table 1 import"), ("H.5", BOARD_CSV, "the board"),
                            ("H.6", KERNEL_PY, "the frozen kernel"), ("H.7", AGE_CSV, "the age lock"),
                            ("H.8", CANON_PP, "01b per-point CSV (read by Paper 8 and verify_gates.py)")):
        h = sha256(path)
        check_true(cid, "SHA-256", "%s: %s..." % (what, SHA[path][:16]), h == SHA[path], "" if h == SHA[path] else "got " + h)
    check_true("H.9", "SHA-256", "01_ board input is byte-identical to the board", sha256(BOARD_COPY) == sha256(BOARD_CSV))
    row = next((r for r in fm if r["relative_path"] == "00_final_board/m33_corbelli_reconstructed_board.csv"), {})
    check_true("H.10", "Data Availability", "file_manifest.csv gives the board's SHA-256", row.get("sha256") == SHA[BOARD_CSV])
    ages = {float(r["t50_gyr"]): float(r["t50_sigma_gyr"]) for r in read_csv(AGE_CSV)}
    check_true("H.11", "Sec 4 para 1", "age lock: t50 6.0 +/- 1.2 and 7.1 +/- 1.3 Gyr", ages == {6.0: 1.2, 7.1: 1.3}, str(ages))

    print("1. Board: Corbelli (2014) Table 1 through the frozen kernel")
    B = load_board()
    R, Vo, e = B["R"], B["Vo"], B["e"]
    src = read_csv(IMPORT_CSV)
    scol = lambda n: np.array([float(r[n]) for r in src])
    reb = kernel.reconstruct(IMPORT_CSV, os.path.join(OUT, "m33_board_rebuilt.csv"))
    check("B.1", "Sec 4 / App A", "rebuilt Vgas vs the board (max |diff| km/s)", np.abs(reb["Vgas_kms"].to_numpy() - B["Vg"]).max(), 0, 1e-9)
    check("B.2", "Sec 4 / App A", "rebuilt Vdisk vs the board (max |diff| km/s)", np.abs(reb["Vdisk_kms"].to_numpy() - B["Vd"]).max(), 0, 1e-9)
    check("B.3", "Sec 4 / App A", "rebuilt Vbar, Sigma_H2, Sigma_gas vs the board (max |diff|)",
          max(np.abs(reb[c].to_numpy() - B["col"](c)).max() for c in ("Vbar_kms", "Sigma_H2_Msun_pc2", "Sigma_gas_total_He_Msun_pc2")), 0, 1e-9)
    check("B.4", "App A", "board R, Vobs, errV, Sigma_HI, Sigma_* = Table 1 (max |diff|)",
          max(np.abs(B["col"](b) - scol(a)).max() for a, b in (("R_kpc", "R_kpc"), ("V_r_kms", "Vobs_kms"), ("sigma_V_kms", "errV_kms"),
              ("Sigma_HI_Msun_pc2", "Sigma_HI_Msun_pc2"), ("Sigma_star_Msun_pc2", "Sigma_star_Msun_pc2"))), 0, 0)
    printed("B.5", "Sec 4 para 3", "radial points", len(R), "58", 0)
    printed("B.6", "Sec 4 para 3", "innermost radius (kpc)", R.min(), "0.24", 0.005)
    printed("B.7", "Sec 4 para 3", "outermost radius (kpc)", R.max(), "22.72", 0.005)
    check("B.8", "App A", "Vbulge = 0 at all radii (max |Vbulge|)", np.abs(B["Vb"]).max(), 0, 0)
    check("B.9", "App A (support)", "Sigma_H2 = 10 exp(-R/2.2), gas = 1.33 (HI + H2) (max |diff|)",
          max(np.abs(B["col"]("Sigma_H2_Msun_pc2") - 10.0 * np.exp(-R / 2.2)).max(),
              np.abs(B["col"]("Sigma_gas_total_He_Msun_pc2") - 1.33 * (B["col"]("Sigma_HI_Msun_pc2") + 10.0 * np.exp(-R / 2.2))).max()), 0, 1e-9)
    with open(KERNEL_PY, encoding="utf-8") as fh:
        ktxt = fh.read()
    kconst = all(s in ktxt for s in ("RMAX_KPC = 23.0", "h_gas = np.full_like(r_grid, 0.5)", "h_star = 0.1 + 0.9 * (r_grid / RMAX_KPC)",
                                     "Sigma_H2 = 10.0 * np.exp(-r_grid / 2.2)", "1.33 * (Sigma_HI + Sigma_H2)"))
    check_true("B.10", "App A (support)", "kernel: gas half-thickness 0.5 kpc; stellar 0.1 to 1.0 kpc at 23 kpc; H2; He 1.33", kconst)
    rg = np.arange(kernel.DR_KPC / 2.0, kernel.RMAX_KPC + 1e-12, kernel.DR_KPC)
    closure = {r["quantity"]: float(r["integrated_value_Msun"]) for r in read_csv(os.path.join(BOARD_DIR, "m33_corbelli_mass_closure_table.csv"))}
    for cid, q, prof in (("B.11", "Stellar mass", kernel.interp_profile(rg, R, scol("Sigma_star_Msun_pc2"))),
                         ("B.12", "H2 mass", 10.0 * np.exp(-rg / 2.2)),
                         ("B.13", "HI profile mass", kernel.interp_profile(rg, R, scol("Sigma_HI_Msun_pc2")))):
        mass = kernel.mass_from_surface_density(rg, prof)
        check(cid, "Sec 4 para 1 (support)", "mass closure: %s (Msun, 3 s.f. as tabulated)" % q,
              float("%.3g" % mass), closure[q], 0)
    printed("B.14", "App A", "distance of the board's radii (01b script's D0_KPC, kpc)", module_constant(CANON_PY, "D0_KPC"), "840", 0.5)

    print("2. Gate, scoring convention and the canonical-gate predictions of v3")
    gb = gbar(B)
    phi = hermes_phi(R, gb)
    g98 = float(np.percentile(gb, 98))
    pp = read_csv(CANON_PP)
    pcol = lambda n: np.array([float(r[n]) for r in pp])
    check("G.1", "Sec 2.1", "log-grid gate = 01b phi_canonical (max |diff|)", np.abs(phi - pcol("phi_canonical")).max(), 0, 1e-12)
    check("G.2", "Sec 2 eq. 3 (support)", "board g98, (km/s)^2/kpc", g98, 2653.949886626496, 1e-9)
    check("G.3", "Sec 2 eq. 2 (support)", "beta at t50 = 6.0 Gyr", beta(6.0, g98), 1.834212, 5e-7)
    check("G.4", "Sec 2 eq. 2 (support)", "beta at t50 = 7.1 Gyr", beta(7.1, g98), 1.698756, 5e-7)
    c_2pi = 299792.458 / (2 * PI) * (KPC_M / 1e3) / (1e9 * 365.25 * 86400.0)
    printed("G.5", "Sec 2", "c/(2 pi) in Gyr (km/s)^2/kpc", c_2pi, "46654", 0.5)
    V = {"H6": v_hermes(R, gb, 6.0), "H7": v_hermes(R, gb, 7.1), "MOND": v_mond(R, gb)}
    # The 01b per-point file was written with beta rounded to 5 dp (1.83421, 1.69876, as its memo
    # quotes); every number here uses the closed form at full precision. The two differ by < 1e-4 km/s.
    diffs = [np.abs(gb - pcol("gbar_kms2_per_kpc")).max() / max(1.0, np.abs(gb).max()),
             np.abs(np.sqrt(gb * R) - pcol("Vbar_kms")).max()]
    dfull = 0.0
    for k, t, tag in (("H6", 6.0, "primary_t60"), ("H7", 7.1, "conservative_t71")):
        V5 = np.sqrt(np.maximum(gb * (1.0 + round(beta(t, g98), 5) * phi), 0.0) * R)
        res = V5 - Vo
        diffs += [np.abs(V5 - pcol("Vhermes_%s_kms" % tag)).max(), np.abs(res - pcol("residual_%s" % tag)).max(),
                  np.abs(res ** 2 / (e ** 2 + SIG2) - pcol("chi2_contribution_%s" % tag)).max()]
        dfull = max(dfull, np.abs(V[k] - pcol("Vhermes_%s_kms" % tag)).max())
    check("G.6", "Sec 4 / Table 2", "01b per-point g_bar, Vbar, V_Hermes, residual, chi2 (beta at 5 dp)", max(diffs), 0, 1e-9)
    check("G.7", "Sec 4 / Table 2", "full-precision beta vs the 01b V_Hermes (max |diff| km/s)", dfull, 0, 1e-4)
    msum = {r["mask"]: r for r in read_csv(CANON_MASKS)}
    names01b = dict(zip([m[0] for m in MASKS], ["Full board", "Main disk R<=8", "Warp R>8", "Main disk R<=10", "Outer tail R>10",
                                                 "Barker break R>=9.1", "Barker S2 R>=11.6", "Putman R>=15"]))
    d01b = max(abs(stats(B, V[k], f(R))["chi2nu"] - float(msum[names01b[s]][c]))
               for s, _, f in MASKS for k, c in (("H6", "chi2nu_primary_t60"), ("H7", "chi2nu_conservative_t71")))
    check("G.8", "Table 2", "01b mask summary, 8 masks x 2 ages at 4 dp (max |diff|)", d01b, 0, 5e-5)
    sup = read_csv(SUPERSEDED)
    printed("G.9", "Sec 2.1", "full board on the superseded gate (01_, t50 6.0)",
            float(next(r for r in sup if float(r["t50_Gyr"]) == 6.0)["chi2_reduced"]), "3.13", 0.005)
    printed("G.10", "Sec 2.1", "full board on the canonical gate (t50 6.0)", stats(B, V["H6"], R > 0)["chi2nu"], "3.22", 0.005)
    printed("G.11", "Sec 2.5", "a0 = 1.2e-10 m/s^2 in (km/s)^2/kpc", A0, "3702.8", 0.05)
    with open(MOND_CFG, encoding="utf-8") as fh:
        a0_cfg = json.load(fh)["a0_kms2_per_kpc"]
    check("G.12", "Sec 2.5", "a0 equals the 02_ MOND run configuration", A0, a0_cfg, 1e-9)
    printed("G.13", "Sec 2.7", "sigma_int = sqrt(386) km/s", np.sqrt(SIG2), "19.6", 0.05)
    share = SIG2 / (e ** 2 + SIG2)
    printed("G.14", "Sec 2.7", "radii where the floor exceeds errV^2 (of 58)", int((SIG2 > e ** 2).sum()), "58", 0)
    printed("G.15", "Sec 2.7", "median floor share of each point's variance (%)", 100 * np.median(share), "98", 0.5)

    print("3. Table 2 (full board and radial masks; 8 printed cells per row)")
    maskf = {s: f for s, _, f in MASKS}
    label = {s: lab for s, lab, _ in MASKS}
    mlabel = dict(MODELS)
    t2rows = []
    for i, (s, k, N, r0, r1, chi_, frac, mean, med, rms) in enumerate(TABLE2, 1):
        st = stats(B, V[k], maskf[s](R))
        got = (st["N"], st["Rmin"], st["Rmax"], st["chi2nu"], st["frac"], st["mean"], st["med"], st["rms"])
        d = max(abs(a - b) for a, b in zip(got, (N, r0, r1, chi_, frac, mean, med, rms)))
        check("T2.%d" % i, "Table 2", "%s, %s (max |diff| of 8 cells)" % (s, mlabel[k][:27]), d, 0, 0.005)
        t2rows.append([label[s], num(st["Rmin"], "%.2f"), num(st["Rmax"], "%.2f"), st["N"], mlabel[k], num(st["chi2nu"]),
                       num(st["frac"]), num(st["mean"]), num(st["med"]), num(st["rms"])])
    write_csv("table2_masks.csv", ["mask", "R_min_kpc", "R_max_kpc", "N", "model", "chi2nu", "chi2_fraction_pct",
                                   "mean_residual_kms", "median_abs_residual_kms", "rms_residual_kms"], t2rows)

    print("4. Table 2b (floor sensitivity; Hermes at t50 = 6.0 Gyr, log-grid gate)")
    t2b, n = [], 0
    for sig2, *cells in TABLE2B:
        h, m = (lambda f: stats(B, V["H6"], f(R), sig2)), (lambda f: stats(B, V["MOND"], f(R), sig2))
        got = (h(lambda r: r > 0)["chi2nu"], h(lambda r: r <= 10)["chi2nu"], h(lambda r: r > 10)["chi2nu"],
               h(lambda r: r > 10)["frac"], m(lambda r: r > 0)["chi2nu"], m(lambda r: r > 10)["chi2nu"])
        for col, g, w in zip(TABLE2B_COLS, got, cells):
            n += 1
            check("T2b.%d" % n, "Table 2b", "sigma_int^2 = %d: %s" % (sig2, col), g, float(w), 0.5 * 10.0 ** -decimals(w))
        t2b.append([sig2] + [num(g) for g in got])
    write_csv("table2b_floor_sweep.csv", ["sigma_int_sq_kms2", "hermes_t6.0_chi2nu", "hermes_t6.0_R_le_10_chi2nu",
                                          "hermes_t6.0_R_gt_10_chi2nu", "hermes_t6.0_chi2_share_R_gt_10_pct",
                                          "mond_chi2nu", "mond_R_gt_10_chi2nu"], t2b)

    print("5. Abstract and Sections 4 and 6 (M33 numbers in the prose)")
    prose_checks(B, emit_check)

    print("6. Table 3 (Vdisk scaling sweep, g98, beta and phi recomputed at each scale)")
    printed("T3.1", "Table 3 rows 1-2", "source-native M*/L at 3.6 um", ML0, "0.683", 0.0005)
    s03 = next(r for r in read_csv(SWEEP_03) if float(r["scale_Vdisk"]) == 1.0)
    check("T3.2", "Table 3 (support)", "same M*/L as the 03_ sweep files", ML0, float(s03["implied_global_ML_W1"]), 1e-12)
    rows = {n: sweep_row(B, n) for n in range(SWEEP_COARSE[0], SWEEP_COARSE[1] + 1, SWEEP_COARSE[2])}
    coarse_min = {t: min(rows, key=lambda n: rows[n][t][0]) for t in AGES}
    for t in AGES:
        for n in range(coarse_min[t] - SWEEP_FINE_HALF, coarse_min[t] + SWEEP_FINE_HALF + 1):
            if n not in rows:
                rows[n] = sweep_row(B, n)
    best = {t: min(rows, key=lambda n: rows[n][t][0]) for t in AGES}
    check_true("T3.3", "Table 3 (support)", "fine-grid minima lie inside their +/-0.005 fine windows, not at an edge (6.0: %.4f, 7.1: %.4f)"
               % (best[6.0] / 1e4, best[7.1] / 1e4), all(abs(best[t] - coarse_min[t]) < SWEEP_FINE_HALF for t in AGES))
    printed("T3.4", "Table 3 row 3", "Vdisk scale of the lowest Hermes chi2/N (t50 6.0)", best[6.0] / 1e4, "1.368", 0.0005)
    printed("T3.5", "Table 3 / Sec 4 para 7", "lowest Hermes chi2/N on the scaling (t50 6.0)", rows[best[6.0]][6.0][0], "2.41", 0.005)
    printed("T3.6", "Sec 4 para 7", "disk mass there, x source-native (scale^2)", (best[6.0] / 1e4) ** 2, "1.87", 0.005)
    printed("T3.7", "Sec 4 para 7", "MOND chi2/N there", rows[best[6.0]]["mond"][0], "0.53", 0.005)
    check("T3.8", "Table 3 (support)", "Vdisk scale of the lowest Hermes chi2/N (t50 7.1)", best[7.1] / 1e4, 1.415, 0.0005)
    check("T3.9", "Table 3 (support)", "lowest Hermes chi2/N on the scaling (t50 7.1)", rows[best[7.1]][7.1][0], 2.48, 0.005)
    b6 = rows[best[6.0]][6.0]
    out_share = 27 * b6[2] / (58 * b6[0])
    check_true("T3.10", "Sec 4 para 7 (support)", "R > 10 still dominates Hermes chi2 at the minimum (share > 50%)", out_share > 0.5,
               "%.1f%%" % (100 * out_share))
    elim = [n for n in sorted(rows) if min(rows[n][6.0][2], rows[n][7.1][2]) <= 1.0]
    mond_there = min(rows[n]["mond"][0] for n in elim) if elim else float("nan")
    check_true("T3.11", "Table 3 caption (support)", "Hermes R > 10 chi2/N <= 1 only where MOND chi2/N > 1 (first at %.4f; MOND >= %.2f)"
               % (elim[0] / 1e4 if elim else float("nan"), mond_there), bool(elim) and mond_there > 1.0)
    t3, dfac, dml = [], 0.0, 0.0
    for i, (sc, fac, ml, t, h, ho, hm, m, mo, mm, note) in enumerate(TABLE3):
        n = int(round(float(sc) * 1e4))
        r = rows[n] if n in rows else sweep_row(B, n)
        got = (r[t][0], r[t][2], r[t][3], r["mond"][0], r["mond"][1], r["mond"][2])
        d = max(abs(a - float(b)) for a, b in zip(got, (h, ho, hm, m, mo, mm)))
        check("T3.%d" % (12 + i), "Table 3 row %d" % (i + 1), "scale %s, t50 %.1f: 6 model cells (max |diff|)" % (sc, t), d, 0, 0.005)
        s = float(sc)
        dfac, dml = max(dfac, abs(s * s - float(fac))), max(dml, abs(ML0 * s * s - float(ml)))
        t3.append([sc, num(s * s, "%.6f"), num(ML0 * s * s, "%.6f"), "%.1f" % t, num(r["g98"]), num(beta(t, r["g98"])),
                   num(got[0]), num(r[t][1]), num(got[1]), num(got[2]), num(got[3]), num(got[4]), num(got[5]), note])
    check("T3.23", "Table 3", "M/L factor column, 11 cells = scale^2 (max |diff|)", dfac, 0, 0.005)
    check("T3.24", "Table 3", "global M*/L column, 11 cells = 0.682647 scale^2 (max |diff|)", dml, 0, 0.0005)
    write_csv("table3_rows_printed.csv", ["scale_Vdisk", "ML_factor", "global_ML_3p6um", "t50_Gyr", "g98", "beta", "hermes_chi2nu",
                                          "hermes_R_le_10_chi2nu", "hermes_R_gt_10_chi2nu", "hermes_R_gt_10_mean_residual_kms",
                                          "mond_chi2nu", "mond_R_gt_10_chi2nu", "mond_R_gt_10_mean_residual_kms", "note"], t3)
    fine = set()
    for t in AGES:
        fine |= set(range(coarse_min[t] - SWEEP_FINE_HALF, coarse_min[t] + SWEEP_FINE_HALF + 1))
    sweep = []
    for n in sorted(rows):
        r, s = rows[n], n / 1e4
        sweep.append(["%.4f" % s, "0.0025" if (n - SWEEP_COARSE[0]) % SWEEP_COARSE[2] == 0 else "0.0001",
                      num(s * s, "%.8f"), num(ML0 * s * s, "%.8f"), num(r["g98"])]
                     + [num(x) for t in AGES for x in r[t]] + [num(x) for x in r["mond"]])
    write_csv("table3_vdisk_sweep.csv", ["scale_Vdisk", "grid_step", "ML_factor", "global_ML_3p6um", "g98"]
              + ["hermes_t%.1f_%s" % (t, c) for t in AGES for c in ("chi2nu", "R_le_10_chi2nu", "R_gt_10_chi2nu", "R_gt_10_mean_residual_kms")]
              + ["mond_chi2nu", "mond_R_gt_10_chi2nu", "mond_R_gt_10_mean_residual_kms"], sweep)

    print("7. Appendix A (Vgas and Vdisk printed to 3 dp)")
    appa = [["%.2f" % r, "%.1f" % v, "%.1f" % s, "%.3f" % g, "%.3f" % d, "%.1f" % b]
            for r, v, s, g, d, b in zip(R, Vo, e, B["Vg"], B["Vd"], B["Vb"])]
    exact = all(float(a[0]) == r and float(a[1]) == v and float(a[2]) == s for a, r, v, s in zip(appa, R, Vo, e))
    check_true("A.1", "App A", "R, Vobs and errV print exactly at 2, 1 and 1 dp (58 rows)", exact)
    write_csv("appendix_a_board_3dp.csv", ["R_kpc", "Vobs_kms", "errV_kms", "Vgas_kms", "Vdisk_kms", "Vbulge_kms"], appa)
    A3 = dict(B, Vg=np.array([float(a[3]) for a in appa]), Vd=np.array([float(a[4]) for a in appa]), Vb=np.array([float(a[5]) for a in appa]))
    g3 = gbar(A3)
    V3 = {"H6": v_hermes(R, g3, 6.0), "H7": v_hermes(R, g3, 7.1), "MOND": v_mond(R, g3)}
    printed("A.2", "App A / Sec 4 para 4", "from the 3-dp table: Hermes chi2/N, t50 6.0", stats(B, V3["H6"], R > 0)["chi2nu"], "3.22", 0.005)
    printed("A.3", "App A / Sec 4 para 4", "from the 3-dp table: Hermes chi2/N, t50 7.1", stats(B, V3["H7"], R > 0)["chi2nu"], "3.39", 0.005)
    printed("A.4", "App A / Sec 4 para 6", "from the 3-dp table: Hermes R <= 10 chi2/N", stats(B, V3["H6"], R <= 10)["chi2nu"], "0.93", 0.005)
    printed("A.5", "App A / Sec 4 para 6", "from the 3-dp table: Hermes R > 10 chi2/N", stats(B, V3["H6"], R > 10)["chi2nu"], "5.86", 0.005)
    printed("A.6", "App A / Sec 4 para 3", "from the 3-dp table: MOND chi2/N", stats(B, V3["MOND"], R > 0)["chi2nu"], "0.12", 0.005)
    printed("A.7", "App A / Table 2", "from the 3-dp table: MOND R > 10 chi2/N", stats(B, V3["MOND"], R > 10)["chi2nu"], "0.23", 0.005)
    A2 = dict(B, Vg=np.round(B["Vg"], 2), Vd=np.round(B["Vd"], 2))
    g2 = gbar(A2)
    check("A.8", "App A note (support)", "v3's 2-dp table gives 3.38, not 3.39 (t50 7.1)", stats(B, v_hermes(R, g2, 7.1), R > 0)["chi2nu"], 3.38, 0.005)
    check("A.9", "App A note (support)", "v3's 2-dp table gives 0.92, not 0.93 (R <= 10)", stats(B, v_hermes(R, g2, 6.0), R <= 10)["chi2nu"], 0.92, 0.005)
    if os.path.isfile(PAPER_TEXT):
        rows_txt = appendix_a_rows(PAPER_TEXT)
        cols = (R, Vo, e, B["Vg"], B["Vd"], B["Vb"])
        ok_cells = sum(1 for i, row in enumerate(rows_txt[:len(R)]) for j, c in enumerate(row[:6])
                       if abs(float(c) - cols[j][i]) <= 0.5 * 10.0 ** -decimals(c) * (1 + 1e-6))
        dp = sorted({decimals(row[3]) for row in rows_txt} | {decimals(row[4]) for row in rows_txt})
        v4 = is_v4_text(PAPER_TEXT)          # decision D-23: the v4 text prints Vgas and Vdisk to 3 dp
        check_true("A.10", "App A", "paper7_main.md (v%d text): %d/%d cells equal the board at their printed precision; Vgas/Vdisk %s dp%s"
                   % (4 if v4 else 3, ok_cells, 6 * len(R), "/".join(map(str, dp)), " (v4 needs 3)" if v4 else ""),
                   ok_cells == 6 * len(R) and (dp == [3] or not v4))
    else:
        check_true("A.10", "App A", "paper7_main.md present", False)
    # What the 3-dp table does and does not reproduce at printed precision
    d3t = max(abs(stats(B, V3[k], maskf[s](R))["chi2nu"] - chi_) for s, k, _, _, _, chi_, *_ in TABLE2)
    for sc, fac, ml, t, h, ho, hm, m, mo, mm, note in TABLE3:
        gs = gbar(A3, float(sc))
        Vh, Vm = v_hermes(R, gs, t), v_mond(R, gs)
        d3t = max(d3t, *(abs(stats(B, VV, f)["chi2nu"] - float(w)) for VV, f, w in
                         ((Vh, R > 0, h), (Vh, R > 10, ho), (Vm, R > 0, m), (Vm, R > 10, mo))))
    check("A.11", "App A note (support)", "from the 3-dp table: 24 Table 2 and 44 Table 3 chi2/N cells (max |diff|)", d3t, 0, 0.005)
    p3 = []
    prose_checks(A3, lambda kind, cid, where, what, got, target, tol=0, detail="": p3.append((kind, cid, got, target, tol)))
    shown = [x for x in p3 if x[0] == "printed"]
    miss = [cid for kind, cid, got, target, tol in p3
            if (kind == "printed" and abs(float(got) - float(target)) > tol * (1 + 1e-6)) or (kind == "true" and not got)]
    check_true("A.12", "App A note (support)", "from the 3-dp table: the prose numbers P.1-P.38 at printed precision except P.6 (%d of %d)"
               % (len(shown) - len(miss), len(shown)), miss == ["P.6"], "misses: " + (", ".join(miss) or "none"))
    for cid, f, want, lab in (("A.13", lambda r: r > 0, 118.51, "full (printed 118.48; Sec 4's 118)"),
                              ("A.14", lambda r: r <= 10, 87.54, "R <= 10 (printed 87.49)")):
        check(cid, "App A note (support)", "from the 3-dp table, no floor: Hermes %s; needs the CSV" % lab,
              stats(B, V3["H6"], f(R), 0.0)["chi2nu"], want, 0.005)

    print("8. The chain-rule gate (Paper 1's published scores): same chi2/N to 2 dp")
    phic = phi_chain(R, gb)
    info("C.1", "Sec 2.1", "max |phi_chain - phi_loggrid| on the board (no printed claim)", np.abs(phic - phi).max())
    d2 = max(abs(stats(B, v_hermes(R, gb, 6.0 if k == "H6" else 7.1, gate=phi_chain), maskf[s](R))["chi2nu"] - chi_)
             for s, k, _, _, _, chi_, *_ in TABLE2 if k != "MOND")
    check("C.2", "Sec 2.1 / Table 2", "16 Hermes chi2/N cells on the chain-rule gate (max |diff|)", d2, 0, 0.005)
    d3 = 0.0
    for sc, fac, ml, t, h, ho, *_ in TABLE3:
        gs = gbar(B, float(sc))
        Vc = v_hermes(R, gs, t, gate=phi_chain)
        d3 = max(d3, abs(stats(B, Vc, R > 0)["chi2nu"] - float(h)), abs(stats(B, Vc, R > 10)["chi2nu"] - float(ho)))
    check("C.3", "Sec 2.1 / Table 3", "22 Hermes chi2/N cells on the chain-rule gate (max |diff|)", d3, 0, 0.005)

    print("9. SPARC context quoted in Section 4 (from shipped files; no download)")
    p1 = read_csv(P1_EXPORT)
    h1 = np.array([float(r["chi2nu_configg"]) for r in p1])
    m8 = np.array([float(r["mond_chi2nu"]) for r in read_csv(P8_SPARC)])
    printed("S.1", "Sec 4 para 5", "Paper 1 published median chi2/N (chain-rule scores)", np.median(h1), "1.312", 0.0005)
    printed("S.2", "Sec 4 para 5", "MOND median chi2/N (a0 3702.813; Paper 8 table)", np.median(m8), "1.143", 0.0005)
    printed("S.3", "Sec 4 para 5", "Hermes galaxies below chi2/N 5 (%, of 133)", 100 * (h1 < 5).mean(), "84", 0.5)
    printed("S.4", "Sec 4 para 5", "MOND galaxies below chi2/N 5 (%, a0 3702.813; Paper 8 table)", 100 * (m8 < 5).mean(), "71", 0.5)

    write_csv("prose_numbers.csv", ["id", "location", "quantity", "value", "printed", "tolerance"],
              [[c, w, q, num(v), p, "%g" % t] for c, w, q, v, p, t in PROSE])

    print("10. Regenerated tables against the shipped files (expected/)")
    want_rows = {"table2_masks.csv": 24, "table2b_floor_sweep.csv": 4, "table3_rows_printed.csv": 11,
                 "appendix_a_board_3dp.csv": 58, "prose_numbers.csv": len(PROSE), "table3_vdisk_sweep.csv": len(sweep)}
    for i, fn in enumerate(sorted(want_rows), 1):
        if not os.path.isfile(os.path.join(EXP, fn)):
            check_true("F.%d" % i, "expected/", "%s present" % fn, False)
            continue
        a, b = read_csv(os.path.join(OUT, fn)), read_csv(os.path.join(EXP, fn))
        same, why = cells_match(a, b)
        check_true("F.%d" % i, "expected/", "%s rebuilt (%d rows)" % (fn, len(b)), same and len(b) == want_rows[fn], why)

    with open(os.path.join(OUT, "checks.csv"), "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, lineterminator="\n")
        w.writerow(["id", "where", "quantity", "got", "expected", "tolerance", "pass"])
        w.writerows(RESULTS + INFO)
    n_fail = sum(1 for r in RESULTS if not r[-1])
    print("\n%d checks, %d failed (plus %d informational); results in %s"
          % (len(RESULTS), n_fail, len(INFO), os.path.relpath(OUT, ROOT)))
    return 1 if n_fail else 0


if __name__ == "__main__":
    sys.exit(main())
