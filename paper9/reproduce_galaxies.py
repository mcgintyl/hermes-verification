#!/usr/bin/env python3
"""Reproduce the galaxy baseline of Paper 9 (Table 1) from the SPARC rotation curves.

    python fetch_sparc.py        (once, to get and verify the data)
    python reproduce_galaxies.py

Writes results/table1_galaxies.csv (the per-galaxy scores) and results/checks_galaxies.txt,
and exits with status 0 only if every check against the paper's values passes.
"""
import csv
import sys
from pathlib import Path

import numpy as np
from scipy.stats import binomtest, wilcoxon

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))
from hermes_clusters import galaxies as gx  # noqa: E402

OUT = ROOT / "results"
EXP = ROOT / "expected"
SPARC = ROOT / "data" / "sparc"
OUT.mkdir(exist_ok=True)
CHECKS = []


def check(name, ok, detail):
    CHECKS.append((bool(ok), name, detail))
    print("  [%s] %s: %s" % ("PASS" if ok else "FAIL", name, detail))


def near(x, target, places):
    return round(float(x), places) == target


manifest = gx.load_manifest()
ages = gx.load_ages()
missing = [g for g, _, _ in ages if gx.rotmod_path(SPARC, g, manifest) is None]
if missing:
    print("The SPARC rotation curves are not in data/sparc (%d missing)." % len(missing))
    print("Run:  python fetch_sparc.py")
    raise SystemExit(2)

names, curves, H, M, per_H, per_M, pooled_H, pooled_M, npts = [], [], [], [], [], [], [], [], 0
for galaxy, t50, g98 in ages:
    d = gx.read_rotmod(gx.rotmod_path(SPARC, galaxy, manifest))
    vh = gx.hermes_velocity(d, t50, g98)
    vm = gx.mond_velocity(d)
    names.append(galaxy)
    curves.append(d)
    H.append(gx.chi2_per_point(d["Vobs"], vh, d["errV"]))
    M.append(gx.chi2_per_point(d["Vobs"], vm, d["errV"]))
    pooled_H.append(np.abs(d["Vobs"] - vh))
    pooled_M.append(np.abs(d["Vobs"] - vm))
    per_H.append(float(np.median(np.abs(d["Vobs"] - vh))))
    per_M.append(float(np.median(np.abs(d["Vobs"] - vm))))
    npts += len(d["R"])
H, M = np.array(H), np.array(M)
pooled_H, pooled_M = np.concatenate(pooled_H), np.concatenate(pooled_M)

print("1. Table 1: the galaxy baseline")
check("133 galaxies and 3,073 data points", len(H) == 133 and npts == 3073,
      "%d galaxies, %d points" % (len(H), npts))
check("median chi2/N", near(np.median(H), 1.311535, 6) and near(np.median(M), 1.142686, 6),
      "Hermes 1 %.6f, MOND %.6f" % (np.median(H), np.median(M)))
check("5% trimmed mean chi2/N", near(gx.trimmed_mean(H), 2.495, 3) and near(gx.trimmed_mean(M), 3.295, 3),
      "Hermes 1 %.4f, MOND %.4f" % (gx.trimmed_mean(H), gx.trimmed_mean(M)))
check("median |dV|, pooled over all points", near(np.median(pooled_H), 20.6, 1) and near(np.median(pooled_M), 25.8, 1),
      "Hermes 1 %.1f km/s, MOND %.1f km/s" % (np.median(pooled_H), np.median(pooled_M)))
check("median |dV|, per galaxy", near(np.median(per_H), 18.9, 1) and near(np.median(per_M), 20.2, 1),
      "Hermes 1 %.1f km/s, MOND %.1f km/s" % (np.median(per_H), np.median(per_M)))
check("catastrophes (chi2/N > 10)", int((H > 10).sum()) == 12 and int((M > 10).sum()) == 15,
      "Hermes 1 %d, MOND %d" % ((H > 10).sum(), (M > 10).sum()))
wins = int((H < M).sum())
check("head-to-head wins", wins == 68 and len(H) - wins == 65, "Hermes 1 %d, MOND %d" % (wins, len(H) - wins))

print("2. The two are statistically indistinguishable (section 3)")
w = wilcoxon(np.log(H), np.log(M))
s = binomtest(wins, len(H), 0.5)
check("paired Wilcoxon on log chi2", near(w.pvalue, 0.9365, 4), "p = %.4f" % w.pvalue)
check("exact sign test on the head-to-head record", near(s.pvalue, 0.8624, 4), "p = %.4f" % s.pvalue)

print("3. Per-galaxy agreement with the published Paper 1 scores")
with open(EXP / "paper1_per_galaxy.csv", encoding="utf-8", newline="") as f:
    ref = {r["galaxy"]: float(r["chi2nu_hermes"]) for r in csv.DictReader(f)}
worst = max((abs(h - ref[n]), n) for n, h in zip(names, H) if n in ref)
check("all 133 galaxies match Paper 1's published scores", len(ref) == 133 and worst[0] < 1e-9,
      "largest difference %.1e (%s)" % (worst[0], worst[1]))


def mond_median(a0):
    return float(np.median([gx.chi2_per_point(d["Vobs"], gx.mond_velocity(d, a0=a0), d["errV"]) for d in curves]))


print("   How the MOND number moves with the conversion of a0 = 1.2e-10 m/s^2:")
print("     %-46s %.6f  <- Table 1" % ("a0 = %.3f, rounded unit constant" % gx.A0, np.median(M)))
print("     %-46s %.6f" % ("a0 = %.3f, exact parsec conversion" % gx.A0_EXACT, mond_median(gx.A0_EXACT)))
print("     %-46s %.6f" % ("a0 = %.1f, the archived Paper 1 export" % gx.A0_ROUNDED_EXPORT,
                           mond_median(gx.A0_ROUNDED_EXPORT)))
print("   The cluster half of this package uses the exact conversion, as its table did.")

with open(OUT / "table1_galaxies.csv", "w", encoding="utf-8", newline="") as f:
    writer = csv.writer(f, lineterminator="\n")
    writer.writerow(["galaxy", "t50_gyr", "g98", "n_points", "chi2N_hermes1", "chi2N_mond",
                     "median_abs_dV_hermes1", "median_abs_dV_mond"])
    for i, (galaxy, t50, g98) in enumerate(ages):
        writer.writerow([galaxy, t50, g98, len(curves[i]["R"]), "%.9f" % H[i], "%.9f" % M[i],
                         "%.4f" % per_H[i], "%.4f" % per_M[i]])
passed = sum(ok for ok, _, _ in CHECKS)
lines = ["%s  %s: %s" % ("PASS" if ok else "FAIL", n, d) for ok, n, d in CHECKS]
lines.append("%d of %d checks passed" % (passed, len(CHECKS)))
(OUT / "checks_galaxies.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
print()
print(lines[-1])
sys.exit(0 if passed == len(CHECKS) else 1)
