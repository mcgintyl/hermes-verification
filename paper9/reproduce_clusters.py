#!/usr/bin/env python3
"""Reproduce the cluster results of Paper 9 of the Hermes series from the published inputs.

    python reproduce_clusters.py

Uses only the files in data/ and the code in hermes_clusters/. Writes
results/table2_seven.csv (Table 2, in the same format as the published table),
results/boards.csv and results/checks.txt, and exits with status 0 only if every
check against the paper's values passes.
"""
import hashlib
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import chi2 as C2

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))
from hermes_clusters import model  # noqa: E402

OUT = ROOT / "results"
EXP = ROOT / "expected"
OUT.mkdir(exist_ok=True)
CHECKS = []
FB_COSMIC = 0.156
STELLAR_BENCHMARK = (0.10, 0.15)  # assumed stellar share of cluster baryons, for the 3-8% estimate


def check(name, ok, detail):
    CHECKS.append((bool(ok), name, detail))
    print("  [%s] %s: %s" % ("PASS" if ok else "FAIL", name, detail))


def near(x, target, places):
    return round(float(x), places) == target


famaey = model.load_famaey()
B = {row[0]: model.build_board(row, famaey) for row in model.SEVEN + model.EXCLUDED}
SEVEN = [row[0] for row in model.SEVEN]
rmax_kpc = pd.read_csv(model.DATA / "rmax_x.csv").set_index("board")["rmax_x_mpc"] * 1000.0

# ---------------------------------------------------------------------------
print("1. Table 2: the seven relaxed clusters")
rows = []
for row in model.SEVEN:
    b = B[row[0]]
    s = model.score(b)
    n = len(b.r)
    rows.append(dict(board=row[6], nodes=n, dof=n - 1, H1_chi2=s["H1_chi2"], H1_chi2_per_dof=s["H1_chi2"] / (n - 1),
                     H1_K=s["H1_K"], MOND_chi2=s["MOND_chi2"], MOND_chi2_per_dof=s["MOND_chi2"] / (n - 1),
                     MOND_K=s["MOND_K"]))
T = pd.DataFrame(rows)
tot = dict(board="TOTAL (seven)", nodes=int(T.nodes.sum()), dof=int(T.dof.sum()), H1_chi2=float(T.H1_chi2.sum()),
           MOND_chi2=float(T.MOND_chi2.sum()))
tot["H1_chi2_per_dof"] = tot["H1_chi2"] / tot["dof"]
tot["MOND_chi2_per_dof"] = tot["MOND_chi2"] / tot["dof"]
tot["H1_K"] = float(np.median(T.H1_K))
tot["MOND_K"] = float(np.median(T.MOND_K))
table = OUT / "table2_seven.csv"
pd.concat([T, pd.DataFrame([tot])], ignore_index=True).to_csv(table, index=False, float_format="%.6f",
                                                                lineterminator="\r\n")
digest = hashlib.sha256(table.read_bytes()).hexdigest()
check("Table 2 is byte-identical to the published table", table.read_bytes() == (EXP / "table2_seven.csv").read_bytes(),
      "SHA-256 " + digest)
ref = pd.read_csv(EXP / "table2_full_precision.csv").set_index("board")
worst = max(abs(r[c] - ref.loc[row[0], c]) for row, r in zip(model.SEVEN, rows)
            for c in ("H1_chi2", "H1_K", "MOND_chi2", "MOND_K"))
check("chi2 and K agree with the full-precision reference", worst < 1e-9, "largest difference %.1e" % worst)
pH, pM = C2.sf(tot["H1_chi2"], tot["dof"]), C2.sf(tot["MOND_chi2"], tot["dof"])
check("Hermes 1 total", near(tot["H1_chi2"], 71.132, 3) and tot["dof"] == 44 and near(tot["H1_chi2_per_dof"], 1.617, 3)
      and near(pH, 0.0059, 4), "%.3f / %d = %.3f, p = %.4f" % (tot["H1_chi2"], tot["dof"], tot["H1_chi2_per_dof"], pH))
check("MOND total", near(tot["MOND_chi2"], 101.787, 3) and near(tot["MOND_chi2_per_dof"], 2.313, 3),
      "%.3f / %d = %.3f, p = %.1e" % (tot["MOND_chi2"], tot["dof"], tot["MOND_chi2_per_dof"], pM))
check("median and mean K", near(tot["H1_K"], 4.24, 2) and near(T.H1_K.mean(), 4.20, 2) and near(tot["MOND_K"], 2.17, 2),
      "Hermes 1 median %.2f, mean %.2f; MOND median %.2f" % (tot["H1_K"], T.H1_K.mean(), tot["MOND_K"]))
check("K ranges", near(T.H1_K.min(), 2.9, 1) and near(T.H1_K.max(), 5.4, 1) and near(T.MOND_K.min(), 1.4, 1)
      and near(T.MOND_K.max(), 2.5, 1),
      "Hermes 1 %.2f-%.2f; MOND %.2f-%.2f" % (T.H1_K.min(), T.H1_K.max(), T.MOND_K.min(), T.MOND_K.max()))
wins = int((T.H1_chi2_per_dof < T.MOND_chi2_per_dof).sum())
check("Hermes 1 wins 5 of 7 boards", wins == 5, "%d of 7 (MOND wins %s)" % (
    wins, ", ".join(T.board[T.H1_chi2_per_dof >= T.MOND_chi2_per_dof])))
check("5 to 9 nodes per board", T.nodes.min() == 5 and T.nodes.max() == 9, "%d to %d" % (T.nodes.min(), T.nodes.max()))
a = B["A611"]
check("A611 worked example", len(a.r) == 8 and round(a.r.min()) == 276 and round(a.r.max()) == 2072
      and round(a.g.min()) == 170 and round(a.g.max()) == 688 and near(T.H1_chi2_per_dof[0], 0.46, 2),
      "%d nodes, %.0f-%.0f kpc, g_bar %.0f-%.0f (km/s)^2/kpc, chi2/dof %.2f" % (
          len(a.r), a.r.min(), a.r.max(), a.g.min(), a.g.max(), T.H1_chi2_per_dof[0]))

# ---------------------------------------------------------------------------
print("2. Excluded clusters quoted in section 4")
for bid, h_t, m_t in (("MACSJ0744", 1.78, 3.14), ("RXJ1347", 4.66, 4.14)):
    b = B[bid]
    s = model.score(b)
    d = len(b.r) - 1
    check(bid, near(s["H1_chi2"] / d, h_t, 2) and near(s["MOND_chi2"] / d, m_t, 2),
          "Hermes 1 %.2f per dof (p %.2g), MOND %.2f per dof (p %.2g), %d nodes" % (
              s["H1_chi2"] / d, C2.sf(s["H1_chi2"], d), s["MOND_chi2"] / d, C2.sf(s["MOND_chi2"], d), d + 1))

# ---------------------------------------------------------------------------
print("3. Age: chi2 does not depend on t50 (section 6)")
res = {}
for t50 in (5.0, 10.0, 13.0):
    sc = [model.score(B[bid], t50) for bid in SEVEN]
    res[t50] = (sum(s["H1_chi2"] for s in sc), float(np.median([s["H1_K"] for s in sc])))
spread = max(v[0] for v in res.values()) - min(v[0] for v in res.values())
check("chi2 identical at t50 = 5, 10 and 13 Gyr", spread < 1e-9,
      "total chi2 %.6f (spread %.1e); median K %.2f, %.2f, %.2f" % (res[10.0][0], spread, res[5.0][1], res[10.0][1],
                                                                     res[13.0][1]))

# ---------------------------------------------------------------------------
print("4. Node admission and X-ray coverage (Appendix A)")
led = pd.concat([B[bid].ledger for bid in SEVEN], ignore_index=True)
n_out = int(led.drop_reason_codes.str.contains("OUTSIDE_BARYONIC_NATIVE_SUPPORT").sum())
n_warn = int(led.drop_reason_codes.str.contains("AUTHOR_WARNED_NODE").sum())
check("node admission", len(led) == 70 and n_out == 16 and n_warn == 3 and int(led.admitted.sum()) == 51,
      "%d raw nodes; %d beyond the baryon model, %d author-warned (RXJ2129); %d scored" % (
          len(led), n_out, n_warn, int(led.admitted.sum())))
beyond = {bid: int(np.sum(B[bid].r > rmax_kpc[bid])) for bid in SEVEN}
n51 = sum(len(B[bid].r) for bid in SEVEN)
check("nodes beyond Rmax_X", sum(beyond.values()) == 33 and n51 == 51, "%d of %d" % (sum(beyond.values()), n51))
inside = {bid: len(B[bid].r) - beyond[bid] for bid in SEVEN}
check("MACSJ0429 has no node on measured gas; MACSJ1115 has one", inside["MACSJ0429"] == 0 and inside["MACSJ1115"] == 1,
      "nodes inside Rmax_X: " + ", ".join("%s %d" % (k.replace("-PARTIAL", ""), v) for k, v in inside.items()))

# ---------------------------------------------------------------------------
print("5. Carrier split inside each board's outermost node radius (section 4)")
fr = []
for bid in SEVEN:
    b = B[bid]
    r, mtot, mg, mbcg, fgal = model.carrier_with_tail(b.spec, famaey)
    assert np.max(np.abs(mtot / b.M - 1)) < 1e-12, bid + ": carrier rebuilt from parts differs"
    i = int(np.argmax(b.r))
    parts = [float(np.interp(np.log(b.r[i]), np.log(r), c)) for c in (np.full_like(mg, mbcg), mg, mg * fgal)]
    fr.append([100.0 * v / sum(parts) for v in parts])
fr = np.array(fr)
check("carrier split", near(fr[:, 1].min(), 92.5, 1) and near(fr[:, 1].max(), 92.7, 1) and near(fr[:, 2].min(), 7.2, 1)
      and near(fr[:, 2].max(), 7.3, 1) and near(fr[:, 0].min(), 0.1, 1) and near(fr[:, 0].max(), 0.3, 1),
      "gas %.1f-%.1f%%, member galaxies %.1f-%.1f%%, BCG + companions %.1f-%.1f%%" % (
          fr[:, 1].min(), fr[:, 1].max(), fr[:, 2].min(), fr[:, 2].max(), fr[:, 0].min(), fr[:, 0].max()))

# ---------------------------------------------------------------------------
print("6. Gas treatment: Mistele et al.'s 1/r^4 alternative (section 5, Appendix A)")
cB = [0.0, 0.0]
for row in model.SEVEN:
    spec = B[row[0]].spec
    RB, MB, *_ = model.carrier_with_tail(spec, famaey, rmax_kpc[row[0]])
    s = model.score(model.assemble(spec, row[6], RB, MB, {}))
    cB[0] += s["H1_chi2"]
    cB[1] += s["MOND_chi2"]
check("1/r^4 tail reverses the ranking", near(cB[0] / 44, 1.963, 3) and near(cB[1] / 44, 1.450, 3),
      "Hermes 1 %.3f, MOND %.3f (chi2/dof, 44 dof)" % (cB[0] / 44, cB[1] / 44))

# ---------------------------------------------------------------------------
print("7. How much the per-cluster amplitude buys (section 5)")
D = {}
for bid in SEVEN:
    b = B[bid]
    D[bid] = dict(t=b.lens - b.mbar, H=b.mbar * model.hermes_q(b.phi, b.g), M=b.mbar * model.mond_q(b.g),
                  P=np.linalg.solve(b.cov, np.eye(len(b.cov))), n=len(b.r))


def chi(bid, k, K):
    r = D[bid]["t"] - K * D[bid][k]
    return float(r @ D[bid]["P"] @ r)


def kfit(bs, k):
    num = sum(float(D[b][k] @ D[b]["P"] @ D[b]["t"]) for b in bs)
    den = sum(float(D[b][k] @ D[b]["P"] @ D[b][k]) for b in bs)
    return max(0.0, num / den)


N = sum(D[b]["n"] for b in SEVEN)
one = {k: sum(chi(b, k, kfit(SEVEN, k)) for b in SEVEN) / (N - 1) for k in "HM"}
check("one K for all clusters (K about 4; Hermes 1 still fits better)",
      near(one["H"], 1.791, 3) and near(one["M"], 2.195, 3) and near(kfit(SEVEN, "H"), 4.0, 1) and one["H"] < one["M"],
      "Hermes 1 %.3f (K %.2f), MOND %.3f (K %.2f)" % (one["H"], kfit(SEVEN, "H"), one["M"], kfit(SEVEN, "M")))
loo = {k: sum(chi(b, k, kfit([o for o in SEVEN if o != b], k)) for b in SEVEN) / N for k in "HM"}
check("each cluster predicted with K from the other six", near(loo["H"], 1.92, 2) and near(loo["M"], 2.22, 2)
      and loo["H"] < loo["M"],
      "Hermes 1 %.3f, MOND %.3f (chi2/dof, %d dof)" % (loo["H"], loo["M"], N))
cal = ["A611", "MACSJ0429", "MACSJ1720"]
test = [b for b in SEVEN if b not in cal]
nt = sum(D[b]["n"] for b in test)
blind = {k: sum(chi(b, k, kfit(cal, k)) for b in test) / nt for k in "HM"}
check("K from A611, MACSJ0429, MACSJ1720 applied blind to the other four",
      near(blind["H"], 3.05, 2) and near(blind["M"], 2.62, 2),
      "Hermes 1 %.2f (K %.2f), MOND %.2f (K %.2f) (chi2/dof, %d dof)" % (blind["H"], kfit(cal, "H"), blind["M"],
                                                                         kfit(cal, "M"), nt))
hp = {}
for k in "HM":
    K = np.array([float(D[b][k] @ D[b]["P"] @ D[b]["t"]) / float(D[b][k] @ D[b]["P"] @ D[b][k]) for b in SEVEN])
    w = np.array([float(D[b][k] @ D[b]["P"] @ D[b][k]) for b in SEVEN])
    Kbar = float(np.sum(w * K) / np.sum(w))
    hp[k] = C2.sf(float(np.sum(w * (K - Kbar) ** 2)), len(K) - 1)
check("does K vary between clusters more than its errors allow?", near(hp["H"], 0.005, 3) and near(hp["M"], 0.24, 2),
      "homogeneity p: Hermes 1 %.4f, MOND %.3f" % (hp["H"], hp["M"]))
k1 = {k: sum(chi(b, k, 1.0) for b in SEVEN) / N for k in "HM"}
print("  (for reference, K = 1, the equations as written: Hermes 1 %.2f, MOND %.2f chi2/dof)" % (k1["H"], k1["M"]))

# ---------------------------------------------------------------------------
print("8. Baryon closure (section 5)")


def k_all(f):
    """Per-cluster K with every baryonic mass scaled by f (the gate rebuilt on the scaled carrier)."""
    out = []
    for bid in SEVEN:
        b = B[bid]
        ph = np.interp(np.log(b.r), np.log(b.R), model.phi_chain(b.R, model.G * f * b.M / b.R ** 2))
        q = model.hermes_q(ph, model.G * f * b.mbar / b.r ** 2)
        t, d, P = b.lens - f * b.mbar, f * b.mbar * q, D[bid]["P"]
        out.append(max(0.0, float(d @ P @ t / float(d @ P @ d))))
    return np.array(out)


m1 = float(np.median(k_all(1.0)))
stellar = float(np.mean(fr[:, 0] + fr[:, 2])) / 100.0
f_lo, f_hi = (1.0 + (s - stellar) / (1.0 - s) for s in STELLAR_BENCHMARK)
d_lo, d_hi = (100.0 * (1.0 - float(np.median(k_all(f))) / m1) for f in (f_lo, f_hi))
check("complete baryonic inventories would reduce K by about 3-8%", round(d_lo) == 3 and round(d_hi) == 8,
      "carrier stellar share %.1f%%; a %.0f-%.0f%% benchmark means f = %.3f-%.3f; median K falls %.1f-%.1f%%" % (
          100 * stellar, 100 * STELLAR_BENCHMARK[0], 100 * STELLAR_BENCHMARK[1], f_lo, f_hi, d_lo, d_hi))
lo, hi = 1.0, 8.0  # beta changes sign near f = 9.6, so the search stays below it
for _ in range(60):
    mid = 0.5 * (lo + hi)
    lo, hi = (mid, hi) if np.median(k_all(mid)) > 1.0 else (lo, mid)
f1 = hi
check("K = 1 would need about 74% of the baryons missing", near(f1, 3.862, 3) and round(100 * (1 - 1 / f1)) == 74,
      "median K = 1 at f = %.3f, i.e. %.1f%% missing" % (f1, 100 * (1 - 1 / f1)))
fb = np.array([f1 * B[bid].mbar[int(np.argmax(B[bid].r))] / B[bid].lens[int(np.argmax(B[bid].r))] for bid in SEVEN])
check("outermost-node baryon fraction at K = 1", near(fb.min(), 0.44, 2) and near(fb.max(), 0.76, 2),
      "%.2f-%.2f, %.1f-%.1f times the cosmic %.3f" % (fb.min(), fb.max(), fb.min() / FB_COSMIC, fb.max() / FB_COSMIC,
                                                       FB_COSMIC))

# ---------------------------------------------------------------------------
det = []
for row in model.SEVEN + model.EXCLUDED:
    b = B[row[0]]
    s = model.score(b)
    d = len(b.r) - 1
    det.append(dict(board=row[6], track=row[4], nodes=len(b.r), r_min_kpc=b.r.min(), r_max_kpc=b.r.max(),
                    g_bar_min=b.g.min(), g_bar_max=b.g.max(), rmax_x_kpc=rmax_kpc[row[0]],
                    nodes_beyond_rmax_x=int(np.sum(b.r > rmax_kpc[row[0]])), gas_model=b.meta["gas_model"],
                    H1_chi2_per_dof=s["H1_chi2"] / d, H1_K=s["H1_K"], MOND_chi2_per_dof=s["MOND_chi2"] / d,
                    MOND_K=s["MOND_K"]))
pd.DataFrame(det).to_csv(OUT / "boards.csv", index=False, float_format="%.4f")
passed = sum(ok for ok, _, _ in CHECKS)
lines = ["%s  %s: %s" % ("PASS" if ok else "FAIL", name, detail) for ok, name, detail in CHECKS]
lines.append("%d of %d checks passed" % (passed, len(CHECKS)))
(OUT / "checks.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
print()
print(lines[-1])
sys.exit(0 if passed == len(CHECKS) else 1)
