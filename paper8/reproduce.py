#!/usr/bin/env python3
"""Rebuild the Paper 8 numbers from public inputs and check them.

    python paper8/reproduce.py --sparc /path/to/rotmod_dir      (or set SPARC_DIR)

Inputs, all public:
  SPARC rotation curves (Lelli, McGaugh & Schombert 2016; Zenodo 16284118), not
    redistributed here; pass the folder of *_rotmod.dat files
  paper1/ages_133.csv                      t50 and g98 for the 133 galaxies
  paper1/Hermes_ConfigG_PerGalaxy_133_Export.csv   r_knee, rho and the Chimera flag
  paper1/hermes_gate_phi.py                the log-grid gate that Paper 8 uses
  paper9/hermes_clusters/gate.py           the chain-rule gate of Paper 1's published scores
  paper7/01b_hermes_result_canonical_gate/ the M33 board and its gate as corrected in Paper 7 v3
  paper8/verify_eta.py                     the operator
  paper8/expected/                         the paper's ten production CSV files (first public here)

Prints one line per check with the place the number appears in the paper, writes
results/ (recomputed per-galaxy table and the three M33 tables), compares them with
expected/, and exits 1 on any FAIL.

Conventions (each one changes a checked value if altered):
  gate      log-grid shear, paper1/hermes_gate_phi.py (the chain-rule gate is checked separately)
  velocity  V = sqrt(max(g_model, 0) R), as in Paper 1
  MOND      simple interpolating function, a0 = 1.2e-10 m/s^2 with the exact parsec
            (3702.813 (km/s)^2/kpc)
  chi2/N    mean over points of (Vobs - V)^2 / (errV^2 + 386)
  outer-5   the last five radii of each rotmod file; residual = V_model - V_obs
  M33       the canonical gate of Paper 7 v3, g98 = 98th percentile of the board's g_bar,
            inner disk R <= 10 kpc at the global age, outer disk R > 10 kpc at the
            outer age; W always uses the global age

Not rebuilt: candidates 1 and 2 of Section 3, whose definitions are not specific enough
to rerun; their table cells are quoted from the original runs on the superseded M33 gate.
"""
import argparse, csv, glob, hashlib, os, re, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(ROOT, "paper1"))
sys.path.insert(0, os.path.join(ROOT, "paper9"))
import verify_eta as ve                                   # noqa: E402
from hermes_gate_phi import hermes_phi                    # noqa: E402  (log-grid)
from hermes_clusters.gate import phi_chain                # noqa: E402  (chain-rule)

PI = np.pi
F = 1.0 / np.sqrt(2.0 * PI)
AK, AU, KDIV, SIG2 = ve.A_KNEE, ve.A_U, ve.K_DIVISOR, ve.SIGMA_INT_SQ
A0 = ve.A0_MOND
EXP = os.path.join(HERE, "expected")
OUT = os.path.join(HERE, "results")
M33_PP = os.path.join(ROOT, "paper7", "01b_hermes_result_canonical_gate",
                      "m33_hermes_per_point_predictions_canonical.csv")
RESULTS = []


def check(cid, where, what, got, want, tol):
    ok = bool(abs(float(got) - float(want)) <= tol)
    RESULTS.append((cid, where, what, float(got), float(want), tol, ok))
    print("  [%s] %-7s %-22s %-52s got %-13.7g want %-11.7g" % ("PASS" if ok else "FAIL", cid, where, what, got, want))


def check_true(cid, where, what, cond, detail=""):
    RESULTS.append((cid, where, what, float(bool(cond)), 1.0, 0, bool(cond)))
    print("  [%s] %-7s %-22s %s %s" % ("PASS" if cond else "FAIL", cid, where, what, detail))


def norm(s):
    s = re.sub(r"[^a-z0-9]", "", s.lower())
    m = re.match(r"^([a-z]+)0*([0-9].*)$", s)
    return m.group(1) + m.group(2) if m else s


def vel(g, R):
    return np.sqrt(np.maximum(g, 0.0) * R)


def chi(Vo, V, e):
    return float(np.mean((Vo - V) ** 2 / (e ** 2 + SIG2)))


def mond(gb):
    x = np.maximum(np.abs(gb) / A0, 1e-30)
    return gb * (1.0 + np.sqrt(1.0 + 4.0 / x)) / 2.0


# ---------------------------------------------------------------- operator variants
def eta_raw(gb, au=AU, use_S=True, use_G=True):
    """1 + S (max(1, sqrt(a_u / Gamma)) - 1) before the pi cap."""
    G = ve.log_compress(gb) if use_G else np.maximum(gb, 0.0)
    with np.errstate(divide="ignore", invalid="ignore"):
        contrast = np.where(G > 0, np.sqrt(au / np.where(G > 0, G, 1.0)), 1.0)
    S = ve.spatial_transition(gb) if use_S else 1.0
    return 1.0 + S * (np.maximum(1.0, contrast) - 1.0)


def g_eta(gb, phi, psi_R, psi_sys, au=AU, k=ve.WEAR_RATE, use_S=True, use_G=True, cap=True, use_W=True):
    eu = eta_raw(gb, au, use_S, use_G)
    if cap:
        eu = np.minimum(PI, eu)
    W = (1.0 - np.exp(-k * psi_sys)) if use_W else 1.0
    return gb * (1 + phi * (PI * np.exp(-psi_R) * (1 + W * (eu - 1)) - F))


def g_base(gb, phi, psi):
    return gb * (1 + phi * (PI * np.exp(-psi) - F))


def g_c4(gb, phi, psi_R, psi_sys):                        # candidate 4: Supplement A formula
    W, G = ve.wear_activation(psi_sys), ve.log_compress(gb)
    with np.errstate(divide="ignore"):
        P = np.exp(-psi_R) * np.sqrt(PI ** 2 + W ** 2 * PI * AK / G)
    return gb * (1 + phi * (np.minimum(PI, np.where(np.isfinite(P), P, PI)) - F))


def g_c5(gb, phi, psi_R, psi_sys, width=1 / PI):          # candidate 5 (the archived run used width 1/pi)
    W, G = ve.wear_activation(psi_sys), ve.log_compress(gb)
    with np.errstate(divide="ignore"):
        M = np.minimum(PI, np.sqrt(1 + AU / G))
    eta = 1 + W * ve.spatial_transition(gb, width=width) * (np.where(np.isfinite(M), M, PI) - 1)
    return gb * (1 + phi * (PI * np.exp(-psi_R) * eta - F))


# ---------------------------------------------------------------- SPARC
def load_sparc(sparc_dir):
    rot = {norm(os.path.basename(f)[:-len("_rotmod.dat")]): f
           for f in glob.glob(os.path.join(sparc_dir, "*_rotmod.dat"))}
    ex = {norm(r["galaxy"]): r for r in csv.DictReader(open(os.path.join(ROOT, "paper1",
                                                       "Hermes_ConfigG_PerGalaxy_133_Export.csv")))}
    gals = []
    for r in csv.DictReader(open(os.path.join(ROOT, "paper1", "ages_133.csv"), encoding="utf-8-sig")):
        k = norm(r["galaxy"])
        if k not in rot:
            sys.exit("missing SPARC file for %s in %s (pass --sparc or set SPARC_DIR)" % (r["galaxy"], sparc_dir))
        a = np.loadtxt(rot[k], comments="#")
        gals.append(dict(name=r["galaxy"], t50=float(r["t50_gyr"]), g98=float(r["g98"]), ex=ex[k],
                         R=a[:, 0], Vo=a[:, 1], e=a[:, 2], Vg=a[:, 3], Vd=a[:, 4], Vb=a[:, 5], SBd=a[:, 6]))
    return gals


def score_sparc(gals):
    out = []
    for G in gals:
        R, Vo, e, Vg, Vd, Vb = G["R"], G["Vo"], G["e"], G["Vg"], G["Vd"], G["Vb"]
        gb = (Vd ** 2 + Vb ** 2 + np.sign(Vg) * Vg ** 2) / R
        phi, phic = hermes_phi(R, gb), phi_chain(R, gb)
        psi = ve.compute_psi(G["t50"], G["g98"])
        gB, gE, gM = g_base(gb, phi, psi), ve.compute_g_model_eta(gb, phi, psi, psi), mond(gb)
        VB, VE, VM = vel(gB, R), vel(gE, R), vel(gM, R)
        # stellar mass-to-light 0.5 (disks) / 0.7 (bulges), g98 recomputed from the rescaled baryons
        gy = (0.5 * Vd ** 2 + 0.7 * Vb ** 2 + np.sign(Vg) * Vg ** 2) / R
        phiy, psiy = hermes_phi(R, gy), ve.compute_psi(G["t50"], np.percentile(gy, 98))
        x = G["ex"]
        rec = dict(galaxy=G["name"], n=len(R), psi=psi, W=float(ve.wear_activation(psi)), t50=G["t50"], g98=G["g98"],
                   base=chi(Vo, VB, e), eta=chi(Vo, VE, e), mond=chi(Vo, VM, e),
                   mond3700=chi(Vo, vel(gb * (1.0 + np.sqrt(1.0 + 4.0 / np.maximum(np.abs(gb) / 3700.0, 1e-30))) / 2.0, R), e),
                   c3=chi(Vo, vel(g_eta(gb, phi, psi, psi, use_W=False), R), e),
                   c4=chi(Vo, vel(g_c4(gb, phi, psi, psi), R), e),
                   c5=chi(Vo, vel(g_c5(gb, phi, psi, psi), R), e),
                   c5w30=chi(Vo, vel(g_c5(gb, phi, psi, psi, width=0.30), R), e),
                   au_e=chi(Vo, vel(g_eta(gb, phi, psi, psi, au=AK / np.e), R), e),
                   k_8pi=chi(Vo, vel(g_eta(gb, phi, psi, psi, k=8 * PI), R), e),
                   k_4pi2=chi(Vo, vel(g_eta(gb, phi, psi, psi, k=4 * PI ** 2), R), e),
                   noS=chi(Vo, vel(g_eta(gb, phi, psi, psi, use_S=False), R), e),
                   noG=chi(Vo, vel(g_eta(gb, phi, psi, psi, use_G=False), R), e),
                   nocap=chi(Vo, vel(g_eta(gb, phi, psi, psi, cap=False), R), e),
                   base_chain=chi(Vo, vel(g_base(gb, phic, psi), R), e),
                   eta_chain=chi(Vo, vel(ve.compute_g_model_eta(gb, phic, psi, psi), R), e),
                   base_y=chi(Vo, vel(g_base(gy, phiy, psiy), R), e),
                   eta_y=chi(Vo, vel(ve.compute_g_model_eta(gy, phiy, psiy, psiy), R), e),
                   mond_y=chi(Vo, vel(mond(gy), R), e),
                   ncap=int((eta_raw(gb) > PI).sum()),
                   mB=float(np.median(np.abs(VB - Vo))), mE=float(np.median(np.abs(VE - Vo))),
                   mM=float(np.median(np.abs(VM - Vo))),
                   o5B=float(np.mean((VB - Vo)[-5:])), o5E=float(np.mean((VE - Vo)[-5:])),
                   o3B=float(np.mean((VB - Vo)[-3:])), o3E=float(np.mean((VE - Vo)[-3:])),
                   rl=float(R[-1] / float(x["r_knee_kpc"])), rho=float(x["rho"]),
                   chim=x["chimera_flag"].strip().lower() == "true",
                   gas=float(np.mean(Vg[-5:] ** 2 / (Vg[-5:] ** 2 + Vd[-5:] ** 2 + Vb[-5:] ** 2))),
                   sb=float(np.median(G["SBd"][G["SBd"] > 0])),
                   vfun=max(float(np.abs(ve.compute_v_model_eta(gb, phi, psi, psi, R) - VE).max()),
                            float(np.abs(ve.compute_g_model_baseline(gb, phi, psi) - gB).max()),
                            float(np.abs(ve.compute_g_mond_simple(gb) - gM).max()),
                            abs(ve.compute_chi2nu(Vo, VE, e) - chi(Vo, VE, e))),
                   _dB=VB - Vo, _dE=VE - Vo, _dM=VM - Vo, _VB=VB, _VE=VE, _VM=VM, _gb=gb)
        out.append(rec)
    return out


# ---------------------------------------------------------------- M33
M33_MASKS = [("full_board_all_R", lambda R: R > 0), ("inner_R_le_10", lambda R: R <= 10),
             ("outer_R_gt_10", lambda R: R > 10), ("far_outer_R_ge_15", lambda R: R >= 15)]
M33_CONFIGS = [(6.0, 6.0, "uniform"), (6.0, 3.0, "outer_3_Gyr"), (6.0, 2.0, "outer_2_Gyr"),
               (6.0, 1.5, "outer_1.5_Gyr"), (6.0, 1.0, "outer_1_Gyr"),
               (7.1, 7.1, "uniform"), (7.1, 3.0, "outer_3_Gyr"), (7.1, 2.0, "outer_2_Gyr"),
               (7.1, 1.5, "outer_1.5_Gyr"), (7.1, 1.0, "outer_1_Gyr")]


def load_m33():
    rows = list(csv.DictReader(open(M33_PP)))
    c = lambda n: np.array([float(r[n]) for r in rows])
    R, Vo, e = c("R_kpc"), c("Vobs_kms"), c("errV_kms")
    gb = (c("Vdisk_kms") ** 2 + c("Vbulge_kms") ** 2 + np.sign(c("Vgas_kms")) * c("Vgas_kms") ** 2) / R
    return dict(R=R, Vo=Vo, e=e, gb=gb, Vbar=np.sqrt(np.maximum(gb * R, 0)), g98=float(np.percentile(gb, 98)),
                phi=hermes_phi(R, gb), phi_file=c("phi_canonical"))


def m33_model(M, t_in, t_out, kind="eta"):
    R, gb, phi = M["R"], M["gb"], M["phi"]
    psi_sys = ve.compute_psi(t_in, M["g98"])
    psi_R = np.where(R <= 10, psi_sys, ve.compute_psi(t_out, M["g98"]))
    if kind == "eta":
        g = ve.compute_g_model_eta(gb, phi, psi_R, psi_sys)
    elif kind == "base":
        g = g_base(gb, phi, psi_R)
    elif kind == "c3":
        g = g_eta(gb, phi, psi_R, psi_sys, use_W=False)
    elif kind == "c4":
        g = g_c4(gb, phi, psi_R, psi_sys)
    elif kind == "c5":
        g = g_c5(gb, phi, psi_R, psi_sys)
    elif kind == "mond":
        g = mond(gb)
    return vel(g, R)


def m33_stats(M, V):
    q = (V - M["Vo"]) ** 2 / (M["e"] ** 2 + SIG2)
    return {name: float(q[f(M["R"])].mean()) for name, f in M33_MASKS}, float(q[M["R"] > 10].sum() / q.sum())


def write_m33_tables(M, outdir):
    """The paper's three M33 tables, in the schema of the production files."""
    R, Vo, e = M["R"], M["Vo"], M["e"]
    den = e ** 2 + SIG2
    summ, longr, rad = [], [], []
    for ti, to, mode in M33_CONFIGS:
        V = m33_model(M, ti, to)
        res = V - Vo
        ci = res ** 2 / den
        for name, f in M33_MASKS:
            m = f(R)
            longr.append(["wear_activated_eta", ti, to, mode, name, int(m.sum()), float(R[m].min()), float(R[m].max()),
                          float(ci[m].mean()), float(ci[m].sum()), float(res[m].mean()),
                          float(np.median(np.abs(res[m]))), float(np.sqrt(np.mean(res[m] ** 2)))])
        summ.append(["wear_activated_eta", ti, to, mode, float(ci.mean()), float(ci[R <= 10].mean()),
                     float(ci[R > 10].mean()), float(ci[R >= 15].mean())])
        for i in range(len(R)):
            vm = round(float(V[i]), 3)
            rad.append([ti, to, mode, R[i], Vo[i], e[i], round(float(M["Vbar"][i]), 3), vm, vm - Vo[i],
                        (vm - Vo[i]) ** 2 / den[i]])
    summ.sort(key=lambda r: (r[1], r[2]))
    files = {
        "paper8_m33_wear_eta_component_age_summary.csv":
            (["model", "inner_t50_Gyr", "outer_t50_Gyr", "age_mode", "full_chi2nu", "inner_R_le_10_chi2nu",
              "outer_R_gt_10_chi2nu", "far_outer_R_ge_15_chi2nu"], summ),
        "paper8_m33_wear_eta_component_age_summary_long.csv":
            (["model", "inner_t50_Gyr", "outer_t50_Gyr", "age_mode", "mask", "N", "R_min_kpc", "R_max_kpc", "chi2nu",
              "chi2_sum", "mean_residual_kms", "median_abs_residual_kms", "rms_residual_kms"], longr),
        "paper8_m33_wear_eta_component_age_radial_profile.csv":
            (["inner_t50_Gyr", "outer_t50_Gyr", "age_mode", "R_kpc", "Vobs_kms", "errV_kms", "Vbar_kms",
              "V_model_wear_eta_kms", "residual_model_minus_obs_kms", "chi2_contribution"], rad),
    }
    for fn, (head, rows) in files.items():
        with open(os.path.join(outdir, fn), "w", newline="") as fh:
            w = csv.writer(fh, lineterminator="\n")
            w.writerow(head)
            w.writerows(rows)
    return list(files)


# ---------------------------------------------------------------- file comparisons
def read_csv(path):
    return list(csv.DictReader(open(path, encoding="utf-8-sig")))


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


def manifest_ok():
    path = os.path.join(HERE, "MANIFEST_SHA256.txt")
    bad = []
    for line in open(path):
        if not line.strip() or line.startswith("#"):
            continue
        h, rel = line.split(None, 1)
        rel = rel.strip()
        if not rel.startswith("expected/"):
            continue
        got = hashlib.sha256(open(os.path.join(HERE, rel), "rb").read()).hexdigest()
        if got != h:
            bad.append(rel)
    return bad


# ---------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--sparc", default=os.environ.get("SPARC_DIR", ""),
                    help="folder of SPARC *_rotmod.dat files (or set SPARC_DIR)")
    sparc_dir = ap.parse_args().sparc
    if not sparc_dir or not glob.glob(os.path.join(sparc_dir, "*_rotmod.dat")):
        sys.exit("no SPARC *_rotmod.dat files in %r; pass --sparc /path/to/rotmod_dir or set SPARC_DIR" % sparc_dir)
    os.makedirs(OUT, exist_ok=True)

    print("0. Shipped data")
    bad = manifest_ok()
    check_true("H.1", "MANIFEST_SHA256", "expected/ files match their SHA-256", not bad, ", ".join(bad))

    G = score_sparc(load_sparc(sparc_dir))
    col = lambda k: np.array([g[k] for g in G])
    B, E, M = col("base"), col("eta"), col("mond")
    imp, cas = (E - B) < -0.5, (E - B) > 0.5
    check("S.0", "Sec 5.1", "galaxies scored", len(G), 133, 0)
    check("S.1", "verify_eta.py", "its V, baseline, MOND and chi2/N functions agree (max |diff|)", col("vfun").max(), 0, 1e-12)

    print("1. Table 5.1 and Section 5.1 (SPARC, unit stellar M/L, log-grid gate)")
    check("T.1", "Table 5.1", "baseline median chi2/N", np.median(B), 1.322593, 5e-7)
    check("T.2", "Table 5.1", "MOND median chi2/N", np.median(M), 1.142696, 5e-7)
    check("T.3", "Table 5.1", "wear-activated eta median chi2/N", np.median(E), 1.188777, 5e-7)
    for i, (lab, x) in enumerate((("baseline", B), ("MOND", M), ("eta", E))):
        check("T.%d" % (4 + 2 * i), "Table 5.1", "%s chi2/N > 5" % lab, (x > 5).sum(), (21, 38, 21)[i], 0)
        check("T.%d" % (5 + 2 * i), "Table 5.1", "%s chi2/N > 10" % lab, (x > 10).sum(), (12, 15, 12)[i], 0)
    check("T.10", "Table 5.1", "MOND worsened > 0.5 vs baseline", (M - B > 0.5).sum(), 44, 0)
    check("T.11", "Table 5.1", "MOND improved > 0.5 vs baseline", (M - B < -0.5).sum(), 44, 0)
    check("T.12", "Table 5.1", "eta casualties (worsened > 0.5)", cas.sum(), 0, 0)
    check("T.13", "Table 5.1", "eta improved > 0.5", imp.sum(), 13, 0)
    check("T.14", "Table 5.1", "eta new chi2/N > 5", ((E > 5) & (B <= 5)).sum(), 0, 0)
    check("T.15", "Table 5.1", "eta new chi2/N > 10", ((E > 10) & (B <= 10)).sum(), 0, 0)
    check("T.16", "Table 5.1", "baseline beats MOND (68/133 = 51.1%)", (B < M).sum(), 68, 0)
    check("T.17", "Table 5.1", "eta beats MOND (77/133 = 57.9%)", (E < M).sum(), 77, 0)
    check("T.18", "Sec 5.1", "per-galaxy median |dV|, eta", np.median(col("mE")), 17.84, 0.005)
    check("T.19", "Sec 5.1", "per-galaxy median |dV|, baseline", np.median(col("mB")), 19.05, 0.005)
    check("T.20", "Sec 5.1", "per-galaxy median |dV|, MOND", np.median(col("mM")), 20.23, 0.005)
    pool = lambda k: np.median(np.abs(np.concatenate([g[k] for g in G])))
    check("T.21", "Sec 5.1", "points pooled", sum(g["n"] for g in G), 3073, 0)
    check("T.22", "Sec 5.1", "pooled median |dV|, eta", pool("_dE"), 18.83, 0.005)
    check("T.23", "Sec 5.1", "pooled median |dV|, baseline", pool("_dB"), 20.67, 0.005)
    check("T.24", "Sec 5.1", "pooled median |dV|, MOND", pool("_dM"), 25.82, 0.005)

    ex1 = {norm(r["galaxy"]): r for r in read_csv(os.path.join(ROOT, "paper1", "Hermes_ConfigG_PerGalaxy_133_Export.csv"))}
    p1m = np.array([float(ex1[norm(g["galaxy"])]["chi2nu_mond"]) for g in G])
    a3700 = np.array([g["mond3700"] for g in G])
    check("P.1", "Sec 5 note", "Paper 1 MOND column = a0 3700 (max |diff|)", np.abs(a3700 - p1m).max(), 0, 1e-12)
    check("P.2", "Sec 5 note", "Paper 1 published MOND median", np.median(p1m), 1.141, 0.0005)

    print("2. Recomputation note: the same comparison on Paper 1's chain-rule gate")
    BC, EC = col("base_chain"), col("eta_chain")
    p1c = np.array([float(ex1[norm(g["galaxy"])]["chi2nu_configg"]) for g in G])
    check("R.0", "Sec 5 note", "chain-rule baseline = Paper 1 chi2nu_configg (max |diff|)", np.abs(BC - p1c).max(), 0, 1e-11)
    check("R.1", "Sec 5 note", "chain-rule baseline median", np.median(BC), 1.311535, 5e-7)
    check("R.2", "Sec 5 note", "chain-rule eta median", np.median(EC), 1.150360, 5e-7)
    check("R.3", "Sec 5 note", "chain-rule eta casualties", (EC - BC > 0.5).sum(), 0, 0)
    check("R.4", "Sec 5 note", "chain-rule eta improved", (EC - BC < -0.5).sum(), 13, 0)
    check("R.5", "Sec 5 note", "chain-rule eta beats MOND", (EC < M).sum(), 76, 0)

    print("3. Stellar mass-to-light 0.5 (disks) / 0.7 (bulges), g98 recomputed")
    BY, EY, MY = col("base_y"), col("eta_y"), col("mond_y")
    check("Y.1", "Sec 5.1", "baseline median", np.median(BY), 1.8968, 5e-5)
    check("Y.2", "Sec 5.1", "eta median", np.median(EY), 1.2009, 5e-5)
    check("Y.3", "Sec 5.1", "MOND median", np.median(MY), 0.4784, 5e-5)
    check("Y.4", "Sec 5.1", "eta casualties", (EY - BY > 0.5).sum(), 0, 0)
    check("Y.5", "Sec 5.1", "eta improved", (EY - BY < -0.5).sum(), 40, 0)
    check("Y.6", "Sec 5.1", "eta chi2/N > 5", (EY > 5).sum(), 35, 0)
    check("Y.7", "Sec 5.1", "MOND chi2/N > 5", (MY > 5).sum(), 6, 0)
    check("Y.8", "Sec 5.1", "MOND chi2/N > 10 (Paper 9: one)", (MY > 10).sum(), 1, 0)

    print("4. Section 3 table and Section 4: candidates, choice of constants, ablations")
    C3, C4, C5, C5w = col("c3"), col("c4"), col("c5"), col("c5w30")
    worst = lambda X: G[int(np.argmax(X - B))]["galaxy"]
    check("C.1", "Table 3 row 3", "candidate 3 (no W) median", np.median(C3), 1.115237, 5e-7)
    check("C.2", "Table 3 row 3", "candidate 3 casualties", (C3 - B > 0.5).sum(), 1, 0)
    check_true("C.3", "Table 3 row 3", "candidate 3 casualty is UGC 05750", worst(C3) == "UGC 05750", worst(C3))
    check("C.4", "Table 3 row 4", "candidate 4 median", np.median(C4), 1.361115, 5e-7)
    check("C.5", "Table 3 row 4", "candidate 4 casualties", (C4 - B > 0.5).sum(), 15, 0)
    check("C.6", "Table 3 row 4", "candidate 4 new chi2/N > 5", ((C4 > 5) & (B <= 5)).sum(), 4, 0)
    check("C.7", "Table 3 row 5", "candidate 5 median (width 1/pi)", np.median(C5), 1.052640, 5e-7)
    check("C.8", "Table 3 row 5", "candidate 5 casualties", (C5 - B > 0.5).sum(), 14, 0)
    check("C.9", "Table 3 row 5", "candidate 5 new chi2/N > 5", ((C5 > 5) & (B <= 5)).sum(), 2, 0)
    check("C.10", "Supp A cand 5", "candidate 5 casualties at width 0.30", (C5w - B > 0.5).sum(), 14, 0)
    check("C.11", "Sec 4 (support)", "eta largest worsening (UGC 05750)", (E - B).max(), 0.370, 0.0005)
    check_true("C.12", "Sec 4", "eta largest worsening is UGC 05750", worst(E) == "UGC 05750", worst(E))
    for cid, k, med, where in (("C.13", "au_e", 1.0790, "a_u = a_k/e, 2 pi^2"), ("C.15", "k_8pi", 1.1540, "a_k/pi, W rate 8 pi")):
        X = col(k)
        check(cid, "Sec 4", where + ": median", np.median(X), med, 5e-5)
        check(cid[:-1] + str(int(cid[-1]) + 1), "Sec 4", where + ": casualties", (X - B > 0.5).sum(), 0, 0)
    X = col("k_4pi2")
    check("C.17", "Sec 4", "a_k/pi, W rate 4 pi^2: casualties", (X - B > 0.5).sum(), 1, 0)
    check_true("C.18", "Sec 4", "limiting galaxy for a_k/e, 8 pi, 4 pi^2 is UGC 05750",
               all(worst(col(k)) == "UGC 05750" for k in ("au_e", "k_8pi", "k_4pi2")))
    NS = col("noS")
    check("A.1", "Sec 6", "operator without S: median", np.median(NS), 1.161950, 5e-6)
    check("A.2", "Sec 6", "operator without S: casualties", (NS - B > 0.5).sum(), 0, 0)
    check("A.3", "Supp A (support)", "largest |change| from removing S", np.abs(NS - E).max(), 0.068, 0.0005)
    check("A.4", "Supp A", "operator without Gamma: casualties", (col("noG") - B > 0.5).sum(), 0, 0)
    check("A.5", "Supp A", "operator without the pi cap: casualties", (col("nocap") - B > 0.5).sum(), 0, 0)
    NC = col("ncap")
    check("A.6", "Supp E", "points where the pi cap binds", NC.sum(), 17, 0)
    check("A.7", "Supp E", "galaxies with a capped point", (NC > 0).sum(), 8, 0)
    lo, hi = 1.0, 200.0                                     # cap threshold with S(g): bisection
    for _ in range(100):
        mid = 0.5 * (lo + hi)
        lo, hi = (mid, hi) if eta_raw(np.array([mid]))[0] > PI else (lo, mid)
    check("A.8", "Supp E", "g_bar below which the cap binds", lo, 49.22, 0.005)
    check("A.9", "Supp E", "Gamma = a_u at g_bar = a_k (e^(1/pi) - 1)", AK * (np.exp(1 / PI) - 1), 594.06, 0.005)
    check("A.10", "Supp A", "S at g_bar = 594.06", ve.spatial_transition(np.array([AK * (np.exp(1 / PI) - 1)]))[0], 0.889, 0.0005)
    g594 = AK * (np.exp(1 / PI) - 1)
    check("A.15", "Supp A", "Gamma below g_bar at the onset (1 - Gamma/g)", 1 - ve.log_compress(np.array([g594]))[0] / g594, 0.151, 0.0005)
    W = col("W")
    check("A.11", "Supp A", "smallest W in the sample", W.min(), 0.224, 0.0005)
    check("A.12", "Supp A (support)", "galaxies with W < 0.5", (W < 0.5).sum(), 13, 0)
    lw = col("t50")[W < 0.5]
    check("A.13", "Supp A (support)", "youngest W < 0.5 galaxy (Gyr)", lw.min(), 3.6, 1e-9)
    check("A.14", "Supp A (support)", "oldest W < 0.5 galaxy (Gyr)", lw.max(), 7.3, 1e-9)

    print("5. Supplement B (outermost five points; residual = V_model - V_obs)")
    o5B, o5E = col("o5B"), col("o5E")
    flip = (o5B < 0) & (o5E > 0)
    check("B.1", "Supp B / Sec 5.1", "median outer-5 residual, baseline", np.median(o5B), -20.51, 0.005)
    check("B.2", "Supp B / Sec 5.1", "median outer-5 residual, eta", np.median(o5E), -10.80, 0.005)
    check("B.3", "Supp B / Sec 5.1", "under-to-over sign flips", flip.sum(), 16, 0)
    check("B.4", "Supp B", "over-to-under sign flips", ((o5B > 0) & (o5E < 0)).sum(), 0, 0)
    check("B.5", "Supp B / Sec 5.1", "flipped galaxies now above +5 km/s", (flip & (o5E > 5)).sum(), 5, 0)
    check_true("B.6", "Supp B", "those five are DDO 161, UGC 06983, UGC 12732, NGC 2366, UGC 05005",
               sorted(g["galaxy"] for g, f in zip(G, flip & (o5E > 5)) if f)
               == sorted(["DDO 161", "UGC 06983", "UGC 12732", "NGC 2366", "UGC 05005"]))
    check("B.7", "Supp B", "galaxies above +5 km/s after correction", (o5E > 5).sum(), 29, 0)
    check("B.8", "Supp B", "galaxies above +5 km/s at baseline", (o5B > 5).sum(), 22, 0)
    check_true("B.9", "Supp B", "every per-galaxy outer-5 shift is >= 0", (o5E - o5B).min() >= -1e-9, "min %.3g" % (o5E - o5B).min())

    print("6. Supplement C (13 improved galaxies vs the other 120)")
    for cid, k, a, b, tol in (("SC.1", "o5B", -41.83, -16.76, 0.005), ("SC.2", "rl", 6.95, 2.44, 0.005),
                              ("SC.3", "gas", 0.518, 0.204, 0.0005), ("SC.4", "sb", 1.1150, 14.9225, 0.0005),
                              ("SC.5", "rho", 0.036, 0.148, 0.0005)):
        x = col(k)
        check(cid + "a", "Supp C", "improved median " + k, np.median(x[imp]), a, tol)
        check(cid + "b", "Supp C", "others median " + k, np.median(x[~imp]), b, tol)
    check("SC.6", "Supp C", "improved with R_last/r_knee > 2", (col("rl")[imp] > 2).sum(), 12, 0)
    check("SC.7", "Supp C", "Chimeras (rho < 0.02) among improved", col("chim")[imp].sum(), 6, 0)
    check("SC.8", "Supp C", "Chimeras among the others", col("chim")[~imp].sum(), 21, 0)
    check_true("SC.9", "Supp C", "Chimera flag equals rho < 0.02 for all 133", np.array_equal(col("chim"), col("rho") < 0.02))
    for cid, name, b0, e0 in (("SC.10", "NGC 6674", 14.54, 13.43), ("SC.11", "NGC 0289", 7.87, 7.28)):
        g = next(g for g in G if g["galaxy"] == name)
        check(cid + "a", "Supp C (support)", name + " baseline chi2/N", g["base"], b0, 0.005)
        check(cid + "b", "Supp C (support)", name + " eta chi2/N", g["eta"], e0, 0.005)

    print("7. M33 (Paper 7 v3 canonical gate)")
    MM = load_m33()
    check("M.0a", "Sec 5.2", "gate equals Paper 7 v3 phi_canonical (max |diff|)", np.abs(MM["phi"] - MM["phi_file"]).max(), 0, 1e-12)
    check("M.0b", "Sec 5.2", "board g98", MM["g98"], 2653.949886626496, 1e-9)
    st, share = m33_stats(MM, m33_model(MM, 6.0, 6.0, "base"))
    check("M.1", "Sec 1", "baseline R <= 10 (Paper 7 v3: 0.93)", st["inner_R_le_10"], 0.9292, 5e-5)
    check("M.2", "Sec 1 / Table 5.2", "baseline R > 10 (Paper 7 v3: 5.86)", st["outer_R_gt_10"], 5.8584, 5e-5)
    check("M.3", "Sec 1", "baseline outer share of chi2 (~85%)", share, 0.8459, 5e-5)
    check("M.4", "Table 5.2", "eta alone R > 10", m33_stats(MM, m33_model(MM, 6.0, 6.0))[0]["outer_R_gt_10"], 2.5972, 5e-5)
    check("M.5", "Table 5.2 / Supp D", "component age alone R > 10 (outer 0 Gyr)",
          m33_stats(MM, m33_model(MM, 6.0, 0.0, "base"))[0]["outer_R_gt_10"], 4.1234, 5e-5)
    check("M.6", "Sec 5.2", "eta + outer 0 Gyr R > 10", m33_stats(MM, m33_model(MM, 6.0, 0.0))[0]["outer_R_gt_10"], 1.0773, 5e-5)
    check("M.7", "Table 3 row 3", "candidate 3 R > 10 (uniform 6.0)", m33_stats(MM, m33_model(MM, 6.0, 6.0, "c3"))[0]["outer_R_gt_10"], 2.5944, 5e-5)
    check("M.8", "Table 3 row 5", "candidate 5 R > 10 (outer 2.0)", m33_stats(MM, m33_model(MM, 6.0, 2.0, "c5"))[0]["outer_R_gt_10"], 1.1036, 5e-5)
    ages = [(ti, to) for ti, to, _ in M33_CONFIGS] + [(6.0, 0.0), (7.1, 0.0)]
    c4o = [m33_stats(MM, m33_model(MM, ti, to, "c4"))[0]["outer_R_gt_10"] for ti, to in ages]
    check("M.9", "Table 3 row 4", "candidate 4 R > 10 (uniform 6.0)", c4o[0], 4.1234, 5e-5)
    check("M.9b", "Table 3 row 4", "candidate 4 has no age axis: R > 10 spread, 12 age pairs", max(c4o) - min(c4o), 0, 1e-12)
    outer = MM["R"] > 10
    under = [int((m33_model(MM, ti, to)[outer] < MM["Vo"][outer]).sum()) for ti, to in ages]
    check("M.10", "Sec 5.2", "outer points under-predicted at every age tested (of 27)", min(under), 27, 0)
    check("M.10b", "Sec 5.2", "outer points (R > 10 kpc)", outer.sum(), 27, 0)
    sm = m33_stats(MM, m33_model(MM, 6.0, 6.0, "mond"))[0]
    check("M.11", "Sec 5.2 / Supp D", "MOND R > 10 on the same board", sm["outer_R_gt_10"], 0.228, 0.0005)
    check("M.12", "Sec 5.2 / Supp D", "MOND full board (Paper 7 v3: 0.12)", sm["full_board_all_R"], 0.12, 0.005)
    want = {  # Supplement D tables: full / inner / outer / far-outer, to the printed 3 dp
        (6.0, 6.0): (1.682, 0.885, 2.597, 3.150), (6.0, 3.0): (1.303, 0.885, 1.782, 2.290),
        (6.0, 2.0): (1.186, 0.885, 1.532, 2.018), (6.0, 1.5): (1.131, 0.885, 1.413, 1.886),
        (6.0, 1.0): (1.077, 0.885, 1.297, 1.756), (7.1, 7.1): (1.844, 0.914, 2.912, 3.474),
        (7.1, 3.0): (1.317, 0.914, 1.780, 2.288), (7.1, 2.0): (1.201, 0.914, 1.531, 2.016),
        (7.1, 1.5): (1.145, 0.914, 1.411, 1.884), (7.1, 1.0): (1.091, 0.914, 1.295, 1.754)}
    n = 0
    for (ti, to), vals in want.items():
        s = m33_stats(MM, m33_model(MM, ti, to))[0]
        got = (s["full_board_all_R"], s["inner_R_le_10"], s["outer_R_gt_10"], s["far_outer_R_ge_15"])
        for lab, gv, wv in zip(("full", "inner", "outer", "far"), got, vals):
            n += 1
            check("D.%d" % n, "Supp D", "t50 %.1f / outer %.1f Gyr, %s" % (ti, to, lab), round(gv, 3), wv, 1e-9)

    print("8. Recomputed tables against the shipped CSV files (expected/)")
    # per-galaxy
    exp = {norm(r["galaxy"]): r for r in read_csv(os.path.join(EXP, "paper8_sparc_133_model_comparison.csv"))}
    pairs = [("baseline_chi2nu", "base"), ("mond_chi2nu", "mond"), ("wear_activated_eta_chi2nu", "eta"),
             ("baseline_median_abs_residual_kms", "mB"), ("mond_median_abs_residual_kms", "mM"),
             ("wear_eta_median_abs_residual_kms", "mE"), ("outer5_signed_residual_baseline_kms", "o5B"),
             ("outer5_signed_residual_wear_eta_kms", "o5E"), ("psi", "psi"), ("wear_activation_W", "W")]
    worst_d = max(abs(float(exp[norm(g["galaxy"])][c]) - g[k]) for g in G for c, k in pairs)
    check("F.1", "Supp F", "per-galaxy table, 10 columns x 133 (max |diff|)", worst_d, 0, 1e-9)
    chi_d = max(abs(float(exp[norm(g["galaxy"])][c]) - g[k]) for g in G
                for c, k in (("baseline_chi2nu", "base"), ("wear_activated_eta_chi2nu", "eta")))
    check("F.0", "Sec 5 note", "log-grid baseline and eta chi2/N vs the per-galaxy table (max |diff|)", chi_d, 0, 1e-11)
    # per-radius
    rows = read_csv(os.path.join(EXP, "paper8_sparc_per_radius_baseline_mond_wear_eta.csv"))
    byg = {}
    for r in rows:
        byg.setdefault(norm(r["galaxy"]), []).append(r)
    wd = 0.0
    for g in G:
        rr = byg[norm(g["galaxy"])]
        for i, r in enumerate(rr):
            wd = max(wd, abs(float(r["V_baseline_kms"]) - g["_VB"][i]), abs(float(r["V_wear_activated_kms"]) - g["_VE"][i]),
                     abs(float(r["V_mond_kms"]) - g["_VM"][i]), abs(float(r["gbar_kms2_per_kpc"]) - g["_gb"][i]))
    check("F.2", "Supp F", "per-radius table, %d rows (max |dV| km/s)" % len(rows), wd, 0, 1e-9)
    # outer-tail audit (tail 3 and 5)
    aud = read_csv(os.path.join(EXP, "paper8_sparc_outer_signed_residual_audit_last3_last5.csv"))
    gi = {norm(g["galaxy"]): g for g in G}
    ad = max(max(abs(float(r["baseline_outer_mean_residual_kms"]) - gi[norm(r["galaxy"])]["o%sB" % r["tail_N"]]),
                 abs(float(r["corrected_outer_mean_residual_kms"]) - gi[norm(r["galaxy"])]["o%sE" % r["tail_N"]])) for r in aud)
    check("F.3", "Supp F", "outer-tail audit, %d rows (max |diff| km/s)" % len(aud), ad, 0, 1e-9)
    imp_set = {norm(g["galaxy"]) for g, f in zip(G, imp) if f}
    flip_set = {norm(g["galaxy"]) for g, f in zip(G, flip) if f}
    check_true("F.4", "Supp F", "improved-galaxy file lists the same 13 galaxies",
               {norm(r["galaxy"]) for r in read_csv(os.path.join(EXP, "paper8_sparc_13_improved_galaxies_characterization.csv"))} == imp_set)
    check_true("F.5", "Supp F", "sign-flip watchlist lists the same 16 galaxies",
               {norm(r["galaxy"]) for r in read_csv(os.path.join(EXP, "paper8_sparc_outer5_sign_flip_watchlist.csv"))} == flip_set)
    it = {r["model"]: r for r in read_csv(os.path.join(EXP, "paper8_iteration_comparison_metrics.csv"))}
    for cid, model, X in (("F.6", "baseline_hermes_recomputed", B), ("F.7", "mond_simple_recomputed", M),
                          ("F.8", "previous_BTFR_safe_positive_eta", C3), ("F.9", "wear_activated_eta_locked", E),
                          ("F.10", "failed_pythagorean_capped", C4), ("F.11", "rejected_phase_transition_FD", C5)):
        r = it[model]
        ok = (abs(float(r["median_chi2nu"]) - np.median(X)) < 1e-9 and int(r["chi2nu_gt5_count"]) == (X > 5).sum()
              and int(r["chi2nu_gt10_count"]) == (X > 10).sum()
              and int(r["worsened_gt0p5_vs_baseline_count"]) == (X - B > 0.5).sum()
              and int(r["improved_gt0p5_vs_baseline_count"]) == (X - B < -0.5).sum()
              and int(r["new_gt5_vs_baseline_count"]) == ((X > 5) & (B <= 5)).sum()
              and int(r["new_gt10_vs_baseline_count"]) == ((X > 10) & (B <= 10)).sum())
        check_true(cid, "Supp F", "iteration table row %s" % model, ok)
    ag = {r["model"]: r for r in read_csv(os.path.join(EXP, "paper8_sparc_aggregate_metrics.csv"))}
    ok = all(abs(float(ag[m]["median_chi2nu"]) - np.median(X)) < 1e-9 for m, X in
             (("baseline_hermes_recomputed", B), ("mond_simple_recomputed", M), ("wear_activated_eta", E)))
    check_true("F.12", "Supp F", "aggregate table medians", ok)
    rows_want = {"paper8_sparc_133_model_comparison.csv": 133, "paper8_sparc_per_radius_baseline_mond_wear_eta.csv": 3073,
                 "paper8_sparc_aggregate_metrics.csv": 3, "paper8_iteration_comparison_metrics.csv": 6,
                 "paper8_m33_wear_eta_component_age_summary.csv": 10, "paper8_m33_wear_eta_component_age_summary_long.csv": 40,
                 "paper8_m33_wear_eta_component_age_radial_profile.csv": 580,
                 "paper8_sparc_outer_signed_residual_audit_last3_last5.csv": 266,
                 "paper8_sparc_13_improved_galaxies_characterization.csv": 13, "paper8_sparc_outer5_sign_flip_watchlist.csv": 16}
    got_rows = {fn: len(read_csv(os.path.join(EXP, fn))) for fn in rows_want}
    check_true("F.16", "Supp F", "row counts of the ten CSV files", got_rows == rows_want,
               "" if got_rows == rows_want else str(got_rows))
    for i, fn in enumerate(write_m33_tables(MM, OUT)):
        same, why = cells_match(read_csv(os.path.join(OUT, fn)), read_csv(os.path.join(EXP, fn)))
        check_true("F.%d" % (13 + i), "Supp F", "M33 table %s rebuilt" % fn.replace("paper8_m33_wear_eta_component_age_", ""), same, why)

    with open(os.path.join(OUT, "paper8_per_galaxy_recomputed.csv"), "w", newline="") as fh:
        keys = [k for k in G[0] if not k.startswith("_")]
        w = csv.writer(fh, lineterminator="\n")
        w.writerow(keys)
        w.writerows([[g[k] for k in keys] for g in G])
    with open(os.path.join(OUT, "checks.csv"), "w", newline="") as fh:
        w = csv.writer(fh, lineterminator="\n")
        w.writerow(["id", "where", "quantity", "got", "expected", "tolerance", "pass"])
        w.writerows(RESULTS)
    n_fail = sum(1 for r in RESULTS if not r[-1])
    print("\n%d checks, %d failed; results in %s" % (len(RESULTS), n_fail, os.path.relpath(OUT, ROOT)))
    return 1 if n_fail else 0


if __name__ == "__main__":
    sys.exit(main())
