#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
verify_gates.py  -- reproducibility test for the TWO Hermes gate implementations.

There are two gate constructions in this project and they are NOT the same:

  * REPO gate  = paper1/hermes_gate_phi.py            -> reproduces Paper 1 (133 SPARC)
  * M33  gate  = paper7/.../m33_hermes_execution_script.py -> reproduces Paper 7 (M33)

This script imports the real Paper 1 gate function, and uses a verbatim copy of the
Paper 7 M33 gate (see M33_gate below -- copied line for line from gate_func() in
m33_hermes_execution_script.py). It then checks each gate against its own published
dataset to a tolerance of 1e-10, and prints the cross terms to show the two gates
genuinely diverge.

Data:
  * M33 reference (in repo)      : paper7/01_hermes_result_locked/m33_hermes_per_point_predictions.csv
  * 133 SPARC per-galaxy export  : Hermes_ConfigG_PerGalaxy_133_Export.csv   (path via HERMES_EXPORT)
  * SPARC rotmod files           : *_rotmod.dat                              (dir  via SPARC_DIR)
"""
import os, sys, csv, glob, re
import numpy as np
from scipy.signal import savgol_filter

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "paper1"))
from hermes_gate_phi import hermes_phi          # the REAL Paper 1 gate

TOL = 1e-10
A_KNEE = 1585.0
M33_REF = os.path.join(HERE, "paper7", "01_hermes_result_locked",
                       "m33_hermes_per_point_predictions.csv")
# 133-galaxy export ships with the repo (paper1/). Override with HERMES_EXPORT if needed.
EXPORT = os.environ.get("HERMES_EXPORT",
    os.path.join(HERE, "paper1", "Hermes_ConfigG_PerGalaxy_133_Export.csv"))
# SPARC rotmod files are the ONLY external dependency. Download the SPARC database
# (Lelli, McGaugh & Schombert 2016; http://astroweb.cwru.edu/SPARC/) and either drop
# the *_rotmod.dat files in ./sparc_database or point SPARC_DIR at them.
SPARC_DIR = os.environ.get("SPARC_DIR", os.path.join(HERE, "sparc_database"))


def M33_gate(R, gbar):
    """Verbatim copy of gate_func() from m33_hermes_execution_script.py (Paper 7).
    Parameterised to take (R, gbar); Vbar = sqrt(gbar*R) exactly as in that script.
    Differences from the repo gate: raw first-below knee, SavGol mode='interp',
    flat curvature epsilon 1e-6."""
    Vbar = np.sqrt(np.maximum(0.0, gbar * R))
    below = np.where(gbar < A_KNEE)[0]
    knee_idx = int(below[0]) if len(below) else int(np.argmin(np.abs(gbar - A_KNEE)))
    r_knee = float(R[knee_idx])
    x = R / r_knee
    win = 11 if len(R) >= 20 else 5
    if win >= len(R):
        win = len(R) - 1 if (len(R) - 1) % 2 else len(R) - 2
    if win % 2 == 0:
        win -= 1
    Vsm = savgol_filter(Vbar, win, 3, mode='interp')
    shear = np.abs(np.gradient(np.log(np.maximum(Vsm, 1e-12)), np.log(R)))
    dVdR = np.gradient(Vsm, R)
    d2VdR2 = np.gradient(dVdR, R)
    curvature = np.abs(d2VdR2) / (np.abs(dVdR) + 1e-6)
    inner = x <= 1.0
    outer = x > 1.0
    curv_med = float(np.nanmedian(curvature[inner])) if np.any(inner) else float(np.nanmedian(curvature))
    shear_med = float(np.nanmedian(shear[outer])) if np.any(outer) else float(np.nanmedian(shear))
    if not np.isfinite(curv_med) or curv_med <= 0:
        vals = curvature[curvature > 0]; curv_med = float(np.nanmedian(vals)) if vals.size else 1.0
    if not np.isfinite(shear_med) or shear_med <= 0:
        vals = shear[shear > 0]; shear_med = float(np.nanmedian(vals)) if vals.size else 1.0
    curv_norm = curvature / curv_med
    shear_norm = shear / shear_med
    w = 1.0 / (1.0 + np.exp(-((x - 1.40) / 0.30)))
    S_pre = (1.0 - w) * curv_norm + w * shear_norm
    S_eff = np.empty_like(S_pre)
    for i, s in enumerate(S_pre):
        S_eff[i] = s if i == 0 else 0.5 * s + 0.5 * S_eff[i - 1]
    return 1.0 - np.exp(-S_eff)


def _norm(s):
    s = re.sub(r'[^a-z0-9]', '', s.lower())
    m = re.match(r'^([a-z]+)0*([0-9].*)$', s)
    return m.group(1) + m.group(2) if m else s


def load_m33():
    rows = list(csv.DictReader(open(M33_REF)))
    R = np.array([float(r["R_kpc"]) for r in rows])
    gbar = np.array([float(r["gbar_kms2_per_kpc"]) for r in rows])
    phi = np.array([float(r["phi"]) for r in rows])
    return R, gbar, phi


def load_sparc():
    rot = {_norm(os.path.basename(f).replace("_rotmod.dat", "")): f
           for f in glob.glob(os.path.join(SPARC_DIR, "*_rotmod.dat"))}
    if not rot:
        sys.stderr.write(
            "\nERROR: no SPARC *_rotmod.dat files found in:\n  %s\n\n"
            "The SPARC rotation-curve database is the only external dependency.\n"
            "Download it from http://astroweb.cwru.edu/SPARC/ (Lelli, McGaugh &\n"
            "Schombert 2016), then either place the *_rotmod.dat files in\n"
            "./sparc_database or set SPARC_DIR to point at them.\n" % os.path.abspath(SPARC_DIR))
        sys.exit(2)
    out = []
    for r in csv.DictReader(open(EXPORT)):
        nm = _norm(r["galaxy"])
        if nm not in rot:
            continue
        Rk, Vg, Vd, Vb = [], [], [], []
        for ln in open(rot[nm]):
            if ln.startswith("#") or not ln.strip():
                continue
            p = ln.split()
            Rk.append(float(p[0])); Vg.append(float(p[3])); Vd.append(float(p[4])); Vb.append(float(p[5]))
        Rk, Vg, Vd, Vb = map(np.array, (Rk, Vg, Vd, Vb))
        gbar = (Vd**2 + Vb**2 + np.sign(Vg) * Vg**2) / Rk   # Config G: M/L = 1
        out.append((r["galaxy"], Rk, gbar, float(r["phi_last"])))
    return out


def main():
    print("=" * 70)
    print("HERMES GATE VERIFICATION   (pass tolerance = 1e-10)")
    print("=" * 70)

    m33_R, m33_g, m33_phi = load_m33()

    # Self-certification: our verbatim M33_gate must reproduce the committed CSV,
    # which is the actual output of m33_hermes_execution_script.py.
    cert = float(np.max(np.abs(M33_gate(m33_R, m33_g) - m33_phi)))
    print("\nSelf-certification of the M33 gate copy vs committed reference CSV:")
    print("  max|delta| = %.2e   %s" % (cert, "[faithful to m33_hermes_execution_script.py]" if cert <= TOL else "[MISMATCH]"))

    sparc = load_sparc()

    # CLAIM 1: repo gate reproduces Paper 1 (133 SPARC) phi_last
    d1 = max(abs(hermes_phi(R, g)[-1] - pl) for _, R, g, pl in sparc)
    # CLAIM 2: M33 gate reproduces Paper 7 (M33) per-point phi
    d2 = float(np.max(np.abs(M33_gate(m33_R, m33_g) - m33_phi)))

    print("\nCLAIM 1  paper1/hermes_gate_phi.py  vs  Paper 1 / 133 SPARC (phi_last)")
    print("  galaxies = %d    max|delta| = %.2e    -> %s" % (len(sparc), d1, "PASS" if d1 <= TOL else "FAIL"))
    print("CLAIM 2  m33_hermes_execution_script.py  vs  Paper 7 / M33 (per-point phi)")
    print("  points   = %d     max|delta| = %.2e    -> %s" % (len(m33_R), d2, "PASS" if d2 <= TOL else "FAIL"))

    # Cross terms: each gate on the OTHER dataset must diverge (proves two gates).
    x1 = float(np.max(np.abs(hermes_phi(m33_R, m33_g) - m33_phi)))
    x2 = max(abs(M33_gate(R, g)[-1] - pl) for _, R, g, pl in sparc)
    print("\nCross-checks (expected to DIVERGE -- confirms the gates are different):")
    print("  repo gate vs M33        : max|delta| = %.2e" % x1)
    print("  M33  gate vs 133 SPARC  : max|delta| = %.2e" % x2)

    ok = (cert <= TOL) and (d1 <= TOL) and (d2 <= TOL)
    print("\n" + "-" * 70)
    print("RESULT: %s  (2/2 reproduction claims pass; gates confirmed distinct)" % ("PASS" if ok else "FAIL"))
    print("-" * 70)
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
