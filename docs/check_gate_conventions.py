#!/usr/bin/env python3
"""Show which gate convention reproduces which published number.

    SPARC_DIR=sparc/sparc_database python docs/check_gate_conventions.py

The gate phi(R) exists in two numerical conventions that differ in ONE line, the
shear. Both are real and both appear in published artifacts:

    log-grid   s = |d lnV / d lnR|   evaluated with np.gradient on log coordinates
    chain-rule s = |(R/V) dV/dR|     evaluated with np.gradient on the native grid

They are the same quantity in the continuum and differ numerically, because
np.gradient on a non-uniform coordinate array is a finite difference whose
truncation error depends on which variable the derivative is taken in.

This script prints, for all 133 galaxies, how each convention compares with the
two published columns of paper1/Hermes_ConfigG_PerGalaxy_133_Export.csv. See
docs/gate_version_history.md ("A third lineage: the Paper 1 scoring
gate") for the discussion. Exits 0; it reports, it does not judge.
"""
import csv
import math
import os
import sys

import numpy as np
from scipy.signal import savgol_filter

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "paper1"))

import verify_hermes as vh                      # noqa: E402  data reader and g_bar
from hermes_gate_phi import hermes_phi          # noqa: E402  the log-grid gate

EXPORT = os.path.join(ROOT, "paper1", "Hermes_ConfigG_PerGalaxy_133_Export.csv")
SPARC_DIR = os.environ.get("SPARC_DIR", os.path.join(ROOT, "sparc_database"))
A_KNEE = 1585.0
KDIV = 46654.0
FLOOR = 1.0 / math.sqrt(2.0 * math.pi)


def phi_chain(R_kpc, g_bar, a_knee=A_KNEE):
    """The same gate as hermes_phi, with the chain-rule shear on the native grid.

    Every other step is identical to paper1/hermes_gate_phi.py. The single
    difference is marked below.
    """
    R = np.asarray(R_kpc, dtype=float).copy()
    g = np.asarray(g_bar, dtype=float).copy()
    order = np.argsort(R)
    R, g = R[order], g[order]
    uniq = np.concatenate(([True], np.diff(R) > 0))
    R, g = R[uniq], g[uniq]
    if len(R) < 3:
        return np.zeros_like(R)
    g = np.where(np.isfinite(g) & (g >= 0), g, 0.0)

    idx = None
    for i in range(len(R) - 1):
        if g[i] >= a_knee > g[i + 1]:
            idx = i
            break
    if idx is not None:
        g1, g2 = g[idx], g[idx + 1]
        t = 0.0 if g2 == g1 else (a_knee - g1) / (g2 - g1)
        r_knee = R[idx] + t * (R[idx + 1] - R[idx])
    else:
        r_knee = R[np.argmin(np.abs(g - a_knee))]
    if (not np.isfinite(r_knee)) or (r_knee <= 0):
        r_knee = np.median(R)

    x = R / r_knee
    V = np.sqrt(np.maximum(g * R, 0.0))
    N = len(R)
    if N < 5:
        V_sm = V
    else:
        window = 5 if N < 20 else 11
        if window > N:
            window = N if (N % 2 == 1) else (N - 1)
        V_sm = savgol_filter(V, window_length=window, polyorder=3, mode="mirror")
    V_sm = np.where(V_sm > 0, V_sm, 1e-12)

    dVdR = np.gradient(V_sm, R, edge_order=1)
    s = np.abs(R / V_sm * dVdR)              # <-- THE ONE DIFFERENCE (chain rule)
    d2VdR2 = np.gradient(dVdR, R, edge_order=1)
    eps_kappa = 1e-4 * np.median(np.abs(dVdR))
    if (not np.isfinite(eps_kappa)) or (eps_kappa <= 0):
        eps_kappa = 1e-12
    kappa = np.abs(d2VdR2) / (np.abs(dVdR) + eps_kappa)

    eps = 1e-12
    inner, outer = x <= 1.0, x > 1.0
    med_k = np.median(kappa[inner]) if np.any(inner) else np.median(kappa)
    med_s = np.median(s[outer]) if np.any(outer) else np.median(s)
    w = 1.0 / (1.0 + np.exp(-(x - 1.40) / 0.30))
    S_pre = (1.0 - w) * (kappa / (med_k + eps)) + w * (s / (med_s + eps))
    S_eff = np.empty_like(S_pre)
    S_eff[0] = S_pre[0]
    for i in range(1, N):
        S_eff[i] = 0.5 * S_pre[i] + 0.5 * S_eff[i - 1]
    return np.clip(1.0 - np.exp(-np.maximum(S_eff, 0.0)), 0.0, 1.0)


def main():
    if not os.path.isdir(SPARC_DIR):
        sys.stderr.write("SPARC rotmod files not found. Set SPARC_DIR; see the README.\n")
        return 2
    with open(EXPORT, encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle))

    chi2 = {"chain-rule": [], "log-grid": []}
    phi = {"chain-rule": [], "log-grid": []}
    medians = {"chain-rule": [], "log-grid": []}
    missing = []
    for row in rows:
        galaxy = row["galaxy"].strip()
        path = vh.find_rotmod(SPARC_DIR, galaxy)
        if path is None:
            missing.append(galaxy)
            continue
        data = vh.read_rotmod(path)
        g_bar = vh.compute_gbar(data)
        beta = math.pi * math.exp(-float(row["t50_gyr"]) * float(row["g98"]) / KDIV) - FLOOR
        for name, gate in (("chain-rule", phi_chain), ("log-grid", hermes_phi)):
            p = gate(data["R"], g_bar)
            V = np.sqrt(np.maximum(g_bar * (1.0 + beta * p), 0.0) * data["R"])
            c = vh.chi2_nu(data["Vobs"], V, data["errV"])
            medians[name].append(c)
            chi2[name].append(abs(c - float(row["chi2nu_configg"])))
            phi[name].append(abs(p[-1] - float(row["phi_last"])))

    n = len(medians["chain-rule"])
    if missing:
        print("not found in %s: %d galaxies" % (SPARC_DIR, len(missing)))
    print()
    print("Published export: paper1/Hermes_ConfigG_PerGalaxy_133_Export.csv   (%d galaxies)" % n)
    print()
    print("%-12s %28s %28s" % ("convention", "vs chi2nu_configg", "vs phi_last"))
    print("-" * 72)
    for name in ("chain-rule", "log-grid"):
        print("%-12s %14.2e  %5d off   %14.2e  %5d off"
              % (name, max(chi2[name]), sum(1 for v in chi2[name] if v > 1e-9),
                 max(phi[name]), sum(1 for v in phi[name] if v > 1e-9)))
    print()
    print("(max |difference|, and how many of the %d galaxies differ by more than 1e-9)" % n)
    print()
    for name in ("chain-rule", "log-grid"):
        print("median chi2/N, %-11s %.6f" % (name + ":", np.median(medians[name])))
    print()
    print("The chi-squared column was produced with the chain-rule convention; the phi")
    print("column with the log-grid convention. Use the one that matches the number you")
    print("are reproducing. See docs/gate_version_history.md.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
