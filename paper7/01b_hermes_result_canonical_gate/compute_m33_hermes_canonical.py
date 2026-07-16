#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
compute_m33_hermes_canonical.py -- the Hermes calculation for M33 (canonical gate).

This is the file that performs the M33 Hermes calculation for Paper 7 v3
(DOI 10.5281/zenodo.21347359): board -> phi -> V_hermes -> chi2nu.
The gate is IMPORTED from paper1/hermes_gate_phi.py, not copied, so there is
exactly one gate implementation in play.

Pipeline (identical to Paper 1's frozen chassis):
  1. Read the frozen Format B board (R, Vobs, errV, Vgas, Vdisk, Vbulge).
  2. Vbar^2 = Vdisk^2 + Vbulge^2 + sign(Vgas)*Vgas^2 ;   g_bar = Vbar^2 / R
  3. phi(R) = canonical gate (paper1/hermes_gate_phi.py)
  4. beta   = pi*exp(-psi) - 1/sqrt(2*pi),  psi = t50*g98 / 46654
              (c/(2*pi) = 46654 Gyr (km/s)^2/kpc; g98 = 98th pct of g_bar)
              ages: t50 = 6.0 Gyr (primary), 7.1 Gyr (conservative)
  5. V_hermes = Vbar * sqrt(1 + beta*phi)
  6. chi2nu  = sum( (V_hermes-Vobs)^2 / (errV^2 + sigma_int^2) ) / N,
              sigma_int^2 = 386 (same scoring convention as Paper 1)

Self-check: phi must match the committed canonical per-point CSV to 1e-12.
Note on V_hermes: the committed CSV's V_hermes columns were generated with beta
rounded to 5 decimals (1.83421 / 1.69876); this script computes beta analytically,
so V_hermes agrees to ~1e-3 km/s rather than machine precision. chi2nu is
unaffected at the reported precision (3.22 / 3.39).
"""
import os, sys, csv
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.normpath(os.path.join(HERE, "..", "..", "paper1")))
from hermes_gate_phi import hermes_phi                     # the single canonical gate

BOARD = os.path.normpath(os.path.join(
    HERE, "..", "01_hermes_result_locked", "m33_corbelli_reconstructed_board_input.csv"))
REF = os.path.join(HERE, "m33_hermes_per_point_predictions_canonical.csv")

C_OVER_2PI = 46654.0       # Gyr (km/s)^2 / kpc
SIGMA_INT_SQ = 386.0       # (km/s)^2
AGES = [("primary", 6.0), ("conservative", 7.1)]


def beta_of(t50, g98):
    return np.pi * np.exp(-t50 * g98 / C_OVER_2PI) - 1.0 / np.sqrt(2.0 * np.pi)


def main():
    rows = list(csv.DictReader(open(BOARD)))
    col = lambda n: np.array([float(r[n]) for r in rows])
    R, Vobs, errV = col("R_kpc"), col("Vobs_kms"), col("errV_kms")
    Vgas, Vdisk, Vbul = col("Vgas_kms"), col("Vdisk_kms"), col("Vbulge_kms")

    Vbar2 = Vdisk**2 + Vbul**2 + np.sign(Vgas) * Vgas**2
    Vbar = np.sqrt(np.maximum(0.0, Vbar2))
    gbar = Vbar2 / R
    g98 = np.percentile(gbar, 98)
    phi = hermes_phi(R, gbar)
    denom = errV**2 + SIGMA_INT_SQ
    N = len(R)

    print("M33 canonical-gate Hermes calculation  (N=%d, g98=%.1f)" % (N, g98))
    print("%-14s %-9s %-10s %-14s %-14s" % ("age", "t50", "beta", "chi2nu full", "chi2nu R>10"))
    for label, t50 in AGES:
        beta = beta_of(t50, g98)
        Vh = Vbar * np.sqrt(np.maximum(0.0, 1.0 + beta * phi))
        ci = (Vh - Vobs)**2 / denom
        outer = R > 10.0
        print("%-14s %-9.1f %-10.5f %-14.3f %-14.3f"
              % (label, t50, beta, ci.sum() / N, ci[outer].sum() / outer.sum()))

    # self-check against the committed canonical per-point CSV
    ref = list(csv.DictReader(open(REF)))
    phi_ref = np.array([float(r["phi_canonical"]) for r in ref])
    dphi = float(np.max(np.abs(phi - phi_ref)))
    ok = dphi <= 1e-12
    print("\nself-check vs committed CSV: max|delta phi| = %.2e  -> %s"
          % (dphi, "PASS" if ok else "FAIL"))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
