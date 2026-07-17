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

Distance sensitivity (--distance / --distance-sensitivity)
----------------------------------------------------------
The board's radii are tabulated at Corbelli et al. (2014)'s adopted distance
D0 = 840 kpc. To test how the fit depends on that distance, this script rescales
the board self-consistently. M33's surface densities (Msun/pc^2) derive from
distance-independent surface brightness / HI column density, so under D -> D':

    R -> R * (D'/D0)                          (physical radii scale with distance)
    Vgas, Vdisk, Vbulge -> * sqrt(D'/D0)      (thin-disk: V^2 ~ Sigma * R, Sigma fixed)
    Vobs, errV -> unchanged                   (kinematic; distance-independent)

This leaves g_bar, g98, beta, and phi very nearly invariant; the fit changes
because V_hermes = Vbar*sqrt(1+beta*phi) scales as sqrt(D') while Vobs does not.
See m33_radius_uncertainty_note.md. (Scaling assumption: standard thin-disk
V^2 ~ Sigma*R; the fixed physical disk thicknesses are a second-order caveat.)

Self-check (default distance only): phi must match the committed canonical
per-point CSV to 1e-12.
"""
import os, sys, csv, argparse
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

# Distance provenance (see m33_radius_uncertainty_note.md)
D0_KPC = 840.0                         # Corbelli et al. 2014 adopted distance
DELTA_MU = 0.07                        # Gieren et al. 2013, mu = 24.62 +/- 0.07 mag
FRAC_DD = (np.log(10) / 5.0) * DELTA_MU  # dD/D = 0.0322


def beta_of(t50, g98):
    return np.pi * np.exp(-t50 * g98 / C_OVER_2PI) - 1.0 / np.sqrt(2.0 * np.pi)


def load_board():
    rows = list(csv.DictReader(open(BOARD)))
    col = lambda n: np.array([float(r[n]) for r in rows])
    return dict(R=col("R_kpc"), Vobs=col("Vobs_kms"), errV=col("errV_kms"),
                Vgas=col("Vgas_kms"), Vdisk=col("Vdisk_kms"), Vbul=col("Vbulge_kms"))


def hermes_at_distance(b, D_kpc):
    """Return (chi2nu_full, chi2nu_outer, phi, beta) per age, with the board
    rescaled self-consistently from D0 to D_kpc."""
    s = D_kpc / D0_KPC
    R = b["R"] * s
    Vgas, Vdisk, Vbul = (b["Vgas"] * np.sqrt(s), b["Vdisk"] * np.sqrt(s), b["Vbul"] * np.sqrt(s))
    Vobs, errV = b["Vobs"], b["errV"]              # kinematic, unchanged
    Vbar2 = Vdisk**2 + Vbul**2 + np.sign(Vgas) * Vgas**2
    Vbar = np.sqrt(np.maximum(0.0, Vbar2))
    gbar = Vbar2 / R
    g98 = np.percentile(gbar, 98)
    phi = hermes_phi(R, gbar)
    denom = errV**2 + SIGMA_INT_SQ
    N = len(R)
    outer = R > 10.0
    out = {}
    for label, t50 in AGES:
        beta = beta_of(t50, g98)
        Vh = Vbar * np.sqrt(np.maximum(0.0, 1.0 + beta * phi))
        ci = (Vh - Vobs)**2 / denom
        out[label] = dict(chi2nu=ci.sum() / N, chi2nu_outer=ci[outer].sum() / outer.sum(),
                          beta=beta, g98=g98)
    return out, phi


def main():
    ap = argparse.ArgumentParser(description="M33 canonical-gate Hermes calculation.")
    ap.add_argument("--distance", type=float, default=D0_KPC,
                    help="assumed M33 distance in kpc (default %.0f, Corbelli 2014)" % D0_KPC)
    ap.add_argument("--distance-sensitivity", action="store_true",
                    help="print chi2nu at D0 and D0 +/- 1sigma (from Gieren 2013, %.1f%%)" % (100 * FRAC_DD))
    args = ap.parse_args()

    b = load_board()

    if args.distance_sensitivity:
        dD = FRAC_DD * D0_KPC
        print("M33 distance sensitivity  (D0 = %.0f kpc, 1sigma = %.1f%% = %.0f kpc; Gieren 2013)"
              % (D0_KPC, 100 * FRAC_DD, dD))
        print("%-10s %-9s | %-22s | %-22s" % ("D (kpc)", "", "primary t50=6.0", "conservative t50=7.1"))
        print("%-10s %-9s | %-10s %-10s | %-10s %-10s"
              % ("", "", "chi2nu", "chi2nu>10", "chi2nu", "chi2nu>10"))
        for D in (D0_KPC - dD, D0_KPC, D0_KPC + dD):
            tag = "(-1sigma)" if D < D0_KPC else ("(+1sigma)" if D > D0_KPC else "(adopted)")
            r, _ = hermes_at_distance(b, D)
            print("%-10.1f %-9s | %-10.3f %-10.3f | %-10.3f %-10.3f"
                  % (D, tag, r["primary"]["chi2nu"], r["primary"]["chi2nu_outer"],
                     r["conservative"]["chi2nu"], r["conservative"]["chi2nu_outer"]))
        print("\nNote: chi2nu moves because V_model ~ sqrt(D) while V_obs is fixed;")
        print("g_bar, g98, beta, phi stay ~invariant. errR_kpc is a coherent scale,")
        print("not independent per-point error -- see m33_radius_uncertainty_note.md.")
        return 0

    r, phi = hermes_at_distance(b, args.distance)
    g98 = r["primary"]["g98"]
    print("M33 canonical-gate Hermes calculation  (N=%d, D=%.1f kpc, g98=%.1f)"
          % (len(b["R"]), args.distance, g98))
    print("%-14s %-9s %-10s %-14s %-14s" % ("age", "t50", "beta", "chi2nu full", "chi2nu R>10"))
    for label, t50 in AGES:
        print("%-14s %-9.1f %-10.5f %-14.3f %-14.3f"
              % (label, t50, r[label]["beta"], r[label]["chi2nu"], r[label]["chi2nu_outer"]))

    # self-check only meaningful at the adopted distance
    if abs(args.distance - D0_KPC) < 1e-9:
        ref = list(csv.DictReader(open(REF)))
        phi_ref = np.array([float(x["phi_canonical"]) for x in ref])
        dphi = float(np.max(np.abs(phi - phi_ref)))
        ok = dphi <= 1e-12
        print("\nself-check vs committed CSV: max|delta phi| = %.2e  -> %s"
              % (dphi, "PASS" if ok else "FAIL"))
        return 0 if ok else 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
