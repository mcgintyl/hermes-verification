#!/usr/bin/env python3
"""
==============================================================================
 SUPERSEDED -- HISTORICAL REPRODUCTION ONLY.  DO NOT USE FOR NEW WORK.
 The canonical gate is  paper1/hermes_gate_phi.py :: hermes_phi()
 See historical/README.md and docs/gate_version_history.md.
==============================================================================

HISTORICAL Tier-1 gate (v1) -- reconstruction of the January 2026 gate that
produced the Stage-1 lock beta ~= 2.807 (rounded 2.81).

This is NOT the current Config G gate. Do not use it for the frozen 133-galaxy
results. It exists only to reproduce the *historical* Stage-1 calibration so the
derivation path is auditable.

DO NOT MIX THE TWO GATES. This gate expects the historical *unsigned* g_bar
(use gbar_unsigned() below). Passing the current signed g_bar into this gate --
or the historical unsigned g_bar into hermes_phi() -- yields a hybrid matching
no configuration, current or historical, and fails silently. Across the 133
galaxies that hybrid moves velocities by a median 5%, up to 20% (NGC 5371).
Use gbar_unsigned() + tier1_gate() together, or neither.

Difference from the current gate (paper1/hermes_gate_phi.py, verify_hermes.py),
see docs/gate_version_history.md:

  quantity              historical (this file)          current Config G gate
  --------------------  ------------------------------  ------------------------
  outer shear           |(R/V) dV/dR|  (chain rule)     |d ln V / d ln R| (log grid)
  curvature stabilizer  1e-12                           1e-4 * median(|dV/dR|)
  zone normalization    arithmetic MEAN                 MEDIAN
  baryonic gas term     unsigned quadrature             sign(Vgas) * Vgas^2

Everything else (knee at a_knee = 1585, Savitzky-Golay window 5/11 poly-3 mirror,
logistic mixer x0 = 1.40 / s = 0.30, EMA lambda = 0.5, gate 1 - exp(-S_eff)) is
identical to the current gate.

Validation targets (reproduced by fit_tier1_beta.py against SPARC / Zenodo
16284118):
  * F583-4 phi matches the Tier-1 spec golden vector to < 1e-3.
  * global Stage-1 beta = 2.806585  (historical lock 2.807 -> 2.81).
"""
import numpy as np
from scipy.signal import savgol_filter

A_KNEE = 1585.0        # (km/s)^2 / kpc
SIGMA_INT2 = 386.0     # (km/s)^2, intrinsic variance added to errV^2

# The seven young, bulgeless galaxies of the Stage-1 calibration set.
TIER1_SEVEN = ["D631-7", "DDO 154", "DDO 168", "ESO 444-G084",
               "F583-4", "UGC 06930", "UGC 07608"]

# SPARC rotmod filename stems for the seven (no radial cuts).
TIER1_STEMS = {
    "D631-7": "D631-7", "DDO 154": "DDO154", "DDO 168": "DDO168",
    "ESO 444-G084": "ESO444-G084", "F583-4": "F583-4",
    "UGC 06930": "UGC06930", "UGC 07608": "UGC07608",
}

# F583-4 golden phi vector from the Tier-1 Pristine Law spec section 6
# (means gate, beta = 2.81). Used as the historical implementation check.
F583_4_GOLDEN_PHI = [0.8365, 0.7314, 0.6577, 0.7410, 0.7476, 0.6931,
                     0.7438, 0.7678, 0.6559, 0.5874, 0.4186, 0.3723]


def gbar_unsigned(data):
    """Baryonic acceleration with UNSIGNED gas quadrature (Tier-1 spec 2.1)."""
    R = data["R"]
    vbar2 = data["Vgas"] ** 2 + data["Vdisk"] ** 2 + data["Vbul"] ** 2
    return np.where(R > 0, vbar2 / R, 0.0)


def tier1_gate(R_kpc, g_bar, a_knee=A_KNEE):
    """Return the HISTORICAL Tier-1 gate phi(R) in [0, 1] at each radius."""
    R = np.asarray(R_kpc, dtype=float).copy()
    g = np.asarray(g_bar, dtype=float).copy()

    order = np.argsort(R)
    R, g = R[order], g[order]
    uniq = np.concatenate(([True], np.diff(R) > 0))
    R, g = R[uniq], g[uniq]
    if len(R) < 3:
        return np.zeros_like(R)
    if np.any(R <= 0):
        raise ValueError("R_kpc must be > 0 everywhere.")
    g = np.where(np.isfinite(g) & (g >= 0), g, 0.0)

    # knee radius (identical to current gate)
    idx = None
    for i in range(len(R) - 1):
        if (g[i] >= a_knee) and (g[i + 1] < a_knee):
            idx = i
            break
    if idx is not None:
        r1, r2, g1, g2 = R[idx], R[idx + 1], g[idx], g[idx + 1]
        t = 0.0 if g2 == g1 else (a_knee - g1) / (g2 - g1)
        r_knee = r1 + t * (r2 - r1)
    else:
        r_knee = R[np.argmin(np.abs(g - a_knee))]
    if (not np.isfinite(r_knee)) or (r_knee <= 0):
        r_knee = np.median(R)
    x = R / r_knee

    # smooth V(R) (identical to current gate)
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

    # HISTORICAL shear: chain-rule form |(R/V) dV/dR|  (not the log-grid gradient)
    dVdR = np.gradient(V_sm, R, edge_order=1)
    s = np.abs((R / V_sm) * dVdR)

    # HISTORICAL curvature: fixed 1e-12 stabilizer (not 1e-4 * median)
    d2VdR2 = np.gradient(dVdR, R, edge_order=1)
    kappa = np.abs(d2VdR2) / (np.abs(dVdR) + 1e-12)

    # HISTORICAL zone normalization: arithmetic MEAN (not median)
    eps = 1e-12
    inner = x <= 1.0
    outer = x > 1.0
    mean_kappa_inner = np.mean(kappa[inner]) if np.any(inner) else np.mean(kappa)
    mean_s_outer = np.mean(s[outer]) if np.any(outer) else np.mean(s)
    kappa_norm = kappa / (mean_kappa_inner + eps)
    s_norm = s / (mean_s_outer + eps)

    # mixer, stress, EMA, gate (identical to current gate)
    w = 1.0 / (1.0 + np.exp(-(x - 1.40) / 0.30))
    S_pre = (1.0 - w) * kappa_norm + w * s_norm
    S_eff = np.empty_like(S_pre)
    S_eff[0] = S_pre[0]
    for i in range(1, N):
        S_eff[i] = 0.5 * S_pre[i] + 0.5 * S_eff[i - 1]
    S_eff = np.maximum(S_eff, 0.0)
    return np.clip(1.0 - np.exp(-S_eff), 0.0, 1.0)
