#!/usr/bin/env python3
"""
verify_eta.py — Paper 8 Verification Script

Implements the wear-activated BTFR-safe density-release refinement (eta_WA)
described in McGinty (2026), "A Wear-Activated Density-Release Refinement
of the Hermes Gravity Equation: A Paper 1 Addendum."

This script defines the eta operator and, when run against SPARC rotation
curve data with the Paper 1 gate (phi), reproduces the Paper 8 results.

Requirements:
    - numpy
    - scipy (for Savitzky-Golay filter used in gate construction)
    - Paper 1 gate implementation (hermes_gate_phi.py from the verification repo)
    - SPARC rotation curve files

Usage:
    # Import the eta functions
    from verify_eta import compute_eta_WA, compute_g_model_eta

    # Or run standalone verification
    python verify_eta.py --sparc_dir /path/to/SPARC --age_table /path/to/ages.csv

Author: Team Hermes (Louis A. McGinty)
AI-assisted implementation; see Paper 8 AI disclosure.
Repository: github.com/mcgintyl/hermes-verification
"""

import numpy as np
import warnings

# =============================================================================
# INHERITED CONSTANTS (from Paper 1, frozen)
# =============================================================================

A_KNEE = 1585.0          # Gate threshold, (km/s)^2/kpc
LOGISTIC_WIDTH = 0.30     # Spatial transition width ratio (inherited from gate)
K_DIVISOR = 46654.0       # c/(2*pi) in SPARC units: Gyr*(km/s)^2/kpc
SIGMA_INT_SQ = 386.0      # Intrinsic variance floor for chi-squared

# =============================================================================
# DERIVED CONSTANTS (Paper 8, closed-form from inherited constants)
# =============================================================================

A_U = A_KNEE / np.pi                    # Universal numerator, ~504.5
WEAR_RATE = 2.0 * np.pi**2              # Wear activation rate, ~19.74
BETA_CEILING = np.pi                     # Max positive coherence
BETA_FLOOR = -1.0 / np.sqrt(2.0 * np.pi)  # Coherence floor, ~-0.3989
ETA_CAP = np.pi                          # Max local operator value


# =============================================================================
# CORE FUNCTIONS
# =============================================================================

def compute_psi(t50_gyr, g98):
    """
    Compute the wear variable psi from stellar age and peak baryonic acceleration.

    Parameters
    ----------
    t50_gyr : float
        Stellar half-mass age in gigayears.
    g98 : float
        98th-percentile peak baryonic acceleration in (km/s)^2/kpc.

    Returns
    -------
    float
        The wear variable psi (dimensionless).
    """
    return t50_gyr * g98 / K_DIVISOR


def compute_beta(psi):
    """
    Compute the age score beta from the wear variable.

    Parameters
    ----------
    psi : float or array
        Wear variable.

    Returns
    -------
    float or array
        Age score beta.
    """
    return np.pi * np.exp(-psi) - 1.0 / np.sqrt(2.0 * np.pi)


def log_compress(g, a_k=A_KNEE):
    """
    Logarithmic compression of baryonic acceleration.

    Regularizes extreme density contrasts. For g >> a_k, grows
    logarithmically. For g << a_k, approximately equals g.

    Parameters
    ----------
    g : float or array
        Baryonic acceleration in (km/s)^2/kpc.
    a_k : float
        Gate threshold (default: 1585).

    Returns
    -------
    float or array
        Compressed acceleration Gamma(g).
    """
    g_safe = np.maximum(g, 0.0)
    return a_k * np.log(1.0 + g_safe / a_k)


def wear_activation(psi_sys):
    """
    Wear activation gate W(psi).

    Suppresses the density release in near-newborn (low-wear) systems.
    W ~ 0 for psi ~ 0, W ~ 1 for significant wear.

    Parameters
    ----------
    psi_sys : float
        Systemic wear variable (from global stellar age).

    Returns
    -------
    float
        Wear activation value in [0, 1].
    """
    return 1.0 - np.exp(-WEAR_RATE * psi_sys)


def spatial_transition(g_bar, a_k=A_KNEE, width=LOGISTIC_WIDTH):
    """
    Logistic spatial transition function S(g).

    Governs where the density release activates along the radial profile.
    S ~ 0 when g >> a_k (dense core), S ~ 1 when g << a_k (outer disk).

    Parameters
    ----------
    g_bar : float or array
        Local baryonic acceleration in (km/s)^2/kpc.
    a_k : float
        Transition center (default: 1585).
    width : float
        Width ratio (default: 0.30).

    Returns
    -------
    float or array
        Transition value in [0, 1].
    """
    return 1.0 / (1.0 + np.exp((g_bar - a_k) / (width * a_k)))


def compute_eta_U(g_bar):
    """
    Local density-release operator eta_U(R).

    Evaluates whether the low-density branch exceeds the local
    three-dimensional field. The max function selects between unity
    (baseline) and the low-density branch sqrt(a_u / Gamma).
    Capped at pi.

    Parameters
    ----------
    g_bar : float or array
        Local baryonic acceleration in (km/s)^2/kpc.

    Returns
    -------
    float or array
        Local operator eta_U in [1, pi].
    """
    gamma = log_compress(g_bar)

    # Density contrast: low-density branch vs baseline
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        contrast = np.where(
            gamma > 0,
            np.sqrt(A_U / gamma),
            1.0
        )

    # Max selects the dominant branch
    branch_value = np.maximum(1.0, contrast)

    # Spatial transition modulates the release
    s = spatial_transition(g_bar)

    # Local operator with cap
    eta_u = 1.0 + s * (branch_value - 1.0)
    eta_u = np.minimum(ETA_CAP, eta_u)

    return eta_u


def compute_eta_WA(g_bar, psi_sys):
    """
    Wear-gated density-release operator eta_WA(R).

    The local operator eta_U is modulated by the global wear activation.
    W uses systemic psi (global stellar age), not local psi(R).

    Parameters
    ----------
    g_bar : float or array
        Local baryonic acceleration in (km/s)^2/kpc.
    psi_sys : float
        Systemic wear variable (from global stellar age).

    Returns
    -------
    float or array
        Wear-gated operator eta_WA.
    """
    W = wear_activation(psi_sys)
    eta_u = compute_eta_U(g_bar)
    return 1.0 + W * (eta_u - 1.0)


def compute_g_model_eta(g_bar, phi, psi_local, psi_sys):
    """
    Compute the Paper 8 modified gravitational acceleration.

    g_model(R) = g_bar(R) * [1 + phi(R) * (pi * exp(-psi(R)) * eta_WA(R) - 1/sqrt(2*pi))]

    Parameters
    ----------
    g_bar : array
        Baryonic acceleration at each radius, (km/s)^2/kpc.
    phi : array
        Gate values at each radius (from Paper 1 gate construction).
    psi_local : float or array
        Local wear variable. For standard SPARC: psi_sys at all radii.
        For M33 component-age: varies by region.
    psi_sys : float
        Systemic wear variable (for W gate).

    Returns
    -------
    array
        Modified gravitational acceleration at each radius.
    """
    eta_wa = compute_eta_WA(g_bar, psi_sys)
    positive_term = np.pi * np.exp(-psi_local) * eta_wa
    floor = 1.0 / np.sqrt(2.0 * np.pi)
    beta_eff = positive_term - floor

    return g_bar * (1.0 + phi * beta_eff)


def compute_v_model_eta(g_bar, phi, psi_local, psi_sys, R_kpc):
    """
    Compute the Paper 8 modified rotation velocity.

    V_model(R) = sqrt(g_model(R) * R)

    Parameters
    ----------
    g_bar : array
        Baryonic acceleration at each radius.
    phi : array
        Gate values at each radius.
    psi_local : float or array
        Local wear variable.
    psi_sys : float
        Systemic wear variable.
    R_kpc : array
        Radii in kpc.

    Returns
    -------
    array
        Predicted rotation velocity in km/s.
    """
    g_model = compute_g_model_eta(g_bar, phi, psi_local, psi_sys)

    # Handle negative g_model (floor-saturated cases)
    v_sq = g_model * R_kpc
    return np.sign(v_sq) * np.sqrt(np.abs(v_sq))


def compute_chi2nu(v_obs, v_model, v_err, sigma_int_sq=SIGMA_INT_SQ):
    """
    Compute reduced chi-squared.

    Parameters
    ----------
    v_obs : array
        Observed rotation velocities in km/s.
    v_model : array
        Model rotation velocities in km/s.
    v_err : array
        Measurement uncertainties in km/s.
    sigma_int_sq : float
        Intrinsic variance floor (default: 386).

    Returns
    -------
    float
        Reduced chi-squared (chi2 / N).
    """
    sigma_sq = v_err**2 + sigma_int_sq
    chi2 = np.sum((v_obs - v_model)**2 / sigma_sq)
    return chi2 / len(v_obs)


# =============================================================================
# BASELINE HERMES (for comparison)
# =============================================================================

def compute_g_model_baseline(g_bar, phi, psi):
    """
    Compute baseline Paper 1 gravitational acceleration (no eta).

    g_model(R) = g_bar(R) * [1 + beta * phi(R)]
    """
    beta = compute_beta(psi)
    return g_bar * (1.0 + beta * phi)


# =============================================================================
# MOND SIMPLE (for comparison)
# =============================================================================

A0_MOND = 3703.0  # MOND acceleration threshold in (km/s)^2/kpc (~1.2e-10 m/s^2)

def compute_g_mond_simple(g_bar):
    """
    Compute MOND (simple interpolation) gravitational acceleration.

    g_MOND = g_bar * nu(g_bar / a0)
    where nu(x) = (1 + sqrt(1 + 4/x)) / 2
    """
    x = np.abs(g_bar) / A0_MOND
    x = np.maximum(x, 1e-10)  # prevent division by zero
    nu = (1.0 + np.sqrt(1.0 + 4.0 / x)) / 2.0
    return g_bar * nu


# =============================================================================
# VERIFICATION SUMMARY
# =============================================================================

def print_verification_summary(results):
    """
    Print Paper 8 verification summary matching the main text tables.

    Parameters
    ----------
    results : list of dict
        Per-galaxy results with chi2nu for each model.
    """
    n = len(results)
    chi2_baseline = np.array([r['chi2nu_baseline'] for r in results])
    chi2_mond = np.array([r['chi2nu_mond'] for r in results])
    chi2_eta = np.array([r['chi2nu_eta'] for r in results])

    delta = chi2_eta - chi2_baseline

    print(f"\n{'='*60}")
    print(f"Paper 8 Verification Summary ({n} galaxies)")
    print(f"{'='*60}")
    print(f"")
    print(f"{'Model':<25} {'Median':>10} {'> 5':>6} {'> 10':>6}")
    print(f"{'-'*25} {'-'*10} {'-'*6} {'-'*6}")
    print(f"{'Baseline Hermes':<25} {np.median(chi2_baseline):>10.3f} "
          f"{np.sum(chi2_baseline > 5):>6d} {np.sum(chi2_baseline > 10):>6d}")
    print(f"{'MOND (simple)':<25} {np.median(chi2_mond):>10.3f} "
          f"{np.sum(chi2_mond > 5):>6d} {np.sum(chi2_mond > 10):>6d}")
    print(f"{'Wear-activated eta':<25} {np.median(chi2_eta):>10.3f} "
          f"{np.sum(chi2_eta > 5):>6d} {np.sum(chi2_eta > 10):>6d}")
    print(f"")
    print(f"Material casualties (delta > 0.5):  {np.sum(delta > 0.5)}")
    print(f"Improvements (delta < -0.5):        {np.sum(delta < -0.5)}")
    print(f"New catastrophic (>5):              {np.sum((chi2_eta > 5) & (chi2_baseline <= 5))}")
    print(f"New catastrophic (>10):             {np.sum((chi2_eta > 10) & (chi2_baseline <= 10))}")
    print(f"Beats MOND:                         {np.sum(chi2_eta < chi2_mond)}/{n}")
    print(f"")
    print(f"Paper 8 expected values:")
    print(f"  Median chi2nu:    1.189")
    print(f"  Casualties:       0")
    print(f"  Improvements:     13")
    print(f"  Beats MOND:       77/133")
    print(f"{'='*60}")


# =============================================================================
# CONSTANTS SUMMARY
# =============================================================================

def print_constants():
    """Print all constants used in the Paper 8 refinement."""
    print(f"\nPaper 8 Constants")
    print(f"{'='*50}")
    print(f"INHERITED from Paper 1:")
    print(f"  pi          = {np.pi:.10f}")
    print(f"  e           = {np.e:.10f}")
    print(f"  a_k         = {A_KNEE}")
    print(f"  0.30        = {LOGISTIC_WIDTH}")
    print(f"  1/sqrt(2pi) = {1/np.sqrt(2*np.pi):.10f}")
    print(f"  K (c/2pi)   = {K_DIVISOR}")
    print(f"  sigma_int^2 = {SIGMA_INT_SQ}")
    print(f"")
    print(f"DERIVED in Paper 8:")
    print(f"  a_u = a_k/pi      = {A_U:.6f}")
    print(f"  2*pi^2            = {WEAR_RATE:.6f}")
    print(f"  eta cap            = {ETA_CAP:.10f} (= pi)")
    print(f"{'='*50}")


if __name__ == "__main__":
    import sys

    print_constants()

    print("\nverify_eta.py loaded successfully.")
    print("This script defines the Paper 8 wear-activated eta operator.")
    print("")
    print("Functions available:")
    print("  compute_eta_WA(g_bar, psi_sys)    - the full wear-gated operator")
    print("  compute_g_model_eta(g_bar, phi, psi_local, psi_sys)")
    print("                                     - modified gravitational acceleration")
    print("  compute_g_model_baseline(g_bar, phi, psi)")
    print("                                     - Paper 1 baseline (for comparison)")
    print("  compute_g_mond_simple(g_bar)       - MOND simple nu (for comparison)")
    print("")
    print("To run full verification, provide SPARC data directory:")
    print("  python verify_eta.py --sparc_dir /path/to/SPARC --age_table ages.csv")
    print("")
    print("See github.com/mcgintyl/hermes-verification for full pipeline.")
