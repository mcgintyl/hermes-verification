#!/usr/bin/env python3
"""
Hermes Equation Verification Tool
==================================
Reproduces the per-galaxy chi-squared scores from:

    McGinty (2026), "The Hermes Equation: Galaxy Rotation Curves
    from Age with No Free Constants"

Usage:
    python verify_hermes.py --sparc <rotmod_dir> --ages <age_csv>

Inputs:
    rotmod_dir : folder of SPARC *_rotmod.dat files (one per galaxy)
    age_csv    : CSV with columns: galaxy, t50_gyr, g98
                 (the Hermes_ConfigG_PerGalaxy_133_Export.csv works)

Outputs:
    Per-galaxy table to stdout with chi-squared for Hermes and MOND.
    Optional --csv <path> writes full results to a CSV file.
"""

import argparse
import csv
import os
import sys
import numpy as np
from scipy.signal import savgol_filter


# ─────────────────────────────────────────────────────────
# Physical / equation constants
# ─────────────────────────────────────────────────────────
A_KNEE_DEFAULT = 1585.0        # (km/s)^2 / kpc
SIGMA_INT_SQ   = 386.0         # intrinsic variance floor  (km/s)^2
A0_MOND        = 1.2e-10       # m/s^2  — MOND acceleration scale
# Unit conversion: 1 (km/s)^2/kpc  =  3.2408e-14 m/s^2
KMS2_KPC_TO_MS2 = 3.2408e-14


# ─────────────────────────────────────────────────────────
# 1.  Gate function φ(R)  — verified code from hermes_gate_phi.py
# ─────────────────────────────────────────────────────────
def hermes_phi(R_kpc, g_bar, a_knee=A_KNEE_DEFAULT):
    """Return gate φ(R) in [0,1] at each radius."""
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

    # Step 1: knee radius
    idx_cross = None
    for i in range(len(R) - 1):
        if (g[i] >= a_knee) and (g[i + 1] < a_knee):
            idx_cross = i
            break

    if idx_cross is not None:
        r1, r2 = R[idx_cross], R[idx_cross + 1]
        g1, g2 = g[idx_cross], g[idx_cross + 1]
        t = 0.0 if g2 == g1 else (a_knee - g1) / (g2 - g1)
        r_knee = r1 + t * (r2 - r1)
    else:
        r_knee = R[np.argmin(np.abs(g - a_knee))]

    if (not np.isfinite(r_knee)) or (r_knee <= 0):
        r_knee = np.median(R)

    # Step 2: normalize
    x = R / r_knee

    # Step 3: smooth V(R)
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

    # Step 4: curvature and shear
    s = np.abs(np.gradient(np.log(V_sm), np.log(R), edge_order=1))

    dVdR = np.gradient(V_sm, R, edge_order=1)
    d2VdR2 = np.gradient(dVdR, R, edge_order=1)
    eps_kappa = 1e-4 * np.median(np.abs(dVdR))
    if (not np.isfinite(eps_kappa)) or (eps_kappa <= 0):
        eps_kappa = 1e-12
    kappa = np.abs(d2VdR2) / (np.abs(dVdR) + eps_kappa)

    eps = 1e-12
    inner = x <= 1.0
    outer = x > 1.0
    med_kappa_inner = np.median(kappa[inner]) if np.any(inner) else np.median(kappa)
    med_s_outer = np.median(s[outer]) if np.any(outer) else np.median(s)
    kappa_norm = kappa / (med_kappa_inner + eps)
    s_norm = s / (med_s_outer + eps)

    # Step 5: logistic mixer
    w = 1.0 / (1.0 + np.exp(-(x - 1.40) / 0.30))
    S_pre = (1.0 - w) * kappa_norm + w * s_norm

    # Step 6: EMA smooth and gate
    lam = 0.5
    S_eff = np.empty_like(S_pre)
    S_eff[0] = S_pre[0]
    for i in range(1, N):
        S_eff[i] = lam * S_pre[i] + (1.0 - lam) * S_eff[i - 1]
    S_eff = np.maximum(S_eff, 0.0)

    phi = 1.0 - np.exp(-S_eff)
    return np.clip(phi, 0.0, 1.0)


# ─────────────────────────────────────────────────────────
# 2.  Age score β
# ─────────────────────────────────────────────────────────
def hermes_beta(t50_gyr, g98_kms2kpc):
    """
    β = π · exp(−2π · t₅₀ · g₉₈ / c) − 1/√(2π)

    t50_gyr      : half-mass age in Gyr
    g98_kms2kpc  : 98th-percentile baryonic acceleration in (km/s)^2/kpc

    The denominator  c / (2π)  in compound units (Gyr·(km/s)^2/kpc) = 46654.
    So  ψ = 2π · t50 · g98 / c  =  t50 · g98 / 46654.
    """
    C_OVER_2PI = 46654.0        # c/(2π) in Gyr·(km/s)^2/kpc
    psi = t50_gyr * g98_kms2kpc / C_OVER_2PI
    beta = np.pi * np.exp(-psi) - 1.0 / np.sqrt(2.0 * np.pi)
    return beta


# ─────────────────────────────────────────────────────────
# 3.  Read a SPARC rotmod file
# ─────────────────────────────────────────────────────────
def read_rotmod(path):
    """
    Return dict with arrays: R, Vobs, errV, Vgas, Vdisk, Vbul
    All in SPARC native units (kpc, km/s).
    """
    R, Vobs, errV, Vgas, Vdisk, Vbul = [], [], [], [], [], []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split()
            if len(parts) < 6:
                continue
            R.append(float(parts[0]))
            Vobs.append(float(parts[1]))
            errV.append(float(parts[2]))
            Vgas.append(float(parts[3]))
            Vdisk.append(float(parts[4]))
            Vbul.append(float(parts[5]))
    return {
        'R':     np.array(R),
        'Vobs':  np.array(Vobs),
        'errV':  np.array(errV),
        'Vgas':  np.array(Vgas),
        'Vdisk': np.array(Vdisk),
        'Vbul':  np.array(Vbul),
    }


# ─────────────────────────────────────────────────────────
# 4.  Baryonic acceleration  g_bar  from SPARC columns
# ─────────────────────────────────────────────────────────
def compute_gbar(data):
    """
    V_bar²(R) = V_disk² + V_bul² + sgn(V_gas)·V_gas²
    g_bar(R)  = V_bar² / R      in (km/s)^2 / kpc
    """
    R = data['R']
    Vdisk = data['Vdisk']
    Vbul  = data['Vbul']
    Vgas  = data['Vgas']

    Vbar_sq = Vdisk**2 + Vbul**2 + np.sign(Vgas) * Vgas**2
    # Guard: where R=0 set gbar=0
    with np.errstate(divide='ignore', invalid='ignore'):
        gbar = np.where(R > 0, Vbar_sq / R, 0.0)
    return gbar


# ─────────────────────────────────────────────────────────
# 5.  Chi-squared (reduced) for a model
# ─────────────────────────────────────────────────────────
def chi2_nu(Vobs, Vmodel, errV, sigma_int_sq=SIGMA_INT_SQ):
    """
    χ²_ν  =  (1/N) Σ (Vobs - Vmodel)² / (errV² + σ_int²)
    """
    sigma_sq = errV**2 + sigma_int_sq
    residuals_sq = (Vobs - Vmodel)**2
    return np.mean(residuals_sq / sigma_sq)


# ─────────────────────────────────────────────────────────
# 6.  Hermes model:  V_model(R)
# ─────────────────────────────────────────────────────────
def hermes_model(data, t50_gyr, g98_kms2kpc):
    """
    g_model(R)  = g_bar(R) · [1 + β · φ(R)]
    V_model(R)  = sqrt( g_model(R) · R )

    Returns V_model array and the computed beta.
    """
    R = data['R']
    gbar = compute_gbar(data)

    # Gate φ(R) from baryonic profile
    phi = hermes_phi(R, gbar)

    # Age score β from age + peak acceleration
    beta = hermes_beta(t50_gyr, g98_kms2kpc)

    # Model acceleration
    g_model = gbar * (1.0 + beta * phi)
    g_model = np.maximum(g_model, 0.0)

    V_model = np.sqrt(g_model * R)
    return V_model, beta, phi


# ─────────────────────────────────────────────────────────
# 7.  MOND model:  simple interpolation function
# ─────────────────────────────────────────────────────────
def mond_model(data):
    """
    g_MOND(R)  = g_bar · ν(g_bar / a₀)
    where  ν(y) = (1 + √(1 + 4/y)) / 2    [simple interpolation function]

    V_MOND(R)  = sqrt( g_MOND(R) · R )
    """
    R = data['R']
    gbar = compute_gbar(data)

    # Convert g_bar to m/s^2 for MOND threshold comparison
    gbar_si = np.abs(gbar) * KMS2_KPC_TO_MS2
    y = gbar_si / A0_MOND
    y = np.maximum(y, 1e-30)

    # Simple interpolation:  ν(y) = (1 + √(1 + 4/y)) / 2
    nu = 0.5 * (1.0 + np.sqrt(1.0 + 4.0 / y))

    g_mond = np.abs(gbar) * nu
    # Restore sign
    g_mond = np.where(gbar >= 0, g_mond, -g_mond)
    g_mond = np.maximum(g_mond, 0.0)

    V_mond = np.sqrt(g_mond * R)
    return V_mond


# ─────────────────────────────────────────────────────────
# 8.  Galaxy-name normalizer  (CSV name → SPARC filename stem)
# ─────────────────────────────────────────────────────────
def normalize_galaxy_name(csv_name):
    """
    Convert CSV galaxy name to the SPARC filename stem.
    'NGC 5371' → 'NGC5371',  'DDO 64' → 'DDO064',  'NGC0100' → 'NGC0100'
    """
    name = csv_name.strip().replace(' ', '')
    # Already matches most cases.  Handle DDO 64 -> DDO064 etc.
    # Split prefix letters from trailing digits/hyphens
    return name


def find_rotmod(sparc_dir, csv_name):
    """Find the rotmod file for a galaxy, trying several name variants."""
    base = normalize_galaxy_name(csv_name)
    candidate = os.path.join(sparc_dir, f"{base}_rotmod.dat")
    if os.path.isfile(candidate):
        return candidate

    # Try zero-padding numeric suffixes  (DDO 64 -> DDO064)
    parts = csv_name.strip().split()
    if len(parts) == 2:
        prefix, num = parts
        try:
            int(num)
            for width in (3, 4, 5):
                padded = prefix + num.zfill(width)
                candidate = os.path.join(sparc_dir, f"{padded}_rotmod.dat")
                if os.path.isfile(candidate):
                    return candidate
        except ValueError:
            pass

    # Case-insensitive fallback
    target = base.lower() + "_rotmod.dat"
    for fn in os.listdir(sparc_dir):
        if fn.lower() == target:
            return os.path.join(sparc_dir, fn)

    return None


# ─────────────────────────────────────────────────────────
# 9.  Main pipeline
# ─────────────────────────────────────────────────────────
def run(sparc_dir, age_csv, output_csv=None):
    # Read age table
    ages = {}
    with open(age_csv, newline='', encoding='utf-8-sig') as f:
        reader = csv.DictReader(f)
        for row in reader:
            galaxy = row['galaxy'].strip()
            t50 = float(row['t50_gyr'])
            g98 = float(row['g98'])
            ages[galaxy] = (t50, g98)

    results = []
    n_found = 0
    n_missing = 0

    for galaxy, (t50, g98) in sorted(ages.items()):
        path = find_rotmod(sparc_dir, galaxy)
        if path is None:
            print(f"  WARNING: no rotmod file for {galaxy}", file=sys.stderr)
            n_missing += 1
            continue

        data = read_rotmod(path)
        if len(data['R']) < 3:
            print(f"  WARNING: too few points for {galaxy}", file=sys.stderr)
            continue

        n_found += 1

        # Hermes prediction
        V_hermes, beta, phi = hermes_model(data, t50, g98)
        chi2_hermes = chi2_nu(data['Vobs'], V_hermes, data['errV'])

        # MOND prediction
        V_mond = mond_model(data)
        chi2_mond = chi2_nu(data['Vobs'], V_mond, data['errV'])

        results.append({
            'galaxy':       galaxy,
            't50_gyr':      t50,
            'g98':          g98,
            'beta_eff':     beta,
            'n_points':     len(data['R']),
            'chi2nu_hermes': chi2_hermes,
            'chi2nu_mond':  chi2_mond,
        })

    # ── Print summary table ──
    print(f"\nHermes Equation Verification  ({n_found} galaxies, {n_missing} missing)")
    print("=" * 95)
    print(f"{'Galaxy':<18s} {'t50':>5s} {'g98':>10s} {'beta':>8s}  {'N':>3s}"
          f"  {'chi2_H':>8s} {'chi2_M':>8s}  {'winner':>6s}")
    print("-" * 95)

    chi2_h_list = []
    chi2_m_list = []
    hermes_wins = 0

    for r in results:
        winner = "Hermes" if r['chi2nu_hermes'] < r['chi2nu_mond'] else "MOND"
        if r['chi2nu_hermes'] < r['chi2nu_mond']:
            hermes_wins += 1
        chi2_h_list.append(r['chi2nu_hermes'])
        chi2_m_list.append(r['chi2nu_mond'])
        print(f"{r['galaxy']:<18s} {r['t50_gyr']:5.1f} {r['g98']:10.1f} {r['beta_eff']:8.4f}"
              f"  {r['n_points']:3d}  {r['chi2nu_hermes']:8.3f} {r['chi2nu_mond']:8.3f}"
              f"  {winner:>6s}")

    print("-" * 95)
    h_arr = np.array(chi2_h_list)
    m_arr = np.array(chi2_m_list)
    print(f"{'Median chi2_nu':<30s}  Hermes: {np.median(h_arr):8.3f}   MOND: {np.median(m_arr):8.3f}")

    # 5% trimmed mean
    def trimmed_mean(x, pct=0.05):
        s = np.sort(x)
        n = len(s)
        lo = int(np.floor(n * pct))
        hi = n - lo
        return np.mean(s[lo:hi])

    print(f"{'Trimmed mean (5%)':<30s}  Hermes: {trimmed_mean(h_arr):8.3f}   MOND: {trimmed_mean(m_arr):8.3f}")
    print(f"{'Hermes wins':<30s}  {hermes_wins} / {len(results)} ({100*hermes_wins/len(results):.1f}%)")

    below5_h = np.sum(h_arr < 5.0)
    below5_m = np.sum(m_arr < 5.0)
    print(f"{'Below chi2=5.0':<30s}  Hermes: {below5_h} ({100*below5_h/len(results):.1f}%)"
          f"   MOND: {below5_m} ({100*below5_m/len(results):.1f}%)")

    # ── Optional CSV output ──
    if output_csv:
        with open(output_csv, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=[
                'galaxy', 't50_gyr', 'g98', 'beta_eff', 'n_points',
                'chi2nu_hermes', 'chi2nu_mond'])
            writer.writeheader()
            writer.writerows(results)
        print(f"\nResults written to {output_csv}")

    return results


# ─────────────────────────────────────────────────────────
if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Verify Hermes equation results against SPARC rotation curves.")
    parser.add_argument('--sparc', required=True,
                        help='Directory containing *_rotmod.dat files')
    parser.add_argument('--ages', required=True,
                        help='CSV with galaxy, t50_gyr, g98 columns')
    parser.add_argument('--csv', default=None,
                        help='Optional output CSV path')
    args = parser.parse_args()

    run(args.sparc, args.ages, args.csv)
