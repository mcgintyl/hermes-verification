#!/usr/bin/env python3
"""
M33 / Corbelli et al. (2014) Format B reconstruction engineering script.

Purpose:
  Build SPARC-style component columns Vgas, Vdisk, Vbulge=0 from public
  Corbelli source-native profiles, without running Hermes or MOND.

Source data:
  Corbelli et al. 2014, A&A 572 A23, online Table 1:
  R_kpc, V_r_kms, sigma_V_kms, Sigma_HI_Msun_pc2, Sigma_star_Msun_pc2.

Kernel:
  Uniform-thickness annular quadrature in the disk plane.
  Gas half-thickness = 0.5 kpc.
  Stellar half-thickness flares linearly from 0.1 kpc at the center to 1.0 kpc at 23 kpc.
  These geometry choices follow the Corbelli et al. description of the source dynamical
  treatment, but this implementation is an engineering reconstruction rather than a
  literal Casertano (1983) code clone.

No model fitting or per-galaxy tuning is performed.
"""

import numpy as np
import pandas as pd
from pathlib import Path

G = 4.30091727003628e-6  # kpc (km/s)^2 / Msun
RMAX_KPC = 23.0
DR_KPC = 0.01
NTHETA = 1440

def load_table(path: str | Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    required = [
        "R_kpc", "V_r_kms", "sigma_V_kms",
        "Sigma_HI_Msun_pc2", "Sigma_star_Msun_pc2"
    ]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns: {missing}")
    return df

def interp_profile(r_grid, R_data, Sigma_data):
    return np.interp(r_grid, R_data, Sigma_data, left=Sigma_data[0], right=Sigma_data[-1])

def compute_disk_v(target_R, r_grid, Sigma_Msun_pc2, h_kpc, ntheta=NTHETA, chunk=128):
    target_R = np.asarray(target_R, float)
    r_grid = np.asarray(r_grid, float)
    Sigma = np.asarray(Sigma_Msun_pc2, float) * 1e6  # Msun/kpc^2
    h = np.asarray(h_kpc, float)
    dr = np.gradient(r_grid)

    theta = np.linspace(0.0, 2.0 * np.pi, ntheta, endpoint=False)
    cos_t = np.cos(theta)
    dtheta = 2.0 * np.pi / ntheta

    out = []
    for R0 in target_R:
        acc = 0.0
        for start in range(0, len(r_grid), chunk):
            end = min(start + chunk, len(r_grid))
            a = r_grid[start:end][:, None]
            H = h[start:end][:, None]
            cos = cos_t[None, :]

            s2 = R0 * R0 + a * a - 2.0 * R0 * a * cos
            numerator = R0 - a * cos
            denom = s2 * np.sqrt(s2 + H * H)
            integrand = np.where(denom > 1e-30, numerator / denom, 0.0)
            angular = integrand.sum(axis=1) * dtheta
            shell = (Sigma[start:end] * r_grid[start:end] * dr[start:end]) * angular
            acc += G * shell.sum()

        V2 = R0 * acc
        out.append(np.sign(V2) * np.sqrt(abs(V2)) if V2 != 0 else 0.0)
    return np.asarray(out)

def mass_from_surface_density(r_grid, Sigma_Msun_pc2):
    dr = np.gradient(r_grid)
    return float(np.sum(2.0 * np.pi * r_grid * dr * Sigma_Msun_pc2 * 1e6))

def reconstruct(table_path: str | Path, output_path: str | Path) -> pd.DataFrame:
    src = load_table(table_path)

    R = src["R_kpc"].to_numpy(float)
    r_grid = np.arange(DR_KPC / 2.0, RMAX_KPC + 1e-12, DR_KPC)

    Sigma_HI = interp_profile(
        r_grid, R, src["Sigma_HI_Msun_pc2"].to_numpy(float)
    )
    Sigma_star = interp_profile(
        r_grid, R, src["Sigma_star_Msun_pc2"].to_numpy(float)
    )
    Sigma_H2 = 10.0 * np.exp(-r_grid / 2.2)
    Sigma_gas_total = 1.33 * (Sigma_HI + Sigma_H2)

    h_gas = np.full_like(r_grid, 0.5)
    h_star = 0.1 + 0.9 * (r_grid / RMAX_KPC)

    Vgas = compute_disk_v(R, r_grid, Sigma_gas_total, h_gas)
    Vdisk = compute_disk_v(R, r_grid, Sigma_star, h_star)
    Vbar = np.sqrt(np.maximum(0.0, Vgas**2 + Vdisk**2))

    out = pd.DataFrame({
        "galaxy": "M33_NGC598",
        "R_kpc": src["R_kpc"],
        "Vobs_kms": src["V_r_kms"],
        "errV_kms": src["sigma_V_kms"],
        "Sigma_HI_Msun_pc2": src["Sigma_HI_Msun_pc2"],
        "Sigma_H2_Msun_pc2": 10.0 * np.exp(-src["R_kpc"].to_numpy(float) / 2.2),
        "Sigma_gas_total_He_Msun_pc2": 1.33 * (
            src["Sigma_HI_Msun_pc2"].to_numpy(float)
            + 10.0 * np.exp(-src["R_kpc"].to_numpy(float) / 2.2)
        ),
        "Sigma_star_Msun_pc2": src["Sigma_star_Msun_pc2"],
        "Vgas_kms": Vgas,
        "Vdisk_kms": Vdisk,
        "Vbulge_kms": 0.0,
        "Vbar_kms": Vbar,
        "kernel": "uniform_thickness_annular_quadrature_source_geometry",
        "source_ml_convention": "Corbelli2014_BVIgi_pixel_SED_stellar_surface_density_not_SPARC_unit_ML",
    })
    out.to_csv(output_path, index=False)
    return out

if __name__ == "__main__":
    base = Path(__file__).resolve().parent
    reconstruct(
        base / "m33_corbelli2014_table1_import.csv",
        base / "m33_corbelli_formatB_reconstructed_board_engineering.csv",
    )
