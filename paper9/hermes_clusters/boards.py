"""Board construction: the baryonic carrier and lensing-node admission.

build_baryonic_carrier and admission_and_covariance, with the constants and helpers they
use, are copied verbatim from the team's sealed board-build script
(run_step2_oos_board_build.py, frozen 2026-08-30). The only change is HERE, which the two
functions use to record file paths in their diagnostics; it points at this package.
"""

from __future__ import annotations

import hashlib
import importlib.util
import math
import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd


HERE = Path(__file__).resolve().parent.parent
G = 4.30091e-6
FLOOR = 1.0 / math.sqrt(2.0 * math.pi)

GAS_CAVEAT_TEXT = (
    "For the 1/r^4 post-Rmax_X extrapolation, source identifies this cluster "
    "as a clear outlier whose gas mass is likely significantly underestimated."
)
NO_CAVEAT_TEXT = "NONE"


@dataclass(frozen=True)
class BoardSpec:
    board: str
    source_name: str
    profile_stem: str
    mistele_stem: str
    track: str
    clean_raw_start: int = 0


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def import_path(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def caveat_fields(board: str) -> dict[str, object]:
    flagged = board in {"MACSJ0429", "MS2137"}
    return {
        "source_gas_underestimate_flag": "YES" if flagged else "NO",
        "source_gas_caveat": GAS_CAVEAT_TEXT if flagged else NO_CAVEAT_TEXT,
        "gas_extrapolation_branch": "BEST_FIT_BETA_EXTENDED",
        "active_baryonic_bias_direction": "NOT_DIRECTLY_ASSIGNED",
        "numerical_correction": "NONE",
        "weight_or_covariance_change": "NONE",
    }


def build_baryonic_carrier(
    spec: BoardSpec,
    cluster_info,
    fraction_interpolator,
    gas_fraction: np.ndarray,
    famaey_root: Path,
) -> tuple[np.ndarray, np.ndarray, dict[str, object]]:
    profile_path = famaey_root / "profiles" / f"results_profile_NOgnfw_all_{spec.profile_stem}.txt"
    values = np.loadtxt(profile_path)
    radius_mpc = np.asarray(values[:, 0], float)
    if len(radius_mpc) != 50 or np.any(~np.isfinite(radius_mpc)) or np.any(np.diff(radius_mpc) <= 0):
        raise RuntimeError(f"{spec.board}: invalid frozen 50-row support")

    bcg_matches = [row for row in cluster_info.cluster_BCG_mass if str(row[0]) == spec.source_name]
    if len(bcg_matches) != 1:
        raise RuntimeError(f"{spec.board}: BCG lookup failure")
    bcg_row = bcg_matches[0]
    cluster_index = int(np.flatnonzero(cluster_info.nameclust == spec.profile_stem)[0])
    gas_index = int(np.flatnonzero(cluster_info.cluster_data[:, 0] == spec.source_name)[0])
    gas = cluster_info.cluster_data[gas_index]
    r200_mpc = float(cluster_info.r200[cluster_index])
    redshift = float(cluster_info.za[cluster_index])
    mbcg = (float(bcg_row[1]) + float(bcg_row[3])) * 1e11
    ne0, r0, alpha, re0, beta0 = (float(gas[index]) for index in range(2, 7))
    ne1, re1, beta1 = (float(gas[index]) for index in range(7, 10))
    gas_model = "double_beta" if ne1 != 0.0 else "single_beta"
    if spec.source_name == "MACS1720":
        ne1 = 0.0
        re1 = 0.0
        beta1 = 0.0
        gas_model = "single_beta_source_directed"

    x = radius_mpc / r200_mpc
    fgal = np.full_like(x, 1.0 / np.max(gas_fraction))
    inner = x < 1.0
    fgal[inner] = 1.0 / fraction_interpolator(x[inner])
    mgas = np.asarray(
        cluster_info.mgas(
            radius_mpc * 1000.0,
            ne0,
            r0,
            alpha,
            re0,
            beta0,
            ne1,
            re1,
            beta1,
        ),
        float,
    )
    mbar = mbcg + mgas * (1.0 + fgal)
    if np.any(~np.isfinite(mbar)) or np.any(mbar <= 0):
        raise RuntimeError(f"{spec.board}: invalid Mbar carrier")
    metadata = {
        "profile_path": str(profile_path.relative_to(HERE)),
        "profile_sha256": sha256(profile_path),
        "native_rows": len(radius_mpc),
        "support_min_mpc": float(radius_mpc[0]),
        "support_max_mpc": float(radius_mpc[-1]),
        "support_min_kpc": float(radius_mpc[0] * 1000.0),
        "support_max_kpc": float(radius_mpc[-1] * 1000.0),
        "u_min": float(np.log(radius_mpc[0] * 1000.0)),
        "u_max": float(np.log(radius_mpc[-1] * 1000.0)),
        "u_span": float(np.log(radius_mpc[-1] / radius_mpc[0])),
        "redshift": redshift,
        "r200_mpc": r200_mpc,
        "mbcg_companions_msun": float(mbcg),
        "gas_model": gas_model,
        "gas_fit_points": int(float(gas[1])),
        "baryonic_formula": "M_BCG+companions + M_gas*(1+f_gal)",
    }
    return radius_mpc * 1000.0, mbar, metadata


def admission_and_covariance(
    spec: BoardSpec,
    radius_native_kpc: np.ndarray,
    mbar_native: np.ndarray,
    mistele_root: Path,
) -> tuple[pd.DataFrame, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, dict[str, object]]:
    mass_path = mistele_root / f"{spec.mistele_stem}.csv"
    corr_path = mistele_root / f"{spec.mistele_stem}-M-corr.csv"
    raw = pd.read_csv(mass_path)
    corr = np.loadtxt(corr_path, delimiter=",")
    if corr.shape != (len(raw), len(raw)):
        raise RuntimeError(f"{spec.board}: correlation dimension mismatch")

    ledger_rows: list[dict[str, object]] = []
    selected: list[int] = []
    support_min = float(radius_native_kpc[0] / 1000.0)
    support_max = float(radius_native_kpc[-1] / 1000.0)
    for index, row in raw.iterrows():
        radius = float(row["r_Mpc"])
        mass = float(row["M_Msun"])
        sigma = float(row["M_stat_err_Msun"])
        density = float(row["ρ_Msun_Mpc3"])
        source_warned = spec.board == "RXJ2129-PARTIAL" and index < spec.clean_raw_start
        reasons: list[str] = []
        if source_warned:
            reasons.append("AUTHOR_WARNED_NODE")
        if not np.isfinite(radius) or radius <= 0:
            reasons.append("MECHANICAL_INVALID_RADIUS")
        if not np.isfinite(mass) or mass <= 0:
            reasons.append("MECHANICAL_NONPOSITIVE_MASS")
        if not np.isfinite(sigma) or sigma <= 0:
            reasons.append("MECHANICAL_NONPOSITIVE_STAT_SIGMA")
        within_support = bool(np.isfinite(radius) and support_min <= radius <= support_max)
        if not within_support:
            reasons.append("OUTSIDE_BARYONIC_NATIVE_SUPPORT")
        admitted = len(reasons) == 0
        if admitted:
            selected.append(int(index))
        base = {
            "board": spec.board,
            "track": spec.track,
            "raw_index_zero_based": int(index),
            "radius_mpc": radius,
            "M_lens_Msun": mass,
            "M_stat_sigma_Msun": sigma,
            "rho_Msun_Mpc3": density,
            "rho_negative_diagnostic": bool(np.isfinite(density) and density < 0),
            "source_warned_node": source_warned,
            "inside_baryonic_native_support": within_support,
            "admitted": admitted,
            "drop_reason_codes": "NONE" if admitted else "|".join(reasons),
            **caveat_fields(spec.board),
        }
        ledger_rows.append(base)

    indices = np.asarray(selected, int)
    if spec.board == "RXJ2129-PARTIAL" and len(indices) < 2:
        status = "DIAGNOSTIC_ONLY" if len(indices) == 1 else "HOLD_NO_SURVIVORS"
    else:
        status = "BUILT" if len(indices) >= 2 else "MECHANICALLY_UNUSABLE"
    if status != "BUILT":
        raise RuntimeError(f"{spec.board}: frozen build did not qualify ({status})")

    radius_kpc = raw.loc[indices, "r_Mpc"].to_numpy(float) * 1000.0
    lens = raw.loc[indices, "M_Msun"].to_numpy(float)
    sigma = raw.loc[indices, "M_stat_err_Msun"].to_numpy(float)
    mbar = np.exp(
        np.interp(np.log(radius_kpc), np.log(radius_native_kpc), np.log(mbar_native))
    )
    corr_selected = corr[np.ix_(indices, indices)]
    covariance = corr_selected * np.outer(sigma, sigma)
    asymmetry = float(np.max(np.abs(covariance - covariance.T)))
    covariance_scale = max(1.0, float(np.max(np.abs(covariance))))
    if asymmetry > 1e-8 * covariance_scale:
        raise RuntimeError(f"{spec.board}: covariance not symmetric")
    eigenvalues = np.linalg.eigvalsh((covariance + covariance.T) / 2.0)
    if np.min(eigenvalues) <= 0:
        raise RuntimeError(f"{spec.board}: covariance not positive definite")
    np.linalg.cholesky((covariance + covariance.T) / 2.0)
    diagnostics = {
        "mass_source_path": str(mass_path.relative_to(HERE)),
        "mass_source_sha256": sha256(mass_path),
        "corr_source_path": str(corr_path.relative_to(HERE)),
        "corr_source_sha256": sha256(corr_path),
        "raw_node_count": len(raw),
        "admitted_raw_indices_zero_based": indices.tolist(),
        "admitted_node_count": len(indices),
        "conditional_dof": len(indices) - 1,
        "covariance_max_abs_asymmetry": asymmetry,
        "covariance_min_eigenvalue_symmetrized": float(np.min(eigenvalues)),
        "covariance_cholesky_spd": True,
        "board_status": status,
    }
    return (
        pd.DataFrame(ledger_rows),
        indices,
        radius_kpc,
        lens,
        sigma,
        mbar,
        covariance,
        diagnostics,
    )
