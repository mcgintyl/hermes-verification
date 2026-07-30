#!/usr/bin/env python3
"""
Reproduce the HISTORICAL Stage-1 beta lock and document the current-gate value.

Runs the single-global-beta Tier-1 fit (unweighted mean of per-galaxy chi2/N
over the seven young galaxies) under two gates:

  * the historical means gate (historical/tier1_gate_v1.py) -> beta = 2.806585
    == the documented Stage-1 lock 2.807 (rounded 2.81);
  * the current shipped Config G gate (paper1/verify_hermes.hermes_phi)
    -> beta = 2.508314  (diagnostic; NOT the historical lock).

Also checks the F583-4 golden phi vector (means gate) from the Tier-1 spec.

The 2.807-vs-2.508 gap is a gate-version change, not a fitting error; see
docs/gate_version_history.md. Usage:

    python historical/fit_tier1_beta.py --sparc /path/to/rotmod_dir
"""
import sys
import os
import argparse

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

import numpy as np
from scipy.optimize import minimize_scalar

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "paper1"))

import importlib.util
_vh_path = os.path.join(os.path.dirname(HERE), "paper1", "verify_hermes.py")
_spec = importlib.util.spec_from_file_location("verify_hermes_hist", _vh_path)
vh = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(vh)

import tier1_gate_v1 as t1

LEGACY_BETA = 2.806585
CURRENT_BETA = 2.508314
TOL_BETA = 1e-3
TOL_GOLDEN = 1e-3


def chi2_nu(data, gbar, phi, beta):
    R = data["R"]
    a = gbar * (1.0 + beta * phi)
    v_model = np.sqrt(np.clip(a * R, 0.0, None))
    denom = data["errV"] ** 2 + t1.SIGMA_INT2
    return np.sum((data["Vobs"] - v_model) ** 2 / denom) / len(R)


def fit_global_beta(boards):
    """boards: list of (data, gbar, phi). Returns beta minimizing mean chi2/N."""
    res = minimize_scalar(
        lambda b: np.mean([chi2_nu(d, g, p, b) for d, g, p in boards]),
        bounds=(-1.0, 15.0), method="bounded",
    )
    return res.x


def main():
    ap = argparse.ArgumentParser(description="Reproduce the historical Stage-1 beta lock.")
    ap.add_argument("--sparc", required=True, help="Directory with *_rotmod.dat files")
    args = ap.parse_args()

    legacy_boards, current_boards = [], []
    for name in t1.TIER1_SEVEN:
        path = vh.find_rotmod(args.sparc, name)
        if path is None:
            print(f"  ERROR: rotmod for {name} not found under {args.sparc}")
            sys.exit(2)
        data = vh.read_rotmod(path)
        # historical: unsigned gas + means gate
        g_hist = t1.gbar_unsigned(data)
        legacy_boards.append((data, g_hist, t1.tier1_gate(data["R"], g_hist)))
        # current: shipped gbar (signed) + shipped medians gate
        g_cur = vh.compute_gbar(data)
        current_boards.append((data, g_cur, vh.hermes_phi(data["R"], g_cur)))

    beta_legacy = fit_global_beta(legacy_boards)
    beta_current = fit_global_beta(current_boards)

    # F583-4 golden vector check (historical means gate)
    f_path = vh.find_rotmod(args.sparc, "F583-4")
    f_data = vh.read_rotmod(f_path)
    f_phi = t1.tier1_gate(f_data["R"], t1.gbar_unsigned(f_data))
    golden_max = float(np.max(np.abs(np.asarray(f_phi) - np.asarray(t1.F583_4_GOLDEN_PHI))))

    print("=" * 70)
    print("HISTORICAL Stage-1 beta reconstruction")
    print("=" * 70)
    print(f"  historical means gate : global beta = {beta_legacy:.6f}   (lock 2.807 -> 2.81)")
    print(f"  current Config G gate : global beta = {beta_current:.6f}   (diagnostic, not the lock)")
    print(f"  F583-4 golden phi     : max|diff| = {golden_max:.2e}")
    print()

    ok = True
    def chk(label, cond):
        nonlocal ok
        print(f"  [{'PASS' if cond else 'FAIL'}] {label}")
        ok = ok and cond

    chk(f"historical beta == {LEGACY_BETA} (+/- {TOL_BETA})", abs(beta_legacy - LEGACY_BETA) < TOL_BETA)
    chk(f"current beta    == {CURRENT_BETA} (+/- {TOL_BETA})", abs(beta_current - CURRENT_BETA) < TOL_BETA)
    chk(f"F583-4 golden phi within {TOL_GOLDEN}", golden_max < TOL_GOLDEN)

    print()
    print("  RESULT:", "PASS" if ok else "FAIL",
          "-- historical gate reproduces the Stage-1 lock; current gate documented.")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
