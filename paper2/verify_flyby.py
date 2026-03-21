#!/usr/bin/env python3
"""
Flyby Anomaly Geometric Score Verification — Paper 2
=====================================================
Independently computes the geometric score S for all 12 documented
Earth flybys (1990-2013) and verifies against the published results
in McGinty (2026).

The score:
    S = |Delta cos delta| x exp(-h/H) x (V0/V_inf)^p

Fixed parameters (no per-mission tuning):
    H  = 2500 km   (altitude scale height)
    V0 = 10 km/s   (reference velocity)
    p  = 1.0       (velocity exponent)

Usage:
    python verify_flyby.py                     # full verification
    python verify_flyby.py --mission NEAR      # single mission
    python verify_flyby.py --csv results.csv   # write CSV output
    python verify_flyby.py --p 0.5             # sensitivity check

Reads: flyby_data.csv (same directory)
"""

import argparse
import csv
import io
import math
import os
import sys

# Force UTF-8 output on Windows
if sys.stdout.encoding != "utf-8":
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
    sys.stderr = io.TextIOWrapper(sys.stderr.buffer, encoding="utf-8", errors="replace")

# =====================================================================
#  Fixed parameters (Section 2.6 of the paper)
# =====================================================================
H_DEFAULT = 2500.0   # km — altitude scale height
V0 = 10.0            # km/s — reference velocity
P_DEFAULT = 1.0      # velocity exponent


# =====================================================================
#  Published results from Section 4.1 (for verification)
# =====================================================================
PUBLISHED_SCORES = {
    "Galileo I":   {"abs_dcos": 0.1486, "exp_term": 0.6811, "vel_term": 1.1173, "S": 0.1131},
    "Galileo II":  {"abs_dcos": 0.1699, "exp_term": 0.8859, "vel_term": 1.1261, "S": 0.1695},
    "NEAR":        {"abs_dcos": 0.6254, "exp_term": 0.8061, "vel_term": 1.4599, "S": 0.7359},
    "Cassini":     {"abs_dcos": 0.0215, "exp_term": 0.6250, "vel_term": 0.6246, "S": 0.0084},
    "MESSENGER":   {"abs_dcos": 0.0044, "exp_term": 0.3911, "vel_term": 2.4655, "S": 0.0042},
    "Rosetta I":   {"abs_dcos": 0.1726, "exp_term": 0.4573, "vel_term": 2.5907, "S": 0.2045},
    "Rosetta II":  {"abs_dcos": 0.0345, "exp_term": 0.1200, "vel_term": 1.0663, "S": 0.0044},
    "Rosetta III": {"abs_dcos": 0.0375, "exp_term": 0.3707, "vel_term": 1.0667, "S": 0.0148},
    "Juno":        {"abs_dcos": 0.1979, "exp_term": 0.7990, "vel_term": 0.9615, "S": 0.1521},
    "EPOXI I":     {"abs_dcos": 0.0349, "exp_term": 0.00198, "vel_term": 2.7778, "S": 0.00019},
    "EPOXI II":    {"abs_dcos": 0.4185, "exp_term": 2.8e-8,  "vel_term": 2.7778, "S": 3.3e-8},
    "EPOXI III":   {"abs_dcos": 0.4593, "exp_term": 5.1e-6,  "vel_term": 2.8571, "S": 6.6e-6},
}


# =====================================================================
#  Core computation
# =====================================================================

def compute_score(delta_in_deg, delta_out_deg, h_km, v_inf_kms,
                  H=H_DEFAULT, p=P_DEFAULT):
    """
    Compute the geometric flyby score S.

    Returns dict with all intermediate values for auditing.
    """
    # Step 1: Geometric asymmetry
    delta_in_rad = math.radians(delta_in_deg)
    delta_out_rad = math.radians(delta_out_deg)
    cos_in = math.cos(delta_in_rad)
    cos_out = math.cos(delta_out_rad)
    delta_cos = cos_in - cos_out
    abs_delta_cos = abs(delta_cos)
    sign = "+" if delta_cos > 0 else ("-" if delta_cos < 0 else "0")

    # Step 2: Altitude attenuation
    exp_term = math.exp(-h_km / H)

    # Step 3: Velocity weighting
    vel_term = (V0 / v_inf_kms) ** p

    # Step 4: Multiply
    S = abs_delta_cos * exp_term * vel_term

    return {
        "cos_in": cos_in,
        "cos_out": cos_out,
        "delta_cos": delta_cos,
        "abs_delta_cos": abs_delta_cos,
        "sign": sign,
        "exp_term": exp_term,
        "vel_term": vel_term,
        "S": S,
    }


def load_flyby_data(csv_path):
    """Load flyby data from CSV."""
    data = []
    with open(csv_path, "r") as f:
        reader = csv.DictReader(f)
        for row in reader:
            data.append({
                "number": int(row["number"]),
                "mission": row["mission"],
                "date": row["date"],
                "h_km": float(row["h_km"]),
                "v_inf_kms": float(row["v_inf_kms"]),
                "delta_in_deg": float(row["delta_in_deg"]),
                "delta_out_deg": float(row["delta_out_deg"]),
                "dv_obs_mms": float(row["dv_obs_mms"]),
                "dv_unc_mms": float(row["dv_unc_mms"]) if row["dv_unc_mms"] else None,
                "source": row["source"],
            })
    return data


def classify_observed(dv_obs, dv_unc, mission=""):
    """
    Classify the observed result following the paper's convention:
    - NULL: no anomaly detected (dv ~ 0 or upper bound only)
    - AMBIGUOUS: measured but consistent with zero at 2-sigma (Cassini)
    - ANOMALY: significant detection (|dv| >> uncertainty)
    """
    if abs(dv_obs) < 0.001:
        return "NULL"
    # MESSENGER: +0.02 +/- 0.01 — paper treats as consistent with zero
    if dv_unc and abs(dv_obs) <= 3 * dv_unc and abs(dv_obs) < 0.1:
        return "NULL"
    # Cassini: -2.0 +/- 1.0 — ambiguous (thruster contamination)
    if dv_unc and abs(dv_obs) < 2.5 * dv_unc:
        return "AMBIGUOUS"
    return "ANOMALY"


def main():
    parser = argparse.ArgumentParser(
        description="Verify flyby anomaly geometric scores (Paper 2)"
    )
    parser.add_argument("--mission", "-m", help="Verify a single mission")
    parser.add_argument("--csv", help="Write results to CSV file")
    parser.add_argument("--p", type=float, default=P_DEFAULT,
                        help=f"Velocity exponent (default: {P_DEFAULT})")
    parser.add_argument("--H", type=float, default=H_DEFAULT,
                        help=f"Scale height in km (default: {H_DEFAULT})")
    args = parser.parse_args()

    # Load data
    script_dir = os.path.dirname(os.path.abspath(__file__))
    csv_path = os.path.join(script_dir, "flyby_data.csv")
    flybys = load_flyby_data(csv_path)

    if args.mission:
        flybys = [f for f in flybys if args.mission.lower() in f["mission"].lower()]
        if not flybys:
            print(f"Mission '{args.mission}' not found.")
            sys.exit(1)

    # Compute scores
    results = []
    for fb in flybys:
        score = compute_score(
            fb["delta_in_deg"], fb["delta_out_deg"],
            fb["h_km"], fb["v_inf_kms"],
            H=args.H, p=args.p,
        )
        obs_class = classify_observed(fb["dv_obs_mms"], fb["dv_unc_mms"])

        # Sign prediction check (for anomaly + ambiguous with nonzero dv)
        sign_match = None
        if obs_class in ("ANOMALY", "AMBIGUOUS") and abs(fb["dv_obs_mms"]) > 0.001:
            obs_sign = "+" if fb["dv_obs_mms"] > 0 else "-"
            sign_match = (score["sign"] == obs_sign)

        # Verify against published
        pub = PUBLISHED_SCORES.get(fb["mission"])
        pub_delta = None
        if pub and args.p == P_DEFAULT and args.H == H_DEFAULT:
            pub_delta = score["S"] - pub["S"]

        results.append({
            **fb,
            **score,
            "obs_class": obs_class,
            "sign_match": sign_match,
            "pub_S": pub["S"] if pub else None,
            "pub_delta": pub_delta,
        })

    # ─── Print results ───
    print(f"Flyby Anomaly Geometric Score Verification  (H={args.H}, V0={V0}, p={args.p})")
    print("=" * 110)
    print(f"{'#':>2s} {'Mission':<12s} {'h(km)':>7s} {'V_inf':>6s} {'d_in':>7s} "
          f"{'d_out':>7s} {'|Dcos|':>7s} {'exp':>8s} {'vel':>7s} "
          f"{'S':>9s} {'sign':>4s} {'obs':>8s} {'class':>9s} {'pub_S':>9s}")
    print("-" * 110)

    for r in results:
        pub_s = f"{r['pub_S']:.4f}" if r['pub_S'] is not None and r['pub_S'] > 0.001 else (
            f"{r['pub_S']:.2e}" if r['pub_S'] is not None else "")
        comp_s = f"{r['S']:.4f}" if r['S'] > 0.001 else f"{r['S']:.2e}"

        sign_flag = ""
        if r["sign_match"] is True:
            sign_flag = " ok"
        elif r["sign_match"] is False:
            sign_flag = " !!"

        print(f"{r['number']:>2d} {r['mission']:<12s} {r['h_km']:>7.0f} {r['v_inf_kms']:>6.2f} "
              f"{r['delta_in_deg']:>7.2f} {r['delta_out_deg']:>7.2f} "
              f"{r['abs_delta_cos']:>7.4f} {r['exp_term']:>8.4f} {r['vel_term']:>7.4f} "
              f"{comp_s:>9s} {r['sign']:>4s} "
              f"{r['dv_obs_mms']:>+8.2f} {r['obs_class']:>9s} {pub_s:>9s}{sign_flag}")

    # ─── Detailed steps for single mission ───
    if args.mission and len(results) == 1:
        r = results[0]
        print()
        print("Detailed computation steps:")
        print(f"  delta_in  = {r['delta_in_deg']}deg -> cos(delta_in)  = {r['cos_in']:.6f}")
        print(f"  delta_out = {r['delta_out_deg']}deg -> cos(delta_out) = {r['cos_out']:.6f}")
        print(f"  Delta cos delta = {r['cos_in']:.6f} - {r['cos_out']:.6f} = {r['delta_cos']:+.6f}")
        print(f"  |Delta cos delta| = {r['abs_delta_cos']:.4f}")
        print(f"  Sign of Delta cos delta: {r['sign']} -> predicted DV sign: {r['sign']}")
        print(f"  exp(-{r['h_km']}/{args.H}) = exp({-r['h_km']/args.H:.4f}) = {r['exp_term']:.6f}")
        print(f"  ({V0}/{r['v_inf_kms']})^{args.p} = {r['vel_term']:.6f}")
        print(f"  S = {r['abs_delta_cos']:.4f} x {r['exp_term']:.6f} x {r['vel_term']:.6f} = {r['S']:.6f}")
        if r['pub_S'] is not None:
            print(f"  Published S = {r['pub_S']}")
            if r['pub_delta'] is not None:
                print(f"  Difference  = {r['pub_delta']:+.6f}")

    # ─── Summary ───
    print()
    print("-" * 110)

    # Classification check
    n_anomaly = sum(1 for r in results if r["obs_class"] == "ANOMALY")
    n_null = sum(1 for r in results if r["obs_class"] == "NULL")
    n_ambig = sum(1 for r in results if r["obs_class"] == "AMBIGUOUS")

    # Classify by score threshold
    # Any threshold between 0.02 and 0.11 works
    threshold = 0.05
    correct = 0
    juno_exception = False
    for r in results:
        predicted = "ANOMALY" if r["S"] >= threshold else "NULL"
        actual = r["obs_class"]
        if actual == "AMBIGUOUS":
            actual = "NULL"  # Cassini consistent with null at 2-sigma
        if r["mission"] == "Juno":
            juno_exception = True
            continue  # Juno is the documented exception
        if predicted == actual:
            correct += 1

    n_classified = len(results) - (1 if juno_exception else 0)
    print(f"Classification: {correct}/{n_classified} correct (excluding Juno)")
    juno_results = [r for r in results if r["mission"] == "Juno"]
    if juno_results:
        print(f"  Juno: S = {juno_results[0]['S']:.4f} "
              f"(high-S group, predicted ANOMALY, observed NULL -> GRACE hypothesis)")

    # Sign predictions
    sign_checks = [(r["mission"], r["sign"], r["dv_obs_mms"])
                   for r in results if r["sign_match"] is not None]
    sign_correct = sum(1 for r in results if r["sign_match"] is True)
    sign_total = len(sign_checks)
    print(f"Sign predictions: {sign_correct}/{sign_total} correct")
    for mission, sign, dv in sign_checks:
        obs_sign = "+" if dv > 0 else "-"
        match = "ok" if sign == obs_sign else "MISMATCH"
        print(f"  {mission:<12s}  predicted {sign}  observed {obs_sign}  [{match}]")

    # Published score verification
    if args.p == P_DEFAULT and args.H == H_DEFAULT:
        print()
        print("Published score verification (S values from paper Table 4.1):")
        max_err = 0
        for r in results:
            if r["pub_S"] is not None and r["pub_delta"] is not None:
                err = abs(r["pub_delta"])
                max_err = max(max_err, err)
                if r["pub_S"] > 0.001:
                    rel = err / r["pub_S"] * 100
                    status = "ok" if rel < 1.0 else f"!  ({rel:.1f}%)"
                    print(f"  {r['mission']:<12s}  computed {r['S']:.4f}  "
                          f"published {r['pub_S']:.4f}  diff {r['pub_delta']:+.4f}  {status}")
                else:
                    print(f"  {r['mission']:<12s}  computed {r['S']:.2e}  "
                          f"published {r['pub_S']:.2e}  diff {r['pub_delta']:+.2e}")
        print(f"  Max absolute error: {max_err:.6f}")

    # CSV output
    if args.csv:
        with open(args.csv, "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(["number", "mission", "date", "h_km", "v_inf_kms",
                             "delta_in_deg", "delta_out_deg",
                             "abs_delta_cos", "exp_term", "vel_term", "S",
                             "sign", "dv_obs_mms", "obs_class",
                             "sign_match", "published_S"])
            for r in results:
                writer.writerow([
                    r["number"], r["mission"], r["date"],
                    r["h_km"], r["v_inf_kms"],
                    r["delta_in_deg"], r["delta_out_deg"],
                    f"{r['abs_delta_cos']:.4f}", f"{r['exp_term']:.6f}",
                    f"{r['vel_term']:.4f}", f"{r['S']:.6f}",
                    r["sign"], r["dv_obs_mms"], r["obs_class"],
                    r["sign_match"] if r["sign_match"] is not None else "",
                    r["pub_S"] if r["pub_S"] is not None else "",
                ])
        print(f"\nResults written to {args.csv}")


if __name__ == "__main__":
    main()
