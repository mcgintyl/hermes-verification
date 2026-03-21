#!/usr/bin/env python3
"""
Chirality Paper Verification Tool — Paper 4
=============================================
Verifies internal consistency of all supplementary data files for:

  McGinty (2026). "The Emergent Plane at Cosmological Scale: A Geometric
  Chirality Prediction, Empirical Exclusion Limits, and a Falsifiable Roadmap."

Checks performed:
  1. Cross-match counts: ALFALFA membership totals are self-consistent
  2. Histogram integrity: bin counts sum to reported N for each group
  3. Summary statistics: N values match across files
  4. Cluster chirality: CW + CCW = N_galaxies, p-values recomputed
  5. Excess variance: chi-squared and p-values verified
  6. Per-void summary: total sources match cross-match counts
  7. Dark candidates: header-only file confirms zero candidates
  8. Paper Table 1: published numbers match CSV data

Usage:
    python verify_chirality.py           # run all checks
    python verify_chirality.py --verbose # show every intermediate check

Reads: all paper4_*.csv files in the same directory.
"""

import argparse
import csv
import io
import math
import os
import sys
from collections import defaultdict

# Force UTF-8 output on Windows
if sys.stdout.encoding != "utf-8":
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
    sys.stderr = io.TextIOWrapper(sys.stderr.buffer, encoding="utf-8", errors="replace")

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))


def load_csv(filename):
    """Load a CSV file from the script directory."""
    path = os.path.join(SCRIPT_DIR, filename)
    with open(path, "r", newline="") as f:
        reader = csv.DictReader(f)
        return list(reader)


# =====================================================================
#  Check 1: Cross-match counts
# =====================================================================

def check_crossmatch_counts(verbose=False):
    """Verify ALFALFA x VAST membership counts are self-consistent."""
    rows = load_csv("paper4_step2_crossmatch_counts.csv")
    counts = {r["metric"]: int(r["value"]) for r in rows}

    checks = []
    n_total = counts["N_total_ALFALFA"]
    n_interior = counts["N_in_deep_void_interior_Rlt0p5"]
    n_body = counts["N_in_deep_void_body_Rlt1"]
    n_shell = counts["N_in_deep_void_shell_0p75to1p25"]
    n_wall = counts["N_wall_not_in_any_void"]

    # Paper claims: 31,502 total ALFALFA sources
    checks.append(("ALFALFA total = 31,502", n_total == 31502, f"got {n_total}"))
    # Paper claims: 39 deep-void-interior sources
    checks.append(("Deep void interior = 39", n_interior == 39, f"got {n_interior}"))
    # Interior must be subset of body
    checks.append(("Interior <= Body", n_interior <= n_body,
                    f"{n_interior} vs {n_body}"))
    # Wall + any-void should account for total (approximately; some in
    # non-deep voids are not counted as wall)
    # The paper notes: "any void" includes all voids, not just deep ones
    checks.append(("Body count = 939", n_body == 939, f"got {n_body}"))
    checks.append(("Shell count = 2,087", n_shell == 2087, f"got {n_shell}"))
    checks.append(("Wall count = 23,374", n_wall == 23374, f"got {n_wall}"))

    if verbose:
        for name, ok, detail in checks:
            print(f"    {'ok' if ok else 'FAIL':>4s}  {name}  ({detail})")

    return checks


# =====================================================================
#  Check 2: Histogram bin counts sum to reported N
# =====================================================================

def check_histograms(verbose=False):
    """Verify histogram bin counts sum to the correct N per group."""
    checks = []

    # Expected N per group from the summary stats file
    stats = load_csv("paper4_hi_property_summary_stats.csv")
    expected_n = {}
    for row in stats:
        group = row["group"]
        expected_n[(group, "logMHI")] = int(row["logMHI_N"])
        expected_n[(group, "W50")] = int(row["W50_N"])
        expected_n[(group, "logMHI_Ms")] = int(row["logMHI_Ms_N"])

    hist_files = {
        "logMHI": "paper4_hist_logMHI.csv",
        "W50": "paper4_hist_W50.csv",
        "logMHI_Ms": "paper4_hist_logMHI_over_Ms.csv",
    }

    for prop, filename in hist_files.items():
        rows = load_csv(filename)
        group_sums = defaultdict(int)
        for row in rows:
            group_sums[row["group"]] += int(row["count"])

        for group, total in group_sums.items():
            key = (group, prop)
            if key in expected_n:
                exp = expected_n[key]
                diff = total - exp
                # Allow off-by-few for bin edge effects (sources outside range)
                ok = abs(diff) <= max(2, int(0.001 * exp))
                exact = total == exp
                label = "ok" if exact else ("~ok (edge)" if ok else "FAIL")
                checks.append((f"hist {prop} [{group}] sum={total} vs N={exp}",
                                ok, f"diff={diff}"))
                if verbose:
                    print(f"    {label:>12s}  hist {prop} [{group}]: "
                          f"sum={total}, expected N={exp}, diff={diff}")
            else:
                if verbose:
                    print(f"    skip  hist {prop} [{group}]: no expected N")

    return checks


# =====================================================================
#  Check 3: Summary stats N values match cross-match counts
# =====================================================================

def check_summary_stats_vs_crossmatch(verbose=False):
    """Verify summary stats N matches cross-match counts where applicable."""
    checks = []

    xm = load_csv("paper4_step2_crossmatch_counts.csv")
    xm_counts = {r["metric"]: int(r["value"]) for r in xm}

    stats = load_csv("paper4_hi_property_summary_stats.csv")
    stats_by_group = {r["group"]: r for r in stats}

    # deep_interior_Rlt0p5 should have N = 39
    if "deep_interior_Rlt0p5" in stats_by_group:
        n_stats = int(stats_by_group["deep_interior_Rlt0p5"]["logMHI_N"])
        n_xm = xm_counts["N_in_deep_void_interior_Rlt0p5"]
        ok = n_stats == n_xm
        checks.append((f"Summary N(deep interior) = crossmatch N", ok,
                        f"stats={n_stats}, xm={n_xm}"))
        if verbose:
            print(f"    {'ok' if ok else 'FAIL':>4s}  Summary deep_interior N={n_stats} "
                  f"vs crossmatch={n_xm}")

    if "deep_body_Rlt1" in stats_by_group:
        n_stats = int(stats_by_group["deep_body_Rlt1"]["logMHI_N"])
        n_xm = xm_counts["N_in_deep_void_body_Rlt1"]
        ok = n_stats == n_xm
        checks.append((f"Summary N(deep body) = crossmatch N", ok,
                        f"stats={n_stats}, xm={n_xm}"))
        if verbose:
            print(f"    {'ok' if ok else 'FAIL':>4s}  Summary deep_body N={n_stats} "
                  f"vs crossmatch={n_xm}")

    return checks


# =====================================================================
#  Check 4: Cluster chirality — CW + CCW = N_galaxies, p-values
# =====================================================================

def check_cluster_chirality(verbose=False):
    """Verify CW+CCW=N, recompute binomial p-values."""
    checks = []
    rows = load_csv("paper4_pure_cluster_member_cuts_summary.csv")

    for row in rows:
        cut = row["cut"]
        n_gal = int(row["N_galaxies"])
        cw = int(row["CW_pos"])
        ccw = int(row["CCW_neg"])

        # CW + CCW must equal N_galaxies
        total = cw + ccw
        ok = total == n_gal
        checks.append((f"CW+CCW={total} vs N={n_gal} [{cut}]", ok, ""))
        if verbose:
            print(f"    {'ok' if ok else 'FAIL':>4s}  {cut}: CW={cw} + CCW={ccw} = {total}, "
                  f"N_galaxies={n_gal}")

        # Verify global binomial p-value direction:
        # p should be > 0.05 for all tests (null result)
        p_binom = float(row["global_binom_p"])
        # The paper's conclusion: no significant result
        # We just verify p > 0.01 (all should be well above)
        p_ok = p_binom > 0.01
        checks.append((f"Binomial p={p_binom:.4f} > 0.01 [{cut}]", p_ok, ""))

        # Verify chi2 df matches N_systems
        n_sys = int(row["N_systems_ge_minN"])
        df = int(row["df"])
        df_ok = df == n_sys
        checks.append((f"df={df} == N_systems={n_sys} [{cut}]", df_ok, ""))
        if verbose and not df_ok:
            print(f"    FAIL  {cut}: df={df} != N_systems={n_sys}")

    return checks


# =====================================================================
#  Check 5: Excess variance — chi2 p-values sensible
# =====================================================================

def check_excess_variance(verbose=False):
    """Verify excess variance table internal consistency."""
    checks = []
    rows = load_csv("paper4_excess_variance_summary_by_richness.csv")

    for row in rows:
        # Skip empty rows (N_clusters = 0)
        n_clusters = int(row["N_clusters"]) if row["N_clusters"] else 0
        if n_clusters == 0:
            continue

        richness = row["Richness_bin"]
        radius = row["Radius"]
        min_n = int(row["Min_N"])

        # chi2 and df should both be present
        if row["chi2"] and row["df"]:
            chi2 = float(row["chi2"])
            df = float(row["df"])

            # df should equal N_clusters
            df_ok = abs(df - n_clusters) < 0.01
            checks.append((f"Variance df={df:.0f} == N={n_clusters} "
                            f"[{radius},{richness},min{min_n}]", df_ok, ""))

            # chi2/df ratio — should be near 1.0 for null (no excess variance)
            # Very small samples (df < 5) can produce extreme ratios by chance
            ratio = chi2 / df if df > 0 else 0
            if df < 5:
                ratio_ok = True  # too few clusters for meaningful ratio
            else:
                ratio_ok = 0.3 < ratio < 3.0  # generous bounds for null
            checks.append((f"chi2/df={ratio:.2f} near 1.0 "
                            f"[{radius},{richness},min{min_n}]", ratio_ok, ""))

            if verbose:
                print(f"    {'ok' if df_ok and ratio_ok else 'WARN':>4s}  "
                      f"{radius} {richness} min{min_n}: "
                      f"N={n_clusters} chi2={chi2:.1f} df={df:.0f} "
                      f"chi2/df={ratio:.2f}")

        # N_p05 should not exceed expected_p05 dramatically
        n_p05 = int(float(row["N_p05"])) if row["N_p05"] else 0
        exp_p05 = float(row["expected_p05"]) if row["expected_p05"] else 0
        if n_clusters >= 10 and exp_p05 > 0:
            # No dramatic excess: N_p05 should be <= 3x expected
            excess_ok = n_p05 <= max(3 * exp_p05, exp_p05 + 10)
            checks.append((f"N(p<0.05)={n_p05} vs expected={exp_p05:.1f} "
                            f"[{radius},{richness},min{min_n}]", excess_ok, ""))

    return checks


# =====================================================================
#  Check 6: Per-void summary — totals match cross-match
# =====================================================================

def check_per_void_summary(verbose=False):
    """Verify per-void source counts sum to cross-match totals."""
    checks = []
    rows = load_csv("paper4_per_void_summary.csv")

    total_interior = 0
    total_body = 0
    n_voids = len(rows)

    for row in rows:
        n_int = int(float(row["N_HI_interior"]))
        n_bod = int(float(row["N_HI_body"]))
        total_interior += n_int
        total_body += n_bod

        # Interior should be <= body for each void
        if n_int > n_bod and n_bod > 0:
            checks.append((f"Void {row['void_id']}: interior <= body",
                            False, f"{n_int} > {n_bod}"))

    # Check against cross-match
    xm = load_csv("paper4_step2_crossmatch_counts.csv")
    xm_counts = {r["metric"]: int(r["value"]) for r in xm}

    # Per-void totals for DEEP voids should match cross-match deep counts
    checks.append((f"N voids in per-void file = {n_voids}", n_voids == 117,
                    f"paper says 117 deepest voids (top 10%)"))

    # Note: per-void sums may not exactly equal cross-match totals because
    # per_void_summary covers only the 117 deepest voids, while cross-match
    # counts cover all deep voids. We check that per-void <= cross-match.
    xm_body = xm_counts["N_in_deep_void_body_Rlt1"]
    xm_interior = xm_counts["N_in_deep_void_interior_Rlt0p5"]
    checks.append((f"Sum(per-void body)={total_body} <= crossmatch {xm_body}",
                    total_body <= xm_body, ""))
    checks.append((f"Sum(per-void interior)={total_interior} <= crossmatch {xm_interior}",
                    total_interior <= xm_interior, ""))

    if verbose:
        print(f"    N voids: {n_voids}")
        print(f"    Sum body: {total_body} (crossmatch: {xm_body})")
        print(f"    Sum interior: {total_interior} (crossmatch: {xm_interior})")

    return checks


# =====================================================================
#  Check 7: Dark candidates — header-only = zero candidates
# =====================================================================

def check_dark_candidates(verbose=False):
    """Verify dark candidates file has header but zero data rows."""
    checks = []
    rows = load_csv("paper4_dark_candidates.csv")
    n = len(rows)
    ok = n == 0
    checks.append((f"Dark candidates: {n} rows (expected 0)", ok, ""))
    if verbose:
        print(f"    {'ok' if ok else 'FAIL':>4s}  Dark candidates: {n} data rows")
    return checks


# =====================================================================
#  Check 8: Paper Table 1 numbers match CSV data
# =====================================================================

def check_paper_table1(verbose=False):
    """Verify numbers from the paper's Table 1 match the CSV files."""
    checks = []
    rows = load_csv("paper4_pure_cluster_member_cuts_summary.csv")

    # Build lookup by cut name
    by_cut = {r["cut"]: r for r in rows}

    # Paper Table 1 entries (manually extracted from the paper):
    # Format: (description, csv_cut_key, expected_n_galaxies, expected_n_systems)
    paper_table1 = [
        ("redMaPPer dz<=0.003 R<=0.5Mpc", "Cut2_dz0.003_R0.5Mpc", 874, 18),
        ("redMaPPer dz<=0.003 R<=1.0Mpc", "Cut2_dz0.003_R1.0Mpc", 1654, 78),
        ("Velocity disp R<=R200", "Cut3_dv<=2sigma_R<=R200c", 2961, 179),
        ("Tempel small groups 5<=N<=20", "Cut4_groups_5<=Ngal<=20", 30971, 1812),
        ("Tempel rich groups N>=30 R<=0.5Mpc", "Cut5_groups_Ngal>=30_Rproj<=0.5Mpc", 2108, 186),
    ]

    for desc, key, exp_gal, exp_sys in paper_table1:
        if key in by_cut:
            row = by_cut[key]
            n_gal = int(row["N_galaxies"])
            n_sys = int(row["N_systems_ge_minN"])
            gal_ok = n_gal == exp_gal
            sys_ok = n_sys == exp_sys
            checks.append((f"Table 1 [{desc}]: N_gal={n_gal} (exp {exp_gal})",
                            gal_ok, ""))
            checks.append((f"Table 1 [{desc}]: N_sys={n_sys} (exp {exp_sys})",
                            sys_ok, ""))
            if verbose:
                print(f"    {'ok' if gal_ok else 'FAIL':>4s}  {desc}: "
                      f"N_gal={n_gal} (exp {exp_gal})")
                print(f"    {'ok' if sys_ok else 'FAIL':>4s}  {desc}: "
                      f"N_sys={n_sys} (exp {exp_sys})")
        else:
            checks.append((f"Table 1 [{desc}]: cut not found in CSV", False, key))
            if verbose:
                print(f"    FAIL  {desc}: cut '{key}' not found")

    # Paper claims: 273,055 total SpArcFiRe galaxies
    # This is the full catalog, not in the cluster CSV — check Tempel total
    tempel_row = by_cut.get("Cut4_groups_5<=Ngal<=20")
    if tempel_row:
        n = int(tempel_row["N_galaxies"])
        checks.append((f"Tempel small groups N={n} matches paper 30,971",
                        n == 30971, ""))

    # Paper: global chirality 49.84% CW, 50.16% CCW
    # Check from Tempel full membership if available
    # (The 273,055 is the full SpArcFiRe catalog, not directly in these CSVs)

    return checks


# =====================================================================
#  Check 9: Binomial p-value recomputation
# =====================================================================

def binomial_p_two_sided(n, k):
    """
    Approximate two-sided binomial p-value for k successes in n trials
    with p=0.5, using the normal approximation.
    """
    if n == 0:
        return 1.0
    expected = n / 2.0
    std = math.sqrt(n) / 2.0
    if std == 0:
        return 1.0
    z = abs(k - expected) / std
    # Two-sided p using error function approximation
    p = math.erfc(z / math.sqrt(2))
    return p


def check_binomial_pvalues(verbose=False):
    """Recompute binomial p-values and compare to published."""
    checks = []
    rows = load_csv("paper4_pure_cluster_member_cuts_summary.csv")

    for row in rows:
        cut = row["cut"]
        n_gal = int(row["N_galaxies"])
        cw = int(row["CW_pos"])
        pub_p = float(row["global_binom_p"])

        # Recompute
        recomp_p = binomial_p_two_sided(n_gal, cw)

        # Allow generous tolerance (normal approx vs exact)
        if pub_p > 0.01:
            ratio = recomp_p / pub_p if pub_p > 0 else float("inf")
            ok = 0.3 < ratio < 3.0
        else:
            ok = recomp_p < 0.05  # both small

        checks.append((f"Binomial p recomp [{cut}]: pub={pub_p:.4f}, "
                        f"recomp={recomp_p:.4f}", ok, ""))
        if verbose:
            print(f"    {'ok' if ok else 'WARN':>4s}  {cut}: "
                  f"CW={cw}/{n_gal}, pub_p={pub_p:.4f}, recomp_p={recomp_p:.4f}")

    return checks


# =====================================================================
#  Main
# =====================================================================

def main():
    parser = argparse.ArgumentParser(
        description="Verify chirality paper supplementary data (Paper 4)"
    )
    parser.add_argument("--verbose", "-v", action="store_true",
                        help="Show every intermediate check")
    args = parser.parse_args()

    all_checks = []

    check_groups = [
        ("Cross-match counts", check_crossmatch_counts),
        ("Histogram integrity", check_histograms),
        ("Summary stats vs cross-match", check_summary_stats_vs_crossmatch),
        ("Cluster chirality (CW+CCW=N, p-values)", check_cluster_chirality),
        ("Excess variance consistency", check_excess_variance),
        ("Per-void summary", check_per_void_summary),
        ("Dark candidates (zero rows)", check_dark_candidates),
        ("Paper Table 1 numbers", check_paper_table1),
        ("Binomial p-value recomputation", check_binomial_pvalues),
    ]

    print("Chirality Paper Data Verification  (Paper 4)")
    print("=" * 75)

    total_pass = 0
    total_fail = 0

    for group_name, check_fn in check_groups:
        print(f"\n  {group_name}")
        print(f"  {'-' * 60}")
        checks = check_fn(verbose=args.verbose)
        all_checks.extend(checks)

        passed = sum(1 for _, ok, _ in checks if ok)
        failed = sum(1 for _, ok, _ in checks if not ok)
        total_pass += passed
        total_fail += failed

        status = "PASS" if failed == 0 else f"FAIL ({failed} issues)"
        print(f"  Result: {passed} passed, {failed} failed  [{status}]")

        if not args.verbose:
            # Show only failures
            for name, ok, detail in checks:
                if not ok:
                    print(f"    FAIL  {name}  {detail}")

    print()
    print("=" * 75)
    print(f"Total: {total_pass} passed, {total_fail} failed")
    if total_fail == 0:
        print("All checks passed.")
    else:
        print(f"{total_fail} check(s) failed. Run with --verbose for details.")

    return 0 if total_fail == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
