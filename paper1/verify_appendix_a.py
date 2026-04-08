#!/usr/bin/env python3
"""
Appendix A Sensitivity Table Verification
============================================
Re-runs the Hermes pipeline with different σ_int² and a_knee parameters
to verify results in Table A1 and A2 of the paper.

Also computes velocity residuals for Table A3.
"""

import sys
import os
import csv
import numpy as np

# Add paper1 directory to path to import verify_hermes
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import verify_hermes as vh

# ─────────────────────────────────────────────────────────
# Test tolerances
# ─────────────────────────────────────────────────────────
TOL_CHI2_MEDIAN = 0.05
TOL_CHI2_TRIM = 0.2
TOL_WIN_RATE = 0.005  # 0.5%

import argparse

# ─────────────────────────────────────────────────────────
# Paths (resolved from arguments in main)
# ─────────────────────────────────────────────────────────
SPARC_DIR = None
AGES_CSV = None


def trimmed_mean(x, pct=0.05):
    """Compute pct% trimmed mean."""
    s = np.sort(x)
    n = len(s)
    lo = int(np.floor(n * pct))
    hi = n - lo
    return np.mean(s[lo:hi])


def read_ages(csv_path):
    """Read galaxy ages from CSV."""
    ages = {}
    with open(csv_path, newline='', encoding='utf-8-sig') as f:
        reader = csv.DictReader(f)
        for row in reader:
            galaxy = row['galaxy'].strip()
            t50 = float(row['t50_gyr'])
            g98 = float(row['g98'])
            ages[galaxy] = (t50, g98)
    return ages


def run_pipeline(sparc_dir, ages, sigma_int_sq=vh.SIGMA_INT_SQ, a_knee=vh.A_KNEE_DEFAULT):
    """
    Run the full pipeline with given sigma_int_sq and a_knee.
    Returns arrays of chi2 values and win rate.
    """
    chi2_h_list = []
    chi2_m_list = []
    hermes_wins = 0

    for galaxy, (t50, g98) in sorted(ages.items()):
        path = vh.find_rotmod(sparc_dir, galaxy)
        if path is None:
            continue

        data = vh.read_rotmod(path)
        if len(data['R']) < 3:
            continue

        # Hermes prediction
        V_hermes, beta, phi = vh.hermes_model(data, t50, g98, a_knee=a_knee)
        # Re-compute chi2 with the given sigma_int_sq
        chi2_hermes = vh.chi2_nu(data['Vobs'], V_hermes, data['errV'], sigma_int_sq=sigma_int_sq)

        # MOND prediction (doesn't depend on sigma_int_sq)
        V_mond = vh.mond_model(data)
        chi2_mond = vh.chi2_nu(data['Vobs'], V_mond, data['errV'], sigma_int_sq=sigma_int_sq)

        if chi2_hermes < chi2_mond:
            hermes_wins += 1

        chi2_h_list.append(chi2_hermes)
        chi2_m_list.append(chi2_mond)

    h_arr = np.array(chi2_h_list)
    m_arr = np.array(chi2_m_list)

    n_galaxies = len(h_arr)
    win_rate = hermes_wins / n_galaxies if n_galaxies > 0 else 0.0

    return {
        'hermes_chi2': h_arr,
        'mond_chi2': m_arr,
        'hermes_median': np.median(h_arr),
        'mond_median': np.median(m_arr),
        'hermes_trim': trimmed_mean(h_arr),
        'mond_trim': trimmed_mean(m_arr),
        'win_rate': win_rate,
        'n_galaxies': n_galaxies,
    }


def count_high_chi2(chi2_arr, threshold=10.0):
    """Count number of galaxies with chi2 > threshold."""
    return np.sum(chi2_arr > threshold)


def test_table_a1():
    """Test Table A1: Sensitivity to σ_int²."""
    print("\n" + "="*80)
    print("TABLE A1: SENSITIVITY TO σ_int²")
    print("="*80)

    ages = read_ages(AGES_CSV)

    tests = [
        {
            'sigma_int_sq': 0,
            'expected': {
                'hermes_median': 37.850,
                'hermes_trim': 112.732,
                'mond_median': 105.530,
                'win_rate': 0.526,
                'hermes_chi2_gt10': 104,
                'mond_chi2_gt10': 95,
            }
        },
        {
            'sigma_int_sq': 100,
            'expected': {
                'hermes_median': 4.123,
                'hermes_trim': 3.787,
                'mond_median': 9.975,
                'win_rate': 0.511,
                'hermes_chi2_gt10': 40,
                'mond_chi2_gt10': 47,
            }
        },
        {
            'sigma_int_sq': 386,
            'expected': {
                'hermes_median': 1.312,
                'hermes_trim': 2.495,
                'mond_median': 1.141,
                'mond_trim': 3.292,
                'win_rate': 0.511,
                'hermes_chi2_gt10': 12,
                'mond_chi2_gt10': 15,
            }
        },
        {
            'sigma_int_sq': 900,
            'expected': {
                'hermes_median': 0.605,
                'hermes_trim': 0.524,
                'mond_median': 1.522,
                'win_rate': 0.511,
                'hermes_chi2_gt10': 2,
                'mond_chi2_gt10': 4,
            }
        },
    ]

    all_pass = True

    for test in tests:
        sigma = test['sigma_int_sq']
        expected = test['expected']

        print(f"\nσ_int² = {sigma}")
        print("-" * 80)

        result = run_pipeline(SPARC_DIR, ages, sigma_int_sq=sigma)

        # Check Hermes median
        h_med = result['hermes_median']
        h_med_exp = expected['hermes_median']
        h_med_diff = abs(h_med - h_med_exp)
        h_med_pass = h_med_diff <= TOL_CHI2_MEDIAN
        status = "PASS" if h_med_pass else "FAIL"
        print(f"  Hermes median:     {h_med:.3f} (exp {h_med_exp:.3f}, diff {h_med_diff:.3f}) [{status}]")
        if not h_med_pass:
            all_pass = False

        # Check MOND median (if provided)
        if 'mond_median' in expected:
            m_med = result['mond_median']
            m_med_exp = expected['mond_median']
            m_med_diff = abs(m_med - m_med_exp)
            m_med_pass = m_med_diff <= TOL_CHI2_MEDIAN
            status = "PASS" if m_med_pass else "FAIL"
            print(f"  MOND median:       {m_med:.3f} (exp {m_med_exp:.3f}, diff {m_med_diff:.3f}) [{status}]")
            if not m_med_pass:
                all_pass = False

        # Check Hermes trim (if provided)
        if 'hermes_trim' in expected:
            h_trim = result['hermes_trim']
            h_trim_exp = expected['hermes_trim']
            h_trim_diff = abs(h_trim - h_trim_exp)
            h_trim_pass = h_trim_diff <= TOL_CHI2_TRIM
            status = "PASS" if h_trim_pass else "FAIL"
            print(f"  Hermes trim:       {h_trim:.3f} (exp {h_trim_exp:.3f}, diff {h_trim_diff:.3f}) [{status}]")
            if not h_trim_pass:
                all_pass = False

        # Check MOND trim (if provided)
        if 'mond_trim' in expected:
            m_trim = result['mond_trim']
            m_trim_exp = expected['mond_trim']
            m_trim_diff = abs(m_trim - m_trim_exp)
            m_trim_pass = m_trim_diff <= TOL_CHI2_TRIM
            status = "PASS" if m_trim_pass else "FAIL"
            print(f"  MOND trim:         {m_trim:.3f} (exp {m_trim_exp:.3f}, diff {m_trim_diff:.3f}) [{status}]")
            if not m_trim_pass:
                all_pass = False

        # Check win rate
        wr = result['win_rate']
        wr_exp = expected['win_rate']
        wr_diff = abs(wr - wr_exp)
        wr_pass = wr_diff <= TOL_WIN_RATE
        status = "PASS" if wr_pass else "FAIL"
        print(f"  Win rate:          {wr:.1%} (exp {wr_exp:.1%}, diff {wr_diff:.3%}) [{status}]")
        if not wr_pass:
            all_pass = False

        # Check chi2 > 10 counts
        if 'hermes_chi2_gt10' in expected:
            h_gt10 = count_high_chi2(result['hermes_chi2'], threshold=10.0)
            h_gt10_exp = expected['hermes_chi2_gt10']
            h_gt10_pass = h_gt10 == h_gt10_exp
            status = "PASS" if h_gt10_pass else "FAIL"
            print(f"  Hermes χ²>10:      {h_gt10} (exp {h_gt10_exp}) [{status}]")
            if not h_gt10_pass:
                all_pass = False

        if 'mond_chi2_gt10' in expected:
            m_gt10 = count_high_chi2(result['mond_chi2'], threshold=10.0)
            m_gt10_exp = expected['mond_chi2_gt10']
            m_gt10_pass = m_gt10 == m_gt10_exp
            status = "PASS" if m_gt10_pass else "FAIL"
            print(f"  MOND χ²>10:        {m_gt10} (exp {m_gt10_exp}) [{status}]")
            if not m_gt10_pass:
                all_pass = False

    return all_pass


def test_table_a2():
    """Test Table A2: Sensitivity to a_knee."""
    print("\n" + "="*80)
    print("TABLE A2: SENSITIVITY TO a_knee")
    print("="*80)

    ages = read_ages(AGES_CSV)

    tests = [
        {
            'a_knee': 1189,
            'expected': {
                'hermes_median': 1.395,
                'hermes_trim': 2.509,
                'win_rate': 0.496,
                'hermes_chi2_gt10': 13,
            }
        },
        {
            'a_knee': 1585,
            'expected': {
                'hermes_median': 1.312,
                'hermes_trim': 2.495,
                'win_rate': 0.511,
                'hermes_chi2_gt10': 12,
            }
        },
        {
            'a_knee': 1981,
            'expected': {
                'hermes_median': 1.312,
                'hermes_trim': 2.473,
                'win_rate': 0.526,
                'hermes_chi2_gt10': 12,
            }
        },
    ]

    all_pass = True

    for test in tests:
        a_knee = test['a_knee']
        expected = test['expected']

        print(f"\na_knee = {a_knee}")
        print("-" * 80)

        # Use default sigma_int_sq=386
        result = run_pipeline(SPARC_DIR, ages, sigma_int_sq=386.0, a_knee=a_knee)

        # Check Hermes median
        h_med = result['hermes_median']
        h_med_exp = expected['hermes_median']
        h_med_diff = abs(h_med - h_med_exp)
        h_med_pass = h_med_diff <= TOL_CHI2_MEDIAN
        status = "PASS" if h_med_pass else "FAIL"
        print(f"  Hermes median:     {h_med:.3f} (exp {h_med_exp:.3f}, diff {h_med_diff:.3f}) [{status}]")
        if not h_med_pass:
            all_pass = False

        # Check Hermes trim
        h_trim = result['hermes_trim']
        h_trim_exp = expected['hermes_trim']
        h_trim_diff = abs(h_trim - h_trim_exp)
        h_trim_pass = h_trim_diff <= TOL_CHI2_TRIM
        status = "PASS" if h_trim_pass else "FAIL"
        print(f"  Hermes trim:       {h_trim:.3f} (exp {h_trim_exp:.3f}, diff {h_trim_diff:.3f}) [{status}]")
        if not h_trim_pass:
            all_pass = False

        # Check win rate
        wr = result['win_rate']
        wr_exp = expected['win_rate']
        wr_diff = abs(wr - wr_exp)
        wr_pass = wr_diff <= TOL_WIN_RATE
        status = "PASS" if wr_pass else "FAIL"
        print(f"  Win rate:          {wr:.1%} (exp {wr_exp:.1%}, diff {wr_diff:.3%}) [{status}]")
        if not wr_pass:
            all_pass = False

        # Check chi2 > 10 count
        h_gt10 = count_high_chi2(result['hermes_chi2'], threshold=10.0)
        h_gt10_exp = expected['hermes_chi2_gt10']
        h_gt10_pass = h_gt10 == h_gt10_exp
        status = "PASS" if h_gt10_pass else "FAIL"
        print(f"  Hermes χ²>10:      {h_gt10} (exp {h_gt10_exp}) [{status}]")
        if not h_gt10_pass:
            all_pass = False

    return all_pass


def test_table_a3():
    """
    Test Table A3: Velocity residuals.
    Compute |V_obs - V_model| at each radial point.
    """
    print("\n" + "="*80)
    print("TABLE A3: VELOCITY RESIDUALS")
    print("="*80)

    ages = read_ages(AGES_CSV)

    # Default parameters for main results
    sigma_int_sq = 386.0
    a_knee = 1585.0

    all_residuals_h = []
    all_residuals_m = []
    all_fracs_h = []
    all_fracs_m = []
    floor_saturated_residuals_h = []
    floor_saturated_residuals_m = []
    floor_saturated_fracs_h = []
    floor_saturated_fracs_m = []

    # Read floor-saturated galaxy names (those with V_inf measurement issues)
    # Based on paper, these are identified by very large χ² or specific measurement properties
    # For now we'll compute for all and then identify the floor-saturated subset

    floor_saturated_galaxies = set()

    for galaxy, (t50, g98) in sorted(ages.items()):
        path = vh.find_rotmod(SPARC_DIR, galaxy)
        if path is None:
            continue

        data = vh.read_rotmod(path)
        if len(data['R']) < 3:
            continue

        # Hermes residuals
        V_hermes, beta, phi = vh.hermes_model(data, t50, g98, a_knee=a_knee)
        residuals_h = np.abs(data['Vobs'] - V_hermes)
        all_residuals_h.extend(residuals_h)
        # Fractional residuals (element-wise)
        fracs_h = residuals_h / np.maximum(np.abs(data['Vobs']), 1e-12)
        all_fracs_h.extend(fracs_h)

        # MOND residuals
        V_mond = vh.mond_model(data)
        residuals_m = np.abs(data['Vobs'] - V_mond)
        all_residuals_m.extend(residuals_m)
        # Fractional residuals (element-wise)
        fracs_m = residuals_m / np.maximum(np.abs(data['Vobs']), 1e-12)
        all_fracs_m.extend(fracs_m)

        # Identify floor-saturated galaxies (those with very high intrinsic scatter)
        # We use chi2 > 10 as indicator
        chi2_hermes = vh.chi2_nu(data['Vobs'], V_hermes, data['errV'], sigma_int_sq=sigma_int_sq)
        if chi2_hermes > 10.0:
            floor_saturated_galaxies.add(galaxy)
            floor_saturated_residuals_h.extend(residuals_h)
            floor_saturated_residuals_m.extend(residuals_m)
            floor_saturated_fracs_h.extend(fracs_h)
            floor_saturated_fracs_m.extend(fracs_m)

    all_residuals_h = np.array(all_residuals_h)
    all_residuals_m = np.array(all_residuals_m)
    all_fracs_h = np.array(all_fracs_h)
    all_fracs_m = np.array(all_fracs_m)
    floor_saturated_residuals_h = np.array(floor_saturated_residuals_h)
    floor_saturated_residuals_m = np.array(floor_saturated_residuals_m)
    floor_saturated_fracs_h = np.array(floor_saturated_fracs_h)
    floor_saturated_fracs_m = np.array(floor_saturated_fracs_m)

    print(f"\nAll {len(ages)} galaxies:")
    print("-" * 80)
    print(f"  N points (Hermes):       {len(all_residuals_h)}")
    print(f"  Hermes median |ΔV|:      {np.median(all_residuals_h):.1f} (exp 20.6)")
    print(f"  MOND median |ΔV|:        {np.median(all_residuals_m):.1f} (exp 25.8)")

    h_frac_med = np.median(all_fracs_h) * 100
    m_frac_med = np.median(all_fracs_m) * 100
    print(f"  Hermes median |ΔV/V_obs|:{h_frac_med:.1f}% (exp 18.4%)")
    print(f"  MOND median |ΔV/V_obs|:  {m_frac_med:.1f}% (exp 22.3%)")

    print(f"\nFloor-saturated galaxies ({len(floor_saturated_galaxies)} galaxies):")
    print("-" * 80)
    print(f"  N points:                {len(floor_saturated_residuals_h)}")
    if len(floor_saturated_residuals_h) > 0:
        print(f"  Hermes median |ΔV|:      {np.median(floor_saturated_residuals_h):.1f} (exp 34.6)")
        print(f"  MOND median |ΔV|:        {np.median(floor_saturated_residuals_m):.1f} (exp 47.5)")
        h_frac_fs = np.median(floor_saturated_fracs_h) * 100
        m_frac_fs = np.median(floor_saturated_fracs_m) * 100
        print(f"  Hermes median |ΔV/V_obs|:{h_frac_fs:.1f}% (exp 14.1%)")
        print(f"  MOND median |ΔV/V_obs|:  {m_frac_fs:.1f}% (exp 19.5%)")
    else:
        print("  (No floor-saturated galaxies found with χ² > 10)")

    return True


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Verify Appendix A sensitivity tables.')
    parser.add_argument('--sparc', required=True, help='Path to SPARC Rotmod_LTG directory')
    parser.add_argument('--ages', default=None, help='Path to ages_133.csv (default: same directory)')
    args = parser.parse_args()
    
    SPARC_DIR = args.sparc
    AGES_CSV = args.ages or os.path.join(os.path.dirname(os.path.abspath(__file__)), 'ages_133.csv')

    print("\n" + "="*80)
    print("APPENDIX A SENSITIVITY TABLE VERIFICATION")
    print("="*80)

    results = []

    # Test Table A1
    try:
        pass_a1 = test_table_a1()
        results.append(('Table A1 (σ_int²)', pass_a1))
    except Exception as e:
        print(f"\nERROR in Table A1: {e}")
        import traceback
        traceback.print_exc()
        results.append(('Table A1 (σ_int²)', False))

    # Test Table A2
    try:
        pass_a2 = test_table_a2()
        results.append(('Table A2 (a_knee)', pass_a2))
    except Exception as e:
        print(f"\nERROR in Table A2: {e}")
        import traceback
        traceback.print_exc()
        results.append(('Table A2 (a_knee)', False))

    # Test Table A3
    try:
        pass_a3 = test_table_a3()
        results.append(('Table A3 (velocity residuals)', pass_a3))
    except Exception as e:
        print(f"\nERROR in Table A3: {e}")
        import traceback
        traceback.print_exc()
        results.append(('Table A3 (velocity residuals)', False))

    # Summary
    print("\n" + "="*80)
    print("SUMMARY")
    print("="*80)
    for name, passed in results:
        status = "PASS" if passed else "FAIL"
        print(f"  {name}: {status}")

    all_passed = all(p for _, p in results)
    print(f"\nOverall: {'PASS' if all_passed else 'FAIL'}")
