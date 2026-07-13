#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
verify_gates.py -- one frozen gate reproduces both published datasets.

The gate function phi(R) is defined once, in paper1/hermes_gate_phi.py, and is
used for every result in the project. This script imports that gate and checks it
against both published datasets to a tolerance of 1e-10:

  CLAIM 1  Paper 1 -- 133 SPARC galaxies (phi_last)
  CLAIM 2  Paper 7 -- M33 canonical-gate result (per-point phi)

(The first M33 run used a separate script with three incidental gate deviations;
it has been superseded -- see paper7/01b_hermes_result_canonical_gate/. This test
uses the canonical, corrected M33 result.)

Data:
  * M33 reference (in repo)      : paper7/01b_hermes_result_canonical_gate/
                                   m33_hermes_per_point_predictions_canonical.csv
  * 133 SPARC per-galaxy export  : paper1/Hermes_ConfigG_PerGalaxy_133_Export.csv
  * SPARC rotmod files           : *_rotmod.dat  (the only external dependency;
                                   set SPARC_DIR or drop them in ./sparc_database)
"""
import os, sys, csv, glob, re
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "paper1"))
from hermes_gate_phi import hermes_phi          # the single frozen gate

TOL = 1e-10
M33_REF = os.path.join(HERE, "paper7", "01b_hermes_result_canonical_gate",
                       "m33_hermes_per_point_predictions_canonical.csv")
EXPORT = os.environ.get("HERMES_EXPORT",
    os.path.join(HERE, "paper1", "Hermes_ConfigG_PerGalaxy_133_Export.csv"))
# SPARC rotmod files are the ONLY external dependency. Download the SPARC database
# (Lelli, McGaugh & Schombert 2016; http://astroweb.cwru.edu/SPARC/) and either drop
# the *_rotmod.dat files in ./sparc_database or point SPARC_DIR at them.
SPARC_DIR = os.environ.get("SPARC_DIR", os.path.join(HERE, "sparc_database"))


def _norm(s):
    s = re.sub(r'[^a-z0-9]', '', s.lower())
    m = re.match(r'^([a-z]+)0*([0-9].*)$', s)
    return m.group(1) + m.group(2) if m else s


def load_m33():
    rows = list(csv.DictReader(open(M33_REF)))
    R = np.array([float(r["R_kpc"]) for r in rows])
    gbar = np.array([float(r["gbar_kms2_per_kpc"]) for r in rows])
    phi = np.array([float(r["phi_canonical"]) for r in rows])
    return R, gbar, phi


def load_sparc():
    rot = {_norm(os.path.basename(f).replace("_rotmod.dat", "")): f
           for f in glob.glob(os.path.join(SPARC_DIR, "*_rotmod.dat"))}
    if not rot:
        sys.stderr.write(
            "\nERROR: no SPARC *_rotmod.dat files found in:\n  %s\n\n"
            "The SPARC rotation-curve database is the only external dependency.\n"
            "Download it from http://astroweb.cwru.edu/SPARC/ (Lelli, McGaugh &\n"
            "Schombert 2016), then either place the *_rotmod.dat files in\n"
            "./sparc_database or set SPARC_DIR to point at them.\n" % os.path.abspath(SPARC_DIR))
        sys.exit(2)
    out = []
    for r in csv.DictReader(open(EXPORT)):
        nm = _norm(r["galaxy"])
        if nm not in rot:
            continue
        Rk, Vg, Vd, Vb = [], [], [], []
        for ln in open(rot[nm]):
            if ln.startswith("#") or not ln.strip():
                continue
            p = ln.split()
            Rk.append(float(p[0])); Vg.append(float(p[3])); Vd.append(float(p[4])); Vb.append(float(p[5]))
        Rk, Vg, Vd, Vb = map(np.array, (Rk, Vg, Vd, Vb))
        gbar = (Vd**2 + Vb**2 + np.sign(Vg) * Vg**2) / Rk   # Config G: M/L = 1
        out.append((r["galaxy"], Rk, gbar, float(r["phi_last"])))
    return out


def main():
    print("=" * 66)
    print("HERMES GATE VERIFICATION   (one frozen gate, tolerance = 1e-10)")
    print("=" * 66)

    m33_R, m33_g, m33_phi = load_m33()
    sparc = load_sparc()

    d1 = max(abs(hermes_phi(R, g)[-1] - pl) for _, R, g, pl in sparc)
    d2 = float(np.max(np.abs(hermes_phi(m33_R, m33_g) - m33_phi)))

    print("\nCLAIM 1  paper1/hermes_gate_phi.py  vs  Paper 1 / 133 SPARC (phi_last)")
    print("  galaxies = %d    max|delta| = %.2e    -> %s" % (len(sparc), d1, "PASS" if d1 <= TOL else "FAIL"))
    print("CLAIM 2  paper1/hermes_gate_phi.py  vs  Paper 7 / M33 (per-point phi)")
    print("  points   = %d     max|delta| = %.2e    -> %s" % (len(m33_R), d2, "PASS" if d2 <= TOL else "FAIL"))

    ok = (d1 <= TOL) and (d2 <= TOL)
    print("\n" + "-" * 66)
    print("RESULT: %s  (one frozen gate reproduces both published datasets)" % ("PASS" if ok else "FAIL"))
    print("-" * 66)
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
