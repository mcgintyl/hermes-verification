"""The galaxy-scale side of Hermes 1: SPARC rotation curves, the model, and MOND.

Units: radius in kpc, velocity in km/s, acceleration in (km/s)^2/kpc.

The rotmod reader and the baryonic acceleration follow the Paper 1 verification tool.
The gate is the same phi_chain used at cluster scale, so one function serves both.
"""
from __future__ import annotations

import csv
from pathlib import Path

import numpy as np

from . import boards
from .gate import phi_chain

G_KNEE = 1585.0                                  # (km/s)^2/kpc, the gate's acceleration knee
KDIV = 46654.0                                   # psi = t50 [Gyr] * g98 [(km/s)^2/kpc] / KDIV
F = boards.FLOOR                                 # 1 / sqrt(2 pi)
SIGMA_INT_SQ = 386.0                             # intrinsic scatter floor, (km/s)^2
# MOND a0 = 1.2e-10 m/s^2, expressed in (km/s)^2/kpc. Two conversions of the same value
# appear in this work and differ by 6.4e-6 relative:
#   A0        = 3702.789, from the rounded unit constant 1 (km/s)^2/kpc = 3.2408e-14 m/s^2.
#               This is what the paper's galaxy table (Table 1) was computed with.
#   A0_EXACT  = 3702.813, from the exact parsec conversion. This is what the cluster half
#               of this package uses, because the cluster table was computed with it.
# The difference moves the galaxy MOND median chi2/N by 1e-5 (1.142686 -> 1.142696) and
# nothing the paper quotes to three decimals. Each half uses the constant its own
# published table used, so both tables reproduce exactly.
A0 = 1.2e-10 / 3.2408e-14
A0_EXACT = 1.2e-10 * 3.0856775814913673e19 / 1.0e6
A0_ROUNDED_EXPORT = 3700.0                       # used by the archived Paper 1 per-galaxy export

DATA = boards.HERE / "data"


def read_rotmod(path):
    """Read one SPARC *_rotmod.dat file. Columns: R, Vobs, errV, Vgas, Vdisk, Vbul."""
    cols = [[], [], [], [], [], []]
    with open(path) as handle:
        for line in handle:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 6:
                continue
            for i in range(6):
                cols[i].append(float(parts[i]))
    keys = ("R", "Vobs", "errV", "Vgas", "Vdisk", "Vbul")
    return {k: np.array(v) for k, v in zip(keys, cols)}


def g_bar(d):
    """Baryonic acceleration: V_bar^2 = V_disk^2 + V_bul^2 + sign(V_gas) V_gas^2, g = V_bar^2 / R."""
    vbar_sq = d["Vdisk"] ** 2 + d["Vbul"] ** 2 + np.sign(d["Vgas"]) * d["Vgas"] ** 2
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.where(d["R"] > 0, vbar_sq / d["R"], 0.0)


def beta(t50_gyr, g98):
    """Hermes 1 wear score: beta = pi exp(-psi) - 1/sqrt(2 pi), psi = t50 g98 / KDIV."""
    return np.pi * np.exp(-t50_gyr * g98 / KDIV) - F


def hermes_velocity(d, t50_gyr, g98):
    """V = sqrt(g_bar [1 + beta phi] R), with phi the Paper 1 chain-rule gate."""
    gb = g_bar(d)
    g = gb * (1.0 + beta(t50_gyr, g98) * phi_chain(d["R"], gb))
    return np.sqrt(np.maximum(g, 0.0) * d["R"])


def mond_velocity(d, a0=A0):
    """MOND simple interpolation: nu(y) = (1 + sqrt(1 + 4/y)) / 2, y = g_bar / a0."""
    gb = g_bar(d)
    y = np.maximum(np.abs(gb) / a0, 1e-30)
    g = np.abs(gb) * (1.0 + np.sqrt(1.0 + 4.0 / y)) / 2.0
    g = np.where(gb >= 0, g, -g)
    return np.sqrt(np.maximum(g, 0.0) * d["R"])


def chi2_per_point(Vobs, Vmodel, errV, sigma_int_sq=SIGMA_INT_SQ):
    """chi^2/N = mean of (Vobs - Vmodel)^2 / (errV^2 + sigma_int^2). Divided by N, not by dof."""
    return float(np.mean((Vobs - Vmodel) ** 2 / (errV ** 2 + sigma_int_sq)))


def load_ages(path=None):
    """The 133 galaxies with their median stellar age t50 and 98th-percentile g_bar."""
    path = DATA / "ages_133.csv" if path is None else path
    with open(path, encoding="utf-8-sig", newline="") as handle:
        return [(r["galaxy"].strip(), float(r["t50_gyr"]), float(r["g98"]))
                for r in csv.DictReader(handle)]


def rotmod_path(sparc_dir, galaxy, manifest=None):
    """Locate a galaxy's rotmod file. The manifest fixes the file name for each galaxy."""
    sparc_dir = Path(sparc_dir)
    if manifest is None:
        manifest = load_manifest()
    if galaxy in manifest:
        candidate = sparc_dir / manifest[galaxy][0]
        if candidate.is_file():
            return candidate
    stem = galaxy.replace(" ", "")
    candidate = sparc_dir / (stem + "_rotmod.dat")
    if candidate.is_file():
        return candidate
    target = (stem + "_rotmod.dat").lower()
    for entry in sparc_dir.iterdir():
        if entry.name.lower() == target:
            return entry
    return None


def load_manifest(path=None):
    """galaxy -> (rotmod file name, sha256, bytes) for the 133 SPARC files this package uses."""
    path = DATA / "sparc_manifest.csv" if path is None else path
    with open(path, encoding="utf-8", newline="") as handle:
        return {r["galaxy"]: (r["rotmod_file"], r["sha256"], int(r["bytes"]))
                for r in csv.DictReader(handle)}


def trimmed_mean(x, fraction=0.05):
    """The paper's 5% trimmed mean: drop floor(0.05 n) values from each end."""
    s = np.sort(np.asarray(x, float))
    k = int(np.floor(fraction * len(s)))
    return float(np.mean(s[k:len(s) - k]))
