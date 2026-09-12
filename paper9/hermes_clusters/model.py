"""Hermes 1 and MOND responses, the GLS amplitude fit, and board assembly.

Units: radius in kpc, mass in solar masses, acceleration in (km/s)^2/kpc.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy import interpolate

from . import boards
from .gate import phi_chain

G = boards.G                                     # kpc (km/s)^2 / Msun
F = boards.FLOOR                                 # 1 / sqrt(2 pi)
KDIV = 46654.0                                   # psi = t50 [Gyr] * g98 [(km/s)^2/kpc] / KDIV (KDIV = c / 2 pi in these units)
A0 = 1.2e-10 * 3.0856775814913673e19 / 1.0e6     # MOND a0 = 1.2e-10 m/s^2, in (km/s)^2/kpc
T50_CLUSTER = 10.0                               # Gyr: an assumption, absorbed exactly by K

DATA = boards.HERE / "data"

# board id, Famaey name, Famaey profile stem, Mistele file stem, track, first clean raw node, table label
SEVEN = [
    ("A611", "A611", "a611", "abell611", "relaxed", 0, "A611"),
    ("MACSJ0429", "MACS0429", "macsj0429", "macsj0429", "relaxed", 0, "MACSJ0429"),
    ("MACSJ1720", "MACS1720", "macsj1720", "macsj1720", "relaxed", 0, "MACSJ1720"),
    ("RXJ2129-PARTIAL", "RXJ2129", "rxj2129", "rxj2129", "relaxed (partial board)", 3, "RXJ2129 (partial)"),
    ("A2261", "A2261", "a2261", "abell2261", "relaxed", 0, "A2261"),
    ("MACSJ1115", "MACS1115", "macsj1115", "macsj1115", "relaxed", 0, "MACSJ1115"),
    ("MACSJ1206", "MACS1206", "macsj1206", "macsj1206", "relaxed", 0, "MACSJ1206"),
]
EXCLUDED = [
    ("MACSJ0744", "MACS0744", "macsj0744", "macsj0744", "excluded: merger shock", 0, "MACSJ0744"),
    ("RXJ1347", "RXJ1347", "rxj1347", "rxj1347", "excluded: disturbed", 0, "RXJ1347"),
]


def load_famaey(root=None):
    """Famaey et al.'s ClusterInfo.py and member-galaxy template."""
    root = DATA / "famaey" if root is None else root
    cinfo = boards.import_path("famaey_cluster_info", root / "ClusterInfo.py")
    x200, gas_fraction = np.loadtxt(root / "fgasMACS.txt", unpack=True)
    finterp = interpolate.interp1d(x200, gas_fraction, kind="cubic")
    return cinfo, finterp, gas_fraction


def hermes_q(phi, g_nodes, t50=T50_CLUSTER):
    """Hermes 1 response q = phi (pi exp(-psi) - 1/sqrt(2 pi)), with g98 taken over the nodes."""
    return phi * (np.pi * np.exp(-t50 * np.percentile(g_nodes, 98) / KDIV) - F)


def mond_q(g_nodes):
    """MOND simple interpolation, nu(y) = (1 + sqrt(1 + 4/y)) / 2 with y = g_bar / a0; q = nu - 1."""
    return (1 + np.sqrt(1 + 4 / (g_nodes / A0))) / 2 - 1.0


def gls_fit(lens, mbar, cov, q):
    """One amplitude K >= 0 by generalised least squares for M_pred = M_bar (1 + K q). Returns (chi2, K)."""
    prec = np.linalg.solve(cov, np.eye(len(cov)))
    t = lens - mbar
    d = mbar * q
    k = max(0.0, float(d @ prec @ t / float(d @ prec @ d)))
    r = t - k * d
    return float(r @ prec @ r), k


@dataclass
class Board:
    spec: boards.BoardSpec
    label: str
    R: np.ndarray        # carrier grid radii (kpc): 50 points, 0.05 Mpc to 1.2 r200
    M: np.ndarray        # baryonic mass on the grid (Msun)
    meta: dict
    ledger: object       # node-admission ledger (pandas DataFrame)
    r: np.ndarray        # admitted node radii (kpc)
    lens: np.ndarray     # lensing M(<r) at the nodes (Msun)
    mbar: np.ndarray     # baryonic M(<r) at the nodes, log-log interpolated from the grid (Msun)
    cov: np.ndarray      # node covariance: correlation matrix x statistical errors
    phi: np.ndarray      # gate at the nodes, built on the grid and interpolated in ln r
    g: np.ndarray        # g_bar at the nodes


def assemble(spec, label, R, M, meta, mistele_root=None):
    """Admit the lensing nodes for a carrier (R, M) and put the gate on them."""
    mistele_root = DATA / "mistele" if mistele_root is None else mistele_root
    ledger, idx, r, lens, sigma, mbar, cov, diag = boards.admission_and_covariance(spec, R, M, mistele_root)
    phi = np.interp(np.log(r), np.log(R), phi_chain(R, G * M / R ** 2))
    g = G * mbar / r ** 2
    return Board(spec, label, R, M, meta, ledger, r, lens, mbar, cov, phi, g)


def build_board(row, famaey, famaey_root=None, mistele_root=None):
    famaey_root = DATA / "famaey" if famaey_root is None else famaey_root
    spec = boards.BoardSpec(*row[:6])
    cinfo, finterp, gas_fraction = famaey
    R, M, meta = boards.build_baryonic_carrier(spec, cinfo, finterp, gas_fraction, famaey_root)
    return assemble(spec, row[6], R, M, meta, mistele_root)


def score(b, t50=T50_CLUSTER):
    """Hermes 1 and MOND, each with one GLS amplitude, on one board."""
    ch, kh = gls_fit(b.lens, b.mbar, b.cov, hermes_q(b.phi, b.g, t50))
    cm, km = gls_fit(b.lens, b.mbar, b.cov, mond_q(b.g))
    return dict(H1_chi2=ch, H1_K=kh, MOND_chi2=cm, MOND_K=km)


def carrier_with_tail(spec, famaey, rmax_kpc=None, famaey_root=None):
    """The carrier rebuilt from its parts, optionally with Mistele et al.'s alternative gas.

    With rmax_kpc=None this reproduces build_baryonic_carrier. With a radius, the gas
    density is continued as a 1/r^4 tail from Rmax_X instead of the best-fit profile.
    Returns radius (kpc), total baryonic mass, gas mass, BCG + companions mass, f_gal.
    """
    famaey_root = DATA / "famaey" if famaey_root is None else famaey_root
    cinfo, finterp, gas_fraction = famaey
    v = np.loadtxt(famaey_root / "profiles" / ("results_profile_NOgnfw_all_%s.txt" % spec.profile_stem))
    r = v[:, 0] * 1000.0
    bcg = [row for row in cinfo.cluster_BCG_mass if str(row[0]) == spec.source_name][0]
    ci = int(np.flatnonzero(cinfo.nameclust == spec.profile_stem)[0])
    gi = int(np.flatnonzero(cinfo.cluster_data[:, 0] == spec.source_name)[0])
    gas = cinfo.cluster_data[gi]
    r200 = float(cinfo.r200[ci])
    mbcg = (float(bcg[1]) + float(bcg[3])) * 1e11
    p = [float(gas[i]) for i in range(2, 10)]
    if spec.source_name == "MACS1720":
        p[5] = p[6] = p[7] = 0.0
    x = (r / 1000.0) / r200
    fgal = np.full_like(x, 1.0 / np.max(gas_fraction))
    inner = x < 1.0
    fgal[inner] = 1.0 / finterp(x[inner])
    mg = np.asarray(cinfo.mgas(r, *p), float)
    if rmax_kpc is not None:
        Rx = float(rmax_kpc)
        mR = float(cinfo.mgas(Rx, *p))
        rhoR = float(cinfo.rho_fit(Rx, *p))
        out = r > Rx
        mg = mg.copy()
        mg[out] = mR + 4.0 * np.pi * rhoR * Rx ** 4 * (1.0 / Rx - 1.0 / r[out])
    return r, mbcg + mg * (1.0 + fgal), mg, mbcg, fgal
