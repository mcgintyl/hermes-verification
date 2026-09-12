"""The Hermes 1 radial gate phi(R): the Paper 1 chain-rule construction, unmodified.

phi_chain(R, g_bar) takes radii in kpc and the baryonic acceleration in (km/s)^2/kpc on
a radial grid (for clusters, the 50-point baryonic carrier grid) and returns phi in [0, 1]:

  1. the knee radius r_k is where g_bar first drops below a_knee = 1585 (km/s)^2/kpc
     (linear interpolation between grid points); x = R / r_k
  2. V = sqrt(g_bar R), smoothed with a Savitzky-Golay filter (window 11 when there are
     at least 20 points, cubic, mirror mode)
  3. shear s_i = |(R_i / V_i)(dV/dR)_i| and curvature
     kappa_i = |(d2V/dR2)_i| / (|(dV/dR)_i| + eps), eps = 1e-4 x median|dV/dR|,
     with both derivatives taken from the smoothed V of step 2 on the native
     radial grid (numpy.gradient, first order at the two end points). Neither is
     computed as a finite difference of ln V against ln R on a log-resampled
     grid: that is a different gate and does not reproduce Paper 1's published
     scores.
  4. kappa is normalised by its median inside the knee (x <= 1), s by its median
     outside (x > 1)
  5. the two are blended with a logistic weight in x (centre 1.40, width 0.30)
  6. an exponential moving average runs outward (weight 0.5), and phi = 1 - exp(-S),
     clipped to [0, 1]

The same function, with the same constants, produces the galaxy-scale results of Paper 1.
"""
import numpy as np
from scipy.signal import savgol_filter

A_KNEE = 1585.0


def phi_chain(R_kpc, g_bar, a_knee=A_KNEE):
    R = np.asarray(R_kpc, float).copy()
    g = np.asarray(g_bar, float).copy()
    o = np.argsort(R)
    R, g = R[o], g[o]
    u = np.concatenate(([True], np.diff(R) > 0))
    R, g = R[u], g[u]
    if len(R) < 3:
        return np.zeros_like(R)
    g = np.where(np.isfinite(g) & (g >= 0), g, 0.0)
    idx = None
    for i in range(len(R) - 1):
        if g[i] >= a_knee > g[i + 1]:
            idx = i
            break
    if idx is not None:
        g1, g2 = g[idx], g[idx + 1]
        t = 0.0 if g2 == g1 else (a_knee - g1) / (g2 - g1)
        rk = R[idx] + t * (R[idx + 1] - R[idx])
    else:
        rk = R[np.argmin(np.abs(g - a_knee))]
    if not np.isfinite(rk) or rk <= 0:
        rk = np.median(R)
    x = R / rk
    V = np.sqrt(np.maximum(g * R, 0))
    N = len(R)
    if N < 5:
        Vs = V
    else:
        w = 5 if N < 20 else 11
        if w > N:
            w = N if N % 2 == 1 else N - 1
        Vs = savgol_filter(V, window_length=w, polyorder=3, mode="mirror")
    Vs = np.where(Vs > 0, Vs, 1e-12)
    dV = np.gradient(Vs, R, edge_order=1)
    s = np.abs(R / Vs * dV)
    d2 = np.gradient(dV, R, edge_order=1)
    ek = 1e-4 * np.median(np.abs(dV))
    if not np.isfinite(ek) or ek <= 0:
        ek = 1e-12
    kap = np.abs(d2) / (np.abs(dV) + ek)
    eps = 1e-12
    inner, outer = x <= 1, x > 1
    mk = np.median(kap[inner]) if inner.any() else np.median(kap)
    ms = np.median(s[outer]) if outer.any() else np.median(s)
    wmix = 1 / (1 + np.exp(-(x - 1.40) / 0.30))
    Sp = (1 - wmix) * (kap / (mk + eps)) + wmix * (s / (ms + eps))
    Se = np.empty_like(Sp)
    Se[0] = Sp[0]
    for i in range(1, N):
        Se[i] = 0.5 * Sp[i] + 0.5 * Se[i - 1]
    return np.clip(1 - np.exp(-np.maximum(Se, 0)), 0, 1)
