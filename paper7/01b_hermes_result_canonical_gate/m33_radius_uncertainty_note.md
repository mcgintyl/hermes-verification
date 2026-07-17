# M33 galactocentric radius (`R_kpc`): derivation and uncertainty

Prepared in response to a reviewer request (V. Neves) for the equation behind
`R_kpc` and a radius-uncertainty column that can be propagated. Citable as-is.

## 1. Where `R_kpc` comes from

The board does **not** re-derive `R_kpc`. Every radius is imported as-published
from **Corbelli et al. (2014), A&A 572, A23, online Table 1**, which tabulates
the M33 rotation curve directly in kiloparsecs. Our reconstruction kernel reads
that column verbatim (`R = R_kpc`); it applies no distance or inclination
transformation of its own.

Corbelli et al. obtained those galactocentric radii from a **tilted-ring model**
of the warped HI disk. In that model a ring's physical galactocentric radius is

$$R_{\text{kpc}} = D \cdot \theta_{\text{ring}}$$

where $\theta_{\text{ring}}$ is the ring's angular galactocentric radius (radians)
in the **deprojected disk plane**, and $D$ is the adopted distance. The
deprojection from sky coordinates to in-plane radius uses each ring's inclination
$i(R)$ and position angle $\text{PA}(R)$. Because M33's outer disk is warped,
Corbelli et al. fit $i$ and PA as **functions of radius** (their Fig. 3), not as
global constants — the inner disk is near $i \approx 54^\circ$, $\text{PA}\approx
22^\circ$, and both drift outward. The inclination is therefore **already folded
into the published `R_kpc`**; we apply no further inclination correction.

## 2. Adopted distance and its uncertainty

$D = 840\ \text{kpc}$. Corbelli et al. (2014) adopt this citing Freedman et al.
(1991) and Gieren et al. (2013), and **state no explicit distance uncertainty**.
We therefore take the uncertainty from the modern Cepheid source they cite,
**Gieren et al. (2013)** (Araucaria Project; ApJ 773, 69; arXiv:1305.4258), which
gives a true distance modulus

$$\mu = 24.62 \pm 0.07\ \text{mag} \quad (\pm 0.03\ \text{stat}, \pm 0.06\ \text{sys})$$

Converting a modulus error to a fractional distance error,
$\delta D / D = (\ln 10 / 5)\,\delta\mu = 0.4605 \times 0.07$:

$$\frac{\delta D}{D} = 0.0322\ (3.2\%), \qquad \delta D \approx 27\ \text{kpc}.$$

## 3. The `errR_kpc` column

$$\text{errR\_kpc} = \frac{\delta D}{D}\,R_{\text{kpc}} = 0.0322 \cdot R_{\text{kpc}}$$

This is a **systematic scale uncertainty, not an independent per-point
measurement error.** For propagation:

- **It is fully correlated across all rings.** If the true distance differs from
  840 kpc, every `R_kpc` scales by the same factor. **Do not combine `errR_kpc`
  values in quadrature as if independent** — a distance shift moves the whole
  curve coherently.
- **It propagates into every distance-scaled quantity.** With $R \propto D$:
  $g_{\text{bar}} = V_{\text{bar}}^2/R \propto D^{-1}$; the age term's
  $g_{98} \propto D^{-1}$; surface densities $\propto D^{-2}$. If you re-scale
  $R$, re-scale those consistently with the same $\delta D/D$, or the board is
  internally inconsistent.
- Velocities (`Vobs`, `Vbar`) are **not** distance-scaled (they come from the
  observed line-of-sight field and the tilted-ring $V/\sin i$ deprojection), so
  `errR_kpc` does not touch them.

## 4. Why there are no independent per-point radius errors

Corbelli's rings sit at **chosen** angular radii; a ring's angular radius has no
meaningful per-ring measurement scatter. The genuine radius uncertainties are:

1. the **global distance scale** — dominant, coherent, captured by `errR_kpc`; and
2. **tilted-ring inclination/PA** uncertainty — radius-dependent and largest in
   the warped outer disk, but **not published by Corbelli as per-ring radius
   errors**.

There is no source for an independent per-point $\sigma_R$ for this dataset.
Manufacturing one would be fabrication, so we do not.

## 5. Flags (adopted, not measured here — verify if you need tighter values)

- **Corbelli (2014) state no distance error.** The 3.2% is *adopted* from Gieren
  et al. (2013), the Cepheid source Corbelli cite. If you prefer a different M33
  distance/uncertainty, rescale `errR_kpc` by your own $\delta D/D$.
- **Inclination/PA vary with radius** (warp); `errR_kpc` captures only the
  distance-scale term, which dominates. A full outer-disk radius error budget
  would additionally require Corbelli's per-ring $i$/PA uncertainties, which are
  not tabulated.
- **`R_kpc` was taken as-published** from Corbelli Table 1, not recomputed. The
  equation in §1 is Corbelli's derivation, given so you can propagate; it is not
  re-executed in this repository.

## References

- Corbelli, E., Thilker, D., Zibetti, S., Giovanardi, C., & Salucci, P. 2014,
  A&A, 572, A23.
- Gieren, W., et al. 2013, ApJ, 773, 69 (arXiv:1305.4258).
- Freedman, W. L., Wilson, C. D., & Madore, B. F. 1991, ApJ, 372, 455.
