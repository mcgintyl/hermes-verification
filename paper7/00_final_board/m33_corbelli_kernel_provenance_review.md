
# Corbelli/M33 Kernel Provenance Review

**Lane:** Corbelli/M33  
**Galaxy:** M33 / NGC 598 / UGC 1117  
**Status locked:** **Format B reconstructed engineering board available; closure checks passed; not yet authorized for Hermes/MOND execution.**  
**Current action:** Kernel/provenance review only.  
**Explicit stop conditions:** No Hermes run. No MOND run. No age anchoring. No pilot-board execution.

---

## 1. Executive decision

The current **uniform-thickness annular quadrature kernel** is acceptable as the **engineering reconstruction kernel** for the M33 Format B candidate board.

It should **not yet be treated as the final frozen physics kernel** until one additional provenance gate is passed:

> Benchmark the annular quadrature output against a literal Corbelli/Casertano-style implementation or against the source paper’s published component curves/figures.

A literal ROTMOD implementation is not the load-bearing requirement here, because Corbelli et al. cite **Casertano (1983)** for the thick-disk force computation, not GIPSY/ROTMOD. The correct final validation target is therefore **Casertano-equivalent thick-disk force**, not ROTMOD by name.

Recommended current label:

> **M33-Corbelli source-native thick-disk reconstruction — engineering kernel accepted; final kernel pending Casertano-equivalence benchmark.**

---

## 2. What the source paper requires

Corbelli et al. (2014) supply the ingredients needed for a source-native reconstruction:

- observed rotation curve,
- velocity uncertainties,
- azimuthally averaged H I surface density,
- modeled stellar surface density,
- molecular gas prescription,
- helium correction,
- no-bulge context,
- source-specific disk-thickness assumptions.

The source baryonic treatment is not SPARC-style unit-M/L photometry. It is an M33-specific stellar mass surface-density construction based on pixel-by-pixel SED / stellar-population modeling. The stellar mass profile therefore carries a radially varying SPS/pixel-SED M/L provenance.

For the force calculation, Corbelli et al. state that the gaseous disk is vertically uniform with half-thickness 0.5 kpc, the stellar disk flares from 100 pc at the center to 1 kpc at the outer disk edge, and disk mass contributions are computed following Casertano (1983), which generalizes the radial-force formula to thick disks.

This matters: the kernel is no longer a background implementation detail. It is the object that converts the public mass surface-density profiles into the board’s \(V_{\rm gas}(R)\) and \(V_{\rm disk}(R)\).

---

## 3. Current kernel

The current engineering kernel uses direct Newtonian annular quadrature.

### Geometry

| Component | Current engineering choice |
|---|---|
| Gas | axisymmetric annuli; vertically uniform density; half-thickness 0.5 kpc |
| Stars | axisymmetric annuli; flaring half-thickness from 0.1 kpc to 1.0 kpc over the disk |
| Bulge | \(V_{\rm bulge}(R)=0\) |
| Vertical profile | uniform slab approximation |
| Force computation | numerical annular quadrature over radius and azimuth |
| Role | deterministic reconstruction kernel, not a fit parameter |

### Numerical settings used in the engineering pass

| Setting | Value |
|---|---:|
| Radial integration step | 0.01 kpc |
| Angular samples | 1440 per annulus |
| Gas half-thickness | 0.5 kpc |
| Stellar half-thickness | linear flare from 0.1 kpc to 1.0 kpc |
| Stellar mass input | source-native SPS / pixel-SED \(\Sigma_\star(R)\) |
| Gas input | \(1.33[\Sigma_{\rm HI}(R)+\Sigma_{\rm H_2}(R)]\) |
| Molecular gas | \(\Sigma_{\rm H_2}(R)=10\exp(-R/2.2)\ M_\odot\,{\rm pc}^{-2}\) |
| \(V_{\rm bulge}\) | 0 |

The kernel is deterministic and physically transparent. It does not tune M/L, does not tune gas mass, and does not fit the rotation curve.

---

## 4. Comparison to Corbelli/Casertano-style treatment

### Agreement

The engineering kernel agrees with the Corbelli source treatment on the load-bearing geometry:

| Item | Corbelli source treatment | Current engineering kernel | Match? |
|---|---|---|---|
| Gas half-thickness | 0.5 kpc | 0.5 kpc | Yes |
| Gas vertical form | vertically uniform | uniform slab | Yes |
| Stellar disk | flaring disk | flaring disk | Yes |
| Stellar half-thickness | 100 pc center to 1 kpc edge | 0.1 to 1.0 kpc | Yes |
| Component source | mass surface-density profiles | same public profiles | Yes |
| Bulge | no bulge / no prominent bar | \(V_{\rm bulge}=0\) | Yes |
| M/L handling | SPS/pixel-SED mass profile | same source-native profile | Yes |
| Purpose | compute disk mass contributions | compute disk mass contributions | Yes |

### Differences

The current kernel is not guaranteed to be bit-equivalent to Corbelli/Casertano.

| Difference | Why it matters |
|---|---|
| Corbelli cites Casertano’s thick-disk radial-force formalism; current kernel uses direct annular quadrature. | These should converge to the same physical force if the density geometry is matched, but numerical details can differ. |
| The current stellar vertical profile is treated as a uniform slab unless further specified. | Corbelli’s line explicitly says the gas disk is vertically uniform; the stellar vertical profile is described as flaring, but not fully specified in the abbreviated public text. |
| Current kernel treats the reconstructed disk as axisymmetric annuli. | M33’s gas disk warps beyond 8 kpc; the board kernel uses radial averages, which is appropriate for a 1D board but not a full 3D warped-potential calculation. |
| The source paper also applies finite-thickness / beam-smearing corrections to the observed H I rotation curve. | Those corrections belong to \(V_{\rm obs}\), not to the baryonic component kernel. The board should not re-apply them. |

The conclusion is that the current kernel is **source-consistent** but not yet **source-identical**.

---

## 5. Comparison to SPARC-style component construction

SPARC-style boards start from already-computed component velocity arrays:

\[
R,\quad V_{\rm obs},\quad errV,\quad V_{\rm gas},\quad V_{\rm disk},\quad V_{\rm bulge}.
\]

In SPARC, \(V_{\rm gas}\) includes the helium correction, while \(V_{\rm disk}\) and \(V_{\rm bulge}\) are tabulated at a 3.6 μm unit-M/L component convention. The Hermes baseline uses those SPARC-provided mass-model velocities with no additional per-galaxy M/L fitting or rescaling.

M33/Corbelli is different:

| Category | SPARC rotmod baseline | M33/Corbelli Format B board |
|---|---|---|
| Component availability | Directly tabulated \(V_{\rm gas}\), \(V_{\rm disk}\), \(V_{\rm bulge}\) | Reconstructed from public surface-density profiles |
| Stellar convention | 3.6 μm component convention, \(M/L=1\) in rotmod table | SPS/pixel-SED \(\Sigma_\star(R)\), radially varying M/L provenance |
| Gas convention | \(V_{\rm gas}\) already includes helium factor | helium factor 1.33 applied during reconstruction |
| Bulge | tabulated if present | \(V_{\rm bulge}=0\) by source morphology |
| Kernel visibility | hidden upstream inside SPARC component construction | explicit and load-bearing in this board |
| Label | SPARC-equivalent rotmod board | source-native thick-disk reconstruction |

Therefore the board must **not** be described as SPARC-equivalent. It may become strict-schema-compatible after reconstruction, but it remains **source-native Corbelli/Casertano provenance**, not SPARC rotmod provenance.

---

## 6. Is a literal ROTMOD/Casertano-equivalent implementation needed?

### ROTMOD

A literal ROTMOD implementation is **not required** for Corbelli/M33 because ROTMOD is not the stated source kernel in Corbelli et al. (2014). ROTMOD was load-bearing in the Namumba papers, not in this Corbelli lane.

### Casertano

A **Casertano-equivalent validation** is required before final model execution.

This does not necessarily mean replacing the current annular quadrature kernel. It means one of the following must happen:

1. implement a literal Casertano-style thick-disk force calculation and compare it to the annular quadrature result;
2. obtain source component arrays or reconstruct source figures closely enough to verify that the current kernel reproduces the published stellar and gas component curves;
3. demonstrate numerical convergence and source-figure agreement to a predeclared tolerance.

If the current kernel passes that benchmark, it can be frozen as the lane kernel. If it fails materially, the kernel must be replaced before any physics comparison.

---

## 7. Minimum final-freeze checklist

Before authorizing Hermes/MOND, the M33 board should pass these checks:

| Gate | Required status |
|---|---|
| Primary table import | Corbelli online Table 1 imported from primary/official source or verified against it |
| Mass closure | already passed at first order: stellar 0.998, H2 1.013; H I excess accepted under source warning |
| Kernel convergence | rerun with finer radial and angular sampling; component velocities stable |
| Casertano equivalence | compare against literal Casertano implementation or source component figures |
| Provenance label | fixed as source-native thick-disk, not SPARC-equivalent |
| M/L convention | fixed as SPS/pixel-SED stellar mass profile, not tunable |
| Gas convention | helium factor 1.33 and H2 prescription preserved |
| Bulge convention | \(V_{\rm bulge}=0\) preserved |
| Model execution | not authorized until PI approval |

Suggested quantitative benchmark for the kernel review:

- component curves stable under doubled angular sampling and halved radial integration step;
- total component masses unchanged within numerical precision;
- \(V_{\rm gas}(R)\), \(V_{\rm disk}(R)\), and \(V_{\rm bar}(R)\) visually and numerically consistent with Corbelli/López Fune component figures wherever those figures show component curves;
- any deviations documented as kernel provenance, not tuned away.

---

## 8. Final board status after kernel review

Current status remains:

> **Corbelli/M33 — Format B reconstructed engineering board available; closure checks passed; not yet authorized for Hermes/MOND execution.**

Kernel decision:

> **Retain the uniform-thickness annular quadrature kernel for engineering reconstruction. Do not freeze it as final until a Casertano-equivalence benchmark is completed.**

Final label if retained:

> **M33-Corbelli source-native thick-disk reconstruction.**

Final warning:

> This lane is a real Format B candidate board, but its baryonic component velocities are reconstructed products. The kernel is therefore part of the data provenance. It must be frozen before physics, not adjusted after physics.

---

## 9. Sources consulted

- Corbelli et al. 2014, *Dynamical signatures of a ΛCDM-halo and the distribution of the baryons in M33*.
- Casertano 1983, thick-disk radial-force formalism referenced by Corbelli et al.
- López Fune, Salucci & Corbelli 2017, M33 radial dark-matter distribution paper, which restates the Corbelli stellar/gas mass-distribution architecture.
- SPARC VizieR ReadMe, J/AJ/152/157.
- Hermes Equation manuscript / Outside SPARC schema materials.
