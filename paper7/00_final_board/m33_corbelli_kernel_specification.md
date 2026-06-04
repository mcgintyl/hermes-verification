# M33-Corbelli Kernel Specification

**Lane label:** M33-Corbelli source-native thick-disk reconstruction  
**Kernel status:** provisionally frozen for this lane after PI/Astra approval and published-figure overlay benchmark.  
**Execution status:** no Hermes run, no MOND run, no gravity comparison.

## Frozen reconstruction kernel

The board uses a deterministic annular quadrature thick-disk force kernel. The kernel is part of the data provenance, not a fit lever.

| Item | Frozen specification |
|---|---|
| Geometry | axisymmetric annular quadrature in disk plane |
| Radial integration step | 0.01 kpc |
| Azimuthal samples | 1440 per annulus |
| Gas vertical geometry | vertically uniform slab |
| Gas half-thickness | 0.5 kpc |
| Stellar vertical geometry | flaring stellar disk |
| Stellar half-thickness | linear flare from 0.1 kpc at center to 1.0 kpc at 23 kpc |
| Gas surface density | Sigma_gas = 1.33 * (Sigma_HI + Sigma_H2) |
| Molecular gas | Sigma_H2 = 10 exp(-R/2.2) Msun pc^-2 |
| Stellar surface density | Corbelli BVIgi / pixel-SED Sigma_star(R) |
| Bulge | V_bulge(R) = 0 |
| Kernel role | source-native reconstruction; not SPARC-equivalent |

## Provenance decision

The target is not ROTMOD by name. The target is a Casertano-equivalent thick-disk force treatment matching Corbelli's stated geometry. The current annular quadrature kernel passed the published Figure 12 top-panel overlay sanity check and is provisionally frozen for the M33-Corbelli lane.

## Non-equivalence statement

This is not a SPARC-equivalent unit-M/L photometric reconstruction. The stellar profile is source-native SPS / pixel-SED stellar mass with radially varying M/L provenance. The M/L convention is part of the data provenance and is not tunable inside Hermes.
