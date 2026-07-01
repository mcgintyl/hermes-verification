# MaNGA Two-Layer Note: Simulation and Assembly-Bias Pass

## Purpose

This pass checks the strongest standard-model counterargument to the controlled MaNGA dynamical residual: ordinary \(\Lambda\)CDM assembly history may already produce an age-organized dynamical residual at fixed stellar mass.

The target question was narrow:

> Has anyone already run a mock MaNGA / DynPop-like analysis through EAGLE, Illustris, IllustrisTNG, or TNG50, constructed a JAM/SPS DML-style quantity, split galaxies into young and old quartiles at fixed stellar mass, and reported the controlled dynamical age residual?

## Result

No exact published comparator was located in this pass.

That does **not** make the residual anomalous. It means the residual is currently unbenchmarked against the proper standard-model simulation test.

## What was found

1. **Assembly bias is a live conventional pathway.**
   Galaxy age is not an isolated stellar-population label in \(\Lambda\)CDM. Halo mass-accretion history, gas loss, quenching, morphology, environment, and star-formation history can remain correlated at fixed mass.

2. **MaNGA-like simulation tools exist.**
   iMaNGA and MaNGIA provide TNG50 / IllustrisTNG-based mock MaNGA products, including MaNGA-like IFU observations and recovered stellar-population / kinematic quantities. These are the right family of tools for a future comparator.

3. **Adjacent MaNGA-simulation comparisons exist.**
   Earlier MaNGA density-slope work compares MaNGA JAM-derived slopes against EAGLE, Illustris, and IllustrisTNG. DynPop VI compares DynPop density slopes against TNG50/TNG100 and reports nontrivial mismatches in simulated mass distributions and dark-matter fractions. SEAGLE compares EAGLE/Illustris/TNG against strong-lensing dark-matter fractions and shows that simulation predictions differ substantially depending on feedback treatment.

4. **No exact DML residual pipeline was found.**
   This pass did not locate a paper that: (a) forward-models simulated galaxies as MaNGA/DynPop-like observations, (b) runs comparable JAM/SPS products, (c) constructs the same DML quantity, (d) uses the same fixed-stellar-mass young/old quartile split, and (e) reports the controlled dynamical \(M/L\) age residual.

## Interpretation for the paper

The paper should not claim:

> \(\Lambda\)CDM cannot explain the controlled residual.

It should also not claim:

> \(\Lambda\)CDM already explains the controlled residual.

The defensible wording is:

> A proper cosmological mock comparison has not been performed here. Adjacent simulation and mock-MaNGA tools exist, and assembly bias remains a plausible standard-model route for the controlled residual. We therefore do not claim the residual is anomalous relative to \(\Lambda\)CDM. The next specialist test is to process iMaNGA/TNG50, EAGLE, or comparable hydrodynamic mock observations through the same two-layer pipeline.

## Draft integration

Integrated into v0.5:

- Abstract: added simulation/assembly-bias limitation paragraph.
- Section 3.4: added simulation and assembly-bias context.
- Section 10.8: expanded cosmological-mock limitation.
- Section 16.2: changed from pending task to completed literature-pass status and future specialist comparator.
- References: added MaNGA/TNG/iMaNGA/DynPop VI/assembly-bias sources.
