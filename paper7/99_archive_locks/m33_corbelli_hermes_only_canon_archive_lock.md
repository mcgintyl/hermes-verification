# Canon Archive Lock — M33 First External Hermes Execution

**Lane:** M33-Corbelli source-native thick-disk reconstruction  
**Galaxy:** M33 / NGC 598  
**Status:** First external Hermes execution complete; Hermes-only package preserved exactly as-is; MOND not run.

## Locked result status

**M33-Corbelli source-native thick-disk reconstruction — first external Hermes execution complete; both bracketed Amber-age lanes pass/non-catastrophic; result elevated; dominant failure mode is coherent outer-disk underprediction.**

## Accepted result summary

| Age lane | t50 | Reduced chi-squared | Classification |
|---|---:|---:|---|
| Primary/source-native mass-weighted | 6.0 ± 1.2 Gyr | 3.126 | PASS / non-catastrophic / elevated |
| Conservative older-weighted bracket | 7.1 ± 1.3 Gyr | 3.296 | PASS / non-catastrophic / elevated |

Both age lanes remain below the current catastrophic threshold of chi-squared_nu = 5. The younger/source-native age performs slightly better. Age sensitivity is modest and does not dominate the result.

## Radial failure structure

The dominant residual structure is a coherent outer-disk underprediction:

| Diagnostic | Accepted value |
|---|---:|
| Outer-disk threshold | R > 10 kpc |
| Fraction of total chi-squared beyond 10 kpc | ~86.3% in both age lanes |
| Outer mean residual | approximately -52 to -53 km/s |
| Interpretation | systematic underprediction in the outer tail |

## Canon interpretation

Recommended paper language:

> Hermes ports to a non-SPARC, source-native reconstructed M33 board at non-catastrophic accuracy under both bracketed age anchors, but exhibits a coherent outer-disk underprediction that must be treated as a real residual structure rather than hidden by the global pass classification.

Do not overstate. Do not call this broad external validation. Do not collapse the two age lanes into one. Do not reinterpret or alter the board after the result.

## Provenance constraints preserved

The execution used the frozen M33-Corbelli source-native board:

- source-native SPS / pixel-SED stellar mass,
- radially varying M/L provenance,
- helium factor 1.33,
- H2 prescription,
- flaring stellar disk,
- gas half-thickness 0.5 kpc,
- Vbul = 0,
- annular thick-disk kernel provisionally frozen after overlay benchmark,
- bracketed Amber age lock carried explicitly.

No M/L, age, gas scaling, kernel parameter, distance, inclination, or uncertainty tuning was performed after execution.

## Archive preservation

The Hermes-only package was copied exactly into the Canon archive location without modification.

| File | SHA256 |
|---|---|
| `m33_first_external_hermes_execution_v2_package.zip` | `81720668c8d76538047b67c7b5f2415f87f9ec869a654830af4c0874fe099fa2` |
| `m33_corbelli_first_external_hermes_execution_HERMES_ONLY_CANON_LOCKED.zip` | `81720668c8d76538047b67c7b5f2415f87f9ec869a654830af4c0874fe099fa2` |

The matching hash confirms the Canon archive copy is byte-identical to the Hermes-only execution package.

## Next gate

MOND comparator authorization may be considered next, using the exact same frozen M33 board and no post-result reconstruction changes.

MOND has not been run for this lane. No MOND result exists in the Canon archive.
