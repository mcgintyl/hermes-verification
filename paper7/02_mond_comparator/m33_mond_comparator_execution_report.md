# M33 MOND Comparator Execution Report

**Lane label:** M33-Corbelli source-native thick-disk reconstruction  
**Execution:** MOND comparator only  
**Hermes package status:** unchanged; canonical Hermes-only package preserved  
**Board edits after Hermes:** none  

## Run configuration

- Board: exact frozen M33-Corbelli reconstructed board used for the Hermes execution.
- Verified input equality against Hermes per-point input columns: max absolute difference is 0.0 for R, Vobs, errV, Vgas, Vdisk, Vbulge, and Vbar.
- MOND interpolation: simple MOND, $\mu(x)=x/(1+x)$, equivalent $\nu(y)=0.5+\sqrt{0.25+1/y}$.
- $a_0$: fixed at $1.2\times10^{-10}\,\mathrm{m\,s^{-2}}$ = 3702.813098 $(\mathrm{km/s})^2/\mathrm{kpc}$.
- $a_0$ fitted: no.
- M/L adjustment: no.
- Gas scaling / helium / H2 / kernel / distance / inclination / uncertainties: unchanged.
- Scoring denominator: $errV^2 + 386$.

## Main result

| Model/lane | N | chi2 | reduced chi2 | median abs dV | mean residual | outer chi2 fraction | classification |
|---|---:|---:|---:|---:|---:|---:|---|
| Hermes, t50=6.0 Gyr | 58 | 181.298 | 3.126 | 29.90 km/s | -27.83 km/s | 86.3% | PASS / non-catastrophic |
| Hermes, t50=7.1 Gyr | 58 | 191.147 | 3.296 | 31.39 km/s | -29.31 km/s | 86.3% | PASS / non-catastrophic |
| MOND fixed simple | 58 | 6.889 | 0.119 | 3.69 km/s | -4.12 km/s | 89.3% | PASS / non-catastrophic |

## Radial structure

| Model/lane | zone | N | reduced chi2 | chi2 fraction | mean residual | median abs dV |
|---|---|---:|---:|---:|---:|---:|
| Hermes primary t50=6.0 | inner_0_3_kpc | 15 | 0.262 | 2.2% | 6.14 km/s | 8.55 km/s |
| Hermes primary t50=6.0 | mid_3_10_kpc | 16 | 1.308 | 11.5% | -19.51 km/s | 24.80 km/s |
| Hermes primary t50=6.0 | outer_gt10_kpc | 27 | 5.794 | 86.3% | -51.62 km/s | 54.83 km/s |
| Hermes conservative t50=7.1 | inner_0_3_kpc | 15 | 0.212 | 1.7% | 4.56 km/s | 6.82 km/s |
| Hermes conservative t50=7.1 | mid_3_10_kpc | 16 | 1.438 | 12.0% | -21.00 km/s | 26.23 km/s |
| Hermes conservative t50=7.1 | outer_gt10_kpc | 27 | 6.109 | 86.3% | -53.05 km/s | 55.87 km/s |
| MOND fixed simple | inner_0_3_kpc | 15 | 0.022 | 4.8% | -0.06 km/s | 2.32 km/s |
| MOND fixed simple | mid_3_10_kpc | 16 | 0.025 | 5.9% | 0.95 km/s | 2.90 km/s |
| MOND fixed simple | outer_gt10_kpc | 27 | 0.228 | 89.3% | -9.39 km/s | 8.55 km/s |

## Comparator readout

MOND cleanly fits the frozen M33 source-native board under the fixed simple-interpolation comparator: reduced chi-square = 0.119, median absolute residual = 3.69 km/s, and mean residual = -4.12 km/s.

The Hermes outer-tail underprediction is not reproduced at the same amplitude by MOND. In the outer disk (R > 10 kpc), Hermes has mean residuals of -51.6 km/s and -53.0 km/s for the two age brackets; MOND has -9.4 km/s. Hermes outer-zone reduced chi-square is 5.794 / 6.109; MOND outer-zone reduced chi-square is 0.228.

The outer disk still contributes most of MOND's total chi-square (89.3%), but the absolute penalty is very small because the total MOND chi-square is only 6.889. This means the M33 outer-tail residual structure is primarily Hermes-specific under the frozen board, not a shared failure of the source-native baryonic reconstruction.

No fitted-MOND branch was run. No MOND M/L adjustment was made. No board reconstruction choices were changed after the Hermes result.
