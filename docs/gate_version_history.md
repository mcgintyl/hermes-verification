# Gate version history

The spatial gate φ(R) was **structurally continuous** across the Hermes
derivation, but **not numerically identical**. Two implementations exist, and
they give different Stage-1 β values. This matters for anyone reproducing the
*derivation path* (not just the frozen model), so both are versioned here.

## The two gates

| Quantity | Historical Tier-1 gate (Jan 2026) | Current Config G gate (published) |
|---|---|---|
| Outer shear | `\|(R/V) · dV/dR\|` (chain rule, linear grid) | `\|d ln V / d ln R\|` (log grid) |
| Curvature stabilizer ε_κ | `1e-12` (fixed) | `1e-4 · median(\|dV/dR\|)` |
| Zone normalization | arithmetic **mean** | **median** |
| Baryonic gas term | unsigned quadrature | `sign(Vgas)·Vgas²` |

Everything else is identical: knee `a_knee = 1585`, Savitzky-Golay window 5/11
poly-3 mirror, logistic mixer `x0 = 1.40 / s = 0.30`, EMA `λ = 0.5`, gate
`φ = 1 − exp(−S_eff)`. (For the seven Stage-1 galaxies the gas is positive, so
signed vs unsigned gas makes no difference to the numbers below.)

## Consequence for the Stage-1 β fit

Fitting **one global β** to minimize the unweighted mean of per-galaxy χ²/N over
the seven young galaxies gives:

| Gate | Global Stage-1 β |
|---|---|
| Historical Tier-1 (means) | **2.806585** → the documented lock 2.807 (rounded 2.81) |
| Current Config G (medians) | **2.508314** |

Intermediate check: keeping the current derivatives but swapping medians → means
gives **2.793**, isolating the normalization change as the dominant effect.

The historical means gate also reproduces the Tier-1 spec's F583-4 golden φ
vector to `< 5e-5`.

## A third lineage: the Paper 1 scoring gate (1.312)

The two gates above are not the whole story. Paper 1's published per-galaxy
score column (`chi2nu_configg`, median **1.3115 ≈ 1.312**) is reproduced by
**neither** of them. It comes from a third combination:

> **chain-rule shear `|(R/V)·dV/dR|`, with every other component at current
> settings** (`ε_κ = 1e-4·median`, median zone normalization, same knee, SavGol,
> mixer and EMA).

Changing only the shear line in the current gate moves the checker from
**0/133 to 133/133** against the published column (max 3.6e-13, median
1.311534725). The two shear forms are algebraically identical —
`d ln V / d ln R = (R/V)(dV/dR)` — but not numerically identical, because
`np.gradient` is applied to different arrays on SPARC's non-uniform radial grid.

Note this is **not** the historical Tier-1 gate: that gate also changes ε_κ,
normalization and gas convention, and tested as a unit it does *not* reproduce
the scores. Only the isolated derivative swap does.

The negative-`g_bar` convention is **score-degenerate** here: signed and clipped
give identical published scores (difference exactly 0.0), because the affected
radii produce zero model velocity either way. The derivative convention is
uniquely identified; the clipping convention is not.

An earlier explanation attributing the 1.312/1.323 gap to NumPy/SciPy version
drift was **refuted by experiment** — the same verifier under the original
production environment (Python 3.13.2 / numpy 2.4.3 / scipy 1.17.1) is
bit-identical, 0.0 across all 133 galaxies.

## What this does and does not affect

- It **does not** change Config G's frozen gate outputs. `verify_gates.py`
  passes — but note carefully what it checks: **gate outputs (`phi_last`), not
  scores.** The current gate reproduces the published `phi_last` column 133/133
  and the published *score* column 0/133. Matching one endpoint per galaxy says
  nothing about the interior φ(R) values the score integrates over. Do not read
  a `verify_gates.py` pass as reproducing Paper 1's scores.
- It **does** mean the *historical* Stage-1 lock (2.807) is reproducible only
  with the historical gate. A global fit with the current public gate correctly
  returns ~2.508 — that is the expected result of the current code, not an error.
- Earlier documents that described the geometry as "reused unchanged" across
  versions were correct in structure but imprecise in these three numerical
  details. This note is the authoritative record of the difference.

## Reproducing both numbers

`historical/tier1_gate_v1.py` implements the historical means gate;
`historical/fit_tier1_beta.py` reproduces both β values and the F583-4 golden
check from public SPARC data (Zenodo 16284118) and is run in CI:

```
python historical/fit_tier1_beta.py --sparc /path/to/rotmod_dir
```

The historical gate is retained **only** for derivation-path auditing. Do not
substitute it for the current gate in any Config G computation.
