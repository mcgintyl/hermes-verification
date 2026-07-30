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

## What this does and does not affect

- It **does not** change Config G or the frozen 133-galaxy results. The current
  gate is the published one and reproduces those outputs to floating-point
  precision (`verify_gates.py`).
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
