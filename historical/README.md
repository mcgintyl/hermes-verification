# ⚠️ SUPERSEDED — historical reproduction only

**Do NOT use anything in this folder for new work.**

These files reconstruct a **January 2026 gate that is no longer in use**. They
exist for one reason: so the historical Stage-1 calibration (β ≈ 2.807) can be
reproduced and audited. They are not part of any published result.

## The canonical gate

> **`paper1/hermes_gate_phi.py` — `hermes_phi()`**

That is the one and only gate for current work. Every published result uses it.
(`paper1/verify_hermes.py` carries a self-contained copy of the same function so
the verifier runs standalone; the two are checked against each other by
`paper1/verify_gate_audit.py` and agree exactly.)

## Why mixing them silently breaks things

The historical gate differs from the canonical gate in four ways — see
[`../docs/gate_version_history.md`](../docs/gate_version_history.md):

| | historical (this folder) | canonical (`paper1/`) |
|---|---|---|
| Outer shear | `\|(R/V)·dV/dR\|` (chain rule) | `\|d ln V / d ln R\|` (log grid) |
| Curvature stabilizer | `1e-12` | `1e-4 · median(\|dV/dR\|)` |
| Zone normalization | arithmetic **mean** | **median** |
| Baryonic gas term | unsigned quadrature | `sign(Vgas)·Vgas²` |

They are **not interchangeable, and they are not mix-and-match.** Feeding the
current signed `g_bar` into `tier1_gate()` — or the historical unsigned `g_bar`
into `hermes_phi()` — produces a hybrid that matches *no* configuration, current
or historical. Measured across the 133-galaxy sample, that particular hybrid
shifts model velocities by a median of **5%**, exceeding 5% in 64 galaxies and
reaching **20%** (NGC 5371). It fails silently: no error, just wrong curves.

If you use this folder at all, use **both** historical pieces together —
`gbar_unsigned()` **and** `tier1_gate()` — or neither.

## Files

| File | Purpose |
|---|---|
| `tier1_gate_v1.py` | Reconstructed Jan-2026 gate. Historical only. |
| `fit_tier1_beta.py` | Reproduces the historical lock (β = 2.806585) and reports the current-gate value (2.508314) for contrast. Run in CI. |

Both modules are loaded by explicit path rather than via `sys.path`, so
`import tier1_gate_v1` does **not** silently become available to other code in
the same process. Please keep it that way.
