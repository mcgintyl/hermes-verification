# Read this before using the phi / V_hermes values in this folder

This folder is the **historical record** of the first M33 Hermes run. Its gate
values (`phi`, `gate_S_eff`, `Vhermes_*`, chi2 columns in
`m33_hermes_per_point_predictions.csv`, and `m33_hermes_execution_script.py`)
were produced by a superseded gate implementation and **do not match the
canonical gate** (`paper1/hermes_gate_phi.py`). Per-point phi differs by up to
0.625.

**For all reporting and reproduction purposes, use
`../01b_hermes_result_canonical_gate/` instead:**

- `m33_hermes_per_point_predictions_canonical.csv` — canonical phi and V_hermes
  (also carries the superseded phi in `phi_superseded_locked` for comparison)
- `compute_m33_hermes_canonical.py` — the script that performs the M33 Hermes
  calculation with the canonical gate
- `m33_canonical_gate_correction_memo.md` — what changed and why
  (full-board chi2nu 3.13 -> 3.22; no conclusion affected; Paper 7 v3,
  DOI 10.5281/zenodo.21347359, reports the canonical numbers)

The **board inputs** in this folder (`m33_corbelli_reconstructed_board_input.csv`:
R, Vobs, errV, Vgas, Vdisk, Vbulge) are unaffected by the gate correction and
remain the frozen Format B board used by both runs.

If you compute phi yourself following the documentation, you will match the
canonical gate — and you will NOT match the phi column in this folder. That is
expected; it is the superseded run.
