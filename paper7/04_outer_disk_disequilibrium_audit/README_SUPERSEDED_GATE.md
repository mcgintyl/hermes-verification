# Read this before using the Hermes values in this folder

This folder is the **historical record** of the outer-disk disequilibrium audit.
Its **Hermes** values were computed with the superseded gate copy of the first M33
run (see `../01_hermes_result_locked/README_SUPERSEDED_GATE.md`): the Hermes rows
of `m33_outer_disk_mask_metrics.csv` (for example a full board of 3.1258), the
Hermes numbers in `m33_outer_disk_disequilibrium_audit_memo.md` and the Hermes
curves in the two plots.

Still current: the **MOND** rows (MOND does not use the gate), the radial masks
themselves, and the literature map in `m33_outer_disk_literature_feature_map.csv`,
which does not depend on any model run.

**For all reporting and reproduction purposes, use:**

- `../expected/table2_masks.csv`: Table 2 on the canonical gate, the same eight
  masks for Hermes at 6.0 and 7.1 Gyr and for MOND, with every residual column
  (full board 3.22, outer R > 10 kpc 5.86 with 84.6% of the chi2);
- `python paper7/reproduce.py`, which rebuilds and checks it.
