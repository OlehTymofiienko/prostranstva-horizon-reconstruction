# Release reproducibility

Do not rerun confirmed stages for routine continuation. Read `CHECKPOINT.json` and the final reports.

| Level | Included | External requirements |
|---|---|---|
| Website and videos | `dist/`, local images, JavaScript and data | Modern browser; WebGL for two 3D viewers |
| B5 from saved B4 | Program, B4 NPZ input, protocol, times, all results | Python, NumPy, SciPy, Matplotlib |
| Final R7B | Program, two QC0 input archives, saved results | Python and dependencies specified by the program |
| B2/B3/B4 from original surfaces | Programs, definitions, numerical results | Original QC0 surface archives listed in reports |
| HF1 from full grid | Frozen functions, criteria, masks, checks, original programs | `BH_HF1_visual_data.zip` with full SXS arrays |
| Yitian from original multipoles | Programs, conventions, derived maps/results | Original author-supplied multipole file |
| New simulation | Saved parameters and provenance where available | Einstein Toolkit/SpEC, initial conditions, computing resources |

This release is not a standalone copy of all simulations. Large source grids and the author-supplied correspondence archive are not duplicated. Historical programs may retain original paths and require separate inputs: this is an archive of the method, not a claim that the entire history runs with one command.

Large HF1/QC0 HTML viewers were split into local data scripts. Frame order and values were checked for exact equality with saved JSON. B5 uses Float32 for viewing; rounding is recorded in `viewer-packaging-checks.json`. Scientific NPZ remains Float64.

B5 requires SciPy ≥1.15 (`sph_harm_y`), NumPy and Matplotlib. Analytical normalisation, rank, conditioning and residual-orthogonality checks are built in. The original run used SciPy 1.17.0; other versions are in `BH_HF4_B5_definition.json`.

Check the publication package without scientific recalculation:

```bash
python scripts/check_release.py
```

`scripts/assemble_release.py --workspace PATH` repackages the original stage workspace. It is not needed to read the site and does not run science. Restoring that workspace requires the external files above.

## English and reading themes

```bash
python scripts/localize_site.py
```

This builds the presentation layer from `docs/en/`, `scripts/reading-ui.json` and `scripts/reading-dynamic.json`. Both languages share scientific arrays, images and videos. Original labels embedded in figures/videos and archived program comments remain in their source language. Theme selection is stored locally; no external translation service or analytics is used. Light, dark, warm and system modes change reading surfaces and plot backgrounds, not scientific colour scales.
