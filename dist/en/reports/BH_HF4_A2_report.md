# HF4-A2 — Atlas of the measured horizon scalar fields

Historical stage: 16 September 2026. English reading edition. This stage continues A1 without repeating spacetime evolution.

**Result:** a reproducible atlas of the mass-curvature and spin-connection scalar maps, with explicit angular conventions and fixed display scales. It is a map on the authors' angular representation, not a uniquely reconstructed three-dimensional horizon surface.

## Data and representations

The same verified Yitian HDF5 file was used (SHA256 `9cdd268e85b783d8a832b72e46945c6b43bb2b3ea0a9e6665e57a5a1b7d56f0f`). Full fields are K=2ΣIℓmYℓm and B=ΣLℓmYℓm. Perturbations δK and δB subtract the 400–500M mean coefficients. The stationary axisymmetric background is retained in the full fields.

The analysis contains 1,501 times from 0 to 150M at 0.1M spacing; the atlas/video uses 301 frames at 0.5M spacing. The H.264 video is 1280×960, 15 fps, approximately 20.0667 seconds. Coefficients and numerical arrays are saved separately from the rendered images.

## Sensitivity to angular truncation

All comparisons use fixed truncations ℓ≤4,6,8.

| Field at t=0 | ℓ≤4 versus ℓ≤8 | ℓ≤6 versus ℓ≤8 |
|---|---:|---:|
| K | 5.20% | 2.02% |
| δK | 8.81% | 3.44% |
| B | 7.61% | 2.74% |
| δB | 36.19% | 13.93% |

Over 5–120M, maximum differences are 5.60% and 0.962% for δK, and 21.23% and 3.55% for δB, respectively. These compare available truncations; they cannot bound unmeasured moments above ℓ=8.

The perturbation norm of δK peaks at 0M, whereas δB peaks near 4.6M. Early evolution is oscillatory and need not decay monotonically.

## Angular rotation and background subtraction

The authors' rotation is exp(−imΩt), with Ω=0.20878380853617037M⁻¹; on the sphere it corresponds to φ→φ−Ωt. The rotation preserves field norms. Relative norm changes are about 3.33×10⁻¹⁶; coefficient phase controls give discrepancies of 3.38×10⁻¹⁶ and 2.78×10⁻¹⁷.

Background subtraction must be transformed consistently: U(C−C̄)=UC−UC̄. Rotating the series first and then subtracting an unrelated late average is not equivalent. That deliberately incorrect procedure gives about 31% error in K and 34% in B at 120M, and errors larger than the perturbation norm near 150M. It is not used for the atlas.

Alternative backgrounds over 300–400, 400–450 and 450–500M change δK by at most 5.44×10⁻⁹, or 0.133% of its norm at 150M, and δB by 2.12×10⁻⁹, or 0.185%. These are processing sensitivities, not total physical uncertainties.

## Rendering controls and interpretation

Real versus complex harmonic reconstruction agrees to 6.42×10⁻¹⁴. An independent JavaScript Legendre implementation agrees with Python to 7.22×10⁻¹⁶. Sixteen mode controls and the minimum DOM/canvas structure were checked. The original A2 report did not claim a complete browser interaction test; later publication checks are separate.

| Map | Fixed colour limits |
|---|---|
| K | −0.30 to 3.10 |
| B | −0.61 to 0.61 |
| δK | −1.66 to 1.66 |
| δB | −0.17 to 0.17 |

The optional asinh display uses parameter 10⁻⁵. Scales do not rescale independently per frame. A narrow sign boundary is not a physical wall. Projection distorts apparent areas; the source coefficients are the quantitative record.

## Scope and next recorded step

A3 was to inventory actual surface geometry, area, metric information and angular transport. Multipoles alone do not determine a unique physical surface, event horizon or interior geometry. A direct comparison with waveforms requires data from the same simulation and consistent time/angular conventions. No new evolution or outgoing correspondence was performed here.

The archived script `BH_HF4_A2.py` accepts the original moments file, the A1 resume-state file, an output directory and optional `--video`. Video generation additionally needs ffmpeg. The NPZ files retain double-precision analysis and atlas coefficients. The original input is supplied separately where it is not part of the public distribution.
