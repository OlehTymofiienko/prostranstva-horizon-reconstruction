# Yitian A2 — Early complex moments and a corrected QC0 reference

Historical stage: 15 September 2026. English reading edition. QC0 remained at 95M, iteration 48640, with both checkpoint files preserved. No evolution or previously completed audit was repeated.

**Result:** mode mixing matters, and a good early fit does not guarantee accurate continuation. A separate audit corrected the interpretation of QC0's final parameters without replacing them with the SpEC values.

## Correction to the earlier prose

The saved R6B program and output use the following QC0 quantities:

| Quantity | QC0 | Yitian / SpEC |
|---|---:|---:|
| Final mass Mf/M | 0.977461203279 | 0.95162 |
| Angular momentum J/M² | 0.648087073319 | — |
| Dimensionless spin χ=J/Mf² | 0.678319499425 | 0.68644 |

The number 0.648087 is J, not χ. R6B's code had used the correct dimensionless spin; the error was in the earlier explanation. A different quoted mass also contained transposed digits. The simulations remain distinct physical configurations.

## Models and validation

For I22 and L32, stored moments are conjugated, their 400–500M mean is subtracted, and exp(−2iΩt) is applied, with Ω=0.20878380853617037M⁻¹. The Kerr spectrum is calculated with qnm 0.4.4 at fixed Mf=0.95162M and χ=0.68644. Only complex amplitudes are fitted by least squares with normalized columns.

| Component | ω, M⁻¹ | α, M⁻¹ |
|---|---:|---:|
| 220 | 0.553471331 | 0.085421831 |
| 221 | 0.541024406 | 0.258312615 |
| 222 | 0.518015073 | 0.436245721 |
| 223 | 0.487436421 | 0.618669512 |
| 320 | 0.792030127 | 0.089004353 |
| 420 | 1.017167729 | 0.091425692 |

The 220/320 beat timescale is about 26.3381M. The spin-weight −2 QNM spectrum is used as a spectral model; the horizon moments themselves use ordinary scalar spherical harmonics.

F2 denotes 220+320; F3 adds 420; O3 adds 221,222,223. Training uses [t0,t0+30]M and continuation uses (t0+30,t0+50]M, for t0/M=0,5,10,15,20,30,40,50,60,70. Reported errors are 100‖Z−Zmodel‖₂/‖Z‖₂ on each whole interval. Validation amplitudes are not refitted. The late background and final parameters use information from the full simulation; this is not a wholly blind forecast, and overlapping windows are not independent experiments.

| Moment | t0/M | Model | Training error | Continuation error |
|---|---:|---|---:|---:|
| I22 | 0 | 220 | 22.3616% | 41.5247% |
| I22 | 0 | F2 | 16.7113% | 33.6142% |
| I22 | 0 | F2O3 | 0.1209% | 2.0975% |
| L32 | 20 | 320 | 86.5950% | 101.8891% |
| L32 | 20 | 220 | 13.1928% | 15.2437% |
| L32 | 20 | F2 | 0.2993% | 1.8319% |
| L32 | 20 | F3 | 0.2792% | 1.8708% |
| I22 | 20 | 220 | 0.3645% | 1.9391% |
| I22 | 20 | F2 | 0.2557% | 2.0712% |
| I22 | 20 | F2O3 | 0.0458% | 1.0563% |

Nineteen of 160 nested comparisons worsened on continuation when the model was expanded. The maximum normalized matrix condition number was about 158. A synthetic three-component signal gives continuation error 4.33×10⁻¹⁶; a consistent change back to the unrotated frame changes errors by at most 4.93×10⁻¹⁶. These check the procedure, not physical accuracy.

## Sensitivity and limits

Controls vary the late mean (350–500, 400–500, 450–500M), training spacing (0.1,0.2,0.4M), and the reported Mf and χ by ±0.000005.

| Case | Mean-window range | Sampling range | Parameter-rounding range |
|---|---|---|---|
| I22, F2O3, t0=0 | 2.09747–2.09747% | 2.09747–2.14075% | 2.09375–2.10119% |
| L32, F2, t0=20 | 1.83185–1.83185% | 1.82852–1.83185% | 1.82450–1.83921% |
| I22, F2, t0=20 | 2.07121–2.07121% | 2.07121–2.07142% | 2.06376–2.07867% |

A sampling check is not spatial-resolution convergence. The AMR target does not define a statistical noise distribution; a residual is not automatically a new structure. Model complexity and conditioning limit component-level interpretation.

## QC0 Kerr reference and the remaining control

For QC0, the earlier approximate reference ω=0.538108542 and α=0.084117941 is replaced by the numerical spectrum ω=0.535340132 and α=0.083450597, in M⁻¹. Comparing the already saved third-period estimates without refitting gives frequency/damping differences of 0.2023%/0.7398% at r=15M and 0.2097%/0.1186% at r=30M.

Old BIC results cannot simply be transferred to a changed reference and objective. A consistent complex fit must be performed, with attention to correlated numerical residuals.

The corrected period is 11.736810M. The third period at r=40M ends near 104.210429M; a continuation to 105M is therefore sufficient. Iteration 53760 was an expectation at this stage, not an already completed run. The original recovery-preparation file had been found; the user did not need to locate it again.

The recorded route was A3 (axisymmetric moments), A4 (m=4), the final QC0 continuation/control, then HF4 geometric reconstruction. Displaying 64 or 128 angular directions does not create independent observers or new numerical evidence. Later reports record completion of the QC0 continuation.

Reproduction: run the archived `BH_YITIAN_A2.py` with `--input` and `--output`. NumPy, SciPy, h5py and Matplotlib are required; qnm is required only with `--recompute-spectrum`. Saved numerical spectra can otherwise be reused.
