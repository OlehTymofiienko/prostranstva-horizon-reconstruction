# Yitian A3 — Axisymmetric common-horizon moments

Historical stage: 15 September 2026, continuing A2. English reading edition. QC0/SpEC evolution and A1/A2 were not repeated.

**Result:** small sets of fundamental Kerr oscillations describe late I20, I40 and L30 with continuation errors of approximately 2.6–3.1% over 100–120M. At early times a successful multimode fit does not ensure successful continuation. An additional nonoscillatory exponential helps in some intermediate windows, but is not identified as a distinct physical mode.

## The stationary background is part of the geometry

Write Zℓ0=Zℓ0∞+δZℓ0, with the background equal to the 400–500M mean. Original moments are retained. For m=0, exp(−imΩt)=1 and the series are real, so these results cannot be attributed to a missing frequency shift mΩ.

| Moment | Late background |
|---|---:|
| I20 | −0.490070977172 |
| I40 | 0.046545899370 |
| L30 | −0.159949839063 |

A1 already checked the Kerr consistency of these dimensionless geometric moments; they are not direct mass/spin measurements.

## Model family and continuation test

Each oscillation contributes exp[−αk(t−t0)] times `[Ak cos(ωk(t−t0))+Bk sin(ωk(t−t0))]`. Two real amplitudes are fitted at fixed Kerr frequency/damping, Mf=0.95162M and χ=0.68644. This is the consistent positive/negative-frequency pair appropriate to a real series. Zero crossings are retained.

| Component | ω, M⁻¹ | α, M⁻¹ |
|---|---:|---:|
| 200 | 0.413206065 | 0.089003429 |
| 201 | 0.390963921 | 0.272523689 |
| 202 | 0.352988107 | 0.470868646 |
| 203 | 0.310236362 | 0.686609563 |
| 300 | 0.656697666 | 0.092902813 |
| 400 | 0.883639654 | 0.094403240 |
| D | 0 | 0.170843662 |

F1=200, F2=200+300 and F3=200+300+400, with 2,4,6 real amplitudes. F1O1/O2/O3 and F3O1/O2/O3 add the first one, two or three overtones 201–203 (two amplitudes each). D is a single real exponential exp[−λ(t−t0)], λ=2α220. D alone and F1D/F2D/F3D have 1,3,5,7 real amplitudes. These are 13 exploratory model families, not a preregistered test.

The algebraic motivation for D is that the squared modulus of a damped m=2 oscillation decays as exp(−2αt), while multiplying m=2 and m=−2 permits m=0. This does not derive a physical coupling law or require a particular horizon moment to contain D.

Training is [t0,t0+30]M; continuation is (t0+30,t0+50]M, with t0/M=0,5,10,15,20,30,40,50,60,70. Errors are 100‖δZ−δZmodel‖₂/‖δZ‖₂ over each whole interval. They are not full-geometry errors, maximum pointwise errors or probabilities. Validation amplitudes are not fitted; background and final parameters still use the full simulation. Overlapping windows are not independent.

## Late and early behaviour

| Moment | Reference model | Training 70–100M | Continuation 100–120M | Absolute validation RMS |
|---|---|---:|---:|---:|
| I20 | F1 | 0.5890% | 2.5939% | 2.419×10⁻⁸ |
| I40 | F3 | 0.9842% | 3.1426% | 8.165×10⁻⁹ |
| L30 | F2 | 0.6100% | 3.0296% | 1.861×10⁻⁸ |

Late I20 is described at the few-percent level by one fundamental; the other two benefit from mixing. This agrees qualitatively with the authors' interpretation, but our continuation error is not their in-window mismatch statistic.

| Moment | Model | Training 0–30M | Continuation 30–50M |
|---|---|---:|---:|
| I20 | F1 | 55.7111% | 855.9300% |
| I20 | F3O3 | 0.4850% | 107.8136% |
| I40 | F3 | 40.1646% | 521.1302% |
| I40 | F3O3 | 0.9051% | 205.4057% |
| L30 | F2 | 52.0373% | 286.7539% |
| L30 | F3O3 | 0.4725% | 42.5852% |

Small training errors do not close the early-time problem. Large relative validation errors partly reflect the much smaller, already damped validation signal. Absolute errors and signal amplitudes remain in the archived JSON. These failures concern the chosen models and procedure, not all possible descriptions within general relativity.

For training at 20–50M and validation at 50–70M, adding D changes errors from 9.7793% to 5.1840% for I20, 30.3344% to 6.4391% for I40, and 10.2802% to 4.9213% for L30. This helps in those windows without uniquely identifying nonlinear relaxation, missing oscillations, changing amplitudes, time-coordinate effects or numerical contributions.

![Continuation errors](BH_YITIAN_A3_prediction.png)

![Example series](BH_YITIAN_A3_examples.png)

## Sensitivity and interpretation

Controls vary the mean window (350–500,400–500,450–500M), training spacing (0.1,0.2,0.4M), Mf/χ rounding by ±0.000005 and training duration (20,30,40M). Damping-compensated residual weights exp[α200(t−t0)] change early F3O3 continuation errors from 107.8136/205.4057/42.5852% to 119.3959/218.3029/47.3645% for I20/I40/L30. These are methodological weights, not a statistical noise model.

| Late reference case | Empirical background | Kerr background at rounded χ | Kerr background at χ from L10 |
|---|---:|---:|---:|
| I20, F1 | 2.5939% | 198.4908% | 2.5943% |
| I40, F3 | 3.1426% | 132.5491% | 3.1430% |
| L30, F2 | 3.0296% | 145.5301% | 3.0304% |

A rounded spin can accurately describe a full nonzero moment yet leave a large error after subtracting nearly equal quantities. These comparisons use the same canonical perturbation norm and reconstruct the full moment before evaluating residuals. The L10 result is internal consistency within the same simulation.

| Model + D | λ=0.8×2α220 | λ=2α220 | λ=1.2×2α220 |
|---|---:|---:|---:|
| I20 | 2.9483% | 5.1840% | 7.5216% |
| I40 | 15.7184% | 6.4391% | 2.7020% |
| L30 | 6.8125% | 4.9213% | 3.7837% |

These 50–70M validation values do not measure λ or define a confidence interval. Establishing a quadratic coupling requires stable amplitudes, a parent-component relation and numerical-resolution controls.

There are 390 primary fits and 558 control cases. Model expansion worsened continuation in 107 of 330 nested comparisons. Maximum normalized matrix condition number: 1732. A known real signal with a background correction gives continuation error 1.76×10⁻¹⁶.

Subsampling is not spatial convergence, late scatter is not early-time uncertainty, and one resolution plus an AMR target cannot supply trustworthy per-moment errors. No significance level or probability of a new structure was calculated.

## Historical checkpoint and reproduction

A3 is complete for these three moments; early dynamics is not declared fully explained. A4 was to examine I44, L54 and I64 with the exp(−4iΩt) rotation, ordinary mixing and separate nonlinear diagnostics. The QC0 95→105M control and subsequent geometric reconstruction remained in the route; later reports record their completion.

Use the archived `BH_YITIAN_A3.py --input /path/common_horizon_moments_data_and_scripts.zip --output repeated`. NumPy, SciPy, h5py and Matplotlib are needed; `--recompute-spectrum` additionally requires qnm 0.4.4. Saved spectra and prior context are included, while the original private input is supplied separately.

References: [Chen et al., sections IV C–D](https://arxiv.org/abs/2208.02965), [qnm documentation](https://qnm.readthedocs.io/en/latest/README.html). D and continuation controls are our exploratory analysis, not a detection claimed in the source paper.
