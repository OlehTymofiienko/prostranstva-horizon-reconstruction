# Yitian A4 — Common-horizon components with m=4

Historical stage: 15 September 2026. English reading edition, continuing A3 from the same verified archive. QC0/SpEC evolution and A1–A3 were not repeated.

**Result:** late I44, L54 and I64 are consistent with sums of Kerr oscillations, with substantial mixing in L54/I64. Added damped templates improve some windows, but a specific 220×220 quadratic coupling is not identified. The earliest dynamics is not closed by the models tested here.

## Definitions and spectrum

The supplied I44 notebook was read without executing it. Processing follows Zℓ4=[Hℓ4*−cℓ4]exp(−4iΩt), where c is the 400–500M mean of H*, Ω=0.20878380853617037M⁻¹, Mf=0.95162M and χ=0.68644. QC0 parameters are not substituted.

| Series | HDF5 family | Column k=ℓ(ℓ+1)+m | Late-mean magnitude | Late RMS scatter |
|---|---|---:|---:|---:|
| I44 | geo_mass/Adaptive | 24 | 1.8966×10⁻¹¹ | 9.0355×10⁻¹² |
| L54 | geo_spin/Adaptive | 34 | 9.2874×10⁻¹² | 6.3778×10⁻¹² |
| I64 | geo_mass/Adaptive | 46 | 3.8389×10⁻¹² | 2.6017×10⁻¹² |

The ±4 relation tests reality of the scalar field, not independent resolution. Unlike axisymmetric moments, these moments vanish in an exactly axisymmetric limit. Late scatter is not treated as the known uncertainty of early data.

Models are sums Ck exp[−i(ωk−iαk)(t−t0)]. Complex amplitudes are fitted simultaneously using normalized-column least squares; frequencies, damping, mass and spin are fixed. Spin-weight −2 QNMs supply a spectral model for moments defined with scalar harmonics.

| Component | ω, M⁻¹ | α, M⁻¹ |
|---|---:|---:|
| 440 | 1.188363807 | 0.089158509 |
| 540 | 1.398835249 | 0.090675201 |
| 640 | 1.608129087 | 0.091728428 |
| 740 | 1.817063866 | 0.092482010 |
| 840 | 2.025912455 | 0.093041537 |
| 441 | 1.183434870 | 0.268044403 |
| 442 | 1.173918243 | 0.448579609 |
| 443 | 1.160455948 | 0.631677823 |
| Q | 1.106942662 | 0.170843662 |

F1…F5 successively contain 440, then 540,640,740,840. O1…O3 add the first one, two or three overtones 441–443. The 17 families are F1–F5; F1O3; F2O1/O2/O3; F4O3; Q; F1Q/F2Q/F3Q/F4Q; F2O3Q; and F4O3Q. Each component has two real amplitudes, so F4O3Q has 16. Reference choices are F2 for I44 and F4 for L54/I64; the latter F4 choice is our control, not an exact mode list specified for those moments in the paper.

Training is [t0,t0+30]M and continuation is (t0+30,t0+50]M, with t0/M=0,5,10,15,20,30,40,50,60,70. All samples are retained. Errors are 100‖Z−Zmodel‖₂/‖Z‖₂ over each interval, not pointwise geometric errors or probabilities. Amplitudes are not refitted on validation, but background/final parameters use the full simulation. Model/window choices are exploratory, with no independent final model-selection test; overlapping windows are not independent.

## Late-time mixing

Validation at 100–120M after training at 70–100M:

| Moment | F1 | F2 | F3 | F4 | F5 |
|---|---:|---:|---:|---:|---:|
| I44 | 3.0130% | 2.7335% | 2.7296% | 2.7313% | 2.7315% |
| L54 | 17.4083% | 3.0451% | 3.0322% | 3.0323% | 3.0322% |
| I64 | 27.1047% | 4.9613% | 3.2283% | 3.2284% | 3.2281% |

Adding 540 markedly helps L54; adding 540/640 helps I64. Further fundamentals provide little improvement in this window, without proving their exact absence.

| Moment / reference | Training error | Validation error | Absolute validation RMS | F4O3Q validation |
|---|---:|---:|---:|---:|
| I44 / F2 | 0.4126% | 2.7335% | 1.5447×10⁻⁷ | 0.5554% |
| L54 / F4 | 0.3329% | 3.0323% | 7.2860×10⁻⁸ | 0.3266% |
| I64 / F4 | 0.3005% | 3.2284% | 2.4725×10⁻⁸ | 0.6422% |

The larger model's success does not establish every fitted component as a physical mode; individual amplitude stability must also be tested.

## Earliest window

| Moment | Model | Training 0–30M | Validation 30–50M |
|---|---|---:|---:|
| I44 | F2 | 60.0159% | 45.8932% |
| I44 | F4O3 | 2.9671% | 53.8131% |
| I44 | F4O3Q | 1.1844% | 45.8100% |
| L54 | F4 | 11.8695% | 54.1114% |
| L54 | F4O3 | 1.3879% | 24.1354% |
| L54 | F4O3Q | 1.1510% | 40.9730% |
| I64 | F4 | 38.3026% | 124.9139% |
| I64 | F4O3 | 3.3185% | 40.1991% |
| I64 | F4O3Q | 2.4180% | 155.9566% |

Q worsens early continuation for L54 and I64 when added to F4O3. Large relative errors refer to an already damped signal; absolute errors and amplitudes are archived. Failure of these models does not establish a failure of general relativity.

![Continuation errors](BH_YITIAN_A4_prediction.png)

![Example series](BH_YITIAN_A4_examples.png)

## What the quadratic template tests

Q=exp[−2i(ω220−iα220)(t−t0)] doubles frequency and damping, and m=2+2 permits m=4. This algebra does not derive a horizon evolution equation or coupling coefficient. Q is a diagnostic template, not a detected nonlinear mode.

For training at 20–50M and validation at 50–70M:

| Moment | Fundamentals | Fundamentals + Q | Fundamentals + squared A2 I22 model |
|---|---:|---:|---:|
| I44 | 19.8076% | 4.3522% | 4.2492% |
| L54 | 22.5039% | 4.9000% | 5.1824% |
| I64 | 28.6917% | 8.3737% | 8.7808% |

The final column uses the saved A2 model 220+320+221+222+223 trained on the same interval, squares its continuation and fits one complex coefficient K. Actual validation I22 data are not used. The parent model has ten additional real amplitudes; this is not an equal-complexity comparison.

| Moment | ωQ−0.03M⁻¹ | ωQ=2ω220 | ωQ+0.03M⁻¹ |
|---|---:|---:|---:|
| I44 | 7.7059% | 4.3522% | 1.0354% |
| L54 | 5.2741% | 4.9000% | 8.8682% |
| I64 | 6.4194% | 8.3737% | 14.5156% |

These diagnostic frequency shifts, at fixed damping, do not uniquely select the doubled frequency and are not uncertainty bounds. K, defined from the Q amplitude divided by the squared 220 amplitude and referred to 20M, is also unstable:

| Moment | Magnitude K at t0=20M | Magnitude K at t0=40M | Magnitude K at t0=70M |
|---|---:|---:|---:|
| I44 | 0.41251 | 0.26072 | 1.75151 |
| L54 | 0.18742 | 0.14097 | 0.70919 |
| I64 | 0.06926 | 0.05837 | 0.20175 |

Its complex phase changes too. At late starts a physical quadratic contribution would have decayed rapidly, increasing sensitivity to other residuals. This procedure does not establish a constant coupling; it does not prove that nonlinear interaction is absent. The normalized 440/Q column inner product over 30M is 0.908403, limiting their separation.

![Diagnostic coupling](BH_YITIAN_A4_coupling.png)

## Controls, checkpoint and reproduction

Controls include late windows 350–500,400–500,450–500M and no subtraction; training spacings 0.1,0.2,0.4M; Mf/χ rounding ±0.000005; training lengths 20,30,40M; damping-compensated residual weights; and Q damping multiplied by 0.8 or 1.2. Varying training length also changes the validation interval. Late offsets have little effect here, while sampling and weighting can alter early errors.

Removing the angular rotation while consistently shifting frequencies to ω−4Ω changes errors by at most 1.471×10⁻¹⁵. Deliberately inconsistent frame/sign controls are retained. A known complex signal with background correction gives 2.828×10⁻¹⁶ error. These are numerical-procedure checks.

There are 510 main fits, 672 sensitivity variants, 48 phase controls and 30 parent-square checks. Expansion failed to improve continuation in 102 of 270 nested comparisons; the maximum normalized primary-matrix condition number is 1014.0. One numerical resolution and an AMR target do not provide per-moment uncertainty, confidence intervals or a probability of a new structure.

A4 is complete for I44, L54 and I64. Early spectral closure and a physical nonlinear coupling remain open. The historical next step was QC0 95→105M and the third period at r=40M, followed by HF4 geometry; the recovery-preparation file had already been found. Subsequent reports record that progress. A scientific movie of measured geometry does not itself prove new physics, reconstruct a global event horizon or determine the black-hole interior.

Use `BH_YITIAN_A4.py --input /path/common_horizon_moments_data_and_scripts.zip --output repeated`. NumPy, SciPy, h5py and Matplotlib are required; qnm 0.4.4 is only needed for `--recompute-spectrum`. `--render-only` regenerates the report/figures from saved JSON. Input SHA256: `9cdd268e85b783d8a832b72e46945c6b43bb2b3ea0a9e6665e57a5a1b7d56f0f`.

Sources: [Chen et al., section IV D](https://arxiv.org/abs/2208.02965), [qnm documentation](https://qnm.readthedocs.io/en/latest/README.html). Continuation tests, Q and the parent-square/specificity diagnostics are our analysis of the supplied series.
