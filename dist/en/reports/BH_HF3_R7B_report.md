# BH-HF3-R7B — final late-ringdown check through 105M

**The planned third-period test at r=40M is complete and compatible with Kerr ringdown at three radii. HF4 can begin.** This completes one control task, not all early dynamics or an independent experimental GR test. The data are numerical QC0, not the directly observed GW150914 signal.

## Technical audit

`BH_HF3_R7A_output.zip` confirms t=105M, it=53760, seven recovery messages from 95M, exit 0 and about 25m18s continuation. Checks covered 125 wave series, the exact 95M seam, full new time grid, three common-horizon QLM slot-2 files and program/parameter hashes. The 136 shared earlier files match R6A byte for byte. R6B propagation/1/r analysis was not repeated. Logs/inventory confirm both checkpoint parts; their gigabyte contents were not included/read here. No evolution was rerun.

## Reference and fitting

The presaved `BH_HF3_R7B_protocol.json` fixes the previously computed numerical Kerr spectrum using Mf=0.977461203278627M, Jf=0.648087073319352M², χ=0.6783194994253878:

ω_K=0.5353401318639867M⁻¹, α_K=0.08345059699511974M⁻¹, T=11.73680980221365M.

The spectrum was not fitted to new waves; see [qnm](https://qnm.readthedocs.io/en/latest/README.html), [Stein 2019](https://arxiv.org/abs/1908.10377). Decimal length provides reproducibility, not physical accuracy.

The third period is age 2T–3T after each local maximum. For 40M: 92.4736196044–104.2104294066M, comprising 47 original samples at 92.50–104.00M, without invented endpoints.

Fit z=C exp[(−α−iω)(t−t₀)]. The fixed model varies complex C; the free model also varies ω,α. Both minimise equally weighted complex RSS=Σ|z−z_model|². C is projected analytically; four numerical starts solve ω,α. This corrects R6B’s separate phase/log-amplitude regressions. Only selected old reference windows were reprocessed with the common criterion.

## Main result

| Radius | ω / M⁻¹ | α / M⁻¹ | ω difference | α difference |
|---|---:|---:|---:|---:|
| 15M | 0.535321949 | 0.082188152 | 0.003% | 1.513% |
| 30M | 0.535714470 | 0.082850136 | 0.070% | 0.720% |
| 40M | 0.535913909 | 0.082331989 | 0.107% | 1.340% |

All meet prior ≤1% frequency/≤5% damping targets. These point-estimate differences are not confidence intervals or total error.

| Radius | Kerr residual | Free residual | Conditional ΔBIC |
|---|---:|---:|---:|
| 15M | 2.387% | 2.355% | −6.524 |
| 30M | 2.594% | 2.585% | −8.407 |
| 40M | 2.661% | 2.632% | −7.029 |

Relative residual=√[RSS/Σ|z|²]. Two extra parameters improve it only slightly; roughly 2.6% remains unexplained. ΔBIC=2N log(RSS_K/RSS_free)−2log(2N), with N complex samples and 2/4 real parameters, is below the old +10 threshold. Correlated numerical residuals lack a calibrated independent-error model: no hypothesis probability or calibrated GR test follows.

## Sensitivity and prediction

Shifts of −0.5,−0.25,0,+0.25,+0.5M at fixed length, plus even/odd sampling, yield frequency differences 0.035–0.118%, damping 0.866–2.080% and Kerr residuals 2.562–2.742%. This checks sampling sensitivity, not resolution convergence. Two analytical signal controls recover known parameters; all four optimisation starts agree without boundary solutions.

Training on the first two thirds of the 40M third period predicts the final 15 samples with errors 2.230% (fixed Kerr) and 5.509% (free). This is one chronological control, not independent observations. Free damping on that shortened fit differs by 5.34%; full-period stability cannot be assumed on any short window.

| Period at 40M | ω difference | α difference | Kerr residual |
|---|---:|---:|---:|
| 1 | 7.486% | 67.613% | 22.878% |
| 2 | 1.717% | 16.444% | 6.351% |
| 3 | 0.107% | 1.340% | 2.661% |

| Outer window | ω | α | ω difference | α difference |
|---|---:|---:|---:|---:|
| 50M, period 2 | 0.526038031 | 0.068636573 | 1.738% | 17.752% |
| 60M, period 1 | 0.494300578 | 0.026086923 | 7.666% | 68.740% |

Signal age after the local peak matters: early disagreement with one late mode is not by itself disagreement with GR. Third-period waveforms agree after radius scaling and one constant phase rotation, without extra time/amplitude fitting:

| Pair | Coherence | Phase-aligned difference |
|---|---:|---:|
| 15–30M | 0.999963815 | 2.521% |
| 15–40M | 0.999969149 | 1.184% |
| 30–40M | 0.999985125 | 1.550% |

Ψ₄²,⁻²≈conjugate(Ψ₄²²) differs by <8×10⁻¹⁶, checking internal imposed symmetry, not merger nature independently.

![R7B checks](BH_HF3_R7B_overview.png)

## Late state and scope

Over 100–105M: Mf=0.9774491229M, Jf=0.6481486193M², χ=0.6784006850. Relative to the fixed reference, changes are about −0.00124% in mass and +0.01197% in χ. The reference intentionally remains fixed; final digits do not measure physical-state certainty. Yitian parameters are not substituted. Its A3/A4 early dynamics/quadratic link remain unconfirmed in full.

`BH_HF3_R7B_THIRD_PERIOD_R15_R30_R40_KERR220_COMPATIBILITY_SUPPORTED_HF4_READY`

No further extension is needed for this task. Conclusions concern one configuration/resolution, finite radii and selected windows: no infinity extrapolation, grid-convergence or independent statistical GR test. Extra dimensions, wormholes and an independent physical “hyperform” are not established.

The historical next step HF4-A1 compares moment definitions, angular basis/orientation/measure, separates scalar maps from metric requirements, checks ℓ≤8 truncation/display scales, then defines admissible scientific visualisations.

Reproduce with `python BH_HF3_R7B.py`; format saved results with `python BH_HF3_R7B_finish.py`. The original archive contains both input ZIPs, programs, protocol, tables/checks/state. NumPy/SciPy/Matplotlib are required; no Einstein Toolkit is launched by these commands.
