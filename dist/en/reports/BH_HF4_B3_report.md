# BH-HF4-B3 — physical area measure and localisation

20 September 2026, following B2. **Broad organisation survives replacing coordinate weights by area reconstructed from the normalised tetrad:** Re Ψ₂ remains longitudinal, Im Ψ₂/|σ| polar. Five late-state QLM area controls pass. Direct same-early-run area comparison is still pending at this historical stage and passes to B4.

## Method and reuse

HF3, Yitian, A1/A2, B1/B2 and evolution were not rerun. B2 fields/masks/weights and axes were reused; <5×10⁻⁹ axis rounding does not change diagnostic-region membership. Inputs: 12 common surfaces from the birth package and five late extracted states from `BH_HF4_A3_inventory.zip`. Late controls do not fill early gaps or establish a continuous movie.

Exported lᵘ,nᵘ,mᵘ are contravariant; mᵘ comes from two physically orthonormal spatial tangents and m⁰=0. Compute Eⁱ_A=∂_A Xⁱ and solve Eⁱ_A mᴬ=mⁱ. A Euclidean pseudoinverse converts components; it is not the physical metric. Then qᴬᴮ=2Re(mᴬconjugate(mᴮ)), q_AB=(qᴬᴮ)⁻¹, dA=√det(q_AB)dθdφ. The saved `qlm_twometric.F90`, `qlm_tetrad1.F90`, `qlm_derivs.F90` confirm conventions.

Interior grid: θᵢ=(i+1/2)π/35, φⱼ=j2π/72. Coordinates/saved shape verify it; see [QLM output documentation](https://einsteintoolkit.org/thornguide/EinsteinAnalysis/QuasiLocalMeasures/documentation.html). Fourth-order centred derivatives use polar reflection plus φ→φ+π and periodic azimuth. Reconstructed late ghost derivatives agree with actual saved ghosts to order 4×10⁻¹⁵.

Positive first-kind Fejér quadrature in cosθ integrates √det(q)/sinθ, with periodic φ summation. Midpoint quadrature is a control on the same grid, not a second simulation.

## Area controls

| t / M | QLM area / M² | Reconstructed / M² | Relative difference |
|---|---:|---:|---:|
| 30 | 41.5843928343 | 41.5843928385 | 9.92×10⁻¹¹ |
| 50 | 41.6642194343 | 41.6642194572 | 5.49×10⁻¹⁰ |
| 70 | 41.6602055236 | 41.6602055584 | 8.35×10⁻¹⁰ |
| 95 | 41.6555859651 | 41.6555859793 | 3.40×10⁻¹⁰ |
| 105 | 41.6527152688 | 41.6527153116 | 1.03×10⁻⁹ |

This validates same-simulation processing, not physical error or convergence; shared systematics remain. Midpoint differs by about 0.025–0.026%; no area-fitting factor was used.

Checks: known anisotropic metric error <10⁻¹⁵; analytical sphere/quadrature moments through degree34 near machine precision; all early metrics finite/positive definite and m⁰=0; early tangency residual ≤6.60×10⁻¹⁵; full-null-tetrad metric agreement order10⁻¹⁵ (internal, not independent); phase rotation preserves density within1.45×10⁻¹⁵. Deliberately using second-order derivatives worsens tangency and shifts early area by0.5–1%; this checks conventions, not simulation convergence.

| Time | AH3 age | Physical area / M² | Coordinate area / M² |
|---|---:|---:|---:|
| 18.625M | 0 | 37.305276 | 4.401483 |
| 18.9375M | 0.3125M | 38.629463 | 4.969929 |
| 20.5625M | 1.9375M | 40.545329 | 5.910045 |

Coordinate areas here use the same derivatives/quadrature with a Euclidean metric, not the earlier triangle estimate. The difference shows why image size is not physical size. Early same-run scalar area comparison remains open here.

## Physical-weight masks

Keep the B2 top≈5% weighted-median-deviation rule and fields/times/axes, but recompute medians/thresholds/integrals with physical weights. Diagnostic |normalised coordinate|≥0.65 remains fixed. Fractions concern the selected mask over12frames.

| Field | Main location | Fraction there | B2/B3 Dice |
|---|---|---:|---:|
| Re Ψ₂ | Longitudinal ends | 94.53–100% | 0.912–0.940 |
| Im Ψ₂ | Poles | 73.03–100% | 0.943–0.961 |
| Magnitude of σ | Poles | 100% | 0.935–0.946 |

Dice uses physical measure. Boundary shifts are expected under changed measure/median. Midpoint masks give Dice1,≥0.99786,≥0.99648 respectively. B2 triangle weights multiplied by local physical/Euclidean density ratio give ≥0.98378,≥0.99832,≥0.99628, separating density from quadrature effects. Percentile90/98 controls are saved, not preregistered detections.

Adjacent physical masks within blocks have Dice0.838–1,0.930–0.997,0.959–0.997 respectively. The1.3125Mgap remains empty. Smooth low-mode patterns can yield correlation/overlap without matter transport or a separate object.

## Limits and continuation

Physical area belongs to this horizon/slice; shape/directions remain coordinate-dependent. m-phase and l→Al,n→n/A preserve metric, but not every diagnostic. |σ| is phase-independent yet depends on null normalisation; spatially varying A may change its mask. Arbitrary tetrad/slicing transformations, second resolution, early-area checks, mask significance and HF1 identity remain untested. QC0, Yitian0389 and HF10305 stay separate; no event horizon or physical Euclidean embedding is reconstructed.

`BH_HF4_B3_fields_3D.html` shows previous surfaces/fields with new masks, fixed scales, optional inner horizons and fixed-time orbit. Caps never enter integrals; gaps are not interpolated. HTML masks match NPZ; static figure inspected, full WebGL unverified then.

Next B4 directly checks early areas and collects ready frames. `BH_HF4_early_geometry_t22_inventory.txt` identifies `/home/alex/simulations/qc0_early_geometry_t22`, final it11264,t22M,exit0. The older `Cactus/qc0-horizon-multipoles-lite-t22-mpi2` path is not the new collection source. Prepared `BH_HF4_B4_collect.py` gathers55AH3VTK/shapes, compactAH1/AH2history and same-run diagnostics, read-only, preserving archives/recording gaps; synthetic tests do not mean it already ran on the user machine.

After collection, compare matching early areas then process new frames only. Status: `BH_HF4_B3_INTRINSIC_MEASURE_RECOVERED_LATE_AREA_CONTROL_PASSED_LOCALIZATION_SUPPORTED_EARLY_AREA_CONTROL_PENDING`. Continue from `BH_HF4_resume_state_B3.json`; the results ZIP contains programs/metrics/masks/viewer/figure and savedB2inputs, with original full archives identified separately.
