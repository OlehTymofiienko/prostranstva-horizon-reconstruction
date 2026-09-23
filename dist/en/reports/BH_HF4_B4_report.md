# HF4-B4 — complete available common-horizon series and early-area validation

21 September 2026. QC0 early run `qc0_early_geometry_t22`.

**The early-area check left open by B3 is complete. All 55 common apparent-horizon frames span 18.625–22M. Broad longitudinal/polar localisation survives the extended interval.** This validates processing and descriptive localisation, not new physics, full tetrad independence or identity with HF1.

## Sources and reuse

`BH_HF4_B4_collection(1).zip` has 818 entries: 816 source files plus collection description/explanation. All 816 sizes/hashes were verified. Archive SHA-256: `34d86aa53bc74f0e1b872761518a145b7612aafbb2f7e735e6993961140936fc`.

| Contents | Count |
|---|---:|
| AH3 VTK surfaces | 55 |
| AH3 original angular-patch shapes | 55 |
| AH1 shapes | 349 |
| AH2 shapes | 349 |
| Diagnostic series | 7 |
| Early-run parameter file | 1 |

The previous 12 VTK files match B3 byte for byte; saved metrics, weights and masks were reused, with 72 arrays checked for exact equality. Thirty old adjacent-mask comparisons were reused. Only 43 new frames and 132 new comparisons were calculated. Evolution, HF3, Yitian and earlier reconstruction were not rerun.

## Run identity and horizon

QLM run ID: `run-qc0-horizon-early-geometry-t22-DESKTOP-GDJ2LUD.localdomain-alex-2026.09.16-10.57.31-1065`. Parameters map AH3 to SphericalSurface 2 and QLM slot 2, with fourth-order derivatives and 35×72 interior nodes plus two ghost points per side.

All 55 frames have valid-data flag 1, calculation-error flag 0, matching `qlm_time[2]`, finite positive area. At first it=9536, `qlm_iteration[2]` is sentinel −1; validation instead uses flags, time, VTK and AHFinderDirect area. It is not used for time derivatives. Later iterations match correctly.

No common horizon is found at 18.5625M; one exists at 18.625M. The latter is the first saved detection, not an exact continuous formation time.

## Physical area

B3 reconstruction is retained: qᴬᴮ=2Re(mᴬ conjugate(mᴮ)), q_AB=(qᴬᴮ)⁻¹, dA=√det(q)dθdφ. Derivatives/polar conditions follow the original run; first-kind Fejér quadrature in cosθ and periodic φ summation use no fitted matching factors.

| Time | Tetrad reconstruction, M² | QLM, M² | AHFinderDirect, M² |
|---|---:|---:|---:|
| 18.625M | 37.3052764802 | 37.3052742283 | 37.31129231 |
| 22M | 41.0176757804 | 41.0176758144 | 41.01880032 |

Maximum relative reconstruction/QLM difference is **6.0364×10⁻⁸**, or **0.00000604%**. AHFinderDirect differs from QLM by 0.00274–0.01613%. This checks same-run area consistency, not physical simulation accuracy: systematic errors may be shared and no second resolution is available.

Area increases at every saved sample, from 37.305276 to 41.017676M², a 9.9514% increase. This discrete apparent-horizon observation is not a separate proof of an area theorem. All new metrics are positive definite; maximum relative tetrad-tangency residual is 7.45×10⁻¹⁵.

## Localisation

The fixed B3 rule selects roughly 5% of physical area by greatest deviation from the weighted median. All fields use the original tetrad. Percentages below refer to the selected region, not the whole horizon. Diagnostic longitudinal/polar regions retain the normalised-coordinate threshold 0.65 and may overlap.

| Field | Main location | Selected-region fraction | Available estimates |
|---|---|---:|---:|
| Re Ψ₂ | Longitudinal ends | 94.53–100% | 53 |
| Im Ψ₂ | Polar regions | 73.03–100% | 55 |
| Magnitude of σ | Polar regions | 100% | 55 |

New longitudinal axes use centres from full AH1/AH2 angular patches. At the 12 overlapping frames they differ from B3 by ≤0.0004841°; all region memberships match. At it=10656 and 10752, individual surfaces are missing, so no longitudinal axis/fraction is invented. AH3 fields, masks, area and polar estimates remain available.

| Field | Adjacent-mask Dice |
|---|---:|
| Re Ψ₂ | 0.838–1 |
| Im Ψ₂ | 0.815–0.997 |
| Magnitude of σ | 0.959–0.998 |

The full series reveals Im Ψ₂ boundary changes missed by sparse sampling: minimum overlap falls to 0.815. Broad polar organisation persists, but the pattern is not stationary. Discrete upper-tail thresholds can change noticeably even for smooth fields.

For the strongest 2%, minimum longitudinal Re Ψ₂ fraction is 86.81% and polar Im Ψ₂ fraction 69.11%. Selecting about 10% by |σ| gives polar fraction ≥91.31%. These are definition-sensitivity checks, not independent detections or preregistered statistics. Matching indices do not prove matter transport or an independent object.

## Viewer and next checkpoint

`BH_HF4_B4_fields_3D.html` offers free rotation, field/time selection, inner horizons and camera orbit. It includes 353 output times over 0–22M; 351 have surfaces, two early times are empty. All 55 AH3 frames are present. AH1/AH2 are absent at 0.375,16.8125,20.8125,21M and gaps are explicit. Individual displays are subsampled, but centre analysis uses full grids. AH3 fields/masks use all 35×72 nodes; grey polar caps only close the image. Scales are fixed, with no temporal interpolation or HF1-mask transfer.

All 165 display masks match numerical results. Checks covered 573 control states, the old B3 gap, empty-frame stops and AH3 without inner surfaces. The static figure was inspected; full WebGL rendering was not verified in that environment. This is early coordinate geometry, not optics, a physical Euclidean embedding or a full movie through 105M.

The then-next stage **B5** freezes fields/masks/weights, measures explained variation at selected harmonic cutoffs and examines residual location. It needs no new evolution. Coordinate harmonics are not invariant geometry; residuals are not automatically new physics. Tetrad, slicing, convergence and HF1 linkage remain separate.

Status: `BH_HF4_B4_COMPLETE_AH3_SERIES_EARLY_AREA_AGREEMENT_CONFIRMED_LOCALIZATION_EXTENDED`. Continue from `BH_HF4_resume_state_B4.json`; `BH_HF4_B4_results.zip` contains programs/results/viewer, required unchanged B3 inputs and B5 protocol. Source archive remains separate by checksum; completed area/mask work need not be repeated.
