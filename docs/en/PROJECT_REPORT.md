# Prostranstva: reconstructing horizon geometry and fields

**Final research release 1.0. 22 September 2026.**

The project progressed from selecting inhomogeneous regions on a numerical horizon to examining their dynamics, reconstructing surfaces and quantifying angular structure. The planned work with the available data, including HF4-B5, is complete. The result is a set of verifiable reconstructions and descriptive findings, not proof of a separate new physical entity.

## Three different simulations

**HF1/HF1R, SXS:BBH:0305, Lev6.** This is the source of the initial finding referred to as a “latent structure” in working discussions. A connected group of nodes is selected by a fixed angular-gradient criterion. Visualisation V2 maps these masks onto the original coordinate surface of common apparent horizon C. A total of 4,698 mask matches across 1,566 snapshots were confirmed without repeating the earlier statistical tests. The film includes R* evolution, a separate camera orbit at fixed time and late control snapshots. The viewer switches between R*, E* and B*.

**QC0, Einstein Toolkit.** Our simulation reaches 105M. The HF3 wave branch and early HF4 geometry branch have explicitly separate time windows. The early run provides individual horizons through 22M and 55 consecutive common-horizon snapshots from 18.625M to 22M. Here the mask selects strong deviations from the field median using physical-area weights. This differs from the HF1 rule; physical identity of the two selections has not been established.

**Yitian, SXS:BBH:0389, SpEC.** Mass and spin multipoles through ℓ=8 provide independent material for studying geometric scalar fields. There are 5,001 samples over 0–500M. Phase transformations and the late background follow the authors’ definitions. Atlas A2 covers 0–150M and shows the effects of angular truncation and reference frame. A spherical field map is not a uniquely reconstructed early 3D shape.

The initial conditions of these three configurations are not assumed identical. The mass, spin or time origin of one branch must not be substituted into another.

## R, E, B and the coloured regions

In HF1, R* refers to the intrinsic scalar curvature of the two-dimensional surface. The electric part E of the Weyl tensor describes the tidal gravitational field; its magnetic part B is related to differential rotation of local inertial directions. These are gravitational quantities, not electromagnetic electric and magnetic fields. This dataset uses their normal components and the normalisation of the original analysis.

All three quantities are defined on the same surface. Different coloured regions arise because the criterion is applied to different fields. Mask colour does not represent matter, temperature or a separate object. A changing mask does not measure matter transport; camera rotation at fixed time does not show physical evolution.

QC0 uses Re Ψ₂, Im Ψ₂ and |σ| in the original tetrad. These quantities are not automatically identified with HF1 notation or Yitian’s multipoles. Tetrad and coordinate conventions matter on a dynamical horizon.

## Ringdown: the final QC0 reference

Final stage R7B corrects the interpretation of J: it must not be treated as the dimensionless spin χ. The current values are:

| Quantity | QC0, R7B |
|---|---:|
| Mf | 0.9774612033 M |
| Jf | 0.6480870733 M² |
| χ = Jf / Mf² | 0.6783194994 |
| Kerr reference frequency | 0.5353401319 M⁻¹ |
| Damping rate | 0.0834505970 M⁻¹ |

In the third period at radii 15, 30 and 40M, a free complex single-mode fit differs from this reference by 0.0034–0.1072% in frequency and 0.72–1.51% in damping. The original comparison targets at these three radii have been met. Earlier periods at more distant radii are not declared to share this asymptotic regime merely because their extraction radius is larger.

BIC comparisons in historical reports remain uncalibrated descriptive measures: there is no independent error model that would convert them into statistical confidence in accepting Kerr. Numerical convergence with resolution was not checked with a separate simulation.

## Geometry and physical area measure

B3/B4 reconstruct the intrinsic metric and physical weights from the saved QC0 tetrad. Areas of all 55 common surfaces were compared with QLM from the same run. The maximum relative difference is 6.0364×10⁻⁸. This is a reconstruction consistency check, not the accuracy of the physical equations or evolution.

Across the saved frames, area increases from 37.30528 to 41.01768 M². The Re Ψ₂ selection is predominantly longitudinal; the Im Ψ₂ and |σ| selections are polar under the operational regions of the original analysis. Boundaries evolve and individual surfaces have gaps; the viewer shows these explicitly.

## Final stage HF4-B5

B5 did not repeat the evolution, area reconstruction or old-mask selection. Saved fields and weights were used. At each time, independent weighted fits used coordinate spherical harmonics with 0≤ℓ≤L, L=2,4,6,8. The largest model has 81 coefficients on 2,520 nodes.

| Field | Fraction of angular variation described at L≤8 | Last frame |
|---|---:|---:|
| Re Ψ₂ | 96.388–99.977% | 99.97664% |
| Im Ψ₂ | 99.074–99.962% | 99.96240% |
| Magnitude of σ | 99.815–99.936% | 99.81524% |

These values show that the main angular pattern admits a compact description. This is fitting to the same data; the quality of a nested model naturally improves as its dimension increases. It is neither an independent prediction of the dynamics nor a unique physical explanation of the structure.

The residual was examined separately. In the first Re Ψ₂ frame, two nodes at θ=90°, φ=40° and 220° account for 83.14% of its squared norm and occupy approximately 0.281% of discrete physical area. The nodes were not excluded. Their origin remains undetermined; a separate numerical check is needed before assigning physical meaning to the feature. The report includes absolute and relative RMS, mask overlap and all times.

B5 coefficients are not Yitian’s geometric multipoles or intrinsic horizon eigenmodes. The |σ| analysis describes the magnitude without complex phase. The squared residual norm is not physical energy.

## Completed scope and open questions

**The reconstruction and analysis release for the available data is complete.** The B5 input, numerical results, programs, checks, visualisations and exact checkpoint are retained. No new heavy evolution is needed to complete this scope.

The following remain beyond the established result:

- A single continuous validated QC0 geometry series from 0 to 105M.
- Convergence of local features with resolution and independence from coordinates and tetrad.
- Proven identity of the selected HF1, QC0 and Yitian regions.
- Complete modelling of Yitian’s early dynamics and unique identification of nonlinear components, not achieved in earlier A3/A4 stages.
- A global event horizon, optical image or complete movie of the numerical spacetime.

Open questions are not hidden on publication or rewritten as negative findings. The project has not established an additional physical “hyperform”; a single B5 test also cannot prove its absence. Positive findings concern numerical geometry, field localisation, its dynamics and descriptive agreement with the accepted ringdown picture.

## Navigating the results

Key final sources within the release:

- HF1: [BH_HF1_V2_report.md](BH_HF1_V2_report.md), original HF1/HF1R/HF1RC protocols.
- QC0 waves: [BH_HF3_R7B_report.md](BH_HF3_R7B_report.md) and its JSON summary.
- QC0 geometry: [BH_HF4_B4_report.md](BH_HF4_B4_report.md), saved NPZ fields and weights.
- QC0 angular structure: [BH_HF4_B5_report.md](BH_HF4_B5_report.md), 660 metric rows, coefficients, residuals and conditioning checks.
- Yitian: A1–A4 and HF4-A1/A2 reports with explicit limits on early reconstruction.

Early reports preserve the research history. The corrected R7B reference and this final report take precedence over earlier statements about mass, spin or “scientific closure” of the entire physical problem. The release date does not imply external peer review.
