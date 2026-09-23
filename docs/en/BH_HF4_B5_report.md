# HF4-B5 — angular structure and residual localisation

## Input and method

Fifty-five QC0 common-horizon frames cover 18.625–22M at 0.0625M spacing, each with 2,520 nodes (35×72). Fields, physical-area weights and original masks come from saved B4 data. Evolution, area, masks, HF3 ringdown and HF1/Yitian were not recalculated.

Input SHA-256: `a5e1c19802ded7f79da5475eb267f1228c6e26208f69c88337346faf617ee8e8`. Preserved protocol SHA-256: `2140113cbd0284e39ed23b30225343160ce150060d6b7c381f2bcb3d96d85569`.

Independent weighted fits at each time use 0≤ℓ≤L, L=2,4,6,8, with 9,25,49,81 real coefficients respectively: 660 field fits, or 220 matrices with three right-hand sides. The basis is Yℓ0, √2 Re Yℓm and √2 Im Yℓm for m>0, with standard unit-sphere normalisation and Condon–Shortley phase, implemented by `scipy.special.sph_harm_y(ell,m,theta,phi)`. Original QC0 polar θ, azimuth φ, axes and indices are preserved.

Minimise Σwᵢ(Fᵢ−Fᴸᵢ)². Define R²=1−Σwᵢ(Fᵢ−Fᴸᵢ)²/Σwᵢ(Fᵢ−⟨F⟩w)². This is explained **weighted angular variation on the fitted nodes**, not energy, probability, simulation accuracy or independent predictive validation. Nested larger spaces cannot worsen the fit. No denominator is degenerate.

Coordinate-sphere harmonics are not computed intrinsic-Laplacian eigenfunctions of the deformed horizon. These coefficients are not Yitian’s geometric multipoles. For |σ| only the real magnitude is expanded, without complex phase or a spin-weighted decomposition.

## Explained variation over 55 frames

| Field | L≤2 | L≤4 | L≤6 | L≤8 |
|---|---:|---:|---:|---:|
| Re Ψ₂ | 88.208–91.337% | 94.434–99.208% | 95.754–99.904% | 96.388–99.977% |
| Im Ψ₂ | 72.330–92.678% | 83.501–98.833% | 95.765–99.829% | 99.074–99.962% |
| Magnitude of σ | 93.902–97.500% | 98.098–99.659% | 99.520–99.854% | 99.815–99.936% |

Ranges cover times, not confidence intervals. L=2 already describes much of the broad pattern, but performance varies with field/time.

| Last frame, t=22M | R² at L≤8 | Residual RMS / original standard deviation |
|---|---:|---:|
| Re Ψ₂ | 99.97664% | 1.528% |
| Im Ψ₂ | 99.96240% | 1.939% |
| Magnitude of σ | 99.81524% | 4.298% |

Relative |σ| residual increases slightly near the end, although absolute residual RMS falls from 0.00577094 to 0.00227902. All fields do not improve monotonically.

## Residual location

The unchanged B4 mask selects roughly 5% of physical area with greatest deviation from the weighted median. The residual mask applies the same rule to r=F−Fᴸ: the weighted 95th percentile of |r−median_w(r)|, including ties. This was fixed in `BH_HF4_B5_definition.json` before fitting.

Weighted Dice overlap and the old-mask fraction of Σwr² were measured. **Σwr² is a numerical norm, not physical energy.**

| Field at L≤8 | Old/residual mask Dice | Squared residual norm inside old mask |
|---|---:|---:|
| Re Ψ₂ | 0.0000–0.2717 | 1.592–83.476% |
| Im Ψ₂ | 0.0000–0.2123 | 1.087–9.268% |
| Magnitude of σ | 0.0000–0.0016 | 1.018–1.848% |

Im Ψ₂ and |σ| residual maxima usually lie outside the old region. Re Ψ₂ varies over time; residuals cannot all be called one persistent localised entity. Zero upper-tail mask overlap does not mean zero residual in the old region.

**First frame:** at t=18.625M, Re Ψ₂ has R²=96.38844%, residual RMS=19.004% of original standard deviation, and 83.476% of squared residual norm in the old mask. Two nodes at θ=90°, φ=40° and 220° contribute **83.144%** of total Σwr² while occupying **0.2807%** of discrete physical area.

High explained variation does not eliminate large local deviations. Dependence on two nodes requires checking original diagnostics/resolution before physical interpretation. Nodes were not removed, smoothed or excluded for refitting. This extra localisation measure is descriptive, without a statistical detection threshold.

`BH_HF4_B5_maps_it9536.png` and `BH_HF4_B5_maps_it11264.png` show field minus mean, approximation minus the same mean, and residual with one scale across each row. Orange contours are B4 masks. These are angular maps, not a new shape reconstruction. All other times are saved in NPZ.

## Numerical checks

All weights are positive; maximum sum-to-one error is 2.22e−16. All 220 matrices have full column rank, with condition numbers 1.405279–3.135913. Analytical Y00, Y10 and addition-theorem error is ≤6.88e−15. Maximum relative normal-equation residual is 6.97e−15; weighted residual mean is consistent with zero. Independent SVD/QR endpoint predictions at L=8 differ by ≤7.55e−15 relative to the original-value scale. Nested fits do not increase squared residual. Input checksum remains unchanged.

These checks validate linear processing, not simulation convergence.

## Outcome and limits

B5 answers the planned question: how much variation selected angular spaces describe and where the residual lies relative to old masks. The main pattern is compactly represented. A nonzero residual, especially early Re Ψ₂, remains; it alone establishes neither new physics nor an independent “latent object”.

Release 1.0 completes reconstruction and descriptive validation of available HF1, QC0/HF3/HF4 and Yitian data. SXS:BBH:0305, QC0 and SXS:BBH:0389 remain separate configurations.

Open limits: no continuous validated QC0 geometry over 0–105M; no second-resolution/coordinate/independent-tetrad test of local residuals; no established origin of the first-frame nodal feature; no proven physical identity of HF1 gradient regions and QC0 deviation regions; no unique full geometry from finite Yitian moments; early Yitian dynamics/nonlinear identification remain open. Displays show apparent horizons or labelled scalar maps, not event horizons or ray-traced optics.

## Reproduction

`BH_HF4_B5.py` performs B5 only; `finish_B5.py` packages saved results. The archive includes B4 input NPZ, geometry table and protocol; README gives commands. Routine continuation needs no rerun.

Tables: `BH_HF4_B5_metrics.csv`, `BH_HF4_B5_coefficients.csv`, `BH_HF4_B5_design_checks.csv`, `BH_HF4_B5_singular_values.csv`, `BH_HF4_B5_residual_node_concentration.csv`. Arrays: `BH_HF4_B5_fits_residuals.npz`. QC0 times/units are retained. [SciPy convention](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.sph_harm_y.html); normalisation was also checked analytically.
