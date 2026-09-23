# HF4-A1 — What can be reconstructed from horizon moments?

Historical stage: 16 September 2026, following HF3-R7B at 105M. English reading edition of the original report. Later stages and the final project report supersede its pending tasks.

**Result:** scalar maps were reconstructed and the late axisymmetric control passed. A finite set of moments does not uniquely determine the full horizon geometry. No spacetime evolution was repeated.

## Inputs and conventions

Yitian's file contains 5,001 times from 0 to 500M and 81 complex coefficients for each of the mass and spin moment families, through ℓ=8. It does not contain the horizon area, pointwise metric, surface grid or waveforms. Its SHA256 is `9cdd268e85b783d8a832b72e46945c6b43bb2b3ea0a9e6665e57a5a1b7d56f0f`.

QC0 is a separate simulation. Its corrected final parameters are Mf=0.977461203278627M, J=0.648087073319352M² and χ=J/Mf²=0.6783194994253878. The number 0.648087 is angular momentum, not dimensionless spin. QC0's 18 real mass/spin diagnostics are not Yitian's 81 complex angular coefficients.

The authors transport the angular basis with a horizon evolution vector X satisfying ℒ_Xε=(2Ṙ/R)ε and ℒ_XYℓm=0. The normalized area element is consequently preserved: dA/R²=dΩ, with R=√(A/4π). Their definitions are Iℓm=(1/4)∮ℛYℓm* dA and Lℓm=(1/2)∮εᵃᵇω_b D_aYℓm* dA. With these conventions the reconstructed scalar maps are:

- K=2Σ IℓmYℓm=R² times Gaussian curvature=R²ℛ/2;
- B=Σ LℓmYℓm=−R² curl(ω)/2.

Stored coefficients are complex-conjugated; column k=ℓ(ℓ+1)+m. This stage uses the original transported angles and retains the stationary coefficients. A frequency-analysis rotation is not silently added to a geometric map. On a dynamical horizon these maps cannot automatically be identified with the real and imaginary parts of Ψ₂. A coloured pattern is not a material substance.

## Numerical reconstruction checks

A 32-point Gauss–Legendre grid and 64 azimuthal points were used at t/M=0,5,10,20,40,80,150,500.

| Check | Maximum discrepancy |
|---|---:|
| Harmonic Gram matrix | 6.44×10⁻¹⁵ |
| Recovered I coefficients | 3.57×10⁻¹⁴ |
| Recovered L coefficients | 4.02×10⁻¹⁵ |
| Imaginary leakage in K | 6.12×10⁻¹⁴ |
| Mean K minus 1 | 1.55×10⁻¹⁰ |
| Mean B minus 0 | 3.62×10⁻¹⁸ |

The monopole identity was verified to 1.11×10⁻¹⁶. These are arithmetic and representation checks, not estimates of physical simulation error.

| t/M | Nonaxisymmetric fraction of K | Relative difference: ℓ≤6 versus ℓ≤8 |
|---:|---:|---:|
| 0 | 56.7741% | 2.01577% |
| 5 | 29.0812% | 0.16116% |
| 10 | 16.9930% | 0.14635% |
| 20 | 6.9340% | 0.06097% |
| 40 | 1.2615% | 0.01445% |
| 80 | 0.04126% | 0.01044% |
| 150 | 0.000524% | 0.01043% |
| 500 | 0.000442% | 0.01043% |

The corresponding nonaxisymmetric B fraction is 10.31% at 0M, 18.56% at 5M and 0.02776% at 80M. It is not monotonic from first appearance.

At t=0, the ℓ≤8 curvature map reaches about −0.261; about 18.3% of its area is negative on the 96×192 evaluation grid. The ℓ≤6 minimum is −0.233. Unknown higher moments prevent a rigorous truncation error bound. Negative Gaussian curvature indicates local saddle behaviour; it does not establish a hole or a change in topology.

## Conditional late axisymmetric geometry

For an axisymmetric metric, with x=cosθ,

`ds²/R² = dx²/f(x) + f(x)dφ²`, and `K=−f''/2`.

Integrating with f(−1)=0 and f′(−1)=2 gives a regular second pole if f(1)=0 and f′(1)=−2. The late 400–500M mean was used, retaining m=0. The omitted nonaxisymmetric fraction was 4.42×10⁻⁶. Tiny monopole/dipole corrections, −3.35×10⁻¹⁴ and +9.45×10⁻¹⁴, enforce regularity only in this conditional model; the original maps are unchanged.

For Kerr, f=(1+b²)(1−x²)/(1+b²x²), where b=χ/(1+√(1−χ²)). The spin inferred from L10 is χ=0.68644097081. This is an internal consistency test using the same simulation.

| Comparison with Kerr | Discrepancy |
|---|---:|
| Maximum absolute f difference | 1.27×10⁻⁷ |
| Maximum relative f difference | 2.30×10⁻⁷ |
| Maximum K difference | 2.54×10⁻⁵ |
| Embedding profile difference | 3.51×10⁻⁷R |

The equatorial and polar embedding radii are approximately 1.0760823R and 0.8390690R. R≈1.76868M follows conditionally from the final Kerr parameters; it was not measured from an area series in this file.

A counterexample demonstrates nonuniqueness. Adding 0.2P10 to K changes the supplied moments through ℓ=8 by at most 1.20×10⁻¹⁵, but changes maximum curvature by 0.2, f by 9.11×10⁻⁴ and the embedding by 0.00248R. The area, regular poles and supplied L coefficients are retained. This is a mathematical construction, not a detected ℓ=10 mode or a complete physical spacetime.

## QC0 inventory and its limits

The R7A archive contains 647 files and passes ZIP CRC checks. Its `sf_radius` series has 161 rows from 95 to 105M with only zero grid indices: these are reductions, not a surface mesh. The 18 real moments are M0…M8 and J0…J8. The area relation A=4πR² agrees to about 6×10⁻¹⁵; M0 differs from the mass diagnostic by about 0.0236%, and J1 from the spin diagnostic by about 0.0000345%.

At the inspected Einstein Toolkit source revision `edd7fe658cb6748cd9ac8efcea314214b87bc043`, QLM multipoles use powers of R and Legendre polynomials of an invariant axial coordinate, rather than the complete Yℓm basis. The exact installed binary revision was not established at this stage. Disabled regular grid output does not imply that recoverable geometry is absent from checkpoints; those were not inspected here.

## Route recorded at this stage

A2 was to build an atlas over 0–150M with a 400–500M background. A physical evolving surface additionally requires area, metric/grid information and the relation to transported angles. The finite scalar moment series alone does not specify that surface. No new user-side simulation was required by A1.

For reproduction, the archived program accepts `--moments`, `--qc0-archive` and `--out`. Original private inputs must be supplied separately where they are excluded from the public release. The recorded environment was Python 3.12.14, NumPy 2.3.5, SciPy 1.17, h5py 3.16 and Matplotlib 3.10.8.

Methodological sources: [Chen et al.](https://arxiv.org/abs/2208.02965), [Ashtekar et al.](https://arxiv.org/abs/gr-qc/0401114), and the Einstein Toolkit source. See the original report and archived outputs for the full numerical records.
