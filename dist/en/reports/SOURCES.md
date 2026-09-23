# Sources and separation of datasets

## Numerical datasets

- **SXS:BBH:0305 / Lev6:** common-horizon surface C, original HF1 branch. Arrays previously obtained by the user and derived masks are used. [SXS catalogue](https://data.black-holes.org/).
- **QC0 / Einstein Toolkit:** the user’s early geometry and late waveform runs. Provenance is tied to parameters, time series and checksums, not just the name QC0. [Einstein Toolkit](https://einsteintoolkit.org/).
- **SXS:BBH:0389 / SpEC:** mass/spin moments provided by Yitian. Mf=0.95162M, χ=0.68644, time starts at first common horizon. [Authors’ paper](https://arxiv.org/abs/2208.02965). Derived results/programs are included; original correspondence and author-supplied data are not redistributed.

Equal initial conditions across datasets have not been established. Their frames are not spliced into one physical process.

## Source checksums

| Source | SHA-256 |
|---|---|
| HF1 `BH_HF1_visual_data.zip` | b47bdb37aa796d281ba292c3908ed334c856360ec1e6e324c0a7406a175b0bb2 |
| QC0 `BH_HF3_R7A_output.zip` | d73d7bce0b88c22022ccc5e571431045370edccbf48ed9b03f86fb261afafa1b |
| QC0 `BH_HF4_B4_collection(1).zip` | 34d86aa53bc74f0e1b872761518a145b7612aafbb2f7e735e6993961140936fc |
| B5 `BH_HF4_B4_fields_geometry_masks.npz` | a5e1c19802ded7f79da5475eb267f1228c6e26208f69c88337346faf617ee8e8 |

Historical manifests retain detailed inventories. `SHA256SUMS.json` covers this release.

## Conventions and software

Electric/magnetic Weyl parts and tidal visualisation: [1208.3034](https://arxiv.org/abs/1208.3034), [1012.4869](https://arxiv.org/abs/1012.4869). Apparent horizons/search: [gr-qc/0306056](https://arxiv.org/abs/gr-qc/0306056). Spherical harmonics: [SciPy sph_harm_y](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.sph_harm_y.html). These establish concept provenance; B5 numbers come from saved project data and report formulas.

Plotly.js 3.3.1 provides interactive viewing, with its licence notice in `dist/assets/plotly.min.js`. Marked 17.0.5 renders Markdown, with its notice in `scripts/vendor/`. NumPy, SciPy and Matplotlib perform analysis/plotting; run versions are in B5 definitions.

Cite project title, release, stage and corresponding data source. These results are not a publication by SXS, Einstein Toolkit or the authors of Yitian’s original dataset.
