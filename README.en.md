# Prostranstva — horizon reconstruction

[Русский](README.md) · **English**

Apparent-horizon surfaces, scalar-field maps, ringdown checks and angular-structure analysis from three separately identified numerical configurations. Research release **1.0**, through **HF4-B5**.

[Website](https://olehtymofiienko.github.io/prostranstva-horizon-reconstruction/en/) · [Final report](docs/en/PROJECT_REPORT.md) · [Reproducibility](docs/en/REPRODUCIBILITY.md) · [Sources](docs/en/SOURCES.md)

## Explore

- **HF1 / SXS:BBH:0305, Lev6:** actual common-horizon coordinate surfaces, with regions selected by the preserved angular-gradient criterion for R*, E* and B*. The viewer contains 235 selected frames; 1,566 snapshots and 4,698 mask matches were checked.
- **QC0 / Einstein Toolkit:** individual and common apparent horizons through 22M; 55 common-horizon frames from 18.625M to 22M, with Re Ψ₂, Im Ψ₂ and |σ|. The separate wave evolution reaches 105M.
- **Yitian / SXS:BBH:0389, SpEC:** scalar angular maps from mass and spin moments through ℓ=8; a 301-frame atlas over 0–150M. These maps do not uniquely determine the early 3D horizon shape.
- **HF4-B5:** 660 weighted angular fits using saved QC0 fields and physical-area weights. At ℓ≤8, the principal angular pattern is described compactly; the first Re Ψ₂ residual is strongly concentrated on two retained nodes.

## Interpretation

The selected orange/yellow regions are operational field masks, not established independent objects. R* is intrinsic surface curvature; E* and B* refer to gravitational Weyl quantities, not electromagnetic fields. All are defined on the same horizon surface.

The project does not establish new physics. It does not reconstruct a global event horizon or a continuous validated QC0 geometry series from 0 to 105M. Resolution, gauge and tetrad robustness of local residuals remain open. The final report preserves these limits and supersedes historical parameter interpretations.

## Reading options

The website provides Russian/English navigation and light, dark, warm and system themes. Theme selection is remembered in your browser. Scientific colour scales are preserved. Original figure/video labels and archived program comments remain in their source language.

## Repository

`dist/` is the self-contained static website; `docs/` holds release documentation; `research/` preserves stage artifacts; `scripts/` contains packaging and validation utilities. Scientific results are retained without rerunning confirmed stages. Large external source grids and the original author-supplied multipole archive are not redistributed; see the reproduction table for each route.

```bash
python scripts/check_release.py
```

Read `CHECKPOINT.json` before continuing research. Scientific release completion applies to the available data and registered HF4-B5 scope, not to every open physical question. It does not imply external peer review.

## Acknowledgements

I am grateful to **Yitian** for kindly providing the time-series data of the
common horizon's mass and spin multipole moments through ℓ=8, together with
plotting scripts and detailed explanations of the data. These materials helped
clarify the definitions and conventions, check the analysis, and study the
geometry and evolution of the common horizon. I sincerely appreciate his help
and the time he devoted to preparing and sharing these data.
