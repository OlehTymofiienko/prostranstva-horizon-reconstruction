# HF1 — selected regions on the original horizon surface

**BH-HF1-VISUAL-V2 complete.** HF1/HF1R masks are placed on original common apparent-horizon C coordinates from SXS:BBH:0305, Lev6, with a one-minute video and interactive 3D scene.

## Data and matching

`BH_HF1_visual_data.zip` contains 1,566 snapshots: all 1,564 HF1R samples over 0–11.9070724M after the first saved common horizon, plus controls near 20/50M. The 10M control is already in the dense series. Total: 2,650,372 nodes. Coordinates, three fields and embedded reference-table checksums were verified. t₀=3685.496267868691M.

Original HF1R functions recover masks through field normalisation, eight nearest direction-sphere neighbours, gradient estimation, top 5% selection and largest group. Only missing node indices and comparison diagnostics were recovered. Permutation/significance/coherence/evolution tests were not repeated; criterion decisions and p-values came from saved results.

**All 4,698 comparisons agree.** Node counts match exactly; maximum centre-direction difference is 3.33×10⁻¹⁶ and gradient-threshold difference 1.78×10⁻¹⁵, consistent with rounding. These are the same components as the original tables.

## Meaning of the display

Surfaces use original `CoordsMeasurementFrame` Cartesian coordinates without displacement, radial exaggeration or smoothing. Connectivity is reconstructed on the direction sphere and transferred to the original coordinates, not via a Cartesian convex hull. All 1,566 meshes are closed: two faces per edge, Euler characteristic two.

Fields/coordinates share snapshot and node indices. Original `CoordCenterInertial` defines directions, without an assumed new frame transformation. Agreement reproduces the original criterion but does not establish coordinate independence.

The video shows R*. **Orange marks the component satisfying the saved criterion; light dots are its original nodes.** If the criterion fails, the largest computed group is grey. Colour is not temperature, energy, matter speed or time.

The boundary is the 1/2 contour of a linearly interpolated binary mask, with selected nodes 1 and others 0. This displays a discrete mask, not a new physical threshold. Faces are split without changing the surface; boundary accuracy is grid-limited.

Both video panels show **one horizon from opposite sides**, not two black holes. Coordinate scale is fixed. This is the common apparent horizon after formation, not a global event horizon or black-hole interior.

## Video timeline

| Time | Action |
|---|---|
| 00:00–00:36 | Surface/region evolve over 0–11.9071M; cameras fixed |
| 00:36–00:48 | Snapshot 2.0004M fixed; cameras orbit |
| 00:48–01:00 | Separate controls near 10,20,50M; gaps not invented |

Playback picks the last saved sample no later than requested time and displays its actual time. No temporal interpolation is used; not every dense sample receives a video frame. Grey does not prove field disappearance. Group motion is not a matter trajectory or proof of an independent object. No new-physics/“hyperform” claim is made.

## Interactive scene

The original `BH_HF1_V2_surface_3D.html` embeds data/library. Mouse/finger rotation, wheel zoom, R*/E*/B*, time, early playback, late jumps, opposite view and fixed-time orbit are available.

The viewer has 235 selected snapshots: uniform time samples, control times and states immediately before/after criterion-decision changes. All masks use original nodes; all 1,566 snapshots remain in stage data. Viewing-copy coordinates are rounded to eight decimals, with maximum error 5×10⁻⁹M. Original arrays/video are not rounded this way. QC0/Yitian data are not mixed into HF1.

## Validation and reproduction

Checks cover masks/meshes, all 1,200 video-frame times, fixed cameras during evolution and fixed time during orbit. Python/JavaScript boundary builders agree on triangles/selected faces. Splitting preserves coordinate area within 2.22×10⁻¹⁶ relatively; this is a technical geometry check, not physical area.

Initial/orbit/late video frames were inspected; control logic was tested separately. Independent browser WebGL rendering was not performed in that environment.

`BH_HF1_V2_results.zip` contains video, HTML, masks/triangles, checks, programs, reference results and checkpoint. The original coordinate archive remains separate; see `REPRODUCE.md`. This completes placing the original finding on its surface. At this historical stage, pending QC0 B5 was not run.
