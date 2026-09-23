# HF4-A3 — Inventory of the available geometric data

Historical stage: 16 September 2026. English reading edition. The inventory requested here was subsequently received and analysed in B1; this report records the earlier checkpoint.

**Status at this stage:** available-data audit complete; local geometry inventory pending. A1/A2 and the numerical evolutions were not repeated.

## What was available

Yitian's series supplies 5,001 times and 81 complex coefficients per moment family. It contains no pointwise metric, area history, surface coordinates or complete transport grid. QC0 R7A contains 647 files with waveforms and horizon diagnostics, but no verified pointwise geometry in that archive.

Its `sf_radius` output has 161 times from 95 to 105M and one value per each of six surface slots. Grid-index columns 6–8 are all zero and coordinate columns 10–12 are all zero. Those reductions do not place a physical surface at the coordinate origin. Metric reductions are not a surface metric grid. Checkpoints were known to exist on the user's computer but had not been inspected here.

Regular HDF5 field output was disabled and its variable list empty; AH surface output and one-/two-dimensional ASCII output were also disabled. Checkpoint writing was enabled. These settings motivate inspecting existing checkpoint contents before proposing another evolution.

## What a checkpoint might preserve

At the inspected Einstein Toolkit revision `edd7fe658cb6748cd9ac8efcea314214b87bc043`, SphericalSurface metadata includes origins, spacing, dimensions, ghost zones and validity. Some QLM metric quantities explicitly disable checkpointing, whereas surface coordinates and tetrad-related quantities may be available. Installed source and actual datasets still require verification.

A conditional route uses a spatial complex dyad m=(u+iv)/√2 with physically orthonormal u,v and vanishing time component. One can solve for the corresponding metric on the tangent span. Four hundred synthetic metric cases gave maximum relative discrepancy 7.47×10⁻¹⁴; a phase-change control gave 5.61×10⁻¹⁴. These are algebraic controls, not tests on a recovered physical grid.

Before interpreting real arrays, one must check time, surface validity, MPI copies, indexing, ghost zones, tangent vectors, positive definiteness and agreement with an independent horizon area.

## External geometry is not interchangeable

An SXS `Horizons.h5` file may contain centres, irreducible masses and spins, allowing A=16πMirr². Some simulations also supply surface dumps. Their existence for SXS:BBH:0389, numerical resolution and correspondence with Yitian's transported angles were not established here. Such data cannot simply be combined with a different simulation or angular convention. No email was sent at this stage.

Coordinate shape, intrinsic geometry and a Euclidean embedding are different objects. Drawing a radial relief from curvature without an established metric reconstruction would misrepresent the data.

## Read-only collector

The prepared collector searches likely QC0 simulation directories and opens HDF5 files read-only. It records metadata and eligible small surface arrays, capped at 2 MiB per array and 96 MiB overall. It does not read the full volume fields, copy whole checkpoints, change the simulation or infer measured time solely from an iteration in a filename.

The collector requires NumPy and h5py for extraction; otherwise it marks the inventory as incomplete. Fourteen synthetic controls covered corrupt files, external HDF5 links, bounds and unchanged inputs. Existing output archives are preserved with timestamps.

The historical user command was to activate `~/gw_env` and run `BH_HF4_A3_collect.py` from Downloads, producing `C:\Users\Admin\BH_HF4_A3_inventory.zip`. This is a record of the completed route, not a request to rerun it now.

| Inventory result | Permitted next step |
|---|---|
| Valid coordinates and sufficient dyad information | Recover the surface metric and compare the area |
| Coordinates alone | Show the coordinate surface with that limitation |
| Only volume/checkpoint data | Plan extraction from existing files before any new evolution |
| Missing early frames | Record the gap explicitly |

The A2 scalar atlas remains valid within its stated scope. The next historical action was to inspect the returned inventory; B1 and the final project report give the subsequent outcome.
