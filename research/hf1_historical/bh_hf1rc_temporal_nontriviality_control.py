#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
BH-HF1R-C — Temporal Nontriviality Control
==========================================

Context
-------
BH-HF1R formally passed on SXS:BBH:0305 / Lev6 common horizon C, but the
native sampling is extremely dense (median cadence ~0.0018 M).  Therefore
adjacent-snapshot Dice≈1 can arise simply because a smooth numerical field
changes negligibly between almost identical times.

HF1R-C asks a stricter question:

    Does the localized morphology remain coherent across physically larger
    fixed lags (0.25, 0.5, 1.0 M), and is that coherence stronger than
    two null models?

Null A: independent random SO(3) rotation of the second component
        (preserves component shape/size, destroys directional alignment).

Null B: distant-time shuffle within the same field
        (preserves actual morphology and coordinate frame, destroys the
        requested local temporal relation).

This stage is still coordinate/gauge dependent.  It is NOT a physical phase
test, NOT a domain-wall detection, NOT QNM-subtracted, and NOT a beyond-GR test.

Inputs
------
1) Frozen HF1R snapshot metrics CSV.  Its candidate_snapshot column is reused;
   HF1R-C does NOT re-run the 199 spatial shuffles.
2) The same SXS horizon HDF5 products used in HF1/HF1R.
3) The already-audited H2R2 reader script.

Frozen lags
-----------
    0.25 M, 0.50 M, 1.00 M

Frozen matching tolerance
-------------------------
For target t+lag, use the nearest morphology-positive snapshot if:
    |(t_j - t_i) - lag| <= max(0.03 M, 0.15*lag)

Native components are mapped to a fixed 4096-point Fibonacci sphere.
Similarity uses one-hop-dilated fixed-grid Dice.

Primary lag-coherence rule, per field/lag
-----------------------------------------
Require all:
    matched_pairs >= 30
    median real dilated Dice >= 0.35
    median centroid step <= 25 deg
    empirical p_vs_SO3 <= 0.01
    real median Dice - SO3-null q99 >= 0.10

A field passes lagged coherence if >=2 of 3 lags pass.
Global lagged coherence requires >=2 of 3 fields.

Time-specificity rule, per field/lag
------------------------------------
Require:
    empirical p_vs_distant_time <= 0.05
    real median Dice - distant-time-null q95 >= 0.05

A field is time-specific if >=2 of 3 lags pass.
Global time-specificity requires >=2 of 3 fields.

Possible classifications
------------------------
BH_HF1RC_LAGGED_COHERENCE_TIME_SPECIFIC
BH_HF1RC_LAGGED_COHERENCE_NOT_TIME_SPECIFIC
BH_HF1RC_NO_LAGGED_COHERENCE

Interpretation
--------------
* TIME_SPECIFIC:
  fixed-lag morphology is stronger than both nulls; proceed to BH-HF2
  QNM/multipole closure.

* NOT_TIME_SPECIFIC:
  morphology survives finite lags and SO(3) randomization, but is not more
  similar at nearby times than at distant times.  This favors a persistent /
  slowly varying coordinate multipolar pattern rather than a distinct transient
  object.  Proceed, if desired, only to a simplified multipole-closure control.

* NO_LAGGED_COHERENCE:
  the dense HF1R pass does not survive finite-lag control.  Stop the transient-
  morphology route.

Outputs
-------
BH_HF1RC_pair_metrics.csv
BH_HF1RC_null_trials.csv
BH_HF1RC_null_summary.csv
BH_HF1RC_field_summary.csv
BH_HF1RC_summary.json
BH_HF1RC_preregistered_protocol.json
BH_HF1RC_lagged_dice.png
BH_HF1RC_centroid_step.png
"""

from __future__ import annotations

import argparse
import importlib.util
import json
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.spatial import cKDTree


STAGE = "BH-HF1R-C"
SID = "SXS:BBH:0305"
LEV = 6
HORIZON = "C"

LAGS_M = [0.25, 0.50, 1.00]
MATCH_TOL_FLOOR_M = 0.03
MATCH_TOL_FRAC = 0.15

KNN_NATIVE = 8
ACTIVE_Q = 0.95

GRID_N = 4096
GRID_KNN = 8

MIN_MATCHED_PAIRS = 30
REAL_DICE_MIN = 0.35
REAL_CENTROID_STEP_MAX_DEG = 25.0
SO3_P_MAX = 0.01
SO3_MARGIN_Q99_MIN = 0.10

TIME_P_MAX = 0.05
TIME_MARGIN_Q95_MIN = 0.05
DISTANT_TIME_MIN_SEPARATION_M = 2.0

MIN_LAGS_PER_FIELD = 2
MIN_FIELDS_GLOBAL = 2

N_NULL_TRIALS = 199
NULL_PAIR_CAP = 120
RNG_SEED = 150914

FIELDS = ["Rstar", "Estar", "Bstar"]


def import_by_path(path, name):
    spec = importlib.util.spec_from_file_location(name, str(path))
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot import {path}")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def robust_z(x):
    x = np.asarray(x, float)
    med = float(np.median(x))
    mad = float(np.median(np.abs(x - med)))
    scale = 1.4826 * mad
    if not np.isfinite(scale) or scale <= 1e-14:
        scale = float(np.std(x))
    if not np.isfinite(scale) or scale <= 1e-14:
        scale = 1.0
    return (x - med) / scale


def unit_directions(coords, center):
    x = np.asarray(coords, float) - np.asarray(center, float)[None, :]
    r = np.linalg.norm(x, axis=1)
    good = np.isfinite(r) & (r > 0)
    if good.sum() < max(20, int(0.9 * len(r))):
        raise RuntimeError("Too many invalid surface radii")
    return x[good] / r[good, None], good


def knn_graph_on_sphere(n, k):
    tree = cKDTree(n)
    _, idx = tree.query(n, k=min(k + 1, len(n)))
    neigh = idx[:, 1:]
    dot = np.einsum("ij,ikj->ik", n, n[neigh])
    dot = np.clip(dot, -1.0, 1.0)
    ang = np.arccos(dot)
    return neigh, ang


def gradient_proxy(z, neigh, ang):
    return np.median(
        np.abs(z[neigh] - z[:, None]) / np.maximum(ang, 1e-8),
        axis=1,
    )


def largest_active_component(active, neigh):
    active = np.asarray(active, bool)
    ids = np.flatnonzero(active)
    if len(ids) == 0:
        return np.array([], dtype=int)

    aset = set(map(int, ids))
    visited = set()
    best = []

    for s0 in ids:
        s0 = int(s0)
        if s0 in visited:
            continue
        stack = [s0]
        visited.add(s0)
        comp = []
        while stack:
            u = stack.pop()
            comp.append(u)
            for v in neigh[u]:
                v = int(v)
                if v in aset and v not in visited:
                    visited.add(v)
                    stack.append(v)
        if len(comp) > len(best):
            best = comp

    return np.asarray(best, dtype=int)


def component_from_field(values, n, neigh, ang):
    z = robust_z(values)
    g = gradient_proxy(z, neigh, ang)
    thr = float(np.quantile(g, ACTIVE_Q))
    active = g >= thr
    idx = largest_active_component(active, neigh)

    if len(idx):
        v = np.sum(n[idx], axis=0)
        vn = np.linalg.norm(v)
        centroid = v / vn if vn > 0 else np.array([np.nan]*3)
    else:
        centroid = np.array([np.nan]*3)

    return idx, centroid


def fibonacci_sphere(n):
    i = np.arange(n, dtype=float)
    phi = (1.0 + np.sqrt(5.0)) / 2.0
    z = 1.0 - 2.0 * (i + 0.5) / n
    theta = 2.0 * np.pi * i / phi
    r = np.sqrt(np.maximum(0.0, 1.0 - z*z))
    return np.column_stack([r*np.cos(theta), r*np.sin(theta), z])


def component_to_grid_mask(component_dirs, grid_tree, grid_n):
    mask = np.zeros(grid_n, dtype=bool)
    if len(component_dirs):
        _, ii = grid_tree.query(component_dirs, k=1)
        mask[np.asarray(ii, dtype=int)] = True
    return mask


def dilate_mask(mask, grid_neigh):
    out = mask.copy()
    ids = np.flatnonzero(mask)
    if len(ids):
        out[np.unique(grid_neigh[ids].ravel())] = True
    return out


def dice(a, b):
    a = np.asarray(a, bool)
    b = np.asarray(b, bool)
    den = int(a.sum()) + int(b.sum())
    if den == 0:
        return 0.0
    return float(2 * np.logical_and(a, b).sum() / den)


def angle_deg(a, b):
    a = np.asarray(a, float)
    b = np.asarray(b, float)
    if not np.all(np.isfinite(a)) or not np.all(np.isfinite(b)):
        return np.nan
    c = float(np.clip(np.dot(a, b), -1.0, 1.0))
    return float(np.degrees(np.arccos(c)))


def random_rotation_matrix(rng):
    # Shoemake uniform random quaternion.
    u1, u2, u3 = rng.random(3)
    q1 = np.sqrt(1-u1) * np.sin(2*np.pi*u2)
    q2 = np.sqrt(1-u1) * np.cos(2*np.pi*u2)
    q3 = np.sqrt(u1) * np.sin(2*np.pi*u3)
    q4 = np.sqrt(u1) * np.cos(2*np.pi*u3)
    x, y, z, w = q1, q2, q3, q4
    return np.array([
        [1-2*(y*y+z*z), 2*(x*y-z*w),   2*(x*z+y*w)],
        [2*(x*y+z*w),   1-2*(x*x+z*z), 2*(y*z-x*w)],
        [2*(x*z-y*w),   2*(y*z+x*w),   1-2*(x*x+y*y)],
    ], dtype=float)


def rotate_raw_mask(raw_mask, grid, grid_tree, rng):
    ids = np.flatnonzero(raw_mask)
    out = np.zeros(len(grid), dtype=bool)
    if len(ids) == 0:
        return out
    R = random_rotation_matrix(rng)
    pts = grid[ids] @ R.T
    _, jj = grid_tree.query(pts, k=1)
    out[np.asarray(jj, dtype=int)] = True
    return out


def catalog_center(h2, hh, horizon, tq):
    t, y = h2.h5_series(hh, horizon, "CoordCenterInertial")
    return h2.interp_cols(t, y[:, :3], [tq])[0]


def mch_at_time(h2, hh, horizon, tq):
    t, y = h2.h5_series(hh, horizon, "ChristodoulouMass")
    vals = y if y.ndim == 1 else y[:, 0]
    return float(np.interp(tq, t, vals))


def field_values(field, R, E, B, mch):
    if field == "Rstar":
        return np.asarray(R, float)
    if field == "Estar":
        return (mch*mch) * np.asarray(E, float)
    if field == "Bstar":
        return (mch*mch) * np.asarray(B, float)
    raise KeyError(field)


def empirical_upper_p(real_value, null_values):
    null_values = np.asarray(null_values, float)
    return float((1 + np.sum(null_values >= real_value - 1e-15)) / (len(null_values) + 1))


def uniformly_cap_records(records, cap):
    if len(records) <= cap:
        return records
    ii = np.linspace(0, len(records)-1, cap).round().astype(int)
    return [records[i] for i in ii]


def build_lag_pairs(field_df, lag):
    g = field_df[field_df["candidate_snapshot"]].sort_values("actual_time_M").copy()
    times = g["actual_time_M"].to_numpy()
    idxs = g["snapshot_index"].astype(int).to_numpy()
    offs = g["offset_M"].to_numpy()

    tol = max(MATCH_TOL_FLOOR_M, MATCH_TOL_FRAC * lag)
    pairs = []
    used = set()

    for pos, (i, t, off) in enumerate(zip(idxs, times, offs)):
        target = t + lag
        jpos = int(np.searchsorted(times, target))
        candidates = []
        for q in [jpos-1, jpos, jpos+1]:
            if 0 <= q < len(times):
                dt_err = abs((times[q] - t) - lag)
                candidates.append((dt_err, q))
        if not candidates:
            continue
        dt_err, q = min(candidates)
        if dt_err > tol:
            continue
        j = int(idxs[q])
        if j == i:
            continue
        key = (i, j)
        if key in used:
            continue
        used.add(key)
        pairs.append({
            "snapshot_i": i,
            "snapshot_j": j,
            "offset_i_M": float(off),
            "offset_j_M": float(offs[q]),
            "actual_dt_M": float(times[q] - t),
            "lag_target_M": float(lag),
            "lag_error_M": float(dt_err),
        })
    return pairs


def make_plots(out, null_summary):
    ns = null_summary.copy()

    fig = plt.figure(figsize=(10.5, 6))
    for field in FIELDS:
        g = ns[ns["field"] == field].sort_values("lag_M")
        plt.plot(g["lag_M"], g["real_median_dilated_dice"], marker="o", label=f"{field}: real")
        plt.plot(g["lag_M"], g["so3_q99"], marker="x", linestyle="--", label=f"{field}: SO(3) q99")
        plt.plot(g["lag_M"], g["time_q95"], marker="s", linestyle=":", label=f"{field}: distant-time q95")
    plt.xlabel("fixed lag [M]")
    plt.ylabel("median one-hop-dilated Dice")
    plt.title("BH-HF1R-C: finite-lag morphology vs null controls")
    plt.ylim(-0.02, 1.02)
    plt.legend(ncol=2, fontsize=8)
    plt.tight_layout()
    fig.savefig(out / "BH_HF1RC_lagged_dice.png", dpi=180)
    plt.close(fig)

    fig = plt.figure(figsize=(9.5, 5.5))
    for field in FIELDS:
        g = ns[ns["field"] == field].sort_values("lag_M")
        plt.plot(g["lag_M"], g["real_median_centroid_step_deg"], marker="o", label=field)
    plt.axhline(REAL_CENTROID_STEP_MAX_DEG, linestyle="--", label="centroid-step gate=25°")
    plt.xlabel("fixed lag [M]")
    plt.ylabel("median centroid great-circle step [deg]")
    plt.title("BH-HF1R-C: finite-lag centroid displacement")
    plt.legend()
    plt.tight_layout()
    fig.savefig(out / "BH_HF1RC_centroid_step.png", dpi=180)
    plt.close(fig)


def run(args):
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    metrics = pd.read_csv(args.snapshot_metrics)
    required_cols = {
        "snapshot_index", "field", "actual_time_M", "offset_M", "candidate_snapshot"
    }
    missing = required_cols - set(metrics.columns)
    if missing:
        raise RuntimeError(f"Snapshot metrics missing columns: {sorted(missing)}")

    metrics["candidate_snapshot"] = metrics["candidate_snapshot"].astype(bool)

    h2 = import_by_path(Path(args.h2_script), "h2r2_frozen")
    rng = np.random.default_rng(RNG_SEED)

    grid = fibonacci_sphere(GRID_N)
    grid_tree = cKDTree(grid)
    grid_neigh, _ = knn_graph_on_sphere(grid, GRID_KNN)

    # Reconstruct only snapshots positive in at least one field.
    positive = metrics[metrics["candidate_snapshot"]].copy()
    pos_by_snapshot = {
        int(i): set(g["field"].tolist())
        for i, g in positive.groupby("snapshot_index")
    }

    components = {}
    unique_snapshots = sorted(pos_by_snapshot.keys())

    print(f"Reconstructing components for {len(unique_snapshots)} positive snapshots...", flush=True)

    with h5py.File(args.horizons, "r") as hh, h5py.File(args.dump, "r") as hd:
        for count, i in enumerate(unique_snapshots, start=1):
            # Use actual time from metrics, common for all fields at same index.
            trow = metrics[metrics["snapshot_index"] == i].iloc[0]
            tsnap = float(trow["actual_time_M"])

            R, exR = h2.read_scalar_snapshot(hd, HORIZON, "DimlessRicciScalar.tdm", int(i))
            E, exE = h2.read_scalar_snapshot(hd, HORIZON, "WeylE_NN.tdm", int(i))
            B, exB = h2.read_scalar_snapshot(hd, HORIZON, "WeylB_NN.tdm", int(i))
            C, exC = h2.read_coord_snapshot(hd, HORIZON, int(i))
            if not (exR == exE == exB == exC):
                raise RuntimeError(f"Extent mismatch at snapshot {i}")

            center = catalog_center(h2, hh, HORIZON, tsnap)
            mch = mch_at_time(h2, hh, HORIZON, tsnap)

            n, good = unit_directions(C, center)
            neigh, ang = knn_graph_on_sphere(n, KNN_NATIVE)

            R = np.asarray(R)[good]
            E = np.asarray(E)[good]
            B = np.asarray(B)[good]

            for field in pos_by_snapshot[i]:
                vals = field_values(field, R, E, B, mch)
                idx, centroid = component_from_field(vals, n, neigh, ang)
                raw_mask = component_to_grid_mask(n[idx], grid_tree, GRID_N)
                dil_mask = dilate_mask(raw_mask, grid_neigh)
                components[(field, int(i))] = {
                    "raw_mask": raw_mask,
                    "dil_mask": dil_mask,
                    "centroid": centroid,
                }

            if count % 100 == 0 or count == len(unique_snapshots):
                print(f"  reconstructed {count}/{len(unique_snapshots)}", flush=True)

    pair_rows = []
    null_trial_rows = []
    null_summary_rows = []

    # Candidate keys per field, for distant-time null.
    field_candidates = {}
    field_time = {}
    for field in FIELDS:
        g = positive[positive["field"] == field].sort_values("actual_time_M")
        field_candidates[field] = g["snapshot_index"].astype(int).tolist()
        field_time[field] = {
            int(r["snapshot_index"]): float(r["actual_time_M"])
            for _, r in g.iterrows()
        }

    for field in FIELDS:
        fdf = metrics[metrics["field"] == field].copy()

        for lag in LAGS_M:
            pairs = build_lag_pairs(fdf, lag)

            # Calculate real metrics for all matched pairs that have reconstructed masks.
            valid_pairs = []
            for p in pairs:
                ki = (field, int(p["snapshot_i"]))
                kj = (field, int(p["snapshot_j"]))
                if ki not in components or kj not in components:
                    continue
                a = components[ki]
                b = components[kj]
                dd = dice(a["dil_mask"], b["dil_mask"])
                cs = angle_deg(a["centroid"], b["centroid"])
                rec = {
                    "field": field,
                    **p,
                    "dilated_dice": dd,
                    "centroid_step_deg": cs,
                }
                pair_rows.append(rec)
                valid_pairs.append(rec)

            if len(valid_pairs) == 0:
                null_summary_rows.append({
                    "field": field,
                    "lag_M": lag,
                    "matched_pairs": 0,
                    "real_median_dilated_dice": np.nan,
                    "real_median_centroid_step_deg": np.nan,
                    "so3_p": np.nan,
                    "so3_q99": np.nan,
                    "so3_margin_vs_q99": np.nan,
                    "time_p": np.nan,
                    "time_q95": np.nan,
                    "time_margin_vs_q95": np.nan,
                    "lag_coherence_pass": False,
                    "time_specific_pass": False,
                })
                continue

            real_med = float(np.median([r["dilated_dice"] for r in valid_pairs]))
            real_step = float(np.median([r["centroid_step_deg"] for r in valid_pairs]))

            sample_pairs = uniformly_cap_records(valid_pairs, NULL_PAIR_CAP)

            so3_trial_medians = []
            time_trial_medians = []

            all_field_ids = field_candidates[field]
            all_times = field_time[field]

            for trial in range(N_NULL_TRIALS):
                so3_vals = []
                time_vals = []

                for p in sample_pairs:
                    ki = (field, int(p["snapshot_i"]))
                    kj = (field, int(p["snapshot_j"]))
                    a = components[ki]
                    b = components[kj]

                    # SO(3) null.
                    rb_raw = rotate_raw_mask(b["raw_mask"], grid, grid_tree, rng)
                    rb_dil = dilate_mask(rb_raw, grid_neigh)
                    so3_vals.append(dice(a["dil_mask"], rb_dil))

                    # Distant-time shuffle null.
                    ti = all_times[int(p["snapshot_i"])]
                    distant = [
                        q for q in all_field_ids
                        if abs(all_times[q] - ti) >= DISTANT_TIME_MIN_SEPARATION_M
                    ]
                    if distant:
                        q = int(rng.choice(distant))
                        cq = components.get((field, q))
                        if cq is not None:
                            time_vals.append(dice(a["dil_mask"], cq["dil_mask"]))

                so3_med = float(np.median(so3_vals)) if so3_vals else np.nan
                time_med = float(np.median(time_vals)) if time_vals else np.nan
                so3_trial_medians.append(so3_med)
                time_trial_medians.append(time_med)
                null_trial_rows.append({
                    "field": field,
                    "lag_M": lag,
                    "trial": trial,
                    "so3_median_dice": so3_med,
                    "distant_time_median_dice": time_med,
                })

            so3_arr = np.asarray(so3_trial_medians, float)
            time_arr = np.asarray(time_trial_medians, float)
            so3_arr = so3_arr[np.isfinite(so3_arr)]
            time_arr = time_arr[np.isfinite(time_arr)]

            so3_q99 = float(np.quantile(so3_arr, 0.99)) if len(so3_arr) else np.nan
            time_q95 = float(np.quantile(time_arr, 0.95)) if len(time_arr) else np.nan
            so3_p = empirical_upper_p(real_med, so3_arr) if len(so3_arr) else np.nan
            time_p = empirical_upper_p(real_med, time_arr) if len(time_arr) else np.nan

            so3_margin = real_med - so3_q99 if np.isfinite(so3_q99) else np.nan
            time_margin = real_med - time_q95 if np.isfinite(time_q95) else np.nan

            lag_pass = bool(
                len(valid_pairs) >= MIN_MATCHED_PAIRS
                and real_med >= REAL_DICE_MIN
                and real_step <= REAL_CENTROID_STEP_MAX_DEG
                and np.isfinite(so3_p) and so3_p <= SO3_P_MAX
                and np.isfinite(so3_margin) and so3_margin >= SO3_MARGIN_Q99_MIN
            )

            time_specific = bool(
                len(valid_pairs) >= MIN_MATCHED_PAIRS
                and np.isfinite(time_p) and time_p <= TIME_P_MAX
                and np.isfinite(time_margin) and time_margin >= TIME_MARGIN_Q95_MIN
            )

            null_summary_rows.append({
                "field": field,
                "lag_M": lag,
                "matched_pairs": int(len(valid_pairs)),
                "real_median_dilated_dice": real_med,
                "real_median_centroid_step_deg": real_step,
                "so3_p": so3_p,
                "so3_q99": so3_q99,
                "so3_margin_vs_q99": so3_margin,
                "time_p": time_p,
                "time_q95": time_q95,
                "time_margin_vs_q95": time_margin,
                "lag_coherence_pass": lag_pass,
                "time_specific_pass": time_specific,
            })

            print(
                f"{field:5s} lag={lag:.2f}M pairs={len(valid_pairs):4d} "
                f"realDice={real_med:.3f} step={real_step:.2f}° "
                f"SO3 p={so3_p:.4f} margin99={so3_margin:.3f} "
                f"time p={time_p:.4f} margin95={time_margin:.3f} "
                f"lagPass={lag_pass} timeSpecific={time_specific}",
                flush=True,
            )

    pair_df = pd.DataFrame(pair_rows)
    null_trials = pd.DataFrame(null_trial_rows)
    null_summary = pd.DataFrame(null_summary_rows)

    field_rows = []
    for field in FIELDS:
        g = null_summary[null_summary["field"] == field]
        n_lag = int(g["lag_coherence_pass"].sum())
        n_time = int(g["time_specific_pass"].sum())
        field_rows.append({
            "field": field,
            "lag_coherence_lags_passed": n_lag,
            "lag_coherence_field_pass": bool(n_lag >= MIN_LAGS_PER_FIELD),
            "time_specific_lags_passed": n_time,
            "time_specific_field_pass": bool(n_time >= MIN_LAGS_PER_FIELD),
        })
    field_summary = pd.DataFrame(field_rows)

    n_fields_lag = int(field_summary["lag_coherence_field_pass"].sum())
    n_fields_time = int(field_summary["time_specific_field_pass"].sum())

    global_lag = bool(n_fields_lag >= MIN_FIELDS_GLOBAL)
    global_time = bool(n_fields_time >= MIN_FIELDS_GLOBAL)

    if not global_lag:
        classification = "BH_HF1RC_NO_LAGGED_COHERENCE"
    elif global_time:
        classification = "BH_HF1RC_LAGGED_COHERENCE_TIME_SPECIFIC"
    else:
        classification = "BH_HF1RC_LAGGED_COHERENCE_NOT_TIME_SPECIFIC"

    pair_df.to_csv(out / "BH_HF1RC_pair_metrics.csv", index=False)
    null_trials.to_csv(out / "BH_HF1RC_null_trials.csv", index=False)
    null_summary.to_csv(out / "BH_HF1RC_null_summary.csv", index=False)
    field_summary.to_csv(out / "BH_HF1RC_field_summary.csv", index=False)

    make_plots(out, null_summary)

    summary = {
        "stage": STAGE,
        "classification": classification,
        "simulation": SID,
        "lev": LEV,
        "horizon": HORIZON,
        "frozen_lags_M": LAGS_M,
        "fields_passing_lagged_coherence": n_fields_lag,
        "fields_passing_time_specificity": n_fields_time,
        "required_fields_global": MIN_FIELDS_GLOBAL,
        "field_summary": field_rows,
        "methodological_status": (
            "coordinate/gauge-dependent finite-lag temporal morphology control; "
            "not intrinsic/covariant, not QNM-subtracted"
        ),
        "interpretation": {
            "BH_HF1RC_LAGGED_COHERENCE_TIME_SPECIFIC": (
                "Finite-lag morphology is coherent and more similar at nearby "
                "times than under both independent SO(3) rotation and distant-time "
                "shuffle nulls. Proceed to BH-HF2 QNM/multipole closure."
            ),
            "BH_HF1RC_LAGGED_COHERENCE_NOT_TIME_SPECIFIC": (
                "Finite-lag morphology survives SO(3) randomization but is not "
                "more similar at requested nearby times than at distant times. "
                "This favors a persistent/slowly varying coordinate multipolar "
                "pattern over a distinct transient object."
            ),
            "BH_HF1RC_NO_LAGGED_COHERENCE": (
                "The dense HF1R pass does not survive fixed-lag control. Stop "
                "the transient-morphology route."
            ),
        }[classification],
        "forbidden_claims": [
            "not a physical phase measurement",
            "not a domain wall detection",
            "not an intrinsic/covariant horizon-gradient measurement",
            "not gauge invariant",
            "not slicing invariant",
            "not numerical-resolution robust yet",
            "not QNM-subtracted",
            "does not establish new physics",
            "does not establish an ether",
            "does not establish AARFI or APK",
        ],
        "next": (
            "BH-HF2 QNM/multipole closure"
            if classification == "BH_HF1RC_LAGGED_COHERENCE_TIME_SPECIFIC"
            else (
                "Simplified multipole/stationarity closure; do not call it a transient object"
                if classification == "BH_HF1RC_LAGGED_COHERENCE_NOT_TIME_SPECIFIC"
                else "Stop this transient-morphology route"
            )
        ),
    }

    (out / "BH_HF1RC_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )

    protocol = {
        "stage": STAGE,
        "status": "FROZEN_BEFORE_LOCAL_EXECUTION",
        "simulation": SID,
        "lev": LEV,
        "horizon": HORIZON,
        "lags_M": LAGS_M,
        "matching_tolerance_M": "max(0.03, 0.15*lag)",
        "fixed_grid_points": GRID_N,
        "native_active_quantile": ACTIVE_Q,
        "primary_lag_rule": {
            "min_matched_pairs": MIN_MATCHED_PAIRS,
            "real_median_dilated_dice_min": REAL_DICE_MIN,
            "real_median_centroid_step_max_deg": REAL_CENTROID_STEP_MAX_DEG,
            "so3_empirical_p_max": SO3_P_MAX,
            "real_minus_so3_q99_min": SO3_MARGIN_Q99_MIN,
            "min_lags_per_field": MIN_LAGS_PER_FIELD,
            "min_fields_global": MIN_FIELDS_GLOBAL,
        },
        "time_specificity_rule": {
            "distant_time_min_separation_M": DISTANT_TIME_MIN_SEPARATION_M,
            "empirical_p_max": TIME_P_MAX,
            "real_minus_distant_time_q95_min": TIME_MARGIN_Q95_MIN,
            "min_lags_per_field": MIN_LAGS_PER_FIELD,
            "min_fields_global": MIN_FIELDS_GLOBAL,
        },
        "nulls": {
            "SO3": "independent random rotation of second component",
            "distant_time": "random positive component from same field >=2M away",
            "trials": N_NULL_TRIALS,
            "pair_cap_per_field_lag_for_null": NULL_PAIR_CAP,
            "rng_seed": RNG_SEED,
        },
        "forbidden_claims": summary["forbidden_claims"],
    }
    (out / "BH_HF1RC_preregistered_protocol.json").write_text(
        json.dumps(protocol, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )

    print("\n" + json.dumps(summary, ensure_ascii=False, indent=2))
    print(f"\nOutputs written to: {out.resolve()}")


def self_test():
    rng = np.random.default_rng(123)
    grid = fibonacci_sphere(1024)
    tree = cKDTree(grid)
    neigh, _ = knn_graph_on_sphere(grid, 8)

    assert np.allclose(np.linalg.norm(grid, axis=1), 1.0, atol=1e-12)

    raw = np.zeros(len(grid), bool)
    raw[:30] = True
    d = dilate_mask(raw, neigh)
    assert dice(d, d) == 1.0

    R = random_rotation_matrix(rng)
    assert np.allclose(R @ R.T, np.eye(3), atol=1e-12)
    assert abs(np.linalg.det(R) - 1.0) < 1e-12

    rr = rotate_raw_mask(raw, grid, tree, rng)
    assert rr.sum() > 0

    assert abs(angle_deg(np.array([1.,0,0]), np.array([0.,1,0])) - 90.0) < 1e-10
    print("BH-HF1R-C self-test: PASS")


def parser():
    p = argparse.ArgumentParser()
    p.add_argument("--self-test", action="store_true")
    p.add_argument(
        "--snapshot-metrics",
        default="/mnt/c/Users/Admin/bh_hf1r_output/BH_HF1R_snapshot_metrics.csv",
    )
    p.add_argument(
        "--h2-script",
        default="/mnt/c/Users/Admin/gw150914_GW_XA4_4E2H2R2_surface_field_recurrence_nonreturn.py",
    )
    p.add_argument(
        "--horizons",
        default="/mnt/c/Users/Admin/Horizons.h5",
    )
    p.add_argument(
        "--dump",
        default="/mnt/c/Users/Admin/HorizonsDump.h5",
    )
    p.add_argument(
        "--out",
        default="/mnt/c/Users/Admin/bh_hf1rc_output",
    )
    return p


if __name__ == "__main__":
    args = parser().parse_args()
    if args.self_test:
        self_test()
    else:
        for path in [args.snapshot_metrics, args.h2_script, args.horizons, args.dump]:
            if not Path(path).exists():
                raise SystemExit(f"Missing required file: {path}")
        run(args)
