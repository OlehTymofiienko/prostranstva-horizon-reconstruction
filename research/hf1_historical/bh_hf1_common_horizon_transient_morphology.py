#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
BH-HF1 — Common-Horizon Transient Morphology Screen
====================================================

Purpose
-------
Exploratory, explicitly gauge-dependent screen for short-lived localized
spatial morphology on the newly formed common apparent horizon C of
SXS:BBH:0305 / Lev6.

This stage does NOT test a "domain wall", does NOT define a physical phase,
and does NOT claim new physics.  It asks only whether the already-available
surface scalar fields contain a localized early-common-horizon morphology
worth subjecting to stronger GR/QNM, resolution, frame, slicing, and
intrinsic-measure controls.

Frozen source fields
--------------------
Rstar  := stored DimlessRicciScalar
Estar  := M_Ch^2 * WeylE_NN
Bstar  := M_Ch^2 * WeylB_NN

The source parser is imported from the previously audited frozen H2R2 script:
  gw150914_GW_XA4_4E2H2R2_surface_field_recurrence_nonreturn.py

Common-horizon times
--------------------
The common-horizon scalar-field birth time is taken directly from the first
AhC/DimlessRicciScalar.tdm index entry.  Frozen target offsets:
  0, 2, 5, 10, 20, 50 M

Morphology proxy
----------------
Because public SXS products used here do not provide the evolving intrinsic
local area measure or full induced 2-metric, HF1 works on an explicitly
coordinate/gauge-dependent angular point-cloud proxy.

For each snapshot:
  * center points by the catalog CoordCenterInertial of horizon C;
  * map each surface point to the unit direction n = (x-center)/|x-center|;
  * robust-standardize each scalar field spatially;
  * build a k-nearest-neighbour graph on the unit sphere;
  * define a local angular-gradient proxy from neighbour differences;
  * threshold the top 5% gradient-score points;
  * measure the largest connected high-gradient component;
  * compare that component fraction with 199 spatial-value shuffles.

This is a screening statistic only.  A candidate may proceed to HF1R/HF2;
a negative result closes this specific morphology-screen definition.

Frozen candidate rule
---------------------
A field is an HF1 morphology candidate only if at least TWO of the six frozen
snapshots satisfy BOTH:
  p_shuffle <= 0.01
  largest_component_fraction_of_active >= 0.35

and at least one positive snapshot occurs at offset <= 10 M.

This rule tests spatial localization, NOT temporal identity of the component.
Temporal tracking is deliberately deferred to HF1R to avoid silently assuming
cross-time surface-point correspondence.

Outputs
-------
BH_HF1_snapshot_metrics.csv
BH_HF1_summary.json
BH_HF1_preregistered_protocol.json
BH_HF1_<field>_<offset>M_surface_proxy.png
BH_HF1_localization_timeline.png
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import math
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from scipy.spatial import cKDTree


STAGE = "BH-HF1"
SID = "SXS:BBH:0305"
LEV = 6
HORIZON = "C"

OFFSETS_M = [0.0, 2.0, 5.0, 10.0, 20.0, 50.0]

KNN = 8
ACTIVE_Q = 0.95
N_SHUFFLE = 199
RNG_SEED = 150914

P_GATE = 0.01
COMPONENT_FRACTION_GATE = 0.35
MIN_POSITIVE_SNAPSHOTS = 2
EARLY_REQUIRED_MAX_OFFSET_M = 10.0

FIELDS = [
    ("Rstar", "DimlessRicciScalar.tdm"),
    ("Estar", "WeylE_NN.tdm"),
    ("Bstar", "WeylB_NN.tdm"),
]


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
    n = x[good] / r[good, None]
    return n, good


def knn_graph_on_sphere(n, k=KNN):
    tree = cKDTree(n)
    d, idx = tree.query(n, k=min(k + 1, len(n)))
    # Remove self
    neigh = idx[:, 1:]
    # Chord distance -> angular distance
    dot = np.einsum("ij,ikj->ik", n, n[neigh])
    dot = np.clip(dot, -1.0, 1.0)
    ang = np.arccos(dot)
    return neigh, ang


def gradient_proxy(z, neigh, ang):
    zi = z[:, None]
    zj = z[neigh]
    g = np.median(np.abs(zj - zi) / np.maximum(ang, 1e-8), axis=1)
    return g


def largest_active_component(active, neigh):
    active = np.asarray(active, bool)
    ids = np.flatnonzero(active)
    if len(ids) == 0:
        return 0, np.array([], dtype=int)

    aset = set(ids.tolist())
    visited = set()
    best = []

    for s in ids:
        if int(s) in visited:
            continue
        stack = [int(s)]
        visited.add(int(s))
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

    return len(best), np.asarray(best, dtype=int)


def morphology_metrics(z, n, neigh, ang, rng):
    g = gradient_proxy(z, neigh, ang)
    thr = float(np.quantile(g, ACTIVE_Q))
    active = g >= thr

    largest_n, largest_idx = largest_active_component(active, neigh)
    active_n = int(active.sum())
    frac = largest_n / max(active_n, 1)

    null_frac = []
    for _ in range(N_SHUFFLE):
        zp = rng.permutation(z)
        gp = gradient_proxy(zp, neigh, ang)
        tp = float(np.quantile(gp, ACTIVE_Q))
        ap = gp >= tp
        ln, _ = largest_active_component(ap, neigh)
        null_frac.append(ln / max(int(ap.sum()), 1))

    null_frac = np.asarray(null_frac, float)
    p = float((1 + np.sum(null_frac >= frac - 1e-15)) / (N_SHUFFLE + 1))

    if len(largest_idx):
        v = np.sum(n[largest_idx], axis=0)
        vn = np.linalg.norm(v)
        centroid = v / vn if vn > 0 else np.array([np.nan, np.nan, np.nan])
    else:
        centroid = np.array([np.nan, np.nan, np.nan])

    return {
        "gradient_q95": thr,
        "gradient_median": float(np.median(g)),
        "gradient_max": float(np.max(g)),
        "active_points": active_n,
        "largest_component_points": int(largest_n),
        "largest_component_fraction_of_active": float(frac),
        "shuffle_p": p,
        "null_component_fraction_median": float(np.median(null_frac)),
        "null_component_fraction_q95": float(np.quantile(null_frac, 0.95)),
        "centroid_x": float(centroid[0]),
        "centroid_y": float(centroid[1]),
        "centroid_z": float(centroid[2]),
        "candidate_snapshot": bool(p <= P_GATE and frac >= COMPONENT_FRACTION_GATE),
        "_g": g,
        "_active": active,
        "_largest_idx": largest_idx,
    }


def catalog_center(h2, hh, horizon, tq):
    t, y = h2.h5_series(hh, horizon, "CoordCenterInertial")
    return h2.interp_cols(t, y[:, :3], [tq])[0]


def mch_at_time(h2, hh, horizon, tq):
    t, y = h2.h5_series(hh, horizon, "ChristodoulouMass")
    # h5_series may already strip time; first numeric data column is mass
    if y.ndim == 1:
        vals = y
    else:
        vals = y[:, 0]
    return float(np.interp(tq, t, vals))


def field_values(raw_name, R, E, B, mch):
    if raw_name == "Rstar":
        return np.asarray(R, float)
    if raw_name == "Estar":
        return (mch * mch) * np.asarray(E, float)
    if raw_name == "Bstar":
        return (mch * mch) * np.asarray(B, float)
    raise KeyError(raw_name)


def plot_proxy(outpath, n, z, metric, title):
    # Mollweide-like longitude/latitude scatter. This is a coordinate-direction
    # visualization only; not an intrinsic horizon map.
    lon = np.arctan2(n[:, 1], n[:, 0])
    lat = np.arcsin(np.clip(n[:, 2], -1.0, 1.0))

    fig = plt.figure(figsize=(10, 5.4))
    ax = fig.add_subplot(111, projection="mollweide")
    sc = ax.scatter(lon, lat, c=z, s=8)
    li = metric["_largest_idx"]
    if len(li):
        ax.scatter(lon[li], lat[li], s=26, facecolors="none", edgecolors="black")
    ax.grid(True, alpha=0.3)
    ax.set_title(title)
    cb = fig.colorbar(sc, ax=ax, shrink=0.8)
    cb.set_label("robust spatial z-score")
    fig.tight_layout()
    fig.savefig(outpath, dpi=180)
    plt.close(fig)


def run(args):
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    h2 = import_by_path(Path(args.h2_script), "h2r2_frozen")

    rng = np.random.default_rng(RNG_SEED)
    rows = []

    with h5py.File(args.horizons, "r") as hh, h5py.File(args.dump, "r") as hd:
        r_index = h2.tdm_index(hd, HORIZON, "DimlessRicciScalar.tdm")
        # H2R2 index semantics: first column is time.
        t_birth = float(np.asarray(r_index)[0, 0])

        print(f"{SID} / Lev{LEV} common-horizon scalar birth: {t_birth:.12f} M")

        for off in OFFSETS_M:
            target = t_birth + off
            i, tsnap, mismatch = h2.nearest_snapshot_idx(r_index, target)

            R, exR = h2.read_scalar_snapshot(hd, HORIZON, "DimlessRicciScalar.tdm", i)
            E, exE = h2.read_scalar_snapshot(hd, HORIZON, "WeylE_NN.tdm", i)
            B, exB = h2.read_scalar_snapshot(hd, HORIZON, "WeylB_NN.tdm", i)
            C, exC = h2.read_coord_snapshot(hd, HORIZON, i)

            if not (exR == exE == exB == exC):
                raise RuntimeError(
                    f"Extent mismatch at offset {off}M: "
                    f"R={exR}, E={exE}, B={exB}, C={exC}"
                )

            center = catalog_center(h2, hh, HORIZON, float(tsnap))
            mch = mch_at_time(h2, hh, HORIZON, float(tsnap))

            n, good = unit_directions(C, center)
            neigh, ang = knn_graph_on_sphere(n, KNN)

            R = np.asarray(R)[good]
            E = np.asarray(E)[good]
            B = np.asarray(B)[good]

            for field, _ in FIELDS:
                f = field_values(field, R, E, B, mch)
                z = robust_z(f)

                m = morphology_metrics(z, n, neigh, ang, rng)

                row = {
                    "stage": STAGE,
                    "sxs_id": SID,
                    "lev": LEV,
                    "horizon": HORIZON,
                    "field": field,
                    "target_offset_M": float(off),
                    "target_time_M": float(target),
                    "actual_time_M": float(tsnap),
                    "time_mismatch_M": float(mismatch),
                    "n_surface_points": int(len(z)),
                    "Mch": float(mch),
                    **{k: v for k, v in m.items() if not k.startswith("_")},
                }
                rows.append(row)

                tag = str(off).replace(".", "p")
                plot_proxy(
                    out / f"BH_HF1_{field}_{tag}M_surface_proxy.png",
                    n,
                    z,
                    m,
                    (
                        f"{SID} Lev{LEV} common horizon C — {field} — "
                        f"t-tC≈{off:g}M\n"
                        "coordinate-direction morphology proxy"
                    ),
                )

                print(
                    f"  {field:5s} +{off:5.1f}M : "
                    f"frac={m['largest_component_fraction_of_active']:.3f}, "
                    f"p={m['shuffle_p']:.4f}, "
                    f"candidate={m['candidate_snapshot']}"
                )

    df = pd.DataFrame(rows)
    df.to_csv(out / "BH_HF1_snapshot_metrics.csv", index=False)

    field_summary = {}
    candidates = []
    for field, _ in FIELDS:
        g = df[df["field"] == field].sort_values("target_offset_M")
        npos = int(g["candidate_snapshot"].sum())
        early = bool(
            ((g["candidate_snapshot"]) &
             (g["target_offset_M"] <= EARLY_REQUIRED_MAX_OFFSET_M)).any()
        )
        is_candidate = bool(npos >= MIN_POSITIVE_SNAPSHOTS and early)
        if is_candidate:
            candidates.append(field)
        field_summary[field] = {
            "positive_snapshot_count": npos,
            "early_positive_at_or_before_10M": early,
            "HF1_morphology_candidate": is_candidate,
            "positive_offsets_M": [
                float(x)
                for x in g.loc[g["candidate_snapshot"], "target_offset_M"].tolist()
            ],
        }

    classification = (
        "BH_HF1_TRANSIENT_LOCALIZED_MORPHOLOGY_CANDIDATE"
        if candidates
        else
        "BH_HF1_NO_TRANSIENT_LOCALIZED_MORPHOLOGY_BY_FROZEN_SCREEN"
    )

    # Timeline figure
    fig = plt.figure(figsize=(10, 5.5))
    for field, _ in FIELDS:
        g = df[df["field"] == field].sort_values("target_offset_M")
        plt.plot(
            g["target_offset_M"],
            g["largest_component_fraction_of_active"],
            marker="o",
            label=field,
        )
    plt.axhline(
        COMPONENT_FRACTION_GATE,
        linestyle="--",
        label=f"component gate={COMPONENT_FRACTION_GATE}",
    )
    plt.xlabel("t - tC [M]")
    plt.ylabel("largest component / active high-gradient points")
    plt.title("BH-HF1 common-horizon localization screen")
    plt.legend()
    plt.tight_layout()
    fig.savefig(out / "BH_HF1_localization_timeline.png", dpi=180)
    plt.close(fig)

    summary = {
        "stage": STAGE,
        "classification": classification,
        "simulation": SID,
        "lev": LEV,
        "horizon": HORIZON,
        "common_horizon_scalar_birth_time_M": float(df["actual_time_M"].min()),
        "frozen_offsets_M": OFFSETS_M,
        "fields": [x[0] for x in FIELDS],
        "field_summary": field_summary,
        "candidate_fields": candidates,
        "methodological_status": (
            "exploratory coordinate/gauge-dependent morphology screen; "
            "not intrinsic/covariant"
        ),
        "interpretation_if_candidate": (
            "At least one stored horizon scalar has spatially localized "
            "high-gradient morphology in the early common-horizon window by "
            "the frozen point-cloud screen. This licenses only a denser "
            "temporal persistence audit (HF1R), not a domain-wall or new-physics claim."
        ),
        "interpretation_if_negative": (
            "No stored field passes the frozen spatial-localization screen. "
            "This specific transient-front route should stop unless a different "
            "physically motivated observable is preregistered independently."
        ),
        "forbidden_claims": [
            "not a physical phase measurement",
            "not a domain wall detection",
            "not an intrinsic horizon-distance or intrinsic-gradient measurement",
            "not gauge invariant",
            "not slicing invariant",
            "not numerical-resolution robust yet",
            "not QNM-subtracted",
            "does not establish new physics",
            "does not establish AARFI or APK",
        ],
        "next_if_candidate": (
            "BH-HF1R dense temporal tracking with frozen component matching; "
            "then BH-HF2 QNM/multipole closure; then resolution/frame/slicing controls."
        ),
        "next_if_negative": "Stop this morphology route.",
    }

    (out / "BH_HF1_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )

    protocol = {
        "stage": STAGE,
        "status": "FROZEN_BEFORE_LOCAL_EXECUTION",
        "simulation": SID,
        "lev": LEV,
        "horizon": HORIZON,
        "offsets_M": OFFSETS_M,
        "fields": [x[0] for x in FIELDS],
        "morphology_proxy": {
            "coordinate_center": "catalog CoordCenterInertial",
            "surface_position": "unit Cartesian direction from horizon center",
            "field_scaling": "robust spatial z-score per snapshot",
            "neighbor_graph": f"{KNN}-nearest-neighbour on unit sphere",
            "gradient_proxy": "median neighbour |dz| / angular separation",
            "active_threshold_quantile": ACTIVE_Q,
            "shuffle_trials": N_SHUFFLE,
        },
        "candidate_rule": {
            "shuffle_p_max": P_GATE,
            "largest_component_fraction_min": COMPONENT_FRACTION_GATE,
            "min_positive_snapshots": MIN_POSITIVE_SNAPSHOTS,
            "must_include_positive_offset_at_or_before_M": EARLY_REQUIRED_MAX_OFFSET_M,
        },
        "forbidden_claims": summary["forbidden_claims"],
    }
    (out / "BH_HF1_preregistered_protocol.json").write_text(
        json.dumps(protocol, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )

    print("\n" + json.dumps(summary, ensure_ascii=False, indent=2))
    print(f"\nOutputs written to: {out.resolve()}")


def self_test():
    rng = np.random.default_rng(123)
    n = rng.normal(size=(500, 3))
    n /= np.linalg.norm(n, axis=1, keepdims=True)
    neigh, ang = knn_graph_on_sphere(n, KNN)

    z = rng.normal(size=500)
    m = morphology_metrics(z, n, neigh, ang, rng)
    assert 0.0 <= m["largest_component_fraction_of_active"] <= 1.0
    assert 0.0 < m["shuffle_p"] <= 1.0
    print("BH-HF1 self-test: PASS")


def parser():
    p = argparse.ArgumentParser()
    p.add_argument("--self-test", action="store_true")
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
        default="/mnt/c/Users/Admin/bh_hf1_output",
    )
    return p


if __name__ == "__main__":
    args = parser().parse_args()
    if args.self_test:
        self_test()
    else:
        for path in [args.h2_script, args.horizons, args.dump]:
            if not Path(path).exists():
                raise SystemExit(f"Missing required file: {path}")
        run(args)
