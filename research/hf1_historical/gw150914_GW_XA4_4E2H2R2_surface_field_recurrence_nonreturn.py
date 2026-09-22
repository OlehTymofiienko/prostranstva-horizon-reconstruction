#!/usr/bin/env python3
"""
GW-X-A4-4E2H2R2 — surface-field phase-recurrence / secular-nonreturn audit
========================================================================

Purpose
-------
H1 established that the public Horizons-only quasi-local mass variables
are below the frozen numerical-resolution floor for one-orbit dynamics.
H2 therefore returns to the seven simulations for which HorizonsDump.h5
is already available: the old SXS:BBH:0305 reference plus the frozen
six-simulation A4 panel.

The question is deliberately narrower than "does spacetime metric change?":

    When the orbital phase returns after 2*pi, does the local horizon
    scalar-field pattern approximately recur, and if it does, is there
    nevertheless a systematic secular drift from orbit to orbit?

Fields
------
Per apparent horizon A/B:
    R* = stored DimlessRicciScalar
    E* = M_Ch^2 * WeylE_NN
    B* = M_Ch^2 * WeylB_NN

These are horizon scalar fields.  H2 does NOT reconstruct the spacetime
metric g_{mu nu}.

Technical repair H2R
--------------------
The original H2 implementation decoded CoordsMeasurementFrame.tdm as
interleaved xyz triples.  That contradicted the already frozen A3-1R
storage convention, which had established blocked_all_x_then_y_then_z
from the old 0305 reference before the independent A4 panel was tested.
H2R changes ONLY this storage decoding.  Scientific states, thresholds,
nulls, panel membership, matched phases, and classification gates are
unchanged.

Critical methodological guard
-----------------------------
Previous A3 work established that the public products do not provide the
evolving intrinsic local area element dA or an intrinsic angular frame.
Therefore H2 does NOT compute or claim true spherical-harmonic multipoles.

Instead it uses the already frozen A3/A4 coordinate-area surrogate:
    F_periodic1 connectivity
    triangle diagonal diag_10_01

and constructs explicitly gauge-dependent "surrogate angular moments" in
one FIXED coordinate frame per simulation.  The frame is frozen at the
first same-phase matched state.

Why a fixed frame?
------------------
The same-phase endpoints are one full orbit apart and should point in the
same coordinate direction.  A half-orbit midpoint should be phase-opposed.
This gives a direct recurrence control:

    d_full  = distance(state_n, state_{n+1})
    d_half1 = distance(state_n, half-phase_n)
    d_half2 = distance(half-phase_n, state_{n+1})

    return_ratio = d_full / mean(d_half1, d_half2)

A value < 1 means the full-orbit same-phase state is closer than the
half-phase control.  This is a recurrence diagnostic, not a new-physics
test.

Shape-only state
----------------
At each snapshot each field is area-surrogate centered and RMS-normalized:

    z_F = (F - <F>_w) / RMS_w(F - <F>_w)

We then compute in the fixed frame:
    3 dipole moments
    5 independent STF quadrupole moments

for each of {R*, E*, B*} and each of horizons {A,B}.

Primary shape dimension = 2 * 3 * 8 = 48.

This removes per-snapshot scalar amplitude and offset from the PRIMARY
state so that the test asks about spatial-pattern recurrence rather than
the already-studied overall amplitude scale.

Quality guards
--------------
H2R2 repairs the SEMANTICS of the quality guards without changing the
scientific recurrence/nonreturn statistics.

The parent H2/H2R incorrectly treated two diagnostics as universal hard
physical identities:

1) coordinate-area centroid == CoordCenterInertial within 0.05 M
2) coordinate-area-surrogate Gauss-Bonnet residual <= 0.05

Neither equality is justified as a hard requirement here:
- CoordsMeasurementFrame.tdm is an inertial Cartesian surface grid, but
  CoordCenterInertial is the catalog horizon center, not a statement that
  the coordinate-area-weighted polygon centroid must coincide with it.
- Gauss-Bonnet is an intrinsic-area statement, while H2 intentionally uses
  a coordinate-area surrogate because local intrinsic dA is unavailable.

Therefore centroid error and GB surrogate residual are now RECORDED
DIAGNOSTICS ONLY.  They do not decide inclusion.

Hard technical guards remain:
1. Frozen reference extractor convention from A4-3R:
       diag_10_01
2. Frozen coordinate storage convention from A3-1R:
       blocked_all_x_then_y_then_z
3. scalar/coordinate extents must match and all required arrays must exist
4. coordinate-area surrogate must be finite and strictly positive
5. nearest HorizonsDump snapshot mismatch:
       <= 1.0 M
6. each matched interval must span approximately one orbit:
       |Delta phi - 2*pi| <= 0.5 rad

Legacy diagnostics retained for transparency:
- |GB surrogate residual| compared with 0.05
- coordinate-area centroid error compared with 0.05 M
These legacy values are not gates in H2R2.

Ordering / secular non-return
-----------------------------
For the 10 same-phase endpoints in each simulation:
- standardize active 48-D shape coordinates within that simulation
- ordered path length is compared with 4999 random state permutations
- |Spearman(PC1, orbit_index)| is compared with 4999 orbit-index permutations

A simulation is "ordered nonreturn positive" if:
    p_path <= 0.05 AND p_pc1 <= 0.05

Panel recurrence
----------------
The six A4 simulations are the independent gate panel.
0305 is a reference diagnostic only.

For each panel simulation define:
    recurrence_effect = mean( mean(d_half1,d_half2) - d_full )

Positive means same-phase full-orbit states are closer than half-phase
controls.  The panel uses an exact one-sided sign-flip test over the six
simulation-level effects.

Recurrence supported if:
    sign-flip p <= 0.05
    >= 4/6 simulations have median return_ratio < 1

Ordered secular non-return replicated if:
    >= 4/6 simulations are ordered-nonreturn positive

Classification
--------------
Both:
    A4_4E2H2R2_PHASE_RECURRENCE_WITH_ORDERED_FIELD_NONRETURN
Recurrence only:
    A4_4E2H2R2_PHASE_RECURRENCE_WITHOUT_ORDERED_FIELD_NONRETURN
Ordered nonreturn only:
    A4_4E2H2R2_ORDERED_FIELD_NONRETURN_WITHOUT_CLEAR_PHASE_RECURRENCE
Neither:
    A4_4E2H2R2_NO_CLEAR_SURFACE_FIELD_ARCHITECTURE
Quality failure:
    A4_4E2H2R2_QUALITY_INCOMPLETE

Forbidden claims
----------------
- no identification of g_{mu nu}
- no proof of a changing physical spacetime metric
- no true intrinsic spherical-harmonic multipoles
- no Tor-APK validation
- no fifth-dimension claim
- no exclusion of ordinary smooth GR or SXS gauge/numerical effects
"""

from __future__ import annotations

import argparse
import itertools
import json
import math
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

try:
    from scipy.stats import spearmanr
except Exception:
    spearmanr = None


EXPECTED_H1 = "A4_4E2H1_QUASILOCAL_DYNAMICS_BELOW_RESOLUTION_FLOOR"
EXPECTED_A43R = "A4_3_REPLICATED_ACROSS_PANEL"
EXPECTED_DIAGONAL = "diag_10_01"

N_MATCHED = 10
N_INTERVALS = 9
MAX_DUMP_TIME_MISMATCH = 1.0
MAX_GB_ABS_RESIDUAL = 0.05
MAX_CENTROID_ERROR_M = 0.05
MAX_ORBIT_PHASE_ERROR = 0.5

PERM_TRIALS = 4999
RNG_SEED = 150914

PANEL = [
    ("SXS:BBH:0318", "v1.2", 4, "close_0305"),
    ("SXS:BBH:1132", "v1.2", 4, "equal_mass_lowspin"),
    ("SXS:BBH:1166", "v1.2", 3, "q2_lowspin"),
    ("SXS:BBH:1220", "v1.2", 4, "q4_lowspin"),
    ("SXS:BBH:0411", "v1.2", 3, "aligned_spin"),
    ("SXS:BBH:0398", "v1.2", 3, "antialigned_spin"),
]

FIELDS = [
    ("Rstar", "DimlessRicciScalar.tdm"),
    ("Estar", "WeylE_NN.tdm"),
    ("Bstar", "WeylB_NN.tdm"),
]


def load_json(path):
    return json.loads(Path(path).read_text(encoding="utf-8"))


def normalize_sid(x):
    s = str(x)
    if "SXS:BBH:" in s:
        i = s.index("SXS:BBH:")
        return s[i:i+12] if len(s) >= i+12 else s[i:]
    digits = "".join(c for c in s if c.isdigit())
    if len(digits) >= 4:
        return "SXS:BBH:" + digits[-4:]
    return s


def infer_time_column(df):
    priority = [
        "time_M", "time", "state_time_M", "phase_locked_time_M",
        "matched_time_M", "source_time_M", "t_M",
    ]
    for c in priority:
        if c in df.columns and np.issubdtype(df[c].dtype, np.number):
            return c
    candidates = [
        c for c in df.columns
        if "time" in c.lower()
        and "lower" not in c.lower()
        and "upper" not in c.lower()
        and np.issubdtype(df[c].dtype, np.number)
    ]
    if len(candidates) == 1:
        return candidates[0]
    if candidates:
        # prefer names ending in _M
        m = [c for c in candidates if c.lower().endswith("_m")]
        return m[0] if m else candidates[0]
    raise RuntimeError(f"Could not infer time column; columns={list(df.columns)}")


def infer_sid_column(df):
    priority = ["sxs_id", "simulation", "sid", "SXS", "sxs"]
    for c in priority:
        if c in df.columns:
            return c
    for c in df.columns:
        if df[c].dtype == object:
            vals = df[c].astype(str)
            if vals.str.contains("SXS:BBH:", regex=False).any():
                return c
    return None


def infer_index_column(df):
    for c in ["state_index", "cycle_index", "orbit_index", "matched_index", "index"]:
        if c in df.columns and np.issubdtype(df[c].dtype, np.number):
            return c
    return None


def load_matched_csv(path, forced_sid=None):
    df = pd.read_csv(path)
    tcol = infer_time_column(df)
    sidcol = infer_sid_column(df)
    idxcol = infer_index_column(df)

    if forced_sid is not None:
        df = df.copy()
        df["__sid"] = forced_sid
    elif sidcol is None:
        raise RuntimeError(
            f"{path}: no SXS simulation column found; columns={list(df.columns)}"
        )
    else:
        df = df.copy()
        df["__sid"] = df[sidcol].map(normalize_sid)

    df["__time"] = pd.to_numeric(df[tcol], errors="coerce")
    if df["__time"].isna().any():
        raise RuntimeError(f"{path}: non-numeric matched-state time values")

    if idxcol is not None:
        df["__index"] = pd.to_numeric(df[idxcol], errors="coerce")
    else:
        df["__index"] = np.nan

    out = {}
    for sid, g in df.groupby("__sid"):
        if g["__index"].notna().all():
            g = g.sort_values("__index")
        else:
            g = g.sort_values("__time")
        times = g["__time"].to_numpy(float)
        if len(times) != N_MATCHED:
            raise RuntimeError(
                f"{path}: {sid} expected {N_MATCHED} matched states, got {len(times)}"
            )
        if not np.all(np.diff(times) > 0):
            raise RuntimeError(f"{path}: {sid} matched times are not strictly increasing")
        out[sid] = times
    return out


def h5_series(h5, horizon, basename):
    path = f"Ah{horizon}.dir/{basename}.dat"
    if path not in h5:
        raise KeyError(path)
    a = np.asarray(h5[path][:], float)
    if a.ndim != 2 or a.shape[1] < 2:
        raise RuntimeError(f"{path}: unexpected shape {a.shape}")
    return a[:, 0], a[:, 1:]


def interp_cols(t, y, tq):
    t = np.asarray(t, float)
    y = np.asarray(y, float)
    tq = np.asarray(tq, float)
    if y.ndim == 1:
        return np.interp(tq, t, y)
    return np.column_stack([np.interp(tq, t, y[:, j]) for j in range(y.shape[1])])


def common_centers_and_phase(h5):
    tA, cA = h5_series(h5, "A", "CoordCenterInertial")
    tB, cB = h5_series(h5, "B", "CoordCenterInertial")
    cA = cA[:, :3]
    cB = cB[:, :3]

    lo = max(tA[0], tB[0])
    hi = min(tA[-1], tB[-1])
    n = min(
        np.sum((tA >= lo) & (tA <= hi)),
        np.sum((tB >= lo) & (tB <= hi)),
    )
    n = max(1000, int(n))
    t = np.linspace(lo, hi, n)
    A = interp_cols(tA, cA, t)
    B = interp_cols(tB, cB, t)
    r = B - A

    phi = np.unwrap(np.arctan2(r[:, 1], r[:, 0]))
    d = np.diff(phi)
    if np.nanmedian(d) < 0:
        phi = -phi
    # minor local coordinate non-monotonicity is removed only for inversion
    phi_mono = np.maximum.accumulate(phi)

    drdt = np.gradient(r, t, axis=0)
    return t, A, B, r, drdt, phi, phi_mono


def fixed_frame(centers, t0):
    t, A, B, r, drdt, phi, phi_mono = centers
    rv = interp_cols(t, r, [t0])[0]
    vv = interp_cols(t, drdt, [t0])[0]
    ex = rv / np.linalg.norm(rv)
    ez = np.cross(rv, vv)
    nez = np.linalg.norm(ez)
    if nez <= 1e-12:
        raise RuntimeError("Could not define frozen orbital-axis frame")
    ez = ez / nez
    ey = np.cross(ez, ex)
    ey = ey / np.linalg.norm(ey)
    # re-orthogonalize ex for numerical hygiene
    ex = np.cross(ey, ez)
    ex = ex / np.linalg.norm(ex)
    return np.vstack([ex, ey, ez])  # rows are basis vectors


def midpoint_times(centers, endpoint_times):
    t, A, B, r, drdt, phi_raw, phi_mono = centers
    out = []
    dphis = []
    for a, b in zip(endpoint_times[:-1], endpoint_times[1:]):
        pa = float(np.interp(a, t, phi_mono))
        pb = float(np.interp(b, t, phi_mono))
        dphi = pb - pa
        dphis.append(dphi)
        if abs(dphi - 2*np.pi) > MAX_ORBIT_PHASE_ERROR:
            raise RuntimeError(
                f"Matched-state interval phase span={dphi:.6f} rad, "
                f"outside 2*pi +/- {MAX_ORBIT_PHASE_ERROR}"
            )
        pm = 0.5 * (pa + pb)
        mask = (t >= a) & (t <= b)
        tt = t[mask]
        pp = phi_mono[mask]
        keep = np.r_[True, np.diff(pp) > 1e-12]
        tt = tt[keep]
        pp = pp[keep]
        if len(tt) < 3:
            raise RuntimeError("Too few monotone phase samples for half-phase inversion")
        tm = float(np.interp(pm, pp, tt))
        out.append(tm)
    return np.asarray(out), np.asarray(dphis)


def tdm_index(h5, horizon, field_tdm):
    p = f"Ah{horizon}.dir/{field_tdm}/index.dat"
    if p not in h5:
        raise KeyError(p)
    return np.asarray(h5[p][:], float)


def tdm_data(h5, horizon, field_tdm):
    p = f"Ah{horizon}.dir/{field_tdm}/data.dat"
    if p not in h5:
        raise KeyError(p)
    return h5[p]


def nearest_snapshot_idx(index, tq):
    times = index[:, 0]
    i = int(np.argmin(np.abs(times - tq)))
    return i, float(times[i]), float(abs(times[i] - tq))


def segment_bounds(index, i):
    end = int(round(index[i, 1]))
    start = 0 if i == 0 else int(round(index[i-1, 1]))
    e0 = int(round(index[i, 2]))
    e1 = int(round(index[i, 3]))
    return start, end, e0, e1


def read_scalar_snapshot(h5, horizon, tdm, i):
    idx = tdm_index(h5, horizon, tdm)
    start, end, e0, e1 = segment_bounds(idx, i)
    data = np.asarray(tdm_data(h5, horizon, tdm)[start:end], float).reshape(-1)
    n = e0 * e1
    if len(data) != n:
        raise RuntimeError(
            f"{horizon}/{tdm} snapshot {i}: scalar segment len={len(data)}, "
            f"extents product={n}"
        )
    return data, (e0, e1)


def read_coord_snapshot(h5, horizon, i):
    """
    TECHNICAL REPAIR H2R.

    CoordsMeasurementFrame.tdm stores the Cartesian components in blocked
    component order for the frozen A3/A4 extractor:

        [x_0,...,x_{N-1}, y_0,...,y_{N-1}, z_0,...,z_{N-1}]

    NOT as interleaved xyz triples.

    This layout was frozen before H2 from the A3-1R coordinate-layout audit.
    It is therefore not selected from the H2 panel results.
    """
    tdm = "CoordsMeasurementFrame.tdm"
    idx = tdm_index(h5, horizon, tdm)
    start, end, e0, e1 = segment_bounds(idx, i)
    data = np.asarray(tdm_data(h5, horizon, tdm)[start:end], float).reshape(-1)
    n = e0 * e1
    if len(data) != 3*n:
        raise RuntimeError(
            f"{horizon}/{tdm} snapshot {i}: coord segment len={len(data)}, "
            f"expected={3*n}"
        )

    # Frozen storage convention: blocked_all_x_then_y_then_z.
    coords = np.column_stack((
        data[0:n],
        data[n:2*n],
        data[2*n:3*n],
    ))
    return coords, (e0, e1)


def surrogate_vertex_weights(coords, extents):
    """
    Frozen A3/A4 convention:
      F_periodic1
      reshape = (Extents[0], Extents[1], 3)
      periodic second logical axis
      diagonal = diag_10_01
    """
    e0, e1 = extents
    g = np.asarray(coords, float).reshape(e0, e1, 3)
    w = np.zeros((e0, e1), float)

    for i in range(e0 - 1):
        for j in range(e1):
            j1 = (j + 1) % e1
            p00 = g[i, j]
            p10 = g[i+1, j]
            p01 = g[i, j1]
            p11 = g[i+1, j1]

            # diag_10_01: triangles (00,10,01) and (11,01,10)
            a1 = 0.5 * np.linalg.norm(np.cross(p10-p00, p01-p00))
            a2 = 0.5 * np.linalg.norm(np.cross(p01-p11, p10-p11))

            for ii, jj in [(i,j), (i+1,j), (i,j1)]:
                w[ii,jj] += a1 / 3.0
            for ii, jj in [(i+1,j1), (i,j1), (i+1,j)]:
                w[ii,jj] += a2 / 3.0

    wf = w.reshape(-1)
    s = wf.sum()
    if not np.isfinite(s) or s <= 0:
        raise RuntimeError("Coordinate-area surrogate produced non-positive total area")
    return wf / s, float(s)


def weighted_mean(x, w):
    return float(np.sum(w*x))


def weighted_rms_centered(x, w, mu=None):
    if mu is None:
        mu = weighted_mean(x, w)
    return float(np.sqrt(np.sum(w*(x-mu)**2)))


def dimensionless_fields(R, E, B, mch):
    m2 = float(mch*mch)
    return {
        "Rstar": np.asarray(R, float),
        "Estar": m2 * np.asarray(E, float),
        "Bstar": m2 * np.asarray(B, float),
    }


def angular_moments(field, coords, center, basis, w):
    d = coords - center[None, :]
    # rows of basis are ex, ey, ez
    xyz = d @ basis.T
    rad = np.linalg.norm(xyz, axis=1)
    good = np.isfinite(rad) & (rad > 1e-12) & np.isfinite(field)
    if good.sum() < 0.95 * len(field):
        raise RuntimeError("Too many invalid surface points for angular moments")
    xyz = xyz[good]
    f = field[good]
    wg = w[good]
    wg = wg / wg.sum()
    n = xyz / np.linalg.norm(xyz, axis=1)[:, None]

    mu = weighted_mean(f, wg)
    sig = weighted_rms_centered(f, wg, mu)
    if not np.isfinite(sig) or sig <= 1e-14:
        raise RuntimeError("Degenerate scalar field RMS in shape normalization")
    z = (f - mu) / sig

    dip = np.sum(wg[:,None] * z[:,None] * n, axis=0)

    Q = np.zeros((3,3), float)
    I = np.eye(3)
    for k in range(len(z)):
        Q += wg[k] * z[k] * (np.outer(n[k], n[k]) - I/3.0)

    # Five independent STF entries; Qzz = -Qxx-Qyy.
    feat = np.array([
        dip[0], dip[1], dip[2],
        Q[0,0], Q[1,1], Q[0,1], Q[0,2], Q[1,2],
    ], float)

    return feat, {
        "weighted_mean": mu,
        "weighted_rms_centered": sig,
    }


def interp_horizon_scalar(h5, horizon, basename, tq):
    t, y = h5_series(h5, horizon, basename)
    return float(np.interp(tq, t, y[:,0]))


def interp_center(h5, horizon, tq):
    t, y = h5_series(h5, horizon, "CoordCenterInertial")
    return interp_cols(t, y[:,:3], [tq])[0]


def snapshot_state(h5, dump, horizon, target_time, basis):
    r_index = tdm_index(dump, horizon, "DimlessRicciScalar.tdm")
    i, tsnap, mismatch = nearest_snapshot_idx(r_index, target_time)
    if mismatch > MAX_DUMP_TIME_MISMATCH:
        raise RuntimeError(
            f"{horizon}: nearest dump snapshot mismatch {mismatch:.6f} M "
            f"> {MAX_DUMP_TIME_MISMATCH}"
        )

    # Structural equality of field grids at this snapshot.
    R, exR = read_scalar_snapshot(dump, horizon, "DimlessRicciScalar.tdm", i)
    E, exE = read_scalar_snapshot(dump, horizon, "WeylE_NN.tdm", i)
    B, exB = read_scalar_snapshot(dump, horizon, "WeylB_NN.tdm", i)
    C, exC = read_coord_snapshot(dump, horizon, i)

    if not (exR == exE == exB == exC):
        raise RuntimeError(
            f"{horizon}: scalar/coordinate extents disagree at snapshot {i}: "
            f"{exR}, {exE}, {exB}, {exC}"
        )

    w, coord_area = surrogate_vertex_weights(C, exR)

    Mirr = interp_horizon_scalar(h5, horizon, "ArealMass", tsnap)
    Mch = interp_horizon_scalar(h5, horizon, "ChristodoulouMass", tsnap)
    center = interp_center(h5, horizon, tsnap)

    fields = dimensionless_fields(R, E, B, Mch)

    expected_R_mean = (Mch*Mch) / (2.0*Mirr*Mirr)
    observed_R_mean = weighted_mean(fields["Rstar"], w)
    gb_resid = observed_R_mean / expected_R_mean - 1.0

    centroid = np.sum(w[:,None] * C, axis=0)
    centroid_error = float(np.linalg.norm(centroid - center))

    # H2R2 guard-semantics repair:
    # These are diagnostics, not hard gates.  Gauss-Bonnet uses intrinsic
    # dA, whereas `w` is only the frozen coordinate-area surrogate; and the
    # catalog coordinate center is not defined as this surrogate centroid.
    gb_legacy_pass = bool(abs(gb_resid) <= MAX_GB_ABS_RESIDUAL)
    centroid_legacy_pass = bool(centroid_error <= MAX_CENTROID_ERROR_M)

    vecs = []
    diag = {
        "target_time_M": float(target_time),
        "snapshot_time_M": float(tsnap),
        "snapshot_time_mismatch_M": float(mismatch),
        "extent0": int(exR[0]),
        "extent1": int(exR[1]),
        "coordinate_area_surrogate": coord_area,
        "M_irr": Mirr,
        "M_ch": Mch,
        "GB_residual": float(gb_resid),
        "GB_legacy_0p05_pass": gb_legacy_pass,
        "centroid_error_M": centroid_error,
        "centroid_legacy_0p05M_pass": centroid_legacy_pass,
    }

    for name in ["Rstar", "Estar", "Bstar"]:
        feat, fd = angular_moments(fields[name], C, center, basis, w)
        vecs.append(feat)
        for j, label in enumerate(
            ["dip_x","dip_y","dip_z","Q_xx","Q_yy","Q_xy","Q_xz","Q_yz"]
        ):
            diag[f"{name}_{label}"] = float(feat[j])
        diag[f"{name}_weighted_mean"] = fd["weighted_mean"]
        diag[f"{name}_weighted_rms_centered"] = fd["weighted_rms_centered"]

    return np.concatenate(vecs), diag


def simulation_paths(data_root, sid, version, lev):
    dname = sid.replace(":", "_")
    base = Path(data_root) / dname / version / f"Lev{int(lev)}"
    return base / "Horizons.h5", base / "HorizonsDump.h5"


def standardized_endpoint_geometry(X):
    X = np.asarray(X, float)
    mu = X.mean(axis=0)
    sd = X.std(axis=0, ddof=0)
    active = sd > 1e-6
    if active.sum() < 2:
        raise RuntimeError("Too few active spatial-moment features after standardization")
    Z = (X[:,active] - mu[active]) / sd[active]
    return Z, active


def path_length(Z):
    return float(np.sum(np.linalg.norm(np.diff(Z, axis=0), axis=1)))


def pca_pc1(Z):
    Zc = Z - Z.mean(axis=0)
    u, s, vt = np.linalg.svd(Zc, full_matrices=False)
    pc1 = Zc @ vt[0]
    frac = float(s[0]**2 / max(np.sum(s*s), np.finfo(float).tiny))
    return pc1, frac


def permutation_order_stats(Z, rng):
    obs_path = path_length(Z)
    pc1, frac = pca_pc1(Z)
    idx = np.arange(len(Z))
    if spearmanr is not None:
        obs_rho = float(abs(spearmanr(pc1, idx).statistic))
    else:
        # ranks are enough here; orbit index ranks are 0..N-1
        r = pd.Series(pc1).rank(method="average").to_numpy(float)
        obs_rho = float(abs(np.corrcoef(r, idx)[0,1]))

    null_path = np.empty(PERM_TRIALS, float)
    null_rho = np.empty(PERM_TRIALS, float)

    for b in range(PERM_TRIALS):
        p = rng.permutation(len(Z))
        null_path[b] = path_length(Z[p])
        if spearmanr is not None:
            null_rho[b] = abs(float(spearmanr(pc1, p).statistic))
        else:
            null_rho[b] = abs(float(np.corrcoef(pd.Series(pc1).rank().to_numpy(), p)[0,1]))

    p_path = (1 + np.sum(null_path <= obs_path)) / (PERM_TRIALS + 1)
    p_rho = (1 + np.sum(null_rho >= obs_rho)) / (PERM_TRIALS + 1)

    return {
        "ordered_path_length": obs_path,
        "path_p_lower": float(p_path),
        "abs_pc1_spearman": obs_rho,
        "pc1_p_upper": float(p_rho),
        "pc1_fraction": frac,
    }


def exact_signflip_p_one_sided_positive(values):
    """
    H1 alternative: mean(values) > 0.
    Enumerate all 2^n sign flips.  n=6 for the independent A4 panel.
    """
    v = np.asarray(values, float)
    obs = float(np.sum(v))
    vals = []
    for signs in itertools.product([-1.0, 1.0], repeat=len(v)):
        vals.append(float(np.sum(v*np.asarray(signs))))
    vals = np.asarray(vals)
    return float(np.mean(vals >= obs - 1e-15))


def raw_distance(a, b):
    a = np.asarray(a, float)
    b = np.asarray(b, float)
    return float(np.sqrt(np.mean((a-b)**2)))


def plot_ratios(intervals, outpath):
    try:
        import matplotlib.pyplot as plt
    except Exception:
        return False

    fig, ax = plt.subplots(figsize=(11,6))
    sims = list(dict.fromkeys(intervals["sxs_id"].tolist()))
    for sid in sims:
        g = intervals[intervals["sxs_id"] == sid]
        ax.plot(g["interval_index"], g["return_ratio"], marker="o", label=sid.split(":")[-1])
    ax.axhline(1.0, linestyle="--", label="same-phase = half-phase control")
    ax.set_xlabel("Matched one-orbit interval index")
    ax.set_ylabel("return ratio = d_full / mean(d_half1,d_half2)")
    ax.set_title("GW-X-A4-4E2H2R2: same-phase recurrence ratio")
    ax.legend(ncol=4, fontsize=8)
    fig.tight_layout()
    fig.savefig(outpath, dpi=180)
    plt.close(fig)
    return True


def plot_panel(per_sim, outpath):
    try:
        import matplotlib.pyplot as plt
    except Exception:
        return False
    g = per_sim[per_sim["role"] == "independent_panel"].copy()
    x = np.arange(len(g))
    fig, ax = plt.subplots(figsize=(10,5))
    ax.bar(x, g["median_return_ratio"].to_numpy(float))
    ax.axhline(1.0, linestyle="--")
    ax.set_xticks(x)
    ax.set_xticklabels([s.split(":")[-1] for s in g["sxs_id"]])
    ax.set_ylabel("Median return ratio")
    ax.set_title("GW-X-A4-4E2H2R2: independent-panel recurrence")
    fig.tight_layout()
    fig.savefig(outpath, dpi=180)
    plt.close(fig)
    return True


def classify(recurrence_supported, ordered_count):
    nonreturn_replicated = ordered_count >= 4
    if recurrence_supported and nonreturn_replicated:
        return "A4_4E2H2R2_PHASE_RECURRENCE_WITH_ORDERED_FIELD_NONRETURN"
    if recurrence_supported:
        return "A4_4E2H2R2_PHASE_RECURRENCE_WITHOUT_ORDERED_FIELD_NONRETURN"
    if nonreturn_replicated:
        return "A4_4E2H2R2_ORDERED_FIELD_NONRETURN_WITHOUT_CLEAR_PHASE_RECURRENCE"
    return "A4_4E2H2R2_NO_CLEAR_SURFACE_FIELD_ARCHITECTURE"


def run(args):
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    h1 = load_json(args.h1_summary)
    if h1.get("classification") != EXPECTED_H1:
        raise RuntimeError(
            f"H1 guard failed: expected {EXPECTED_H1}, "
            f"got {h1.get('classification')}"
        )

    a43r = load_json(args.a4_3r_summary)
    if a43r.get("decision", {}).get("classification") != EXPECTED_A43R:
        raise RuntimeError(
            "A4-3R guard failed: independent dump-equipped panel is not frozen/ready"
        )
    ref_guard = a43r.get("guards", {}).get("reference_extractor_guard", {})
    if not ref_guard.get("passed", False):
        raise RuntimeError("A4-3R reference extractor guard did not pass")
    if ref_guard.get("selected_diagonal") != EXPECTED_DIAGONAL:
        raise RuntimeError(
            f"Frozen triangulation mismatch: expected {EXPECTED_DIAGONAL}, "
            f"got {ref_guard.get('selected_diagonal')}"
        )

    matched_panel = load_matched_csv(args.matched10)
    matched_ref = load_matched_csv(args.reference_states, forced_sid="SXS:BBH:0305")

    expected_ids = {x[0] for x in PANEL}
    if set(matched_panel.keys()) != expected_ids:
        raise RuntimeError(
            f"Matched panel IDs mismatch.\nexpected={sorted(expected_ids)}\n"
            f"got={sorted(matched_panel.keys())}"
        )

    sims = []
    for sid, version, lev, stratum in PANEL:
        hp, dp = simulation_paths(args.data_root, sid, version, lev)
        sims.append({
            "sxs_id": sid,
            "version": version,
            "lev": lev,
            "stratum": stratum,
            "role": "independent_panel",
            "horizons": hp,
            "dump": dp,
            "times": matched_panel[sid],
        })

    sims.append({
        "sxs_id": "SXS:BBH:0305",
        "version": "reference_local",
        "lev": 6,
        "stratum": "reference_0305",
        "role": "reference_only",
        "horizons": Path(args.reference_horizons),
        "dump": Path(args.reference_dump),
        "times": matched_ref["SXS:BBH:0305"],
    })

    for s in sims:
        if not s["horizons"].exists():
            raise FileNotFoundError(s["horizons"])
        if not s["dump"].exists():
            raise FileNotFoundError(s["dump"])

    all_state_rows = []
    interval_rows = []
    sim_rows = []
    errors = []
    rng = np.random.default_rng(RNG_SEED)

    for k, sim in enumerate(sims, 1):
        sid = sim["sxs_id"]
        print(f"[{k}/7] {sid}: extracting same-phase endpoints and half-phase controls")
        try:
            with h5py.File(sim["horizons"], "r") as hh, h5py.File(sim["dump"], "r") as hd:
                centers = common_centers_and_phase(hh)
                basis = fixed_frame(centers, sim["times"][0])
                mids, dphis = midpoint_times(centers, sim["times"])

                endpoint_vecs = []
                midpoint_vecs = []

                for n, tq in enumerate(sim["times"]):
                    parts = []
                    diag_join = {
                        "sxs_id": sid,
                        "role": sim["role"],
                        "stratum": sim["stratum"],
                        "kind": "same_phase",
                        "state_index": n,
                        "interval_index": np.nan,
                        "requested_time_M": float(tq),
                    }
                    for H in ["A","B"]:
                        v, d = snapshot_state(hh, hd, H, float(tq), basis)
                        parts.append(v)
                        for key, val in d.items():
                            diag_join[f"{H}_{key}"] = val
                    V = np.concatenate(parts)
                    endpoint_vecs.append(V)
                    for j, val in enumerate(V):
                        diag_join[f"shape_{j:02d}"] = float(val)
                    all_state_rows.append(diag_join)

                for n, tq in enumerate(mids):
                    parts = []
                    diag_join = {
                        "sxs_id": sid,
                        "role": sim["role"],
                        "stratum": sim["stratum"],
                        "kind": "half_phase",
                        "state_index": np.nan,
                        "interval_index": n,
                        "requested_time_M": float(tq),
                    }
                    for H in ["A","B"]:
                        v, d = snapshot_state(hh, hd, H, float(tq), basis)
                        parts.append(v)
                        for key, val in d.items():
                            diag_join[f"{H}_{key}"] = val
                    V = np.concatenate(parts)
                    midpoint_vecs.append(V)
                    for j, val in enumerate(V):
                        diag_join[f"shape_{j:02d}"] = float(val)
                    all_state_rows.append(diag_join)

                endpoint_vecs = np.asarray(endpoint_vecs)
                midpoint_vecs = np.asarray(midpoint_vecs)

                ratios = []
                effects = []
                for n in range(N_INTERVALS):
                    d_full = raw_distance(endpoint_vecs[n], endpoint_vecs[n+1])
                    d_h1 = raw_distance(endpoint_vecs[n], midpoint_vecs[n])
                    d_h2 = raw_distance(midpoint_vecs[n], endpoint_vecs[n+1])
                    half = 0.5*(d_h1+d_h2)
                    rr = d_full / max(half, np.finfo(float).tiny)
                    eff = half - d_full
                    ratios.append(rr)
                    effects.append(eff)
                    interval_rows.append({
                        "sxs_id": sid,
                        "role": sim["role"],
                        "stratum": sim["stratum"],
                        "interval_index": n,
                        "delta_phase_rad": float(dphis[n]),
                        "d_full_same_phase": d_full,
                        "d_half1": d_h1,
                        "d_half2": d_h2,
                        "d_half_mean": half,
                        "return_ratio": rr,
                        "recurrence_effect": eff,
                    })

                Z, active = standardized_endpoint_geometry(endpoint_vecs)
                order = permutation_order_stats(Z, rng)
                ordered_positive = (
                    order["path_p_lower"] <= 0.05
                    and order["pc1_p_upper"] <= 0.05
                )

                sim_rows.append({
                    "sxs_id": sid,
                    "role": sim["role"],
                    "stratum": sim["stratum"],
                    "median_return_ratio": float(np.median(ratios)),
                    "mean_return_ratio": float(np.mean(ratios)),
                    "mean_recurrence_effect": float(np.mean(effects)),
                    "fraction_intervals_return_ratio_lt_1": float(np.mean(np.asarray(ratios) < 1.0)),
                    "active_shape_features": int(active.sum()),
                    **order,
                    "ordered_nonreturn_positive": bool(ordered_positive),
                })

        except Exception as e:
            errors.append({
                "sxs_id": sid,
                "role": sim["role"],
                "error": repr(e),
            })
            print(f"  ERROR: {e}")

    errors_df = pd.DataFrame(errors)
    errors_path = out / "GW150914_GW_XA4_4E2H2R2_errors.csv"
    errors_df.to_csv(errors_path, index=False)

    if errors:
        summary = {
            "status": "GW-X-A4-4E2H2R2 surface-field phase-recurrence / secular-nonreturn audit",
            "classification": "A4_4E2H2R2_QUALITY_INCOMPLETE",
            "completed_simulations": len(sim_rows),
            "expected_simulations": 7,
            "errors": len(errors),
            "errors_csv": str(errors_path),
        }
        sp = out / "GW150914_GW_XA4_4E2H2R2_summary.json"
        sp.write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
        print(json.dumps(summary, ensure_ascii=False, indent=2))
        return

    states = pd.DataFrame(all_state_rows)
    intervals = pd.DataFrame(interval_rows)
    per_sim = pd.DataFrame(sim_rows)

    states_path = out / "GW150914_GW_XA4_4E2H2R2_surface_field_states.csv"
    intervals_path = out / "GW150914_GW_XA4_4E2H2R2_interval_recurrence.csv"
    per_sim_path = out / "GW150914_GW_XA4_4E2H2R2_per_simulation.csv"
    states.to_csv(states_path, index=False)
    intervals.to_csv(intervals_path, index=False)
    per_sim.to_csv(per_sim_path, index=False)

    panel_sim = per_sim[per_sim["role"] == "independent_panel"].copy()
    effects = panel_sim["mean_recurrence_effect"].to_numpy(float)
    signflip_p = exact_signflip_p_one_sided_positive(effects)
    n_ratio_lt1 = int(np.sum(panel_sim["median_return_ratio"] < 1.0))
    recurrence_supported = (signflip_p <= 0.05 and n_ratio_lt1 >= 4)

    ordered_count = int(np.sum(panel_sim["ordered_nonreturn_positive"]))
    classification = classify(recurrence_supported, ordered_count)

    ratio_plot = out / "GW150914_GW_XA4_4E2H2R2_return_ratio.png"
    panel_plot = out / "GW150914_GW_XA4_4E2H2R2_panel_recurrence.png"
    plot_ratios(intervals, ratio_plot)
    plot_panel(per_sim, panel_plot)

    gb_cols = [c for c in states.columns if c.endswith("_GB_residual")]
    cent_cols = [c for c in states.columns if c.endswith("_centroid_error_M")]
    mismatch_cols = [c for c in states.columns if c.endswith("_snapshot_time_mismatch_M")]

    summary = {
        "status": "GW-X-A4-4E2H2R2 surface-field phase-recurrence / secular-nonreturn audit",
        "classification": classification,
        "trigger": {
            "H1": h1.get("classification"),
            "reason": "Horizons-only one-orbit mass dynamics were below the frozen H0 resolution floor",
        },
        "scope": {
            "independent_panel_simulations": 6,
            "reference_simulations": 1,
            "matched_same_phase_states_per_simulation": 10,
            "half_phase_controls_per_simulation": 9,
            "surface_fields": ["Rstar", "Estar", "Bstar"],
            "primary_shape_dimension": 48,
            "true_intrinsic_multipoles_computed": False,
            "physical_spacetime_metric_measured": False,
            "Tor_APK_tested": False,
        },
        "frozen_surface_method": {
            "measure": "coordinate-area surrogate; F_periodic1; fixed diag_10_01",
            "angular_frame": "single coordinate frame frozen at first matched state per simulation",
            "shape_normalization": "per snapshot, per horizon, per field: subtract surrogate weighted mean and divide by surrogate weighted centered RMS",
            "moments_per_field": [
                "dip_x","dip_y","dip_z",
                "Q_xx","Q_yy","Q_xy","Q_xz","Q_yz"
            ],
            "reason_not_spherical_harmonics": "public products do not provide evolving intrinsic local dA and an intrinsic angular frame",
        },
        "quality": {
            "legacy_reference_GB_abs_residual": MAX_GB_ABS_RESIDUAL,
            "GB_residual_is_hard_gate": False,
            "max_observed_GB_abs_residual": float(np.nanmax(np.abs(states[gb_cols].to_numpy(float)))) if gb_cols else None,
            "legacy_reference_centroid_error_M": MAX_CENTROID_ERROR_M,
            "centroid_error_is_hard_gate": False,
            "max_observed_centroid_error_M": float(np.nanmax(states[cent_cols].to_numpy(float))) if cent_cols else None,
            "max_allowed_dump_time_mismatch_M": MAX_DUMP_TIME_MISMATCH,
            "max_observed_dump_time_mismatch_M": float(np.nanmax(states[mismatch_cols].to_numpy(float))) if mismatch_cols else None,
            "max_allowed_orbit_phase_error_rad": MAX_ORBIT_PHASE_ERROR,
            "max_observed_orbit_phase_error_rad": float(np.nanmax(np.abs(intervals["delta_phase_rad"].to_numpy(float)-2*np.pi))),
            "errors": 0,
            "hard_quality_gates": [
                "frozen diag_10_01 extractor convention",
                "blocked_all_x_then_y_then_z coordinate storage convention",
                "required datasets and matching extents",
                "finite positive coordinate-area surrogate",
                "dump snapshot mismatch <= 1.0 M",
                "one-orbit phase error <= 0.5 rad"
            ],
            "diagnostic_only": [
                "coordinate-area centroid error vs CoordCenterInertial",
                "coordinate-area-surrogate Gauss-Bonnet residual"
            ],
        },
        "recurrence_test": {
            "definition": "same-phase full-orbit distance vs mean of the two half-phase distances",
            "return_ratio": "d_full / mean(d_half1,d_half2)",
            "panel_simulation_level_effect": "mean(mean(d_half1,d_half2)-d_full)",
            "exact_signflip_p_one_sided": signflip_p,
            "panel_simulations_with_median_return_ratio_lt_1": n_ratio_lt1,
            "required_count": 4,
            "supported": bool(recurrence_supported),
            "per_simulation_median_return_ratio": {
                r["sxs_id"]: float(r["median_return_ratio"])
                for _, r in panel_sim.iterrows()
            },
        },
        "ordered_nonreturn_test": {
            "permutation_trials_per_simulation": PERM_TRIALS,
            "seed": RNG_SEED,
            "positive_rule": "path_p_lower<=0.05 AND pc1_p_upper<=0.05",
            "independent_panel_positive_count": ordered_count,
            "required_for_replication": 4,
            "replicated": bool(ordered_count >= 4),
            "per_simulation": {
                r["sxs_id"]: {
                    "ordered_nonreturn_positive": bool(r["ordered_nonreturn_positive"]),
                    "path_p_lower": float(r["path_p_lower"]),
                    "pc1_p_upper": float(r["pc1_p_upper"]),
                    "abs_pc1_spearman": float(r["abs_pc1_spearman"]),
                    "pc1_fraction": float(r["pc1_fraction"]),
                }
                for _, r in panel_sim.iterrows()
            },
        },
        "interpretation": {
            "if_both": "the phase-locked spatial field pattern shows recurrence relative to half-phase controls, while the same-phase sequence still drifts in an ordered secular way",
            "if_recurrence_only": "a phase-linked recurrent pattern is resolved, but ordered same-phase secular drift is not replicated across the independent panel",
            "if_nonreturn_only": "ordered same-phase field evolution is replicated, but a distinct full-orbit recurrence advantage over half-phase controls is not clear",
            "if_neither": "this surrogate spatial-moment representation does not show a clear recurrence/nonreturn architecture",
        },
        "forbidden_claims": [
            "does not identify g_mu_nu",
            "does not prove a changing physical spacetime metric",
            "surrogate angular moments are not true intrinsic spherical-harmonic multipoles",
            "does not validate Tor-APK",
            "does not prove a fifth physical dimension",
            "does not rule out ordinary smooth GR",
            "does not eliminate SXS coordinate-gauge or numerical effects",
        ],
        "outputs": {
            "states_csv": str(states_path),
            "interval_recurrence_csv": str(intervals_path),
            "per_simulation_csv": str(per_sim_path),
            "return_ratio_png": str(ratio_plot),
            "panel_recurrence_png": str(panel_plot),
            "errors_csv": str(errors_path),
        }
    }

    sp = out / "GW150914_GW_XA4_4E2H2R2_summary.json"
    sp.write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")

    print(json.dumps({
        "classification": classification,
        "recurrence_supported": recurrence_supported,
        "recurrence_signflip_p": signflip_p,
        "panel_median_return_ratio_lt1_count": n_ratio_lt1,
        "ordered_nonreturn_positive_count": ordered_count,
        "output": str(out),
    }, ensure_ascii=False, indent=2))


def self_test():
    # 1) triangulation on a synthetic cylinder-like grid
    e0, e1 = 5, 8
    z = np.linspace(-1,1,e0)
    th = np.linspace(0,2*np.pi,e1,endpoint=False)
    coords = []
    for zz in z:
        for tt in th:
            coords.append([np.cos(tt), np.sin(tt), zz])
    coords = np.asarray(coords)
    w, A = surrogate_vertex_weights(coords, (e0,e1))
    assert np.isfinite(A) and A > 0
    assert abs(w.sum()-1) < 1e-12

    # 2) angular moments finite
    basis = np.eye(3)
    center = np.zeros(3)
    f = coords[:,0] + 0.2*coords[:,2]
    feat, d = angular_moments(f, coords, center, basis, w)
    assert feat.shape == (8,)
    assert np.all(np.isfinite(feat))

    # 3) exact sign-flip p is valid
    p = exact_signflip_p_one_sided_positive([1,2,3,4,5,6])
    assert 0 < p <= 0.05

    # 4) ordering statistic works
    X = np.column_stack([
        np.linspace(0,1,10),
        np.linspace(0,1,10)**2,
        np.linspace(0,1,10)**3,
    ])
    Z, active = standardized_endpoint_geometry(X)
    rng = np.random.default_rng(1)
    q = permutation_order_stats(Z, rng)
    assert np.isfinite(q["ordered_path_length"])
    assert np.isfinite(q["abs_pc1_spearman"])

    print("self-test: PASS")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--h1-summary", type=Path)
    ap.add_argument("--a4-3r-summary", type=Path)
    ap.add_argument("--matched10", type=Path)
    ap.add_argument("--reference-states", type=Path)
    ap.add_argument("--data-root", type=Path)
    ap.add_argument("--reference-horizons", type=Path)
    ap.add_argument("--reference-dump", type=Path)
    ap.add_argument("--out", type=Path)
    ap.add_argument("--self-test", action="store_true")
    args = ap.parse_args()

    if args.self_test:
        self_test()
    else:
        for name in [
            "h1_summary","a4_3r_summary","matched10","reference_states",
            "data_root","reference_horizons","reference_dump","out"
        ]:
            if getattr(args, name) is None:
                ap.error(f"--{name.replace('_','-')} is required")
        run(args)
