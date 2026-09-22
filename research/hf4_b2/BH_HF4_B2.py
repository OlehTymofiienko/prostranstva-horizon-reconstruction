#!/usr/bin/env python3
"""HF4-B2: descriptive localization and temporal comparisons, after B1.

No evolution, ringdown fits or B1 shape/statistics recalculation. Reads existing
VTK surfaces. Masks use weighted deviations, not the HF1 gradient algorithm.
"""
from pathlib import Path
import argparse
import csv
import hashlib
import json
import re

import numpy as np

HERE = Path(__file__).resolve().parent
CONFIG = {
    "stage": "BH-HF4-B2", "date_UTC": "2026-09-20",
    "purpose": "Descriptive localization on QC0 common apparent horizon; not a new detection test",
    "fields": ["repsi2", "impsi2", "abs_npsigma"],
    "threshold_quantiles": [0.90, 0.95, 0.98], "primary_display_quantile": 0.95,
    "deviation": "absolute deviation from weighted median of the current snapshot",
    "weights": "Euclidean coordinate triangle area, lumped to vertices; not proper horizon area",
    "region_definition": "abs(normalized coordinate)>=0.65 along each local axis; diagnostic regions can overlap",
    "local_axes": "u from coordinate area centroids AH1 to AH2, w=simulation z, v=w cross u; u projected to xy",
    "time_comparison": "same stored angular-grid indices, average normalized coordinate-area weights",
    "maximum_contiguous_gap_M": 0.0625001,
    "HF1_mask_reused": False, "Yitian_field_reused": False,
    "B1_metrics_recomputed": False,
    "confirmatory_or_preregistered": False,
}


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def read_vtk(path):
    text = path.read_text()
    m = re.search(r"POINTS\s+(\d+)\s+\w+\s*\n", text)
    n = int(m[1]); end = text.index("POLYGONS", m.end())
    points = np.fromstring(text[m.end():end], sep=" ").reshape(n, 3)
    m = re.search(r"POLYGONS\s+(\d+)\s+(\d+)\s*\n", text)
    end = text.index("POINT_DATA", m.end())
    cells = np.fromstring(text[m.end():end], sep=" ", dtype=np.int64).reshape(int(m[1]), 5)
    assert np.all(cells[:, 0] == 4)
    quads = cells[:, 1:]
    tris = np.concatenate((quads[:, [0, 1, 2]], quads[:, [0, 2, 3]]))
    fields = {}
    headers = list(re.finditer(r"SCALARS\s+(\S+)\s+\w+\s+1\s*\nLOOKUP_TABLE\s+default\s*\n", text))
    for i, h in enumerate(headers):
        end = headers[i+1].start() if i+1 < len(headers) else len(text)
        fields[h[1]] = np.fromstring(text[h.end():end], sep=" ")
        assert fields[h[1]].shape == (n,)
    assert np.all(np.isfinite(points)) and np.all((tris >= 0) & (tris < n))
    a, b, c = points[tris[:, 0]], points[tris[:, 1]], points[tris[:, 2]]
    triangle_area = np.linalg.norm(np.cross(b-a, c-a), axis=1)/2
    weights = np.zeros(n)
    for k in range(3):
        np.add.at(weights, tris[:, k], triangle_area/3)
    assert np.all(weights > 0)
    fields["abs_npsigma"] = np.hypot(fields["renpsigma"], fields["imnpsigma"])
    return dict(points=points, tris=tris, quads=quads, fields=fields, weights=weights,
                centroid=np.average(points, axis=0, weights=weights), sha256=sha(path), filename=path.name)


def wquantile(values, weights, q):
    idx = np.argsort(values, kind="stable")
    cumulative = np.cumsum(weights[idx])
    return float(values[idx[min(np.searchsorted(cumulative, q*cumulative[-1]), len(idx)-1)]])


def connected(mask, edges, weights):
    remaining = set(np.flatnonzero(mask).tolist()); groups = []
    while remaining:
        todo = [min(remaining)]; remaining.remove(todo[0]); component = []
        while todo:
            u = todo.pop(); component.append(u)
            for v in edges[u]:
                if v in remaining:
                    remaining.remove(v); todo.append(v)
        groups.append(component)
    return sorted(groups, key=lambda g: -weights[g].sum())


def write_csv(path, rows):
    if not rows:
        return
    with path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)


def analyze(input_dir, b1_csv, out):
    out.mkdir(parents=True, exist_ok=True)
    (out/"BH_HF4_B2_analysis_definition.json").write_text(json.dumps(CONFIG, indent=2)+"\n")
    # B1 numbers are used as frozen inputs to the viewer/report, never refitted.
    with b1_csv.open() as f:
        b1 = {int(r["iteration"]): r for r in csv.DictReader(f)}
    files = sorted(input_dir.glob("surface0[123]_*.vtk"))
    by_time = {}
    for path in files:
        hn, it = map(int, re.match(r"surface(\d+)_(\d+)", path.stem).groups())
        by_time.setdefault(it, {})[hn] = read_vtk(path)
    common_times = [t for t in sorted(by_time) if 3 in by_time[t]]
    birth = min(common_times)
    assert birth == 9536 and set(common_times) == set(b1)
    mask_rows, pair_rows, audits, viewer = [], [], [], []
    saved = {}; previous = {}; primary_masks = {}
    reference_tris = None
    for it, surfaces in sorted(by_time.items()):
        frame = dict(iteration=it, t_M=it/512, tau_M=(it-birth)/512, objects=[])
        for hn, s in surfaces.items():
            frame["objects"].append(dict(horizon=hn, xyz=np.round(s["points"], 8).ravel().tolist()))
            if reference_tris is None:
                reference_tris = s["tris"]
            assert np.array_equal(s["tris"], reference_tris)
            assert s["points"].shape == (2520, 3)
            audits.append(dict(file=s["filename"], sha256=s["sha256"], points=2520, triangles=len(s["tris"]),
                               finite_fields=all(np.isfinite(s["fields"][f]).all() for f in CONFIG["fields"])))
        if 3 not in surfaces:
            viewer.append(frame); continue
        s = surfaces[3]; P = s["points"]; W = s["weights"]; W = W/W.sum()
        dv = surfaces[2]["centroid"]-surfaces[1]["centroid"]
        dv[2] = 0; u = dv/np.linalg.norm(dv); w = np.array([0., 0., 1.]); v = np.cross(w, u)
        basis = np.stack([u, v, w], axis=1)
        local = (P-s["centroid"])@basis
        scale = np.max(np.abs(local), axis=0)
        normalized = local/scale
        edges = [set() for _ in P]
        for tri in s["tris"]:
            for a, b in ((tri[0],tri[1]), (tri[1],tri[2]), (tri[2],tri[0])):
                edges[a].add(int(b)); edges[b].add(int(a))
        frame.update(fields={}, masks={}, B1_span_ratio=float(b1[it]["AH3_coordinate_span_max_over_min"]),
                     axes=np.round(basis, 8).tolist())
        saved[f"points_{it}"] = P
        saved[f"weights_{it}"] = W
        for f in CONFIG["fields"]:
            F = s["fields"][f]
            assert np.isfinite(F).all()
            median = wquantile(F, W, .5); D = np.abs(F-median)
            frame["fields"][f] = np.round(F, 8).tolist()
            saved[f"field_{f}_{it}"] = F
            for q in CONFIG["threshold_quantiles"]:
                cutoff = wquantile(D, W, q)
                mask = D >= cutoff; groups = connected(mask, edges, W)
                mass = W[mask].sum()
                second = np.average(normalized[mask]**2, weights=W[mask], axis=0)
                row = dict(iteration=it, t_M=it/512, tau_M=(it-birth)/512, field=f, quantile=q,
                    median=median, deviation_cutoff=cutoff, selected_vertices=int(mask.sum()),
                    coordinate_area_fraction=float(mass), components=len(groups),
                    largest_component_fraction=float(W[groups[0]].sum()/mass),
                    long_axis_region_fraction=float(W[mask & (np.abs(normalized[:,0])>=.65)].sum()/mass),
                    transverse_region_fraction=float(W[mask & (np.abs(normalized[:,1])>=.65)].sum()/mass),
                    polar_region_fraction=float(W[mask & (np.abs(normalized[:,2])>=.65)].sum()/mass),
                    u_rms=float(np.sqrt(second[0])), v_rms=float(np.sqrt(second[1])), w_rms=float(np.sqrt(second[2])))
                mask_rows.append(row)
                if q == .95:
                    frame["masks"][f] = np.flatnonzero(mask).tolist()
                    saved[f"mask_{f}_{it}"] = mask
                    primary_masks[(it,f)] = mask
                    if f in previous:
                        pit, pF, pD, pM, pW = previous[f]
                        avgW = (pW+W)/2
                        def corr(x,y):
                            x=x-np.sum(x*avgW); y=y-np.sum(y*avgW)
                            return float(np.sum(avgW*x*y)/np.sqrt(np.sum(avgW*x*x)*np.sum(avgW*y*y)))
                        pair_rows.append(dict(field=f, iteration_from=pit, iteration_to=it, gap_M=(it-pit)/512,
                            contiguous=(it-pit)/512<=CONFIG["maximum_contiguous_gap_M"],
                            weighted_mask_Dice=float(2*avgW[pM&mask].sum()/(avgW[pM].sum()+avgW[mask].sum())),
                            signed_field_correlation=corr(pF,F), deviation_correlation=corr(pD,D)))
                    previous[f] = (it, F, D, mask, W)
        viewer.append(frame)
    saved["triangles"] = reference_tris
    np.savez_compressed(out/"BH_HF4_B2_fields_and_masks.npz", **saved)
    write_csv(out/"BH_HF4_B2_localization.csv", mask_rows)
    write_csv(out/"BH_HF4_B2_temporal_pairs.csv", pair_rows)
    summary = dict(stage="BH-HF4-B2", classification="QC0_SURFACE_FIELD_LOCALIZATION_AND_TWO_SHORT_TIME_BLOCKS_CHARACTERIZED",
                   definition=CONFIG, input_files=audits, common_iterations=common_times,
                   frozen_B1_metrics_sha256=sha(b1_csv), all_inputs_finite=all(x["finite_fields"] for x in audits),
                   same_mesh_connectivity=True, primary_field_summaries={},
                   gap_between_blocks_M=(10368-9696)/512,
                   physical_metric_reconstructed=False, HF1_structure_identified_in_QC0=False,
                   observed_common_time_interval_M=[min(common_times)/512,max(common_times)/512])
    for f in CONFIG["fields"]:
        rows=[r for r in mask_rows if r["field"]==f and r["quantile"]==.95]
        pairs=[r for r in pair_rows if r["field"]==f and r["contiguous"]]
        summary["primary_field_summaries"][f] = dict(
            long_axis_fraction_range=[min(r["long_axis_region_fraction"] for r in rows),max(r["long_axis_region_fraction"] for r in rows)],
            polar_fraction_range=[min(r["polar_region_fraction"] for r in rows),max(r["polar_region_fraction"] for r in rows)],
            mask_Dice_range=[min(r["weighted_mask_Dice"] for r in pairs),max(r["weighted_mask_Dice"] for r in pairs)],
            field_correlation_range=[min(r["signed_field_correlation"] for r in pairs),max(r["signed_field_correlation"] for r in pairs)],
            selected_area_range=[min(r["coordinate_area_fraction"] for r in rows),max(r["coordinate_area_fraction"] for r in rows)])
    (out/"BH_HF4_B2_summary.json").write_text(json.dumps(summary,indent=2,allow_nan=False)+"\n")
    ranges={}
    for f in CONFIG["fields"]:
        allF=np.concatenate([s[3]["fields"][f] for s in by_time.values() if 3 in s])
        mx=float(np.max(np.abs(allF)))
        ranges[f]=[0,mx] if f=="abs_npsigma" else [-mx,mx]
    db=dict(frames=viewer, triangles=reference_tris.tolist(), fields=CONFIG["fields"], ranges=ranges,
            nt=35,np=72,birth_iteration=birth)
    (out/"BH_HF4_B2_view_data.json").write_text(json.dumps(db,separators=(",",":"),allow_nan=False))
    print(json.dumps(summary["primary_field_summaries"],indent=2))
    return summary


if __name__ == "__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--input",type=Path,default=HERE.parent/"hf4_latest/full")
    p.add_argument("--b1",type=Path,default=HERE.parent/"hf4_latest/provenance/BH_HF4_B1_metrics.csv")
    p.add_argument("--output",type=Path,default=HERE/"results")
    a=p.parse_args();analyze(a.input,a.b1,a.output)
