#!/usr/bin/env python3
"""HF4-A3: read existing QC0 geometry; never run or modify the simulation.

Default: python BH_HF4_A3_collect.py
Only standard library is required for inventory. numpy+h5py enable small
surface-array extraction. Large 3-D evolution datasets are metadata-only.
"""
from __future__ import annotations

import argparse
import collections
import datetime as dt
import hashlib
import json
import math
import os
from pathlib import Path
import re
import sys
import zipfile

VERSION = "HF4-A3-1.0"
MIB = 1024 * 1024
LIMITS = dict(files=50000, depth=7, checkpoints=80, hdf_objects=60000,
              dataset_bytes=2*MIB, array_bytes=96*MIB, copied_bytes=24*MIB,
              source_bytes=2*MIB, attr_bytes=16384)
CP_RE = re.compile(r"checkpoint.*?\.it_(\d+)\.file_(\d+)\.h5$", re.I)
SURFACE_THORNS = ("sphericalsurface::", "quasilocalmeasures::", "ahfinderdirect::")
Q_NAMES = set(("qlm_shape qlm_shape_p qlm_shape_p_p qlm_x qlm_y qlm_z "
               "qlm_x_p qlm_y_p qlm_z_p qlm_x_p_p qlm_y_p_p qlm_z_p_p "
               "qlm_qtt qlm_qtp qlm_qpp qlm_rsc qlm_time qlm_time_p qlm_time_p_p "
               "qlm_area qlm_radius qlm_radius_p qlm_radius_p_p qlm_mass qlm_spin "
               "qlm_calc_error qlm_have_valid_data qlm_have_valid_data_p "
               "qlm_have_valid_data_p_p qlm_iteration qlm_timederiv_order "
               "qlm_ntheta qlm_nphi qlm_nghoststheta qlm_nghostsphi").split())
Q_NAMES.update("qlm_"+prefix+str(i) for prefix in ("m", "l", "n") for i in range(4))
SOURCE_THORNS = {
    "CactusNumerical/SphericalSurface": ("interface.ccl", "param.ccl", "schedule.ccl", "src/setup.cc", "src/setup.c"),
    "EinsteinAnalysis/AHFinderDirect": ("interface.ccl", "param.ccl", "schedule.ccl"),
    "EinsteinAnalysis/QuasiLocalMeasures": ("interface.ccl", "param.ccl", "schedule.ccl", "src/qlm_tetrad1.F90", "src/qlm_gram_schmidt.F90", "src/qlm_twometric.F90", "src/qlm_derivs.F90", "src/qlm_set_coordinates.F90", "src/qlm_set_surface.F90", "src/qlm_boundary.F90"),
    "EinsteinBase/ADMBase": ("interface.ccl", "param.ccl"),
    "McLachlan/ML_BSSN": ("interface.ccl", "param.ccl"),
}


def sha(b):
    return hashlib.sha256(b).hexdigest()


def dump(path, obj):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, ensure_ascii=False, indent=2, allow_nan=False)+"\n", encoding="utf-8")


def jsonable(v):
    if hasattr(v, "tolist"):
        return jsonable(v.tolist())
    if isinstance(v, bytes):
        return v.decode("utf-8", errors="replace")
    if isinstance(v, complex):
        return {"real": jsonable(v.real), "imag": jsonable(v.imag)}
    if isinstance(v, float) and not math.isfinite(v):
        return str(v)
    if isinstance(v, (tuple, list)):
        return [jsonable(x) for x in v]
    if isinstance(v, dict):
        return {str(k): jsonable(x) for k, x in v.items()}
    if v is None or isinstance(v, (str, int, float, bool)):
        return v
    return str(v)


def stat_record(path):
    s = path.stat()
    return dict(bytes=s.st_size, mtime_ns=s.st_mtime_ns)


def variable_name(name):
    return re.split(r"[\s\[]", name.rsplit("/", 1)[-1].lower(), maxsplit=1)[0]


def is_surface_array(name):
    var = variable_name(name)
    if var.startswith("sphericalsurface::sf_"):
        return True
    if var.startswith("ahfinderdirect::ah_"):
        return not var.endswith("::ahmask")
    if var.startswith("quasilocalmeasures::"):
        q = var.split("::", 1)[1]
        return q in Q_NAMES or q.startswith(("qlm_origin_", "qlm_delta_"))
    return False


def small_attrs(obj):
    result = {}
    for i, key in enumerate(obj.attrs):
        if i >= 64:
            result["__remaining_attributes_omitted__"] = True
            break
        try:
            aid = obj.attrs.get_id(key)
            dtype, shape = aid.dtype, aid.shape
            size = math.prod(shape) * dtype.itemsize
            # Variable-length and reference attributes can hide unbounded reads.
            if size > LIMITS["attr_bytes"] or dtype.hasobject:
                result[key] = {"omitted": "large or variable-length attribute", "shape": list(shape), "dtype": str(dtype)}
            else:
                result[key] = jsonable(obj.attrs[key])
        except Exception as e:
            result[key] = {"error": str(e)}
    return result


def inspect_hdf(path, file_id, out, budget, np, h5py):
    """Metadata for all datasets; payload only for allowlisted small arrays."""
    record = dict(id=file_id, path=str(path), before=stat_record(path), datasets=[], groups=[],
                  arrays=[], large_evolution_data_read=False, errors=[])
    arrays = {}
    count = 0
    try:
        with h5py.File(path, "r") as hf:
            record["root_attrs"] = small_attrs(hf)
            def visit(name, obj):
                nonlocal count
                count += 1
                if count > LIMITS["hdf_objects"]:
                    record["truncated"] = True
                    return "limit"
                if isinstance(obj, h5py.Group):
                    if len(obj.attrs):
                        record["groups"].append(dict(name=name, attrs=small_attrs(obj)))
                    return None
                if not isinstance(obj, h5py.Dataset):
                    return None
                size = math.prod(obj.shape)*obj.dtype.itemsize if obj.shape is not None else 0
                row = dict(name=name, variable=variable_name(name), shape=list(obj.shape) if obj.shape is not None else None,
                           dtype=str(obj.dtype), logical_bytes=size, attrs=small_attrs(obj))
                record["datasets"].append(row)
                if not is_surface_array(name):
                    return None
                if obj.shape is None or obj.dtype.hasobject or obj.dtype.kind not in "biufcV":
                    row["payload_status"] = "unsupported_or_empty_dtype"
                    return None
                if obj.is_virtual or obj.external:
                    row["payload_status"] = "external_or_virtual_skipped"
                    return None
                if size > LIMITS["dataset_bytes"] or budget["arrays"] + size > LIMITS["array_bytes"]:
                    row["payload_status"] = "size_limit"
                    return None
                try:
                    data = np.asarray(obj[()])
                    key = "a%05d" % len(arrays)
                    arrays[key] = data
                    budget["arrays"] += data.nbytes
                    info = dict(key=key, dataset=name, sha256_raw_c_order=sha(data.tobytes(order="C")), bytes=data.nbytes)
                    if data.dtype.kind in "biufc":
                        info["finite_count"] = int(np.isfinite(data).sum())
                        info["element_count"] = int(data.size)
                        if data.size <= 32:
                            info["small_values"] = jsonable(data)
                    record["arrays"].append(info)
                    row["payload_status"] = "extracted"
                except Exception as e:
                    row["payload_status"] = "read_error"
                    row["error"] = str(e)
                return None
            hf.visititems(visit)  # HDF5 hard links only; do not follow external links.
    except Exception as e:
        record["errors"].append(str(e))
    try:
        record["after"] = stat_record(path)
        record["source_unchanged_by_size_mtime"] = record["before"] == record["after"]
    except OSError as e:
        record["source_unchanged_by_size_mtime"] = False
        record["errors"].append(str(e))
    record["usable_snapshot_candidate"] = bool(arrays) and not record["errors"] and record["source_unchanged_by_size_mtime"] and not record.get("truncated", False)
    if arrays:
        dest = out / "arrays" / (file_id+".npz")
        dest.parent.mkdir(exist_ok=True)
        np.savez_compressed(dest, **arrays)
        record["npz"] = dest.relative_to(out).as_posix()
        record["npz_sha256"] = sha(dest.read_bytes())
    dump(out / "hdf_metadata" / (file_id+".json"), record)
    return {k: v for k, v in record.items() if k not in ("datasets", "groups", "arrays", "root_attrs")} | {
        "dataset_count": len(record["datasets"]), "extracted_array_count": len(arrays),
        "extracted_variables": sorted({variable_name(x["dataset"]) for x in record["arrays"]}),
        "metadata_file": "hdf_metadata/"+file_id+".json"}


def scan_qc0(root, notes):
    """Only QC0 project directories; bounded traversal, no directory symlinks."""
    if not root.is_dir():
        notes.append("Каталог расчётов отсутствует: "+str(root))
        return []
    bases = [root] if root.name.lower().startswith("qc0") else sorted(
        p for p in root.iterdir() if p.name.lower().startswith("qc0") and p.is_dir())
    files = []
    for base in bases:
        for directory, dirs, names in os.walk(base, followlinks=False):
            rel = Path(directory).relative_to(base)
            dirs[:] = sorted(d for d in dirs if not d.startswith(".") and not (Path(directory)/d).is_symlink())
            if len(rel.parts) >= LIMITS["depth"]:
                if dirs:
                    notes.append("Ограничение глубины обхода: "+str(directory))
                dirs[:] = []
            for name in sorted(names):
                p = Path(directory)/name
                if not name.startswith(".") and p.is_file() and not p.is_symlink():
                    files.append(p)
                if len(files) >= LIMITS["files"]:
                    notes.append("Достигнут предел числа файлов; инвентаризация неполная.")
                    return files
    return files


def copy_small(path, out, label, budget, records, limit=None):
    try:
        before = stat_record(path)
        if before["bytes"] > (limit or LIMITS["source_bytes"]) or budget["copied"]+before["bytes"] > LIMITS["copied_bytes"]:
            records.append(dict(path=str(path), status="size_limit", bytes=before["bytes"]))
            return
        raw = path.read_bytes()
        dst = out / "project_files" / label
        dst.parent.mkdir(parents=True, exist_ok=True)
        dst.write_bytes(raw)
        budget["copied"] += len(raw)
        records.append(dict(path=str(path), output=dst.relative_to(out).as_posix(), sha256=sha(raw),
                            source_unchanged_by_size_mtime=(before == stat_record(path))))
    except OSError as e:
        records.append(dict(path=str(path), error=str(e)))


def collect(args):
    root = args.simulation_root.expanduser().absolute()
    cactus = args.cactus_root.expanduser().absolute()
    output_root = args.output_root.expanduser().absolute()
    # Results must not be written inside the simulation or source directories.
    for protected in (root, cactus):
        if output_root.resolve() == protected.resolve() or protected.resolve() in output_root.resolve().parents:
            raise ValueError("Каталог результата должен находиться вне каталогов расчёта и Cactus.")
    output_root.mkdir(parents=True, exist_ok=True)
    name = "BH_HF4_A3_inventory"
    if (output_root/name).exists() or (output_root/(name+".zip")).exists():
        name += "_"+dt.datetime.now(dt.timezone.utc).strftime("%Y%m%dT%H%M%S_%fZ")
    out = output_root/name
    out.mkdir()
    notes, copied, budget = [], [], {"arrays": 0, "copied": 0}
    try:
        import numpy as np
        import h5py
        hdf_support = dict(available=True, numpy=np.__version__, h5py=h5py.__version__)
    except ImportError as e:
        np = h5py = None
        hdf_support = dict(available=False, reason=str(e))
        notes.append("h5py/numpy недоступны: собран список файлов; массивы HDF5 не извлечены.")
    print("HF4-A3: читаю только существующие файлы QC0. Эволюция не запускается.", flush=True)
    files = scan_qc0(root, notes)
    checkpoints = [p for p in files if CP_RE.search(p.name)]
    geometry_files = [p for p in files if p not in checkpoints and (
        p.suffix.lower() in (".vtk", ".vtp", ".vtu") or re.match(r"h\.t\d+\.ah\d+\.(gp|h5)$", p.name, re.I)
        or (p.suffix.lower() in (".h5", ".hdf5") and any(x in p.name.lower() for x in ("sf_", "qlm_", "horizon", "ah_radius"))))]
    inventory = []
    for p in checkpoints+geometry_files:
        row = dict(path=str(p), **stat_record(p))
        m = CP_RE.search(p.name)
        if m:
            row.update(kind="checkpoint", iteration_from_filename=int(m[1]), part_from_filename=int(m[2]))
        else:
            row["kind"] = "possible_surface_output"
        inventory.append(row)
    print("Найдено частей checkpoint:", len(checkpoints), "; отдельных геометрических файлов:", len(geometry_files), flush=True)
    # Read every retained checkpoint up to the explicit cap, earlier iterations first.
    inspected = []
    selected = sorted(checkpoints, key=lambda p: (int(CP_RE.search(p.name)[1]), str(p)))
    selected += [p for p in geometry_files if p.suffix.lower() in (".h5", ".hdf5")]
    if len(selected) > LIMITS["checkpoints"]:
        notes.append("Достигнут предел чтения HDF5; остальные файлы перечислены без анализа.")
    if h5py is not None:
        for i, p in enumerate(selected[:LIMITS["checkpoints"]]):
            print("HDF5", i+1, "/", min(len(selected), LIMITS["checkpoints"]), ":", p.parent.name, p.name, flush=True)
            inspected.append(inspect_hdf(p, "hdf_%03d" % i, out, budget, np, h5py))
    # Preserve actual project configuration and small surface files, not raw checkpoints.
    extras = [p for p in files if p.suffix.lower() == ".par"]
    for base in (args.windows_root, args.windows_root/"Downloads"):
        if base.is_dir():
            extras.extend(p for p in base.iterdir() if p.is_file() and p.name.lower().startswith("qc0") and p.suffix.lower() == ".par")
    for i, p in enumerate(sorted(set(extras))):
        copy_small(p, out, "parameters/%03d_" % i+p.name, budget, copied)
    useful_asc = [p for p in files if (
        p.name.lower().startswith(("quasilocalmeasures-qlm_scalars.", "quasilocalmeasures-qlm_state.",
                                   "quasilocalmeasures-qlm_grid_", "sphericalsurface-sf_")) and p.name.endswith("..asc"))
        or (p.name.lower().startswith("bh_diagnostics.") and p.suffix.lower() == ".gp")]
    for i, p in enumerate(useful_asc):
        copy_small(p, out, "diagnostics/%04d_" % i+p.name, budget, copied, limit=512*1024)
    for i, p in enumerate(geometry_files):
        if p.suffix.lower() not in (".h5", ".hdf5"):
            copy_small(p, out, "surface_output/%04d_" % i+p.name, budget, copied)
    for thorn, names in SOURCE_THORNS.items():
        for filename in names:
            p = cactus/"arrangements"/thorn/filename
            if p.is_file():
                copy_small(p, out, "source/"+thorn+"/"+filename, budget, copied)
            else:
                copied.append(dict(path=str(p), status="not_found"))
    for relative in ("configs/sim/ThornList", "configs/sim/config-data/cctk_Config.h"):
        p = cactus/relative
        if p.is_file():
            copy_small(p, out, "build/"+p.name, budget, copied)
    cp_sets = collections.defaultdict(list)
    for row in inventory:
        if row["kind"] == "checkpoint":
            cp_sets[(str(Path(row["path"]).parent), row["iteration_from_filename"])].append(row["part_from_filename"])
    sets = [dict(directory=k[0], iteration=k[1], parts=sorted(v),
                 matches_previously_used_two_MPI_parts=sorted(v)==[0, 1]) for k, v in cp_sets.items()]
    summary = dict(stage="BH-HF4-A3", collector_version=VERSION,
        created_utc=dt.datetime.now(dt.timezone.utc).isoformat(), python=sys.version,
        simulation_root=str(root), cactus_root=str(cactus), hdf_support=hdf_support,
        limits=LIMITS, scanned_files=len(files), checkpoint_parts=len(checkpoints), checkpoint_sets=sets,
        inventory=inventory, inspected_hdf=inspected, project_files=copied,
        extracted_uncompressed_bytes=budget["arrays"], copied_bytes=budget["copied"], notes=notes,
        simulation_started=False, source_files_opened_for_write=False,
        checkpoint_payloads_in_bundle=False, large_evolution_arrays_read=False,
        interpretation="Наличие массивов ещё не подтверждает геометрическую пригодность: проверить slot, время QLM, флаги, сетку, дубликаты MPI и нормировку тетрады.")
    summary["status"] = ("HDF_INVENTORY_READY_FOR_REVIEW" if hdf_support["available"] and inspected
                         else "FILE_INVENTORY_ONLY_REVIEW_REQUIRED")
    dump(out/"BH_HF4_A3_inventory.json", summary)
    overview = ["HF4-A3 — инвентаризация существующих геометрических данных", "",
        "Статус: "+summary["status"], "Checkpoint-частей: "+str(len(checkpoints)),
        "HDF5-файлов проверено: "+str(len(inspected)),
        "Извлечено небольших массивов: "+str(sum(x["extracted_array_count"] for x in inspected)),
        "Эволюция не запускалась. Исходные файлы открывались только для чтения.",
        "Большие 3D-массивы не считывались и контрольные точки целиком не копировались.",
        "Извлечённые массивы требуют научной проверки перед реконструкцией.", ""] + notes
    (out/"READ_ME.txt").write_text("\n".join(overview)+"\n", encoding="utf-8")
    program = Path(__file__).resolve()
    if program.is_file():
        (out/"BH_HF4_A3_collect.py").write_bytes(program.read_bytes())
    manifest = {p.relative_to(out).as_posix(): sha(p.read_bytes()) for p in sorted(out.rglob("*")) if p.is_file()}
    dump(out/"MANIFEST_SHA256.json", manifest)
    archive = output_root/(name+".zip")
    with zipfile.ZipFile(archive, "x", zipfile.ZIP_DEFLATED, compresslevel=6) as z:
        for p in sorted(out.rglob("*")):
            if p.is_file():
                z.write(p, arcname=p.relative_to(out).as_posix())
    print("\nГотово. Пришлите только этот ZIP:\n"+str(archive), flush=True)
    if not hdf_support["available"]:
        print("HDF5 не прочитан: в окружении отсутствует numpy или h5py. Это отмечено в отчёте.", flush=True)
    return archive


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--simulation-root", type=Path, default=Path.home()/"simulations")
    ap.add_argument("--cactus-root", type=Path, default=Path.home()/"EinsteinToolkit_2026_05/Cactus")
    ap.add_argument("--windows-root", type=Path, default=Path("/mnt/c/Users/Admin"))
    ap.add_argument("--output-root", type=Path, default=Path("/mnt/c/Users/Admin") if Path("/mnt/c/Users/Admin").is_dir() else Path.home())
    args = ap.parse_args()
    try:
        collect(args)
    except Exception as e:
        print("Ошибка сбора: "+str(e), file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
