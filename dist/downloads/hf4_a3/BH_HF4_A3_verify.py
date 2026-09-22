#!/usr/bin/env python3
"""Synthetic verification only. No QC0/SpEC scientific results are recomputed."""
import argparse
import contextlib
import hashlib
import importlib.util
import io
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import zipfile

import h5py
import numpy as np

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("collector", HERE/"BH_HF4_A3_collect.py")
c = importlib.util.module_from_spec(spec)
spec.loader.exec_module(c)


def hashes(root):
    return {str(p.relative_to(root)): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in root.rglob("*") if p.is_file() and not p.is_symlink()}


def collector_checks():
    checks = {}
    with tempfile.TemporaryDirectory(prefix="hf4_a3_verify_") as td:
        base = Path(td)
        sims, cactus, win, outputs = [base/x for x in ("simulations", "Cactus", "Windows", "output")]
        project = sims/"qc0_lite_t105_recover"
        project.mkdir(parents=True)
        win.mkdir()
        (sims/"unrelated").mkdir()
        (sims/"unrelated/private.txt").write_text("Must not be collected")
        (project/"run.par").write_text("AHFinderDirect::output_h_every = 0\n")
        src = cactus/"arrangements/CactusNumerical/SphericalSurface/interface.ccl"
        src.parent.mkdir(parents=True)
        src.write_text("fixture source, not actual toolkit code\n")
        expected = np.arange(30, dtype=float).reshape(6, 5) + 1
        for part in (0, 1):
            with h5py.File(project/f"checkpoint.chkpt.it_53760.file_{part}.h5", "w") as h:
                h.attrs["time"] = 105.
                h.attrs["iteration"] = 53760
                d = h.create_dataset("SphericalSurface::sf_radius[2] it=53760 tl=0", data=expected)
                d.attrs["origin"] = np.array([0., 0.])
                h.create_dataset("QuasiLocalMeasures::qlm_m1[2] it=53760 tl=0", data=expected*(1+2j))
                h.create_dataset("QuasiLocalMeasures::qlm_m2[2] it=53760 tl=0", data=np.ones((3, 2), dtype=[("real", "f8"), ("imag", "f8")]))
                h.create_dataset("QuasiLocalMeasures::qlm_time[2] it=53760 tl=0", data=105.)
                h.create_dataset("QuasiLocalMeasures::qlm_have_valid_data[2] it=53760 tl=0", data=1)
                # Logical size 1 GiB, unallocated on disk; payload reads are forbidden below.
                h.create_dataset("ML_BSSN::gt11 it=53760 tl=0 rl=0", shape=(512, 512, 512), dtype="f8")
                h.create_dataset("QuasiLocalMeasures::qlm_m3[2] it=53760 tl=0", shape=(10,), dtype="f8", external=[("absent.bin", 0, 80)])
                h["foreign"] = h5py.ExternalLink("absent.h5", "/")
        (project/"checkpoint.chkpt.it_99999.file_0.h5").write_bytes(b"truncated file fixture")
        before = hashes(sims) | {"source": src.read_bytes().hex()}
        args = argparse.Namespace(simulation_root=sims, cactus_root=cactus, windows_root=win, output_root=outputs)
        reads = []
        original = h5py.Dataset.__getitem__
        def guarded(obj, key):
            if not c.is_surface_array(obj.name):
                raise AssertionError("Forbidden bulk dataset read: "+obj.name)
            reads.append(obj.name)
            return original(obj, key)
        h5py.Dataset.__getitem__ = guarded
        try:
            with contextlib.redirect_stdout(io.StringIO()):
                archive = c.collect(args)
        finally:
            h5py.Dataset.__getitem__ = original
        report_dir = archive.with_suffix("")
        summary = json.loads((report_dir/"BH_HF4_A3_inventory.json").read_text())
        meta = json.loads((report_dir/summary["inspected_hdf"][0]["metadata_file"]).read_text())
        arrays = np.load(report_dir/meta["npz"], allow_pickle=False)
        item = next(x for x in meta["arrays"] if "sf_radius" in x["dataset"])
        checks["surface_payload_exact"] = bool(np.array_equal(arrays[item["key"]], expected))
        checks["complex_and_compound_payload_preserved"] = any(arrays[k].dtype.kind == "c" for k in arrays) and any(arrays[k].dtype.fields for k in arrays)
        checks["inputs_byte_identical"] = before == (hashes(sims) | {"source": src.read_bytes().hex()})
        checks["no_bulk_dataset_getitem"] = all(c.is_surface_array(x) for x in reads) and bool(reads)
        checks["logical_gigabyte_metadata_only"] = any(x["logical_bytes"] == 1024**3 and "payload_status" not in x for x in meta["datasets"])
        checks["external_dataset_skipped"] = any(x.get("payload_status") == "external_or_virtual_skipped" for x in meta["datasets"])
        checks["truncated_checkpoint_reported"] = bool(summary["inspected_hdf"][-1]["errors"]) and not summary["inspected_hdf"][-1]["usable_snapshot_candidate"]
        checks["MPI_parts_recorded"] = summary["checkpoint_sets"][0]["matches_previously_used_two_MPI_parts"]
        with zipfile.ZipFile(archive) as z:
            checks["no_raw_checkpoint_in_zip"] = not any("checkpoint.chkpt" in x for x in z.namelist())
            checks["unrelated_project_not_collected"] = not any(b"Must not be collected" in z.read(n) for n in z.namelist())
        old_hash = hashlib.sha256(archive.read_bytes()).hexdigest()
        old_limit = c.LIMITS["dataset_bytes"]
        c.LIMITS["dataset_bytes"] = 32
        try:
            with contextlib.redirect_stdout(io.StringIO()):
                archive2 = c.collect(args)
        finally:
            c.LIMITS["dataset_bytes"] = old_limit
        s2 = json.loads((archive2.with_suffix("")/"BH_HF4_A3_inventory.json").read_text())
        m2 = json.loads((archive2.with_suffix("")/s2["inspected_hdf"][0]["metadata_file"]).read_text())
        checks["large_surface_array_bound"] = any(x.get("payload_status") == "size_limit" for x in m2["datasets"])
        checks["previous_output_preserved"] = archive2 != archive and old_hash == hashlib.sha256(archive.read_bytes()).hexdigest()
        run = subprocess.run([sys.executable, "-S", str(HERE/"BH_HF4_A3_collect.py"),
            "--simulation-root", str(sims), "--cactus-root", str(cactus),
            "--windows-root", str(win), "--output-root", str(base/"no_h5py")], capture_output=True, text=True)
        no_h5 = json.loads((base/"no_h5py/BH_HF4_A3_inventory/BH_HF4_A3_inventory.json").read_text())
        checks["stdlib_fallback_honest"] = run.returncode == 0 and not no_h5["hdf_support"]["available"] and no_h5["status"] == "FILE_INVENTORY_ONLY_REVIEW_REQUIRED"
        args.output_root = project/"forbidden"
        try:
            c.collect(args)
            checks["output_inside_simulation_rejected"] = False
        except ValueError:
            checks["output_inside_simulation_rejected"] = not args.output_root.exists()
    assert all(checks.values()), checks
    return checks


def dyad_identity_checks():
    """Independent algebraic control in generic positive spatial metrics.

    E is the tangent-coordinate Jacobian (3x2). u,v are orthonormal with
    respect to gamma, as in the QLM Gram-Schmidt construction. E T = [u v].
    Therefore q^{-1}=T T^T. This is NOT q=E^T E (Euclidean rendering metric).
    """
    rng = np.random.default_rng(4303)
    errors, phases, projections = [], [], []
    for unused in range(400):
        E = rng.normal(size=(3, 2))
        G = rng.normal(size=(3, 3))
        gamma = G.T@G + np.eye(3)*.6
        known = E.T@gamma@E
        u = E[:, 1].copy()
        u /= np.sqrt(u@gamma@u)
        v = E[:, 0] - u*(u@gamma@E[:, 0])
        v /= np.sqrt(v@gamma@v)
        m = (u+1j*v)/np.sqrt(2.)
        def reconstruct(mm):
            ma = np.linalg.lstsq(E, mm, rcond=None)[0]
            inverse = 2*np.real(ma[:, None]*ma.conj()[None, :])
            return np.linalg.inv(inverse), np.linalg.norm(E@ma-mm)/np.linalg.norm(mm)
        q, projection = reconstruct(m)
        qrot, unused = reconstruct(m*np.exp(1j*rng.uniform(-np.pi, np.pi)))
        errors.append(np.linalg.norm(q-known)/np.linalg.norm(known))
        phases.append(np.linalg.norm(qrot-q)/np.linalg.norm(q))
        projections.append(projection)
    result = dict(cases=400, synthetic_only=True, max_relative_metric_error=max(errors),
                  max_phase_rotation_error=max(phases), max_tangent_projection_residual=max(projections))
    assert max(errors) < 1e-10 and max(phases) < 1e-10 and max(projections) < 1e-10, result
    return result


if __name__ == "__main__":
    r = dict(collector=collector_checks(), conditional_metric_identity=dyad_identity_checks())
    (HERE/"results").mkdir(exist_ok=True)
    (HERE/"results/BH_HF4_A3_verification.json").write_text(json.dumps(r, indent=2)+"\n")
    print(json.dumps(r, indent=2))
