#!/usr/bin/env python3
"""Inventory already supplied files; no ringdown or multipole refits."""
from pathlib import Path
import hashlib
import json
import re
import zipfile
import h5py

BASE = Path(__file__).resolve().parent.parent
OUT = BASE/"hf4_a3/results"


def digest(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def main():
    par = BASE/"hf3_r7b/input/provenance/qc0-horizon-multipoles-lite-t105-recover.par"
    txt = par.read_text()
    keys = ["IOHDF5::out_every", "IOHDF5::out_vars", "AHFinderDirect::output_h_every", "IOASCII::out0D_every",
            "IOASCII::out1D_every", "IOASCII::out2D_every", "SphericalSurface::maxntheta", "SphericalSurface::maxnphi",
            "IO::checkpoint_on_terminate", "QuasiLocalMeasures::output_vtk_every"]
    params = {k: re.findall(r"^\s*"+re.escape(k)+r"\s*=\s*([^\n]+)", txt, re.M) for k in keys}
    h5 = BASE/"hf4_a1/input/Moments_500M_L30.h5"
    datasets = []
    with h5py.File(h5, "r") as f:
        def visit(n, d):
            if isinstance(d, h5py.Dataset):
                datasets.append(dict(name=n, shape=list(d.shape), dtype=str(d.dtype)))
        f.visititems(visit)
    arc = BASE/"upload/BH_HF3_R7A_output.zip"
    with zipfile.ZipFile(arc) as z:
        members = [n for n in z.namelist() if not n.endswith("/")]
    geom = [n for n in members if re.search(r"\.(h5|hdf5|vtk|vtp|vtu)$", n, re.I) or re.search(r"/h\.t\d+\.ah\d+\.gp$", n)]
    p = BASE/"hf3_r7b/input/t95_t105/sphericalsurface-sf_radius..asc"
    lines = p.read_text().splitlines()
    rows = [x.split() for x in lines if x.strip() and not x.startswith("#")]
    # Confirm the ASCII header before interpreting metadata columns.
    assert any("6:ix 7:iy 8:iz" in x and "10:x 11:y 12:z" in x for x in lines)
    evidence = dict(path=str(p), rows=len(rows), columns=sorted(set(map(len, rows))),
        column_format="1:it 2:tl 3:rl 4:c 5:ml 6:ix 7:iy 8:iz 9:time 10:x 11:y 12:z 13..18:sf_radius[0..5]",
        distinct_ijk=sorted(set(tuple(x[5:8]) for x in rows)), distinct_xyz=sorted(set(tuple(x[9:12]) for x in rows)),
        first_time=float(rows[0][8]), last_time=float(rows[-1][8]),
        interpretation="One angular grid sample per surface, not a full surface mesh")
    r = dict(date_utc="2026-09-16", task="HF4-A3 geometry inventory; no A1/A2/R7B physics recomputed",
        parameters_path=str(par), parameters_sha256=digest(par), parameters=params,
        Yitian_HDF5_path=str(h5), Yitian_HDF5_sha256=digest(h5), Yitian_datasets=datasets,
        QC0_archive_sha256=digest(arc), QC0_archive_members=len(members), QC0_archive_full_geometry_candidates=geom,
        QC0_existing_radius_series=evidence, reference_source_commit="edd7fe658cb6748cd9ac8efcea314214b87bc043",
        reference_is_not_user_local_binary_verification=True,
        QLM_surface_and_tetrad_storage_declared=True, QLM_twometric_explicit_checkpoint_no=True,
        actual_local_checkpoints_examined=False,
        source_checksums={str(p.relative_to(BASE/"hf4_a3")): digest(p) for p in sorted((BASE/"hf4_a3/reference").glob("*")) if p.is_file()})
    OUT.mkdir(exist_ok=True)
    (OUT/"BH_HF4_A3_available_data_audit.json").write_text(json.dumps(r, ensure_ascii=False, indent=2)+"\n")
    print(json.dumps(dict(archive_members=len(members), full_geometry_candidates=geom, radius_series=evidence), ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
