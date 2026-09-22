#!/usr/bin/env python3
"""HF4-A1: dimensionless horizon scalar maps and a reconstruction audit.

No Einstein Toolkit evolution is run. Inputs are read-only.
Use original geometric moments, including stationary m=0 terms.
Requirements: Python >=3.10, numpy, scipy>=1.15, h5py, matplotlib.
"""
from __future__ import annotations

import argparse
import hashlib
import io
import json
import platform
import re
import zipfile
from pathlib import Path

import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
import numpy as np
from numpy.polynomial import Legendre, Polynomial
import scipy
from scipy.integrate import cumulative_trapezoid
from scipy.special import sph_harm_y

ELL = np.array([l for l in range(9) for m in range(-l, l+1)])
EM = np.array([m for l in range(9) for m in range(-l, l+1)])
EXPECTED_HDF_SHA = "9cdd268e85b783d8a832b72e46945c6b43bb2b3ea0a9e6665e57a5a1b7d56f0f"
EXPECTED_QC0_SHA = "d73d7bce0b88c22022ccc5e571431045370edccbf48ed9b03f86fb261afafa1b"
SNAPSHOTS = [0, 5, 10, 20, 40, 80, 150, 500]


def sha(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1048576), b""):
            h.update(b)
    return h.hexdigest()


def json_write(path, obj):
    path.write_text(json.dumps(obj, ensure_ascii=False, indent=2, allow_nan=False)+"\n")


def basis(theta, phi):
    """Columns ordered k=l(l+1)+m; theta is colatitude, phi azimuth."""
    theta, phi = np.broadcast_arrays(theta, phi)
    return np.column_stack([sph_harm_y(l, m, theta.ravel(), phi.ravel())
                            for l, m in zip(ELL, EM)])


def quadrature(n=32, nphi=64):
    x, wx = np.polynomial.legendre.leggauss(n)
    phi = np.arange(nphi)*2*np.pi/nphi
    y = basis(np.arccos(x)[:, None], phi[None, :])
    w = np.repeat(wx, nphi)*2*np.pi/nphi
    return x, phi, w, y


def pair_residual(a):
    return max(float(np.max(np.abs(a[:, l*(l+1)-m]-(-1)**m*a[:, l*(l+1)+m].conj())))
               for l in range(9) for m in range(l+1))


def coefficient_metrics(t, mass, spin):
    # Parseval norms on the normalized area sphere, not on an unknown embedding.
    data = {"time_M": t}
    for name, a in [("K", 2*mass), ("B", spin)]:
        p = np.abs(a)**2
        norm = np.sqrt(p.sum(axis=1))
        data[name+"_L2"] = norm
        data[name+"_nonaxisymmetric_L2"] = np.sqrt(p[:, EM != 0].sum(axis=1))
        data[name+"_nonaxisymmetric_fraction"] = data[name+"_nonaxisymmetric_L2"]/norm
        for cutoff in [2, 4, 6]:
            data[f"{name}_relative_difference_L8_L{cutoff}"] = np.sqrt(p[:, ELL > cutoff].sum(axis=1))/norm
        data[name+"_shell_L2"] = np.column_stack([np.sqrt(p[:, ELL == l].sum(axis=1)) for l in range(9)])
    return data


def kerr_fields(x, chi):
    b = chi/(1+np.sqrt(1-chi**2))
    d = 1+b*b*x*x
    f = (1+b*b)*(1-x*x)/d
    k = (1+b*b)**2*(1-3*b*b*x*x)/d**3
    B = 0.5*(1+b*b)**2*(3*b*x-b**3*x**3)/d**3
    return f, k, B


def f_from_K(k):
    # f''=-2K, f(-1)=0, f'(-1)=2, evaluated in a Legendre basis.
    j = k.integ(lbnd=-1)
    jj = j.integ(lbnd=-1)
    return Legendre([2, 2])-2*jj


def profile_from_f(f, n=1601):
    th = np.linspace(0, np.pi, n)
    x = np.cos(th)
    p = f.convert(kind=Polynomial)
    h, rem = divmod(p, Polynomial([1, 0, -1]))
    assert np.max(np.abs(rem.coef)) < 2e-12, "Pole conditions not satisfied"
    hh = h(x)
    rad = 1-f.deriv()(x)**2/4
    assert min(hh) > 0 and min(rad) > -2e-12, "No real surface-of-revolution embedding"
    rho = np.sin(th)*np.sqrt(hh)
    z = -cumulative_trapezoid(np.sqrt(np.maximum(rad, 0)/hh), th, initial=0)
    z -= z[n//2]
    return th, rho, z, hh, rad


def late_geometry(mass, spin, t):
    late = t >= 400
    mi, sl = mass[late].mean(axis=0), spin[late].mean(axis=0)
    b = float(sl[2].real/np.sqrt(3*np.pi))
    chi_fit = 2*b/(1+b*b)
    coeff = np.array([2*mi[l*(l+1)].real*np.sqrt((2*l+1)/(4*np.pi)) for l in range(9)])
    raw_coeff = coeff.copy()
    # Only this conditional axisymmetric metric enforces exact regularity.
    # Scalar maps elsewhere retain all original data, with no edits.
    coeff[0], coeff[1] = 1, 0
    kc = Legendre(coeff)
    f = f_from_K(kc)
    th, rho, z, h, rad = profile_from_f(f)
    x = np.cos(th)
    fk, kk, bk = kerr_fields(x, chi_fit)
    bfit = chi_fit/(1+np.sqrt(1-chi_fit**2))
    hk = (1+bfit*bfit)/(1+bfit*bfit*x*x)
    # Analytic Kerr f' and the same Euclidean embedding quadrature.
    fpk = -2*x*(1+bfit*bfit)**2/(1+bfit*bfit*x*x)**2
    zk = -cumulative_trapezoid(np.sqrt(np.maximum(1-fpk*fpk/4, 0)/hk), th, initial=0)
    zk -= zk[len(zk)//2]
    rhok = np.sin(th)*np.sqrt(hk)
    # A constructive counterexample to unique inversion from finite multipoles.
    epsilon = 0.2
    kc_alt = kc+epsilon*Legendre.basis(10)
    f_alt = f_from_K(kc_alt)
    _, rho_alt, z_alt, h_alt, rad_alt = profile_from_f(f_alt)
    qx, qw = np.polynomial.legendre.leggauss(64)
    differences = [np.pi*np.sqrt((2*l+1)/(4*np.pi))*np.dot(qw, (kc_alt(qx)-kc(qx))*Legendre.basis(l)(qx)) for l in range(9)]
    def poles(g):
        return {"f_minus1": float(g(-1)), "f_plus1": float(g(1)),
                "fp_minus1": float(g.deriv()(-1)), "fp_plus1": float(g.deriv()(1))}
    nonaxis_late = float(np.linalg.norm(mi[EM != 0])/np.linalg.norm(mi))
    summary = {
        "window_M": [400, 500], "sample_count": int(late.sum()),
        "chi_from_L10": float(chi_fit), "chi_email": 0.68644,
        "axisymmetric_approximation_only": True,
        "late_mean_mass_nonaxisymmetric_L2_fraction": nonaxis_late,
        "K_legendre_l0_regularization": float(coeff[0]-raw_coeff[0]),
        "K_legendre_l1_regularization": float(coeff[1]-raw_coeff[1]),
        "poles": poles(f),
        "min_f_divided_by_1_minus_x2": float(h.min()),
        "max_absolute_f_difference_Kerr_fitted_spin": float(np.max(np.abs(f(x)-fk))),
        "max_relative_f_difference_including_pole_limits": float(np.max(np.abs(h/hk-1))),
        "max_absolute_K_difference_Kerr_fitted_spin": float(np.max(np.abs(kc(x)-kk))),
        "max_profile_displacement_over_R": float(np.max(np.hypot(rho-rhok, z-zk))),
        "equatorial_rho_over_R": float(rho[len(rho)//2]),
        "polar_z_over_R": float(z[0]),
        "Euclidean_rotation_embedding_condition_pass": bool(np.min(rad) > -2e-12),
        "counterexample": {
            "label": "Mathematical alternative, not a measured mode or a simulated history",
            "compared_set": "Axisymmetric control: same I_l0 for l<=8 and all m!=0 identically zero, not a replacement of all measured early-time coefficients",
            "added_K": "0.2 P_10(x)", "same_area": True,
            "max_abs_change_I_l0_l_le_8": float(max(np.abs(differences))),
            "poles": poles(f_alt),
            "min_f_divided_by_1_minus_x2": float(h_alt.min()),
            "min_K_on_grid": float(kc_alt(x).min()),
            "max_K_change": float(np.max(np.abs(kc_alt(x)-kc(x)))),
            "max_f_change": float(np.max(np.abs(f_alt(x)-f(x)))),
            "max_profile_displacement_over_R": float(np.max(np.hypot(rho_alt-rho, z_alt-z))),
            "Euclidean_rotation_embedding_condition_pass": bool(np.min(rad_alt) > -2e-12),
            "higher_modes_constrained_by_data": False,
            "same_spin_moments_possible": "Keep the same curl two-form; both metrics have the same area form."
        },
        "absolute_radius_at_late_time_only_under_Kerr_and_email_mass": float(np.sqrt(2*0.95162**2*(1+np.sqrt(1-0.68644**2)))),
        "absolute_radius_time_series_present": False,
        "Kerr_comparison_is_internal_consistency_not_resolution_error": True
    }
    arrays = {"theta": th, "x": x, "K_L8": kc(x), "K_Kerr": kk,
              "f_L8": f(x), "f_Kerr": fk, "rho_over_R_L8": rho, "z_over_R_L8": z,
              "rho_over_R_Kerr": rhok, "z_over_R_Kerr": zk,
              "K_alternative": kc_alt(x), "f_alternative": f_alt(x),
              "rho_alternative": rho_alt, "z_alternative": z_alt}
    return summary, arrays


def read_ascii(text):
    header = next(line for line in text.splitlines() if line.startswith("# data columns:"))
    columns = {name: int(col)-1 for col, name in re.findall(r"(\d+):([^\s]+)", header)}
    return np.loadtxt(io.StringIO(text)), columns


def qc0_inventory(zip_path):
    with zipfile.ZipFile(zip_path) as z:
        names = z.namelist()
        sf, sh = read_ascii(z.read("t95_t105/sphericalsurface-sf_radius..asc").decode())
        mp, mh = read_ascii(z.read("t95_t105/quasilocalmeasures-qlm_multipole_moments..asc").decode())
        sc, ch = read_ascii(z.read("t95_t105/quasilocalmeasures-qlm_scalars..asc").decode())
        assert np.array_equal(mp[:, 8], sc[:, 8])
        m0 = mp[:, mh["qlm_mp_m0[2]"]]; m = sc[:, ch["qlm_mass[2]"]]
        j1 = mp[:, mh["qlm_mp_j1[2]"]]; j = sc[:, ch["qlm_spin[2]"]]
        area = sc[:, ch["qlm_area[2]"]]; r = sc[:, ch["qlm_radius[2]"]]
        par = z.read("provenance/qc0-horizon-multipoles-lite-t105-recover.par").decode()
        settings = {}
        for key in ["AHFinderDirect::output_h_every", "IOASCII::out1D_every", "IOASCII::out2D_every", "IOHDF5::out_every", "IOHDF5::out_vars"]:
            match = re.search(r"^\s*"+re.escape(key)+r"\s*=\s*([^\n]+)", par, re.M|re.I)
            settings[key] = match.group(1).strip() if match else "not found"
        return {"archive_sha256": sha(zip_path), "member_count": len(names),
            "archive_CRC_pass": z.testzip() is None,
            "common_horizon_slot": 2, "sf_radius_rows": len(sf), "sf_radius_columns": sf.shape[1],
            "sf_radius_unique_saved_grid_indices": np.unique(sf[:, 5:8], axis=0).astype(int).tolist(),
            "sf_radius_is_full_surface": False,
            "qlm_slot2_multipole_names": [s for s in mh if s.endswith("[2]")],
            "qlm_slot2_multipoles_are_81_complex_modes": False,
            "max_relative_M0_minus_qlm_mass": float(np.max(np.abs(m0/m-1))),
            "max_relative_J1_minus_qlm_spin": float(np.max(np.abs(j1/j-1))),
            "max_relative_area_minus_4piR2": float(np.max(np.abs(area/(4*np.pi*r*r)-1))),
            "time_window_checked_M": [float(sc[0, 8]), float(sc[-1, 8])],
            "qlm_two_metric_files": [n for n in names if n.startswith("t95_t105/") and "qlm_twometric" in n],
            "output_settings": settings,
            "pointwise_two_metric_in_returned_archive": False,
            "ET_reference_commit": "edd7fe658cb6748cd9ac8efcea314214b87bc043",
            "ET_reference_source": "https://bitbucket.org/einsteintoolkit/einsteinanalysis/src/edd7fe658cb6748cd9ac8efcea314214b87bc043/QuasiLocalMeasures/src/qlm_multipoles.F90",
            "exact_user_binary_commit_known": False,
            "QC0_and_Yitian_are_not_interchangeable_datasets": True}


def plot_maps(out, t, mass, spin):
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10})
    lon = np.linspace(-np.pi, np.pi, 181)
    lat = np.linspace(-np.pi/2, np.pi/2, 91)
    y = basis(np.pi/2-lat[:, None], lon[None, :])
    times = [0, 5, 20, 80]
    idx = [int(np.argmin(abs(t-v))) for v in times]
    maps = [(2*mass[idx]@y.T).real.reshape(4, len(lat), len(lon)),
            (spin[idx]@y.T).real.reshape(4, len(lat), len(lon))]
    fig, axes = plt.subplots(2, 4, figsize=(15.2, 7.0), subplot_kw={"projection": "mollweide"})
    fig.subplots_adjust(left=.035, right=.965, top=.85, bottom=.20, wspace=.18, hspace=.31)
    fig.suptitle("HF4-A1 · Восстановленные скалярные поля общего горизонта", fontsize=18, y=.975)
    fig.text(.5, .918, "Yitian / SpEC · исходные моменты до ℓ = 8 · авторские переносимые углы", ha="center", fontsize=12)
    fig.text(.027, .872, r"Кривизна  $\mathcal{K}=R^2K_G$", ha="left", fontsize=12, weight="bold")
    fig.text(.027, .523, r"Поле вращения  $\mathcal{B}=-R^2\,\mathrm{curl}\,\omega/2$", ha="left", fontsize=12, weight="bold")
    for row, field in enumerate(maps):
        lo, hi = float(field.min()), float(field.max())
        if row == 1:
            hi = max(abs(lo), abs(hi)); lo = -hi
        for col, tm in enumerate(times):
            ax = axes[row, col]
            pc = ax.pcolormesh(lon, lat, field[col], shading="auto", rasterized=True,
                              cmap="coolwarm" if row == 1 else "viridis", norm=Normalize(lo, hi))
            ax.grid(alpha=.22, color="#67738a", linewidth=.5)
            ax.set_xticklabels([]); ax.set_yticklabels([])
            ax.set_title(f"t = {tm} M", fontsize=11, pad=8)
        cb = fig.colorbar(pc, ax=list(axes[row, :]), orientation="vertical", fraction=.013, pad=.014)
        cb.ax.tick_params(labelsize=9)
    fig.text(.5, .103, "Цветовые шкалы одинаковы во всех кадрах каждой строки; поздние постоянные моменты сохранены.", ha="center", fontsize=11)
    fig.text(.5, .062, "Овалы — карты угловых координат. Их контур не изображает физическую форму горизонта.", ha="center", fontsize=11, color="#475569")
    fig.savefig(out/"BH_HF4_A1_scalar_maps.png", dpi=170)
    plt.close(fig)
    return {"times_M": times, "longitude": lon, "latitude": lat, "K": maps[0], "B": maps[1]}


def plot_checks(out, t, metrics, late):
    plt.rcParams.update({"axes.spines.top": False, "axes.spines.right": False, "axes.grid": True, "grid.alpha": .2})
    fig, ax = plt.subplots(2, 2, figsize=(12.6, 9.2), layout="constrained")
    fig.suptitle("HF4-A1 · Проверки реконструкции и её границы", fontsize=17)
    use = t <= 80
    for cutoff, color in [(2, "#c0841a"), (4, "#e36c45"), (6, "#2858a5")]:
        ax[0, 0].semilogy(t[use], 100*metrics[f"K_relative_difference_L8_L{cutoff}"][use], color=color, label=f"ℓ ≤ {cutoff} против ℓ ≤ 8")
    ax[0, 0].set(title="Влияние известных высших гармоник", xlabel="Время SpEC, M", ylabel="Различие карт кривизны по L², %")
    ax[0, 0].legend(fontsize=9)
    for name, label, color in [("K", "Кривизна", "#2858a5"), ("B", "Поле вращения", "#a13e71")]:
        ax[0, 1].semilogy(t[use], 100*metrics[name+"_nonaxisymmetric_fraction"][use], label=label, color=color)
    ax[0, 1].set(title="Доля неосесимметричных компонентов", xlabel="Время SpEC, M", ylabel="Доля нормы поля по L², %")
    ax[0, 1].legend(fontsize=9)
    for suffix, label, style, color in [("L8", "Реконструкция ℓ ≤ 8", "-", "#2858a5"), ("Kerr", "Керр, χ из L₁₀", "--", "#e36c45")]:
        rho, z = late["rho_over_R_"+suffix], late["z_over_R_"+suffix]
        ax[1, 0].plot(np.r_[rho, -rho[::-1]], np.r_[z, z[::-1]], style, color=color, label=label)
    ax[1, 0].set(title="Поздний осесимметричный контроль", xlabel="ρ / R", ylabel="z / R", aspect="equal")
    ax[1, 0].legend(fontsize=9, loc="lower center")
    ax[1, 1].plot(late["x"], late["K_L8"], color="#2858a5", label="Из известных ℓ ≤ 8")
    ax[1, 1].plot(late["x"], late["K_alternative"], color="#a13e71", label="Добавлен гипотетический 0,2 P₁₀")
    ax[1, 1].set(title="Разные кривизны с теми же ℓ ≤ 8", xlabel="x = cos θ", ylabel=r"Безразмерная кривизна $\mathcal{K}$")
    ax[1, 1].legend(fontsize=9)
    fig.savefig(out/"BH_HF4_A1_reconstruction_checks.png", dpi=170)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--moments", type=Path, required=True)
    ap.add_argument("--qc0-archive", type=Path, required=True)
    ap.add_argument("--out", type=Path, default=Path("results"))
    args = ap.parse_args()
    out = args.out; out.mkdir(parents=True, exist_ok=True)
    assert sha(args.moments) == EXPECTED_HDF_SHA, "Unexpected source HDF5"
    assert sha(args.qc0_archive) == EXPECTED_QC0_SHA, "Unexpected QC0 archive"
    with h5py.File(args.moments) as f:
        t = f["time_adaptive"][:]
        mass = f["geo_mass/Adaptive"][:].conj()
        spin = f["geo_spin/Adaptive"][:].conj()
        datasets = {}
        f.visititems(lambda name, obj: datasets.update({name: {"shape": list(obj.shape), "dtype": str(obj.dtype)}}) if isinstance(obj, h5py.Dataset) else None)
    assert t.shape == (5001,) and mass.shape == spin.shape == (5001, 81)
    assert np.all(np.isfinite(mass)) and np.all(np.isfinite(spin)) and np.all(np.diff(t)>0)
    idx = np.array([int(np.argmin(abs(t-v))) for v in SNAPSHOTS])
    x, phi, w, y = quadrature()
    gram = y.conj().T@(w[:, None]*y)
    kc, bc = 2*mass[idx]@y.T, spin[idx]@y.T
    # Use real maps, so the test also checks the physical reality convention.
    back_i = (kc.real*w)@y.conj()/2
    back_l = (bc.real*w)@y.conj()
    checks = {
        "quadrature": "32 Gauss-Legendre nodes in cos(theta), 64 uniform azimuths",
        "basis_Gram_max_abs_error": float(np.max(abs(gram-np.eye(81)))),
        "mass_reality_pair_max_abs_residual": pair_residual(mass),
        "spin_reality_pair_max_abs_residual": pair_residual(spin),
        "K_max_imaginary_part_snapshots": float(np.max(abs(kc.imag))),
        "B_max_imaginary_part_snapshots": float(np.max(abs(bc.imag))),
        "mass_roundtrip_max_abs_error": float(np.max(abs(back_i-mass[idx]))),
        "spin_roundtrip_max_abs_error": float(np.max(abs(back_l-spin[idx]))),
        "K_area_mean_max_abs_difference_from_one": float(np.max(abs(2*mass[:, 0].real/np.sqrt(4*np.pi)-1))),
        "B_area_mean_max_abs": float(np.max(abs(spin[:, 0]/np.sqrt(4*np.pi)))),
        "sphere_monopole_K_error": float(np.max(abs(2*np.sqrt(np.pi)*y[:, 0]-1))),
        "roundtrip_is_numerical_consistency_not_physical_accuracy": True,
        "higher_ell_error_bound_available": False,
        "new_evolution_run": False,
    }
    assert checks["basis_Gram_max_abs_error"] < 1e-12
    assert checks["mass_roundtrip_max_abs_error"] < 1e-11
    assert checks["spin_roundtrip_max_abs_error"] < 1e-11
    assert checks["K_area_mean_max_abs_difference_from_one"] < 1e-8
    metrics = coefficient_metrics(t, mass, spin)
    late_summary, late_arrays = late_geometry(mass, spin, t)
    inventory = qc0_inventory(args.qc0_archive)
    # Finer quadrature is used for approximate extrema/area fractions, not error bars.
    _, _, wf, yf = quadrature(96, 192)
    kf, bf = (2*mass[idx]@yf.T).real, (spin[idx]@yf.T).real
    k6 = (2*mass[idx, :49]@yf[:, :49].T).real
    snapshots = []
    for j, i in enumerate(idx):
        s = {"time_M": float(t[i]), "K_min_sampled": float(kf[j].min()), "K_max_sampled": float(kf[j].max()),
             "B_min_sampled": float(bf[j].min()), "B_max_sampled": float(bf[j].max()),
             "K_negative_area_fraction_quadrature": float(np.dot(wf, kf[j]<0)/(4*np.pi)),
             "K_L6_min_sampled": float(k6[j].min())}
        s.update({key: float(value[i]) for key, value in metrics.items() if key != "time_M" and value.ndim == 1})
        snapshots.append(s)
    maps = plot_maps(out, t, mass, spin)
    plot_checks(out, t, metrics, late_arrays)
    np.savez_compressed(out/"BH_HF4_A1_scalar_snapshots.npz", **maps)
    np.savez_compressed(out/"BH_HF4_A1_norms_and_late_profile.npz", **metrics, **{"late_"+k: v for k,v in late_arrays.items()})
    classification = "BH_HF4_A1_SCALAR_MAPS_RECONSTRUCTED_LATE_AXISYMMETRIC_CONTROL_PASSED_FINITE_MULTIPOLES_DO_NOT_FIX_FULL_GEOMETRY"
    summary = {"stage": "BH-HF4-A1", "date_UTC": "2026-09-16", "classification": classification,
        "sources": {"Yitian_HDF5_sha256": sha(args.moments), "QC0_R7A_sha256": sha(args.qc0_archive)},
        "environment": {"python": platform.python_version(), "numpy": np.__version__, "scipy": scipy.__version__, "h5py": h5py.__version__, "matplotlib": matplotlib.__version__},
        "source_datasets": datasets,
        "conventions": {"index": "k=l(l+1)+m", "complex_conjugation_of_raw_HDF5": True,
            "harmonics": "orthonormal complex Y_lm, Condon-Shortley, exp(+im phi)",
            "frame": "original author-transported angles, no added time-dependent phase rotation",
            "area_form": "dA/R(t)^2 = sin(theta) dtheta dphi",
            "K_definition": "R^2 K_G = 2 sum I_lm Y_lm",
            "B_definition": "-R^2 curl(omega)/2 = sum L_lm Y_lm",
            "background_subtracted_for_maps": False,
            "synthesis_cutoff": 8,
            "early_Psi2_identity_assumed": False},
        "checks": checks, "snapshots": snapshots, "late_geometry": late_summary, "QC0_inventory": inventory,
        "next_step": "HF4-A2: time-dependent scalar atlas, fixed scales, original and late-subtracted maps, cutoffs 4/6/8, frame audit. No inferred radial displacement or unique full geometry.",
        "not_established": ["Metric reconstruction of the early non-axisymmetric horizon", "Error bound from missing ell>8", "Absolute area history of Yitian horizon", "Complete merger spacetime", "Extra dimensions or new physics"]}
    json_write(out/"BH_HF4_A1_summary.json", summary)
    print(json.dumps({"classification": classification, "checks": checks, "late_geometry": late_summary,
                      "snapshot_t0": snapshots[0], "QC0_M0_mass_difference": inventory["max_relative_M0_minus_qlm_mass"]}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
