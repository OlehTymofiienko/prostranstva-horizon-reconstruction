#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
BH-HF3-R6B — поздняя многорадиусная проверка до 95M.

Einstein Toolkit НЕ запускается.

Основная идея:
- проверить исходящее распространение на r=15,30,40,50,60M;
- проверить масштабирование амплитуды примерно как 1/r;
- сравнить одинаковый "возраст кольцевания" после локального максимума
  Psi4_(2,2), а не одинаковое координатное время;
- проверить переход от первого к второму керровскому периоду сразу на
  r=15,30,40M;
- проверить третий период на r=15,30M, где он уже доступен;
- проверить, перестаёт ли жёсткая Kerr-(2,2,0) модель быть существенно
  хуже свободной одно-модовой модели на третьем периоде.

Важно:
этот этап является формализацией уже просмотренных данных до 95M,
поэтому это проверка устойчивости и внутренней согласованности,
а не независимая слепая предрегистрация.
"""

from __future__ import annotations

import argparse
import csv
import io
import json
import math
import re
import zipfile
from pathlib import Path

import numpy as np

VERSION = "BH-HF3-R6B-2026-09-06"

DEFAULT_BUNDLE = "/mnt/c/Users/Admin/BH_HF3_R6A_audit_bundle.zip"
DEFAULT_OUTDIR = "/mnt/c/Users/Admin/BH_HF3_R6B_output"

SEGS = ("t0_t22", "t22_t30", "t30_t50", "t50_t70", "t70_t95")
RADII_ALL = [15.0, 30.0, 40.0, 50.0, 60.0]
RADII_TWO_PERIOD = [15.0, 30.0, 40.0]
RADII_THREE_PERIOD = [15.0, 30.0]

SLOT = 2
LATE_REMNANT_WINDOW = (90.0, 95.0)

BCW = {
    "f1": 1.5251,
    "f2": -1.1568,
    "f3": 0.1292,
    "q1": 0.7000,
    "q2": 1.4187,
    "q3": -0.4990,
}

GATES = {
    "adjacent_peak_delay_kerr_rstar_relative_error_max": 0.05,
    "r_times_peak_amplitude_relative_spread_max": 0.10,

    "period0_frequency_relative_spread_max": 0.03,
    "period0_damping_relative_spread_max": 0.20,

    "period1_frequency_relative_spread_max": 0.01,
    "period1_damping_relative_spread_max": 0.05,

    "period1_frequency_error_vs_kerr_max": 0.03,
    "period1_damping_error_vs_kerr_max": 0.20,
    "period1_phase_r2_min": 0.995,
    "period1_logamp_r2_min": 0.98,

    "period2_frequency_error_vs_kerr_max": 0.01,
    "period2_damping_error_vs_kerr_max": 0.05,
    "period2_phase_r2_min": 0.995,
    "period2_logamp_r2_min": 0.99,
    "period2_fixed_kerr_delta_bic_max": 10.0,

    "plus_minus_frequency_relative_difference_max": 0.005,

    "complex_shape_coherence_min": 0.995,
    "complex_shape_nrmse_max": 0.05,
}


def member_text(zf, names, suffix):
    hits = [n for n in names if n.endswith(suffix)]
    if len(hits) != 1:
        raise RuntimeError(
            f"Ожидался один файл с окончанием {suffix}, найдено {len(hits)}"
        )
    return zf.read(hits[0]).decode("utf-8", errors="replace")


def arr_from_text(text):
    a = np.loadtxt(io.StringIO(text), comments="#")
    if a.ndim == 1:
        a = a.reshape(1, -1)
    return a


def columns(text):
    for line in text.splitlines():
        if line.startswith("# data columns:"):
            return {
                name: int(i) - 1
                for i, name in re.findall(r"(\d+):([^\s]+)", line)
            }
    raise RuntimeError("Не найден заголовок '# data columns:'")


def merge_mode(zf, names, ell, m, radius):
    tt, zz = [], []
    last = -np.inf

    for seg in SEGS:
        text = member_text(
            zf, names,
            f"/{seg}/mp_Psi4_l{ell}_m{m}_r{radius:.2f}.asc"
        )
        a = arr_from_text(text)

        for row in a:
            ti = float(row[0])
            if ti > last + 1e-12:
                tt.append(ti)
                zz.append(complex(row[1], row[2]))
                last = ti

    return np.asarray(tt), np.asarray(zz)


def merge_scalar(zf, names, variable):
    tt, yy = [], []
    last = -np.inf

    for seg in SEGS:
        text = member_text(
            zf, names,
            f"/{seg}/quasilocalmeasures-qlm_scalars..asc"
        )
        a = arr_from_text(text)
        mp = columns(text)

        for row in a:
            ti = float(row[8])
            if ti > last + 1e-12:
                tt.append(ti)
                yy.append(float(row[mp[variable]]))
                last = ti

    return np.asarray(tt), np.asarray(yy)


def linfit(x, y):
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    X = np.column_stack([np.ones_like(x), x])
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    pred = X @ beta
    ssr = float(np.sum((y - pred) ** 2))
    sst = float(np.sum((y - np.mean(y)) ** 2))
    r2 = float(1.0 - ssr / sst) if sst > 0 else float("nan")
    return beta, r2


def complex_fit(t, z, start, end):
    mask = (t >= start - 1e-12) & (t <= end + 1e-12)
    tt = t[mask]
    zz = z[mask]

    if len(tt) < 5:
        return None

    phase = np.unwrap(np.angle(zz))
    amp = np.abs(zz)

    bph, r2ph = linfit(tt, phase)
    blog, r2amp = linfit(tt, np.log(np.maximum(amp, 1e-300)))

    return {
        "start_M": float(start),
        "end_M": float(end),
        "samples": int(len(tt)),
        "omega_signed_rad_per_M": float(bph[1]),
        "omega_abs_rad_per_M": float(abs(bph[1])),
        "alpha_per_M": float(-blog[1]),
        "phase_r2": float(r2ph),
        "log_amplitude_r2": float(r2amp),
        "phase_cycles": float(
            abs(phase[-1] - phase[0]) / (2.0 * math.pi)
        ),
        "amplitude_start": float(amp[0]),
        "amplitude_end": float(amp[-1]),
    }


def kerr_rstar(r, M, chi):
    a = chi * M
    root = math.sqrt(max(M * M - a * a, 0.0))
    rp = M + root
    rm = M - root

    return (
        r
        + (2.0 * M * rp / (rp - rm)) * math.log(abs(r - rp))
        - (2.0 * M * rm / (rp - rm)) * math.log(abs(r - rm))
    )


def best_complex_amplitude(tau, z, alpha, omega_signed):
    g = np.exp((-alpha + 1j * omega_signed) * tau)
    c = np.vdot(g, z) / np.vdot(g, g)
    pred = c * g
    rss = float(np.sum(np.abs(z - pred) ** 2))
    return c, rss


def bic(rss, k, n_real):
    rss = max(float(rss), 1e-300)
    return float(
        n_real * math.log(rss / n_real) + k * math.log(n_real)
    )


def relative_spread(values):
    a = np.asarray(values, float)
    med = float(np.median(a))
    return float(
        (np.max(a) - np.min(a)) / max(abs(med), 1e-300)
    )


def interp_complex(tnew, t, z):
    re_ = np.interp(tnew, t, np.real(z))
    im_ = np.interp(tnew, t, np.imag(z))
    return re_ + 1j * im_


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bundle", default=DEFAULT_BUNDLE)
    ap.add_argument("--outdir", default=DEFAULT_OUTDIR)
    args = ap.parse_args()

    bundle = Path(args.bundle)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    if not bundle.is_file():
        raise FileNotFoundError(bundle)

    with zipfile.ZipFile(bundle, "r") as zf:
        names = zf.namelist()

        # ----------------------------------------------------------
        # 1. Поздние параметры остаточной чёрной дыры.
        # ----------------------------------------------------------
        tM, mass = merge_scalar(zf, names, f"qlm_mass[{SLOT}]")
        tJ, spin = merge_scalar(zf, names, f"qlm_spin[{SLOT}]")

        if not np.array_equal(tM, tJ):
            raise RuntimeError("Сетки времени M и J не совпадают")

        late = (
            (tM >= LATE_REMNANT_WINDOW[0] - 1e-12)
            & (tM <= LATE_REMNANT_WINDOW[1] + 1e-12)
        )

        Mf = float(np.nanmedian(mass[late]))
        Jf = float(np.nanmedian(spin[late]))
        chi = float(Jf / (Mf * Mf))

        Mf_omega = (
            BCW["f1"] + BCW["f2"] * (1.0 - chi) ** BCW["f3"]
        )
        Q = (
            BCW["q1"] + BCW["q2"] * (1.0 - chi) ** BCW["q3"]
        )

        omega_kerr = float(Mf_omega / Mf)
        alpha_kerr = float(omega_kerr / (2.0 * Q))
        period_kerr = float(2.0 * math.pi / omega_kerr)

        # ----------------------------------------------------------
        # 2. Моды +/-2, максимумы и амплитуды на пяти радиусах.
        # ----------------------------------------------------------
        modes = {}
        peak_rows = []

        for r in RADII_ALL:
            tp, zp = merge_mode(zf, names, 2, +2, r)
            tm, zm = merge_mode(zf, names, 2, -2, r)

            if not np.array_equal(tp, tm):
                raise RuntimeError(
                    f"Сетки времени m=+2 и m=-2 не совпадают на r={r}"
                )

            sel = (tp >= 30.0) & (tp <= 95.0)
            inds = np.where(sel)[0]
            ip = inds[int(np.argmax(np.abs(zp[sel])))]

            peak_t = float(tp[ip])
            peak_amp = float(abs(zp[ip]))

            modes[r] = {
                "t": tp,
                "zp": zp,
                "zm": zm,
                "peak_t": peak_t,
                "peak_amp": peak_amp,
            }

            peak_rows.append({
                "radius_M": r,
                "peak_time_M": peak_t,
                "peak_abs_Psi4_22": peak_amp,
                "r_times_peak_abs_Psi4_22": float(r * peak_amp),
                "kerr_periods_after_peak_to_95":
                    float((95.0 - peak_t) / period_kerr),
            })

        # ----------------------------------------------------------
        # 3. Приращения времени максимума против Kerr r*.
        # ----------------------------------------------------------
        propagation_rows = []

        for r1, r2 in zip(RADII_ALL[:-1], RADII_ALL[1:]):
            obs = modes[r2]["peak_t"] - modes[r1]["peak_t"]
            pred = (
                kerr_rstar(r2, Mf, chi)
                - kerr_rstar(r1, Mf, chi)
            )
            rel = float(abs(obs - pred) / abs(pred))

            propagation_rows.append({
                "r1_M": r1,
                "r2_M": r2,
                "observed_peak_delay_increment_M": float(obs),
                "kerr_rstar_increment_M": float(pred),
                "relative_difference": rel,
            })

        propagation_pass = all(
            x["relative_difference"]
            <= GATES["adjacent_peak_delay_kerr_rstar_relative_error_max"]
            for x in propagation_rows
        )

        scaled_peak_values = [
            x["r_times_peak_abs_Psi4_22"] for x in peak_rows
        ]
        scaled_peak_spread = relative_spread(scaled_peak_values)

        inverse_radius_pass = bool(
            scaled_peak_spread
            <= GATES["r_times_peak_amplitude_relative_spread_max"]
        )

        # ----------------------------------------------------------
        # 4. Одинаковый возраст кольцевания.
        #
        # period_index=0: первый ожидаемый Kerr-период после максимума.
        # period_index=1: второй.
        # period_index=2: третий.
        # ----------------------------------------------------------
        age_fit_rows = []

        max_period_by_radius = {
            15.0: 3,
            30.0: 3,
            40.0: 2,
            50.0: 1,
            60.0: 0,
        }

        for r in RADII_ALL:
            t = modes[r]["t"]
            peak = modes[r]["peak_t"]

            for k in range(max_period_by_radius[r]):
                start = peak + k * period_kerr
                end = peak + (k + 1) * period_kerr

                if end > 95.0 + 1e-12:
                    continue

                for m, z in ((+2, modes[r]["zp"]), (-2, modes[r]["zm"])):
                    fit = complex_fit(t, z, start, end)
                    if fit is None:
                        continue

                    fit.update({
                        "radius_M": r,
                        "m": m,
                        "period_index": k,
                        "ringdown_age_start_M": float(k * period_kerr),
                        "ringdown_age_end_M": float((k + 1) * period_kerr),
                        "omega_relative_error_vs_kerr": float(
                            abs(fit["omega_abs_rad_per_M"] - omega_kerr)
                            / omega_kerr
                        ),
                        "alpha_relative_error_vs_kerr": float(
                            abs(fit["alpha_per_M"] - alpha_kerr)
                            / alpha_kerr
                        ),
                    })

                    age_fit_rows.append(fit)

        plus_rows = [x for x in age_fit_rows if x["m"] == +2]
        minus_rows = [x for x in age_fit_rows if x["m"] == -2]

        def get_fit(r, k, m=+2):
            rows = [
                x for x in age_fit_rows
                if x["radius_M"] == r
                and x["period_index"] == k
                and x["m"] == m
            ]
            if len(rows) != 1:
                raise RuntimeError(
                    f"Не найден однозначный fit: r={r}, k={k}, m={m}"
                )
            return rows[0]

        # ----------------------------------------------------------
        # 5. Симметрия +/-m.
        # ----------------------------------------------------------
        symmetry_rows = []

        for p in plus_rows:
            m = get_fit(p["radius_M"], p["period_index"], -2)
            rel = float(
                abs(
                    p["omega_abs_rad_per_M"]
                    - m["omega_abs_rad_per_M"]
                )
                / max(
                    0.5 * (
                        p["omega_abs_rad_per_M"]
                        + m["omega_abs_rad_per_M"]
                    ),
                    1e-300,
                )
            )
            symmetry_rows.append({
                "radius_M": p["radius_M"],
                "period_index": p["period_index"],
                "plus_minus_frequency_relative_difference": rel,
            })

        symmetry_pass = all(
            x["plus_minus_frequency_relative_difference"]
            <= GATES["plus_minus_frequency_relative_difference_max"]
            for x in symmetry_rows
        )

        # ----------------------------------------------------------
        # 6. Межрадиусная согласованность периодов.
        # ----------------------------------------------------------
        period0_radii = [15.0, 30.0, 40.0, 50.0]
        p0 = [get_fit(r, 0) for r in period0_radii]

        p0_omega_spread = relative_spread(
            [x["omega_abs_rad_per_M"] for x in p0]
        )
        p0_alpha_spread = relative_spread(
            [x["alpha_per_M"] for x in p0]
        )

        period0_consistency_pass = bool(
            p0_omega_spread
            <= GATES["period0_frequency_relative_spread_max"]
            and p0_alpha_spread
            <= GATES["period0_damping_relative_spread_max"]
        )

        p1 = [get_fit(r, 1) for r in RADII_TWO_PERIOD]

        p1_omega_spread = relative_spread(
            [x["omega_abs_rad_per_M"] for x in p1]
        )
        p1_alpha_spread = relative_spread(
            [x["alpha_per_M"] for x in p1]
        )

        period1_consistency_pass = bool(
            p1_omega_spread
            <= GATES["period1_frequency_relative_spread_max"]
            and p1_alpha_spread
            <= GATES["period1_damping_relative_spread_max"]
            and all(
                x["omega_relative_error_vs_kerr"]
                <= GATES["period1_frequency_error_vs_kerr_max"]
                for x in p1
            )
            and all(
                x["alpha_relative_error_vs_kerr"]
                <= GATES["period1_damping_error_vs_kerr_max"]
                for x in p1
            )
            and all(
                x["phase_r2"] >= GATES["period1_phase_r2_min"]
                for x in p1
            )
            and all(
                x["log_amplitude_r2"]
                >= GATES["period1_logamp_r2_min"]
                for x in p1
            )
        )

        # ----------------------------------------------------------
        # 7. Сходимость от первого периода ко второму на r=15,30,40.
        # ----------------------------------------------------------
        convergence_rows = []

        for r in RADII_TWO_PERIOD:
            a = get_fit(r, 0)
            b = get_fit(r, 1)

            row = {
                "radius_M": r,
                "period0_frequency_error": a["omega_relative_error_vs_kerr"],
                "period1_frequency_error": b["omega_relative_error_vs_kerr"],
                "frequency_error_decreased": bool(
                    b["omega_relative_error_vs_kerr"]
                    < a["omega_relative_error_vs_kerr"]
                ),
                "period0_damping_error": a["alpha_relative_error_vs_kerr"],
                "period1_damping_error": b["alpha_relative_error_vs_kerr"],
                "damping_error_decreased": bool(
                    b["alpha_relative_error_vs_kerr"]
                    < a["alpha_relative_error_vs_kerr"]
                ),
            }
            convergence_rows.append(row)

        convergence_pass = all(
            x["frequency_error_decreased"]
            and x["damping_error_decreased"]
            for x in convergence_rows
        )

        # ----------------------------------------------------------
        # 8. Третий период на r=15,30 и жёсткий Kerr против свободной моды.
        # ----------------------------------------------------------
        third_period_rows = []

        for r in RADII_THREE_PERIOD:
            fit = get_fit(r, 2)
            t = modes[r]["t"]
            z = modes[r]["zp"]

            mask = (
                (t >= fit["start_M"] - 1e-12)
                & (t <= fit["end_M"] + 1e-12)
            )
            tt = t[mask]
            zz = z[mask]
            tau = tt - fit["start_M"]

            _, rss_kerr = best_complex_amplitude(
                tau, zz, alpha_kerr, -omega_kerr
            )
            _, rss_free = best_complex_amplitude(
                tau, zz,
                fit["alpha_per_M"],
                fit["omega_signed_rad_per_M"],
            )

            bic_kerr = bic(rss_kerr, 2, 2 * len(tt))
            bic_free = bic(rss_free, 4, 2 * len(tt))
            delta_bic = float(bic_kerr - bic_free)

            row = dict(fit)
            row.update({
                "rss_fixed_kerr": rss_kerr,
                "rss_free_single_mode": rss_free,
                "delta_BIC_fixed_kerr_minus_free": delta_bic,
                "fixed_kerr_not_strongly_disfavored":
                    bool(delta_bic <= GATES["period2_fixed_kerr_delta_bic_max"]),
            })
            third_period_rows.append(row)

        period2_consistency_pass = bool(
            all(
                x["omega_relative_error_vs_kerr"]
                <= GATES["period2_frequency_error_vs_kerr_max"]
                for x in third_period_rows
            )
            and all(
                x["alpha_relative_error_vs_kerr"]
                <= GATES["period2_damping_error_vs_kerr_max"]
                for x in third_period_rows
            )
            and all(
                x["phase_r2"] >= GATES["period2_phase_r2_min"]
                for x in third_period_rows
            )
            and all(
                x["log_amplitude_r2"]
                >= GATES["period2_logamp_r2_min"]
                for x in third_period_rows
            )
            and all(
                x["fixed_kerr_not_strongly_disfavored"]
                for x in third_period_rows
            )
        )

        # ----------------------------------------------------------
        # 9. Комплексная форма волны на втором периоде:
        # возраст кольцевания одинаков, амплитуда умножена на r.
        # Для r=30 и 40 разрешается только постоянный комплексный множитель
        # относительно r=15 (постоянная фаза/нормировка).
        # ----------------------------------------------------------
        age_grid = np.arange(
            period_kerr,
            2.0 * period_kerr + 1e-12,
            0.25,
        )

        tref = modes[15.0]["peak_t"] + age_grid
        zref = interp_complex(
            tref, modes[15.0]["t"], modes[15.0]["zp"]
        )
        yref = 15.0 * zref

        shape_rows = []

        for r in [30.0, 40.0]:
            tr = modes[r]["peak_t"] + age_grid
            zr = interp_complex(
                tr, modes[r]["t"], modes[r]["zp"]
            )
            yr = r * zr

            c = np.vdot(yr, yref) / np.vdot(yr, yr)
            pred = c * yr

            denom = math.sqrt(
                float(np.mean(np.abs(yref) ** 2))
            )
            nrmse = float(
                math.sqrt(float(np.mean(np.abs(pred - yref) ** 2)))
                / max(denom, 1e-300)
            )

            coherence = float(
                abs(np.vdot(yref, yr))
                / math.sqrt(
                    float(np.vdot(yref, yref).real)
                    * float(np.vdot(yr, yr).real)
                )
            )

            shape_rows.append({
                "reference_radius_M": 15.0,
                "radius_M": r,
                "period_index": 1,
                "best_complex_factor_real": float(np.real(c)),
                "best_complex_factor_imag": float(np.imag(c)),
                "complex_coherence": coherence,
                "normalized_rmse": nrmse,
            })

        shape_pass = all(
            x["complex_coherence"]
            >= GATES["complex_shape_coherence_min"]
            and x["normalized_rmse"]
            <= GATES["complex_shape_nrmse_max"]
            for x in shape_rows
        )

        # ----------------------------------------------------------
        # 10. Общая классификация.
        # ----------------------------------------------------------
        main_pass = all([
            propagation_pass,
            inverse_radius_pass,
            symmetry_pass,
            period0_consistency_pass,
            period1_consistency_pass,
            convergence_pass,
            period2_consistency_pass,
            shape_pass,
        ])

        if main_pass:
            classification = (
                "BH_HF3_R6B_FIVE_RADIUS_OUTGOING_PROPAGATION_SUPPORTED_"
                "INVERSE_RADIUS_SCALING_SUPPORTED_"
                "RINGDOWN_AGE_COLLAPSE_SUPPORTED_"
                "KERR220_ASYMPTOTIC_CLOSURE_SUPPORTED_R15_R30_"
                "R40_THIRD_PERIOD_PENDING"
            )
        else:
            classification = (
                "BH_HF3_R6B_MULTIRADIUS_LATE_CLOSURE_REVIEW_REQUIRED"
            )

        summary = {
            "stage": "BH-HF3-R6B",
            "version": VERSION,
            "classification": classification,
            "main_pass": main_pass,

            "remnant": {
                "late_window_M": list(LATE_REMNANT_WINDOW),
                "Mf": Mf,
                "Jf": Jf,
                "chi": chi,
            },

            "kerr_220_reference": {
                "omega_R_rad_per_M": omega_kerr,
                "alpha_per_M": alpha_kerr,
                "period_M": period_kerr,
                "Q": Q,
            },

            "five_radius_peak_propagation": {
                "rows": propagation_rows,
                "pass": propagation_pass,
            },

            "inverse_radius_peak_scaling": {
                "r_times_peak_amplitude_relative_spread":
                    scaled_peak_spread,
                "pass": inverse_radius_pass,
            },

            "period0_consistency": {
                "radii_M": period0_radii,
                "frequency_relative_spread": p0_omega_spread,
                "damping_relative_spread": p0_alpha_spread,
                "pass": period0_consistency_pass,
            },

            "period1_consistency": {
                "radii_M": RADII_TWO_PERIOD,
                "frequency_relative_spread": p1_omega_spread,
                "damping_relative_spread": p1_alpha_spread,
                "pass": period1_consistency_pass,
            },

            "period0_to_period1_convergence": {
                "rows": convergence_rows,
                "pass": convergence_pass,
            },

            "period2_r15_r30": {
                "rows": third_period_rows,
                "pass": period2_consistency_pass,
            },

            "plus_minus_symmetry": {
                "pass": symmetry_pass,
            },

            "second_period_complex_shape": {
                "rows": shape_rows,
                "pass": shape_pass,
            },

            "gates": GATES,

            "scientific_boundary": (
                "Результат является сильной внутренней многорадиусной "
                "проверкой стандартного исходящего кольцевания Kerr/ОТО. "
                "Он не является независимым численным воспроизведением, "
                "поскольку все радиусы получены в одной симуляции, и не "
                "является доказательством новой физики."
            ),

            "next_question": (
                "Нужен ли ещё расчёт примерно до 105M для третьего полного "
                "периода на r=40M, либо текущего многорадиусного замыкания "
                "достаточно для завершения HF3 и перехода к HF4."
            ),
        }

        (outdir / "BH_HF3_R6B_summary.json").write_text(
            json.dumps(summary, indent=2, ensure_ascii=False),
            encoding="utf-8",
        )

        def write_csv(name, rows):
            if not rows:
                return
            fields = sorted(set().union(*(x.keys() for x in rows)))
            with (outdir / name).open(
                "w", newline="", encoding="utf-8"
            ) as f:
                w = csv.DictWriter(f, fieldnames=fields)
                w.writeheader()
                w.writerows(rows)

        write_csv("BH_HF3_R6B_peaks.csv", peak_rows)
        write_csv("BH_HF3_R6B_peak_propagation.csv", propagation_rows)
        write_csv("BH_HF3_R6B_age_matched_fits.csv", age_fit_rows)
        write_csv("BH_HF3_R6B_plus_minus_symmetry.csv", symmetry_rows)
        write_csv("BH_HF3_R6B_convergence.csv", convergence_rows)
        write_csv("BH_HF3_R6B_third_period_kerr.csv", third_period_rows)
        write_csv("BH_HF3_R6B_complex_shape.csv", shape_rows)

        report = [
            "BH-HF3-R6B: поздняя многорадиусная проверка",
            "=" * 76,
            f"Классификация: {classification}",
            "",
            "Поздний остаток:",
            f"  Mf = {Mf:.9f}",
            f"  Jf = {Jf:.9f}",
            f"  chi = {chi:.9f}",
            "",
            "Kerr (2,2,0):",
            f"  omega = {omega_kerr:.9f} rad/M",
            f"  alpha = {alpha_kerr:.9f} 1/M",
            f"  период = {period_kerr:.6f} M",
            "",
            "Максимумы на радиусах:",
        ]

        for x in peak_rows:
            report.append(
                f"  r={x['radius_M']:.0f}M: "
                f"t_peak={x['peak_time_M']:.3f}M, "
                f"|Psi4|={x['peak_abs_Psi4_22']:.9e}, "
                f"r|Psi4|={x['r_times_peak_abs_Psi4_22']:.9e}"
            )

        report += [
            "",
            "Приращения времени максимума против Kerr r*:",
        ]

        for x in propagation_rows:
            report.append(
                f"  {x['r1_M']:.0f}->{x['r2_M']:.0f}M: "
                f"наблюдение={x['observed_peak_delay_increment_M']:.3f}M, "
                f"Kerr={x['kerr_rstar_increment_M']:.3f}M, "
                f"ошибка={100*x['relative_difference']:.3f}%"
            )

        report += [
            "",
            f"Разброс r|Psi4_peak| по пяти радиусам: "
            f"{100*scaled_peak_spread:.3f}%",
            "",
            "Первый период после максимума:",
            f"  разброс omega={100*p0_omega_spread:.3f}%",
            f"  разброс alpha={100*p0_alpha_spread:.3f}%",
            "",
            "Второй период после максимума, r=15,30,40M:",
            f"  разброс omega={100*p1_omega_spread:.3f}%",
            f"  разброс alpha={100*p1_alpha_spread:.3f}%",
        ]

        for x in p1:
            report.append(
                f"  r={x['radius_M']:.0f}M: "
                f"|omega|={x['omega_abs_rad_per_M']:.9f}, "
                f"ошибка omega={100*x['omega_relative_error_vs_kerr']:.3f}%, "
                f"alpha={x['alpha_per_M']:.9f}, "
                f"ошибка alpha={100*x['alpha_relative_error_vs_kerr']:.3f}%"
            )

        report += [
            "",
            "Третий период, r=15,30M:",
        ]

        for x in third_period_rows:
            report.append(
                f"  r={x['radius_M']:.0f}M: "
                f"|omega|={x['omega_abs_rad_per_M']:.9f}, "
                f"ошибка omega={100*x['omega_relative_error_vs_kerr']:.3f}%, "
                f"alpha={x['alpha_per_M']:.9f}, "
                f"ошибка alpha={100*x['alpha_relative_error_vs_kerr']:.3f}%, "
                f"ΔBIC={x['delta_BIC_fixed_kerr_minus_free']:.3f}"
            )

        report += [
            "",
            "Комплексная форма второго периода относительно r=15M:",
        ]

        for x in shape_rows:
            report.append(
                f"  r={x['radius_M']:.0f}M: "
                f"когерентность={x['complex_coherence']:.9f}, "
                f"нормированная ошибка={x['normalized_rmse']:.6f}"
            )

        report += [
            "",
            f"Пятирадиусное распространение: {propagation_pass}",
            f"Масштабирование 1/r: {inverse_radius_pass}",
            f"Симметрия +/-m: {symmetry_pass}",
            f"Согласованность первого периода: {period0_consistency_pass}",
            f"Согласованность второго периода: {period1_consistency_pass}",
            f"Сходимость первый -> второй период: {convergence_pass}",
            f"Третий период r=15,30M совместим с Kerr: "
            f"{period2_consistency_pass}",
            f"Комплексная форма второго периода согласована: {shape_pass}",
            "",
            "Граница вывода:",
            f"  {summary['scientific_boundary']}",
            "",
            "Следующий вопрос:",
            f"  {summary['next_question']}",
        ]

        (outdir / "BH_HF3_R6B_report.txt").write_text(
            "\n".join(report) + "\n",
            encoding="utf-8",
        )

        outzip = outdir.parent / (outdir.name + ".zip")
        if outzip.exists():
            outzip.unlink()

        with zipfile.ZipFile(
            outzip, "w",
            compression=zipfile.ZIP_DEFLATED,
            compresslevel=9,
        ) as zout:
            for p in sorted(outdir.iterdir()):
                if p.is_file():
                    zout.write(p, arcname=p.name)

    print("=" * 88)
    print("BH-HF3-R6B ЗАВЕРШЁН")
    print("=" * 88)
    print("Классификация:", classification)
    print()
    print(f"Mf={Mf:.9f}, Jf={Jf:.9f}, chi={chi:.9f}")
    print(
        f"Kerr: omega={omega_kerr:.9f}, "
        f"alpha={alpha_kerr:.9f}, период={period_kerr:.6f}M"
    )
    print()
    for x in peak_rows:
        print(
            f"r={x['radius_M']:.0f}M: "
            f"peak={x['peak_time_M']:.3f}M, "
            f"r|Psi4_peak|={x['r_times_peak_abs_Psi4_22']:.9e}"
        )
    print()
    for x in propagation_rows:
        print(
            f"{x['r1_M']:.0f}->{x['r2_M']:.0f}M: "
            f"Δt_peak={x['observed_peak_delay_increment_M']:.3f}M, "
            f"Kerr r*={x['kerr_rstar_increment_M']:.3f}M, "
            f"ошибка={100*x['relative_difference']:.3f}%"
        )
    print()
    print(
        "Разброс r|Psi4_peak|:",
        f"{100*scaled_peak_spread:.3f}%"
    )
    print(
        "Второй период: разброс omega:",
        f"{100*p1_omega_spread:.3f}%"
    )
    print(
        "Второй период: разброс alpha:",
        f"{100*p1_alpha_spread:.3f}%"
    )
    print()
    for x in third_period_rows:
        print(
            f"Третий период r={x['radius_M']:.0f}M: "
            f"ошибка omega={100*x['omega_relative_error_vs_kerr']:.3f}%, "
            f"ошибка alpha={100*x['alpha_relative_error_vs_kerr']:.3f}%, "
            f"ΔBIC={x['delta_BIC_fixed_kerr_minus_free']:.3f}"
        )
    print()
    print("Пятирадиусное распространение:", propagation_pass)
    print("Масштабирование 1/r:", inverse_radius_pass)
    print("Согласованность второго периода:", period1_consistency_pass)
    print("Сходимость к Kerr:", convergence_pass)
    print("Третий период r=15,30M:", period2_consistency_pass)
    print("Комплексная форма второго периода:", shape_pass)
    print()
    print("ZIP создан:")
    print(outzip)
    print()
    print("Пришлите этот один ZIP в ChatGPT.")
    print("Новый расчёт Einstein Toolkit пока не запускайте.")


if __name__ == "__main__":
    main()
