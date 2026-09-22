#!/usr/bin/env python3
"""HF3/Yitian A2: фиксированные керровские моды и отложенный отрезок времени.

Не запускает SpEC или Einstein Toolkit. Не исполняет присланные блокноты.
По умолчанию использует сохранённую таблицу частот; --recompute-spectrum
пересчитывает её с qnm 0.4.4. См. отчёт об ограничениях интерпретации.
"""
from __future__ import annotations

import argparse
import hashlib
import io
import json
import platform
import shutil
import zipfile
from pathlib import Path

import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import scipy

VERSION = "BH-YITIAN-A2-2026-09-15-v1"
EXPECTED_H5_SHA256 = "9cdd268e85b783d8a832b72e46945c6b43bb2b3ea0a9e6665e57a5a1b7d56f0f"
MF, CHI = 0.95162, 0.68644
OMEGA = 0.20878380853617037
STARTS = [0, 5, 10, 15, 20, 30, 40, 50, 60, 70]
SPECTRAL_MODES = [(2, 0), (2, 1), (2, 2), (2, 3), (3, 0), (4, 0)]
MODELS = {
    "F220": ["220"],
    "F320": ["320"],
    "F2": ["220", "320"],
    "F3": ["220", "320", "420"],
    "F2_O1": ["220", "320", "221"],
    "F2_O2": ["220", "320", "221", "222"],
    "F2_O3": ["220", "320", "221", "222", "223"],
    "F3_O1": ["220", "320", "420", "221"],
    "F3_O2": ["220", "320", "420", "221", "222"],
    "F3_O3": ["220", "320", "420", "221", "222", "223"],
}
PAPER = "https://arxiv.org/pdf/2208.02965"
QNM_DOC = "https://qnm.readthedocs.io/en/latest/README.html"


def dump(path, obj):
    Path(path).write_text(json.dumps(obj, ensure_ascii=False, indent=2,
                                    allow_nan=False) + "\n", encoding="utf-8")


def cj(z):
    return {"real": float(np.real(z)), "imag": float(np.imag(z))}


def compute_spectrum():
    import qnm
    # Настройки qnm по умолчанию; более строгий допуск не означает
    # автоматически более точный корень и здесь не требуется.
    sequences = {}
    for ell, n in SPECTRAL_MODES:
        seq = qnm.spinsequence.KerrSpinSeq(s=-2, l=ell, m=2, n=n, a_max=0.70)
        seq.do_find_sequence()
        sequences[f"{ell}2{n}"] = seq

    def spectrum(mass, chi):
        result = {}
        for key, seq in sequences.items():
            w, _, _ = seq(a=chi)
            cf_err, terms = seq.solver.get_cf_err()
            result[key] = {"omega": float(w.real / mass),
                           "alpha": float(-w.imag / mass),
                           "cf_error_indicator": float(cf_err),
                           "cf_terms": int(terms)}
        return result

    central = spectrum(MF, CHI)
    variations = {
        "mass_minus_half_last_digit": spectrum(MF - 5e-6, CHI),
        "mass_plus_half_last_digit": spectrum(MF + 5e-6, CHI),
        "chi_minus_half_last_digit": spectrum(MF, CHI - 5e-6),
        "chi_plus_half_last_digit": spectrum(MF, CHI + 5e-6),
    }
    # Сверка опубликованных округлённых частот и времени затухания.
    table = {"220": (.5535, 11.707, 5e-5, 5e-4),
             "221": (.5410, 3.8713, 5e-5, 5e-5),
             "222": (.5180, 2.2923, 5e-5, 5e-5),
             "320": (.7920, 11.235, 5e-5, 5e-4),
             "420": (1.0172, 10.938, 5e-5, 5e-4)}
    table_checks = {}
    for k, (w, tau, dw, dtau) in table.items():
        table_checks[k] = bool(abs(central[k]["omega"] - w) <= dw and
                               abs(1 / central[k]["alpha"] - tau) <= dtau)
    if not all(table_checks.values()):
        raise ValueError("Численный спектр расходится с округлённой таблицей статьи")
    wq, _, _ = sequences["220"](a=0.6783194994253878)
    massq = 0.977461203278627
    return {"method": "qnm, решение радиального и углового уравнений Керра",
            "qnm_version": qnm.__version__, "s": -2, "a_max": .70,
            "l_max": 20, "root_tol": float(np.sqrt(np.finfo(float).eps)),
            "cf_tol": 1e-10, "mass_M": MF, "chi": CHI,
            "frequency_units": "M^-1; M — исходная полная масса ADM",
            "central": central, "rounding_variations": variations,
            "rounded_paper_table_checks": table_checks,
            "QC0_220": {"mass_M": massq, "chi": .6783194994253878,
                        "omega": float(wq.real / massq),
                        "alpha": float(-wq.imag / massq)},
            "sources": [PAPER, QNM_DOC],
            "precision_note": "Округление Mf и chi ограничивает физическую точность; cf_error_indicator не является полной погрешностью частоты."}


def read_moments(path):
    raw = path.read_bytes()
    input_sha = hashlib.sha256(raw).hexdigest()
    if zipfile.is_zipfile(io.BytesIO(raw)):
        with zipfile.ZipFile(io.BytesIO(raw)) as z:
            candidates = [n for n in z.namelist() if n.endswith('.h5') and
                          not n.startswith('__MACOSX/')]
            if len(candidates) != 1:
                raise ValueError("В архиве должен быть один файл .h5")
            raw = z.read(candidates[0])
    h5sha = hashlib.sha256(raw).hexdigest()
    if h5sha != EXPECTED_H5_SHA256:
        raise ValueError("Этот этап продолжает проверенный A1: контрольная сумма входа другая")
    with h5py.File(io.BytesIO(raw), 'r') as f:
        t = f['time_adaptive'][()]
        values = {"I22": f['geo_mass/Adaptive'][()][:, 8].conj(),
                  "L32": f['geo_spin/Adaptive'][()][:, 14].conj()}
    if t.shape != (5001,) or not np.all(np.diff(t) > 0):
        raise ValueError("Неожиданная временная сетка")
    return t, values, {"input_name": path.name, "input_sha256": input_sha,
                       "hdf5_sha256": h5sha}


def freq_array(table, keys):
    return np.array([table[k]['omega'] - 1j * table[k]['alpha'] for k in keys])


def matrix(t, t0, frequencies):
    return np.exp(-1j * np.outer(t - t0, frequencies))


def complex_lstsq(a, y):
    norms = np.linalg.norm(a, axis=0)
    scaled = a / norms
    cs, _, rank, singular = np.linalg.lstsq(scaled, y, rcond=None)
    if rank != a.shape[1]:
        raise ValueError("Линейная модель потеряла ранг")
    return cs / norms, float(singular[0] / singular[-1])


def err(y, pred):
    return float(np.linalg.norm(y - pred) / np.linalg.norm(y))


def fit_case(t, z, table, model, t0, stride=1):
    keys = MODELS[model]
    w = freq_array(table, keys)
    train = np.flatnonzero((t >= t0 - 1e-9) & (t <= t0 + 30 + 1e-9))[::stride]
    valid = np.flatnonzero((t > t0 + 30 + 1e-9) & (t <= t0 + 50 + 1e-9))
    a = matrix(t[train], t0, w)
    c, cond = complex_lstsq(a, z[train])
    p, v = a @ c, matrix(t[valid], t0, w) @ c
    at20 = c * np.exp(-1j * w * (20 - t0))
    r = {"model": model, "modes": keys, "start_M": t0,
         "train_end_M": t0 + 30, "validation_start_M": float(t[valid[0]]),
         "validation_end_M": t0 + 50, "training_samples": len(train),
         "validation_samples": len(valid), "stride": stride,
         "fit_relative_L2": err(z[train], p),
         "validation_relative_L2": err(z[valid], v),
         "validation_absolute_rms": float(np.sqrt(np.mean(abs(z[valid] - v) ** 2))),
         "scaled_condition_number": cond,
         "coefficients_at_start": {k: cj(value) for k, value in zip(keys, c)},
         "coefficients_at_reference_20M": {k: cj(value) for k, value in zip(keys, at20)}}
    return r, (train, valid, p, v)


def analyze(t, raw, spec):
    table = spec['central']
    offsets = {k: value[t >= 400].mean() for k, value in raw.items()}
    series = {k: (value - offsets[k]) * np.exp(-2j * OMEGA * t)
              for k, value in raw.items()}
    rows = []
    for name, z in series.items():
        for start in STARTS:
            for model in MODELS:
                row, _ = fit_case(t, z, table, model, start)
                row['moment'] = name
                rows.append(row)

    controls = []
    for name in series:
        for model in ['F220', 'F2', 'F3', 'F2_O3', 'F3_O3']:
            for start in [0, 20, 50]:
                for tail in [350, 400, 450]:
                    z = (raw[name] - raw[name][t >= tail].mean()) * np.exp(-2j * OMEGA * t)
                    row, _ = fit_case(t, z, table, model, start)
                    controls.append({"moment": name, "model": model, "start_M": start,
                                     "control": "baseline", "value": tail,
                                     "fit_relative_L2": row['fit_relative_L2'],
                                     "validation_relative_L2": row['validation_relative_L2']})
                for stride in [1, 2, 4]:
                    row, _ = fit_case(t, series[name], table, model, start, stride)
                    controls.append({"moment": name, "model": model, "start_M": start,
                                     "control": "training_stride", "value": stride,
                                     "fit_relative_L2": row['fit_relative_L2'],
                                     "validation_relative_L2": row['validation_relative_L2']})
                for label, changed in spec['rounding_variations'].items():
                    row, _ = fit_case(t, series[name], changed, model, start)
                    controls.append({"moment": name, "model": model, "start_M": start,
                                     "control": "reported_parameter_rounding", "value": label,
                                     "fit_relative_L2": row['fit_relative_L2'],
                                     "validation_relative_L2": row['validation_relative_L2']})

    # Это алгебраическая проверка обработки, а не независимая проверка ОТО.
    shifted = {k: dict(v, omega=v['omega'] - 2 * OMEGA) for k, v in table.items()}
    frame_errors = []
    for name in series:
        for model in MODELS:
            a, _ = fit_case(t, series[name], table, model, 20)
            b, _ = fit_case(t, raw[name] - offsets[name], shifted, model, 20)
            frame_errors += [abs(a[k] - b[k]) for k in ['fit_relative_L2', 'validation_relative_L2']]
    # Известная сумма экспонент должна восстанавливаться и вне окна подгонки.
    w = freq_array(table, MODELS['F3'])
    synthetic = matrix(t, 20, w) @ np.array([.3 + .4j, -.07 + .1j, .015 - .003j])
    sr, _ = fit_case(t, synthetic, table, 'F3', 20)
    checks = {"frame_equivalence_max_error": float(max(frame_errors)),
              "synthetic_fit_relative_L2": sr['fit_relative_L2'],
              "synthetic_prediction_relative_L2": sr['validation_relative_L2']}
    if max(checks.values()) > 1e-10:
        raise ValueError(f"Не прошла проверка вычислительного метода: {checks}")

    floor = {}
    for name in raw:
        noise = raw[name][t >= 400] - offsets[name]
        floor[name] = {"late_offset": cj(offsets[name]),
                       "late_offset_abs": float(abs(offsets[name])),
                       "late_fluctuation_rms": float(np.sqrt(np.mean(abs(noise) ** 2))),
                       "note": "Позднее рассеяние не задаёт погрешность раннего сигнала; второй численной точности нет."}

    def pick(name, model, start):
        return next(r for r in rows if r['moment'] == name and r['model'] == model and r['start_M'] == start)

    examples = [pick('I22', model, 0) for model in ['F220', 'F2', 'F2_O3']]
    examples += [pick('L32', model, 20) for model in ['F320', 'F220', 'F2', 'F3']]
    examples += [pick('I22', model, 20) for model in ['F220', 'F2', 'F2_O3']]
    # Отдельно сохраняем все случаи ухудшения проверки при добавлении моды.
    pairs = [('F220', 'F2'), ('F2', 'F3'), ('F2', 'F2_O1'),
             ('F2_O1', 'F2_O2'), ('F2_O2', 'F2_O3'),
             ('F3', 'F3_O1'), ('F3_O1', 'F3_O2'), ('F3_O2', 'F3_O3')]
    not_improved = []
    for name in series:
        for start in STARTS:
            for base, more in pairs:
                a, b = pick(name, base, start), pick(name, more, start)
                if b['validation_relative_L2'] > a['validation_relative_L2']:
                    not_improved.append({"moment": name, "start_M": start,
                                         "base": base, "expanded": more,
                                         "base_validation_relative_L2": a['validation_relative_L2'],
                                         "expanded_validation_relative_L2": b['validation_relative_L2']})
    summary = {"version": VERSION,
               "status": "BH_YITIAN_A2_KNOWN_MODE_MIXING_REPRODUCED_PREDICTION_NOT_EXACT",
               "scope": "I22 и L32 одного численного расчёта; не вся геометрия и не наблюдения GW150914",
               "rotation_Omega": OMEGA, "mass_M": MF, "chi": CHI,
               "training_duration_M": 30, "validation_duration_M": 20,
               "start_times_M": STARTS, "model_sets": MODELS,
               "fixed_frequency_fit": True,
               "validation_is_independent_simulation": False,
               "validation_note": "Коэффициенты не подгоняются на отложенном участке, но Mf, chi, Omega и поздний уровень заданы по всей симуляции/письму; это не слепой перспективный прогноз.",
               "metric": "E=||z-z_model||_2/||z||_2; равные веса комплексных отсчётов; не поточечная ошибка и не статистическая значимость",
               "checks": checks, "floor_diagnostics": floor,
               "beat_period_220_320_M": float(2 * np.pi / (table['320']['omega'] - table['220']['omega'])),
               "examples": examples,
               "expanded_model_worse_on_validation_count": len(not_improved),
               "nested_model_comparison_count": len(series) * len(STARTS) * len(pairs),
               "maximum_scaled_condition_number": max(r['scaled_condition_number'] for r in rows),
               "limitations": [
                   "Частоты фиксированы параметрами конечной ЧД; амплитуды каждой модели подгоняются отдельно.",
                   "Это повторение и расширение проверок известного механизма из статьи авторов, не открытие смешения мод.",
                   "Успех нескольких мод не означает физическое измерение каждой обертонной амплитуды.",
                   "Остаток включает возможные нелинейные, координатные и численные эффекты; новой физики он сам по себе не доказывает.",
                   "Один уровень разрешения; AMR 5e-8 не используется как ошибка каждого момента.",
                   "Смена шага анализа не является проверкой сходимости пространственной сетки.",
                   "Отложенные отрезки перекрываются у разных начальных времён и не являются независимыми опытами.",
                   "Данные одной симуляции ОТО не дают независимого свидетельства физики за пределами ОТО.",
                   "Другие моменты, особенно m=0 и старшие гармоники, в A2 ещё не исследованы."]}
    return summary, rows, controls, not_improved, series


def qc0_state(prior, spec):
    d = json.loads(prior.read_text(encoding='utf-8'))
    remnant = d['remnant']
    ref = spec['QC0_220']
    if abs(remnant['chi'] - ref['chi']) > 1e-12 or abs(remnant['Mf'] - ref['mass_M']) > 1e-12:
        raise ValueError("Параметры прежнего R6B отличаются от сохранённой спектральной сверки")
    rows = []
    for r in d['period2_r15_r30']['rows']:
        rows.append({"radius_M": r['radius_M'], "unchanged_fit_start_M": r['start_M'],
                     "unchanged_fit_end_M": r['end_M'],
                     "omega_from_existing_fit": r['omega_abs_rad_per_M'],
                     "alpha_from_existing_fit": r['alpha_per_M'],
                     "relative_omega_difference_numerical_Kerr": abs(r['omega_abs_rad_per_M'] - ref['omega']) / ref['omega'],
                     "relative_alpha_difference_numerical_Kerr": abs(r['alpha_per_M'] - ref['alpha']) / ref['alpha']})
    period = 2 * np.pi / ref['omega']
    return {"source_summary_sha256": hashlib.sha256(prior.read_bytes()).hexdigest(),
            "stored_classification_preserved": d['classification'],
            "remnant_from_R6B": remnant, "old_approximate_reference": d['kerr_220_reference'],
            "numerical_Kerr220_reference": ref, "existing_fits_compared_without_refitting": rows,
            "reference_omega_relative_shift": (ref['omega'] - d['kerr_220_reference']['omega_R_rad_per_M']) / d['kerr_220_reference']['omega_R_rad_per_M'],
            "reference_alpha_relative_shift": (ref['alpha'] - d['kerr_220_reference']['alpha_per_M']) / d['kerr_220_reference']['alpha_per_M'],
            "period_M": float(period), "r40_three_period_end_M": float(69 + 3 * period),
            "checkpoint_state_from_prior_verified_dialog": {"t_M": 95, "iteration": 48640, "exit_status": 0},
            "planned_final_time_M": 105, "expected_iteration_if_step_unchanged": 53760,
            "no_new_evolution_run": True,
            "next_step": "BH-YITIAN-A3: осесимметричные I20, I40, L30; затем отдельный контроль m=4. Работать с уже предоставленными рядами.",
            "QC0_completion_pending": "Восстановление 95→105M из it=48640, аудит шва 95M и третий период r40 с численным спектральным ориентиром.",
            "HF4_gate": "Установить соответствие коэффициентов угловому базису/мере; различать скалярную карту на вспомогательной сфере и физическую метрику горизонта.",
            "corrections": [
                "В старом текстовом пересказе Jf=0.648087073319352 было ошибочно названо безразмерным спином; правильный chi=Jf/Mf^2.",
                "Правильное Mf в файле 0.977461203278627; число 0.977641203 в пересказе было опечаткой.",
                "Ориентир R6B — приближённая формула BCW, а не точное численное значение квазинормальной частоты.",
                "В старом сравнении BIC свободные частота и затухание брались из отдельных регрессий фазы и лог-модуля, а не из минимума того же комплексного RSS. Для строгого сравнения моделей это надо исправить.",
                "Коррелированные остатки и отсутствие модели численной ошибки не позволяют трактовать старые BIC как калиброванную статистическую проверку ОТО.",
                "R6A/R6B не пересчитывались: прочитаны готовый итог и исходный код; архив исходных временных рядов здесь не анализировался."]}


def make_figures(out, t, series, spec, rows):
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10,
                         'axes.spines.top': False, 'axes.spines.right': False,
                         'axes.grid': True, 'grid.alpha': .18})
    colors = {'F220': '#27649a', 'F2': '#db7932', 'F3': '#996aaf',
              'F2_O3': '#2d8a71', 'F3_O3': '#2d8a71'}
    fig, axs = plt.subplots(1, 2, figsize=(12.4, 4.8), layout='constrained')
    for ax, name in zip(axs, ['I22', 'L32']):
        models = ['F220', 'F2', 'F2_O3'] if name == 'I22' else ['F220', 'F2', 'F3_O3']
        labels = {'F220': 'Одна основная мода', 'F2': 'Две основные моды',
                  'F2_O3': 'Две основные + три обертона',
                  'F3_O3': 'Три основные + три обертона'}
        for model in models:
            rr = [r for r in rows if r['moment'] == name and r['model'] == model]
            ax.semilogy([r['start_M'] for r in rr],
                        [100 * r['validation_relative_L2'] for r in rr],
                        'o-', ms=4, color=colors[model], label=labels[model])
        ax.set(xlabel=r'Начало подгонки $t_0/M$', ylabel='Ошибка на отложенном отрезке, %',
               title=rf'${name[0]}_{{{name[1:]}}}$', ylim=(.3, 150))
        ax.legend(fontsize=8, loc='upper right')
    fig.suptitle('Проверка предсказания: дополнительные моды помогают не всегда',
                 fontsize=14, fontweight='bold')
    fig.supxlabel('Подгонка на 30M; проверка на следующих 20M. Частоты фиксированы.\n'
                  'Коэффициенты на проверочном отрезке не подбираются; поздний уровень задан заранее.',
                  fontsize=9)
    fig.savefig(out / 'BH_YITIAN_A2_prediction.png', dpi=155, bbox_inches='tight')
    plt.close(fig)

    fig, axs = plt.subplots(2, 2, figsize=(12.4, 7.6), layout='constrained')
    for col, (name, start, selected) in enumerate([('I22', 0, ['F220', 'F2_O3']),
                                                  ('L32', 20, ['F220', 'F2'])]):
        z = series[name]
        use = (t >= start) & (t <= start + 50)
        for model in selected:
            rr, _ = fit_case(t, z, spec['central'], model, start)
            c = np.array([
                complex(rr['coefficients_at_start'][k]['real'], rr['coefficients_at_start'][k]['imag'])
                for k in MODELS[model]])
            pred = matrix(t[use], start, freq_array(spec['central'], MODELS[model])) @ c
            label = 'Одна мода' if model == 'F220' else ('С обертонами' if model == 'F2_O3' else 'Две моды')
            axs[0, col].plot(t[use], pred.real, color=colors[model], lw=1.4, label=label)
            axs[1, col].semilogy(t[use], abs(z[use] - pred), color=colors[model], label=label)
        axs[0, col].plot(t[use][::5], z[use].real[::5], '.', ms=2.4, color='#28303a', label='Данные')
        for ax in axs[:, col]:
            ax.axvspan(start + 30, start + 50, color='#d8e3ed', alpha=.55)
            ax.set_xlabel(r'Время $t/M$')
        axs[0, col].set(title=rf'${name[0]}_{{{name[1:]}}}$: подгонка {start}–{start+30}M',
                        ylabel='Вещественная часть момента')
        axs[1, col].set(ylabel='Модуль комплексного остатка')
        axs[0, col].legend(fontsize=8)
    fig.suptitle('Подгонка и продолжение за её пределы', fontsize=15, fontweight='bold')
    fig.supxlabel('Затенённая область — отложенные данные. Абсолютный остаток не является оценкой новой физики.', fontsize=9)
    fig.savefig(out / 'BH_YITIAN_A2_examples.png', dpi=155, bbox_inches='tight')
    plt.close(fig)


def write_report(out, summary, rows, controls, state, spec):
    def p(x):
        return f'{100*x:.4f}'
    lines = [
        '# Продолжение HF3 с данными Yitian — этап A2', '',
        'Дата: 15 сентября 2026 года. Расчёты эволюции SpEC и Einstein Toolkit не запускались.', '',
        '**Получен новый воспроизводимый результат:** смешение обычных керровских мод существенно улучшает описание L32 и его продолжение на соседний отрезок времени. Для раннего I22 обертоны сильно улучшают подгонку, однако погрешность продолжения остаётся заметно больше. Полное объяснение всех ранних моментов и дополнительная гиперформа этим не установлены.', '',
        '## Подтверждённая точка остановки', '',
        'QC0 завершён до t=95M, it=48640, штатное завершение и две контрольные части зафиксированы в прежнем диалоге. Сейчас прочитаны сохранённые BH_HF3_R6B_output.zip и BH_HF3_R6B.py. Файлы контрольной точки на компьютере пользователя заново не проверялись. Первичный A1 не пересчитывался: исходный HDF5 опознан по контрольной сумме.', '',
        '| Величина | QC0: из файла R6B | Yitian: из письма |',
        '| --- | ---: | ---: |',
        f"| Mf/M | {state['remnant_from_R6B']['Mf']:.12f} | {MF:.5f} |",
        f"| Jf/M² | {state['remnant_from_R6B']['Jf']:.12f} | Не сопоставлялось напрямую |",
        f"| χf=Jf/Mf² | {state['remnant_from_R6B']['chi']:.12f} | {CHI:.5f} |", '',
        'В прежнем ответе ассистент ошибочно назвал Jf безразмерным спином и переставил цифры в Mf. Программа R6B использовала правильное χf. Начальные условия и нормировки двух симуляций нельзя считать тождественными; параметры между ними не подменяются.', '',
        '## Что сделано в A2', '',
        'Выбраны два момента одного набора: массовый I22 и спиновый L32. Модельные наборы взяты из обсуждения смешения мод авторами и расширены контрольными вложенными наборами. Их физическая мотивировка известна заранее, но A2 является исследовательской проверкой, а не заранее зарегистрированным статистическим испытанием.', '',
        'После сопряжения столбца и вычитания среднего за 400–500M применён авторский множитель exp(−2iΩt), Ω=0.20878380853617037/M. Затем подгонялись только комплексные амплитуды:', '',
        '$$z_{\\mathrm{model}}(t)=\\sum_k C_k\\exp[-(\\alpha_k+i\\omega_k)(t-t_0)].$$', '',
        'Частоты и затухания вычислены qnm 0.4.4 при Mf=0.95162M, χ=0.68644. Амплитуды определяются комплексным методом наименьших квадратов с равными весами отсчётов; перед решением нормируются столбцы матрицы. Частоты, масса и спин под момент не подбирались. Керровские частоты спинового веса −2 применены согласно модели статьи; это не утверждение, что сами скалярные геометрические моменты имеют спиновый вес −2.', '',
        '| Мода (ℓ,m,n) | ω, M⁻¹ | α, M⁻¹ |', '| --- | ---: | ---: |']
    for k, v in spec['central'].items():
        lines.append(f"| ({','.join(k)}) | {v['omega']:.9f} | {v['alpha']:.9f} |")
    lines += ['', 'Обертон — более быстро затухающее колебание того же семейства. Индекс сферического момента, например L32, не заставляет его состоять только из одной моды (3,2,0). Дополнительные моды могут возникать из проекции угловой структуры возмущения на выбранные гармоники.', '',
              f"Период биений основных мод (2,2,0) и (3,2,0): {summary['beat_period_220_320_M']:.4f}M.", '',
              'Подгонка: [t0,t0+30]M. Проверка: (t0+30,t0+50]M, без повторной подгонки амплитуд. t0=0,5,10,15,20,30,40,50,60,70M. Различные окна перекрываются; их нельзя считать независимыми экспериментами. Поздний постоянный уровень и параметры конечной ЧД известны по другим частям этой же симуляции, поэтому это не полностью слепое предсказание будущего.', '',
              '$$E=\\frac{\\left[\\sum_j|z_j-z_{\\mathrm{model},j}|^2\\right]^{1/2}}{\\left[\\sum_j|z_j|^2\\right]^{1/2}}.$$', '',
              'Ниже показано 100E. Это относительная комплексная норма ошибки на отрезке; она не является максимальной поточечной ошибкой или вероятностью ошибки.', '',
              '| Момент | Подгонка, M | Проверка, M | Набор мод | Ошибка подгонки, % | Ошибка проверки, % |',
              '| --- | --- | --- | --- | ---: | ---: |']
    for r in summary['examples']:
        lines.append(f"| {r['moment']} | {r['start_M']}–{r['train_end_M']} | ({r['train_end_M']}, {r['validation_end_M']}] | {' + '.join(r['modes'])} | {p(r['fit_relative_L2'])} | {p(r['validation_relative_L2'])} |")
    lines += ['', '**Содержательный вывод по L32:** две фундаментальные моды существенно улучшают как подгонку, так и проверку на следующем отрезке по сравнению с одной. Это воспроизведение известного механизма смешения, а не новое открытие этого механизма.', '',
              '**Содержательный вывод по I22:** на самом раннем окне добавление трёх обертонов резко снижает ошибку подгонки, но продолжение на 30–50M менее точно. На окне 20–50M добавление второй фундаментальной моды даже немного ухудшает проверку 50–70M. Поэтому формулировка «больше мод — лучшее физическое объяснение» здесь недостаточна.', '',
              f"Во всей сетке {summary['expanded_model_worse_on_validation_count']} из {summary['nested_model_comparison_count']} вложенных сравнений дали ухудшение проверки при расширении модели. Это учёт отрицательных результатов, не частота статистической ошибки. Полные результаты и амплитуды сохранены в JSON.", '',
              '![Проверка на отложенном отрезке](BH_YITIAN_A2_prediction.png)', '',
              '![Примеры подгонки и проверки](BH_YITIAN_A2_examples.png)', '',
              '## Устойчивость и вычислительные проверки', '',
              'Проверены: среднее по 350–500, 400–500 и 450–500M; прореживание обучающих отсчётов в 2 и 4 раза; отдельные изменения Mf и χ на половину последнего указанного разряда (±0.000005). Последнее — чувствительность к округлению, не полная погрешность массы или спина. Проверочные отсчёты при прореживании остаются прежними.', '',
              '| Пример | Изменение позднего среднего: ошибка проверки, % | Шаг обучения 0.1–0.4M: ошибка проверки, % | Округление Mf, χ: ошибка проверки, % |',
              '| --- | ---: | ---: | ---: |']
    for name, model, start in [('I22', 'F2_O3', 0), ('L32', 'F2', 20), ('I22', 'F2', 20)]:
        ranges = []
        for control in ['baseline', 'training_stride', 'reported_parameter_rounding']:
            vals = [r['validation_relative_L2'] for r in controls if r['moment'] == name and r['model'] == model and r['start_M'] == start and r['control'] == control]
            ranges.append(f'{100*min(vals):.5f}–{100*max(vals):.5f}')
        lines.append(f"| {name}, {model}, t0={start}M | {' | '.join(ranges)} |")
    lines += ['', f"Максимальное различие относительных ошибок при согласованном пересчёте в неповёрнутой системе: {summary['checks']['frame_equivalence_max_error']:.3g}. Это проверка алгебры, а не независимое физическое подтверждение.", '',
              f"Синтетическая известная сумма трёх мод восстановлена с ошибкой продолжения {summary['checks']['synthetic_prediction_relative_L2']:.3g}. Максимальное число обусловленности нормированной матрицы во всех основных подгонках: {summary['maximum_scaled_condition_number']:.3g}.", '',
              'Поздние флуктуации и вычитаемые уровни сохранены отдельно. Они не заменяют сравнение двух разрешений. Значение AMR из письма 5×10⁻⁸ не использовано как погрешность каждого мультиполя. Сам по себе остаток после конкретного набора мод не выделяет новую структуру: возможны дополнительные обычные моды, изменение их эффективных амплитуд, нелинейность, выбор времени/координат и численная ошибка.', '',
              '## Уточнение керровского ориентира QC0', '',
              'В исходном R6B применялась приближённая формула BCW. По его собственным Mf и χ численный расчёт даёт:', '',
              '| Ориентир | ω, M⁻¹ | α, M⁻¹ |', '| --- | ---: | ---: |',
              f"| R6B, приближённая формула | {state['old_approximate_reference']['omega_R_rad_per_M']:.9f} | {state['old_approximate_reference']['alpha_per_M']:.9f} |",
              f"| Численное решение Керра | {state['numerical_Kerr220_reference']['omega']:.9f} | {state['numerical_Kerr220_reference']['alpha']:.9f} |", '',
              'Для ранее сохранённых оценок третьего периода, без повторения подгонок или смены старых окон:', '',
              '| Радиус | Различие частоты с численным Керром, % | Различие затухания, % |', '| --- | ---: | ---: |']
    for r in state['existing_fits_compared_without_refitting']:
        lines.append(f"| {r['radius_M']:g}M | {p(r['relative_omega_difference_numerical_Kerr'])} | {p(r['relative_alpha_difference_numerical_Kerr'])} |")
    lines += ['', 'Уточнение ориентира улучшает это описательное согласие. Старые ΔBIC автоматически не переносятся: их надо заново оценивать на исходных рядах с одной и той же целевой функцией. В старой программе свободные частота и затухание извлекались из отдельных регрессий фазы и лог-модуля, а не оптимизировались по тому же комплексному остатку, который входил в BIC. Коррелированные численные остатки также не задают готовую статистическую модель.', '',
              f"При численном периоде T={state['period_M']:.6f}M конец третьего периода на r=40M ожидается около {state['r40_three_period_end_M']:.6f}M, если использовать прежний максимум 69M. Граница 105M остаётся достаточной. При прежнем шаге 1/512M ожидаемая конечная итерация 53760; это ожидание, а не выполненный расчёт.", '',
              '## Обновлённый маршрут', '',
              '| Шаг | Что делаем | Условие завершения и перехода |', '| --- | --- | --- |',
              '| A1 — завершён | Целостность данных Yitian, определения, поздние моменты | Сохранённый отчёт A1; заново не пересчитываем |',
              '| A2 — завершён | I22, L32: смешение, обертоны, отложенные отрезки | Настоящий отчёт; успехи и ограничения сохранены |',
              '| A3 — следующий шаг | Осесимметричные I20, I40, L30 и затем m=4: проверить поздний фон и ранние остатки | Сопоставить модели, окна, изменение коэффициентов и доступный численный уровень; не выбирать результат только по малому остатку |',
              '| QC0, завершение HF3 | Подготовить восстановление 95→105M из it=48640; сохранить сетку и диагностику; проверить шов 95M и третий период r=40M | Частотный ориентир — собственные Mf,χ и численное решение; старые данные сохраняются; сравнение моделей по одному комплексному остатку |',
              '| HF4-A | Установить, что именно можно восстановить из моментов и определения углового базиса | Отдельно карта скалярного поля и собственная метрика; для последней нужны площадь, базис/мера или qAB |',
              '| HF4-B/C | Проверенная двумерная геометрия, затем обоснованное трёхмерное представление и фильм | Контроль углового разрешения и ошибок; видимые деформации не выдаются за внутреннее устройство ЧД |',
              '| HF5 | Связь горизонта с исходящей волной одного и того же расчёта | Нужны соответствующие волновые ряды и согласование времени/конвенций |', '',
              'Порядок ближайших действий: A3 по имеющимся данным, затем окончательный контроль QC0 перед заявлением о завершении его HF3. Найденный prepare_qc0_t95_recovery.py сохраняется как основа следующего восстановления; повторно просить пользователя искать его не требуется.', '',
              'Ранее обсуждавшиеся 64/128 направлений остаются возможными сетками визуализации, но не независимыми наблюдателями. Новые направления, рассчитанные из тех же коэффициентов, сами по себе не дают независимой проверки реконструкции. Подбор числа направлений должен учитывать ℓ≤8, симметрии, угловую меру и конкретную задачу.', '',
              'Под общим горизонтом этого набора понимается поверхность, найденная в численной симуляции как общий кажущийся горизонт, а не восстановленный глобальный горизонт событий. Мультиполи на нём не задают однозначно внутреннюю геометрию чёрной дыры. Исследование «гиперформы» остаётся гипотезой о согласованной динамической структуре; вывод о дополнительных измерениях отсутствует.', '',
              '## Повторение', '',
              'Распаковать пакет. Исходный архив Yitian сохраняется отдельно; его не требуется скачивать или пересчитывать снова. Пример команды из папки пакета:', '',
              '```bash', 'python BH_YITIAN_A2.py --input /путь/common_horizon_moments_data_and_scripts.zip --output repeated', '```', '',
              'Нужны numpy, scipy, h5py и matplotlib. qnm требуется только при флаге --recompute-spectrum. Пакет содержит вычисленную таблицу спектра и её происхождение. Изменять среду Windows пользователя для чтения результатов не нужно.', '',
              'Источники методики: [Chen и соавторы, разделы IV A–D](https://arxiv.org/pdf/2208.02965); [документация qnm](https://qnm.readthedocs.io/en/latest/README.html). Числа настоящего отчёта получены из предоставленного файла и отдельно выполненного расчёта спектра, а не переписаны из статьи.', '']
    (out / 'BH_YITIAN_A2_report.md').write_text('\n'.join(lines), encoding='utf-8')


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    root = Path(__file__).resolve().parent
    ap.add_argument('--input', type=Path, required=True)
    ap.add_argument('--output', type=Path, default=Path('BH_YITIAN_A2_results'))
    ap.add_argument('--spectrum', type=Path, default=root / 'BH_YITIAN_A2_spectrum.json')
    ap.add_argument('--r6b-summary', type=Path, default=root / 'BH_HF3_R6B_summary.json')
    ap.add_argument('--recompute-spectrum', action='store_true')
    args = ap.parse_args()
    out = args.output.resolve()
    out.mkdir(parents=True, exist_ok=True)
    spec = compute_spectrum() if args.recompute_spectrum else json.loads(args.spectrum.read_text())
    if spec['mass_M'] != MF or spec['chi'] != CHI:
        raise ValueError('Не та таблица частот')
    t, raw, provenance = read_moments(args.input)
    summary, rows, controls, negative, series = analyze(t, raw, spec)
    state = qc0_state(args.r6b_summary, spec)
    provenance.update({'version': VERSION, 'python': platform.python_version(),
                       'numpy': np.__version__, 'scipy': scipy.__version__,
                       'h5py': h5py.__version__, 'matplotlib': matplotlib.__version__,
                       'program_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                       'qnm_spectrum_was_recomputed': args.recompute_spectrum})
    for name, obj in [('BH_YITIAN_A2_summary.json', summary),
                      ('BH_YITIAN_A2_all_fits.json', rows),
                      ('BH_YITIAN_A2_sensitivity.json', controls),
                      ('BH_YITIAN_A2_negative_controls.json', negative),
                      ('BH_YITIAN_A2_spectrum.json', spec),
                      ('BH_HF3_resume_state.json', state),
                      ('BH_YITIAN_A2_provenance.json', provenance)]:
        dump(out / name, obj)
    make_figures(out, t, series, spec, rows)
    write_report(out, summary, rows, controls, state, spec)
    if Path(__file__).resolve() != out / Path(__file__).name:
        shutil.copy2(__file__, out / Path(__file__).name)
    if args.r6b_summary.resolve() != out / args.r6b_summary.name:
        shutil.copy2(args.r6b_summary, out / args.r6b_summary.name)
    (out / 'requirements.txt').write_text('numpy\nscipy\nh5py\nmatplotlib\n# Только для --recompute-spectrum:\n# qnm==0.4.4\n', encoding='utf-8')
    print(json.dumps({'status': summary['status'], 'fits': len(rows),
                      'sensitivity_cases': len(controls), 'checks': summary['checks'],
                      'output': str(out)}, ensure_ascii=False, indent=2))


if __name__ == '__main__':
    main()
