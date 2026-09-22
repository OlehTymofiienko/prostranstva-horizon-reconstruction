#!/usr/bin/env python3
"""A3: осесимметричные моменты, постоянный фон и проверка продолжения.

Не запускает эволюцию и не повторяет этапы A1/A2.
Вход: уже проверенный архив Yitian или его HDF5.
По умолчанию спектр берётся из прилагаемого JSON; --recompute-spectrum
пересчитывает только новые частоты m=0, используя qnm 0.4.4.
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
import numpy as np
import scipy
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

VERSION = 'BH-YITIAN-A3-2026-09-15-v1'
H5_SHA = '9cdd268e85b783d8a832b72e46945c6b43bb2b3ea0a9e6665e57a5a1b7d56f0f'
MF, CHI = .95162, .68644
STARTS = [0, 5, 10, 15, 20, 30, 40, 50, 60, 70]
NATIVE = {'I20': 'F1', 'I40': 'F3', 'L30': 'F2'}
MODELS = {
    'F1': ['200'], 'F2': ['200', '300'], 'F3': ['200', '300', '400'],
    'F1O1': ['200', '201'], 'F1O2': ['200', '201', '202'],
    'F1O3': ['200', '201', '202', '203'],
    'F3O1': ['200', '300', '400', '201'],
    'F3O2': ['200', '300', '400', '201', '202'],
    'F3O3': ['200', '300', '400', '201', '202', '203'],
    'D': ['D'], 'F1D': ['200', 'D'], 'F2D': ['200', '300', 'D'],
    'F3D': ['200', '300', '400', 'D'],
}


def dump(p, obj):
    Path(p).write_text(json.dumps(obj, ensure_ascii=False, indent=2, allow_nan=False) + '\n', encoding='utf-8')


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def compute_spectrum(context):
    import qnm
    seqs = {}
    for ell, n in [(2, 0), (2, 1), (2, 2), (2, 3), (3, 0), (4, 0)]:
        seq = qnm.spinsequence.KerrSpinSeq(s=-2, l=ell, m=0, n=n, a_max=.70)
        seq.do_find_sequence()
        seqs[f'{ell}0{n}'] = seq

    def compute(mass, chi, a220):
        rows = {}
        for key, seq in seqs.items():
            w, _, _ = seq(a=chi)
            cf, terms = seq.solver.get_cf_err()
            rows[key] = {'omega': float(w.real / mass), 'alpha': float(-w.imag / mass),
                         'cf_indicator': float(cf), 'cf_terms': int(terms)}
        rows['D'] = {'omega': 0., 'alpha': 2 * a220,
                     'interpretation': 'Диагностическое затухание, НЕ линейная квазинормальная мода'}
        return rows

    a2 = context['A2_spectrum']
    central = compute(MF, CHI, a2['central']['220']['alpha'])
    variations = {}
    for key, mass, chi in [('mass_minus_half_last_digit', MF - 5e-6, CHI),
                          ('mass_plus_half_last_digit', MF + 5e-6, CHI),
                          ('chi_minus_half_last_digit', MF, CHI - 5e-6),
                          ('chi_plus_half_last_digit', MF, CHI + 5e-6)]:
        variations[key] = compute(mass, chi, a2['rounding_variations'][key]['220']['alpha'])
    check = abs(central['200']['omega'] - .4132) <= 5e-5 and abs(1 / central['200']['alpha'] - 11.236) <= 5e-4
    if not check:
        raise ValueError('Спектр 200 не согласуется с округлённой таблицей статьи')
    return {'qnm_version': qnm.__version__, 'mass_M': MF, 'chi': CHI, 's': -2,
            'units': 'M^-1, M — исходная полная масса ADM', 'a_max': .70, 'l_max': 20,
            'root_tol': float(np.sqrt(np.finfo(float).eps)), 'cf_tol': 1e-10,
            'central': central, 'rounding_variations': variations, 'paper_200_check': bool(check),
            'sources': ['https://arxiv.org/pdf/2208.02965', 'https://qnm.readthedocs.io/en/latest/README.html'],
            'D_source': 'Удвоенное затухание 220 из сохранённого A2; не подгонялось к m=0.'}


def read_data(path):
    raw = path.read_bytes()
    src_sha = hashlib.sha256(raw).hexdigest()
    if zipfile.is_zipfile(io.BytesIO(raw)):
        with zipfile.ZipFile(io.BytesIO(raw)) as z:
            names = [n for n in z.namelist() if n.endswith('.h5') and not n.startswith('__MACOSX/')]
            if len(names) != 1:
                raise ValueError('Ожидался один HDF5')
            raw = z.read(names[0])
    if hashlib.sha256(raw).hexdigest() != H5_SHA:
        raise ValueError('Контрольная сумма не соответствует продолжению A1/A2')
    with h5py.File(io.BytesIO(raw), 'r') as f:
        t = f['time_adaptive'][()]
        fields = {'I20': f['geo_mass/Adaptive'][:, 6], 'I40': f['geo_mass/Adaptive'][:, 20],
                  'L30': f['geo_spin/Adaptive'][:, 12]}
    imag = {k: float(np.max(abs(z.imag))) for k, z in fields.items()}
    if max(imag.values()) > 1e-12:
        raise ValueError('Неожиданная мнимая часть m=0')
    return t, {k: z.real for k, z in fields.items()}, {'input_sha256': src_sha, 'hdf5_sha256': H5_SHA,
                                                        'max_imaginary_parts': imag}


def design(tau, model, table):
    columns, labels = [], []
    for key in MODELS[model]:
        w, a = table[key]['omega'], table[key]['alpha']
        env = np.exp(-a * tau)
        if key == 'D':
            columns.append(env)
            labels.append('D')
        else:
            columns.extend([env * np.cos(w * tau), env * np.sin(w * tau)])
            labels.extend([key + '_cos', key + '_sin'])
    return np.array(columns).T, labels


def fit(t, canonical_delta, table, model, start, duration=30, stride=1, shift=0., weighting='uniform'):
    # canonical_delta = Z - c_400. Если фон изменён на shift,
    # подгоняем Z-(c_400+shift), затем восстанавливаем полный уровень.
    train = np.flatnonzero((t >= start - 1e-9) & (t <= start + duration + 1e-9))[::stride]
    valid = np.flatnonzero((t > start + duration + 1e-9) & (t <= start + duration + 20 + 1e-9))
    if len(valid) != 200:
        raise ValueError('Неполное проверочное окно')
    tau = t[train] - start
    a, labels = design(tau, model, table)
    b, _ = design(t[valid] - start, model, table)
    weights = np.ones(len(train)) if weighting == 'uniform' else np.exp(table['200']['alpha'] * tau)
    wa = weights[:, None] * a
    norms = np.linalg.norm(wa, axis=0)
    cc, _, rank, singular = np.linalg.lstsq(wa / norms, weights * (canonical_delta[train] - shift), rcond=None)
    if rank != a.shape[1]:
        raise ValueError('Потеря ранга модели')
    c = cc / norms
    pred, forecast = a @ c + shift, b @ c + shift
    yt, yv = canonical_delta[train], canonical_delta[valid]
    row = {'model': model, 'components': MODELS[model], 'start_M': start, 'train_end_M': start + duration,
           'validation_start_M': float(t[valid[0]]), 'validation_end_M': start + duration + 20,
           'training_samples': len(train), 'validation_samples': len(valid), 'training_stride': stride,
           'training_weighting': weighting, 'baseline_shift': float(shift), 'real_parameter_count': len(c),
           'fit_relative_L2': float(np.linalg.norm(yt - pred) / np.linalg.norm(yt)),
           'validation_relative_L2': float(np.linalg.norm(yv - forecast) / np.linalg.norm(yv)),
           'validation_absolute_rms': float(np.sqrt(np.mean((yv - forecast) ** 2))),
           'validation_signal_rms': float(np.sqrt(np.mean(yv ** 2))),
           'scaled_condition_number': float(singular[0] / singular[-1]),
           'coefficients_at_start': dict(zip(labels, c.tolist()))}
    return row, (train, valid, pred, forecast)


def analyze(t, raw, spec, context):
    table = spec['central']
    baseline = {k: float(z[t >= 400].mean()) for k, z in raw.items()}
    delta = {k: z - baseline[k] for k, z in raw.items()}
    rows, controls, baseline_checks = [], [], {}
    for name, z in delta.items():
        reference = context['A1_backgrounds'][name]
        if abs(baseline[name] - reference['late_mean']) > 1e-12:
            raise ValueError('Поздний фон не совпадает с сохранённым A1')
        baseline_checks[name] = {'late_mean': baseline[name],
                                'late_fluctuation_rms': float(np.sqrt(np.mean(z[t >= 400] ** 2))),
                                'rounded_Kerr_background_shift': reference['kerr_reported_chi'] - baseline[name],
                                'L10_inferred_Kerr_background_shift': reference['kerr_L10_inferred_chi'] - baseline[name]}
        for start in STARTS:
            for model in MODELS:
                row, _ = fit(t, z, table, model, start)
                row['moment'] = name
                rows.append(row)

        # Контроли меняют по одному условию. Все относительные ошибки
        # нормированы одной и той же канонической разностью Z-c_400.
        for start in [0, 20, 60, 70]:
            for model in [NATIVE[name], 'F3O3', NATIVE[name] + 'D']:
                jobs = []
                for lo in [350, 400, 450]:
                    jobs.append(('late_baseline', str(lo), table, {'shift': float(raw[name][t >= lo].mean() - baseline[name])}))
                for label in ['kerr_reported_chi', 'kerr_L10_inferred_chi']:
                    jobs.append(('theoretical_baseline', label, table, {'shift': reference[label] - baseline[name]}))
                for stride in [1, 2, 4]:
                    jobs.append(('training_stride', str(stride), table, {'stride': stride}))
                for label, changed in spec['rounding_variations'].items():
                    jobs.append(('spectrum_parameter_rounding', label, changed, {}))
                jobs.append(('training_weights', 'compensated_decay', table, {'weighting': 'compensated_decay'}))
                if start < 70:
                    for duration in [20, 40]:
                        jobs.append(('training_duration', str(duration), table, {'duration': duration}))
                if 'D' in MODELS[model]:
                    for mult in [.8, 1., 1.2]:
                        changed = dict(table, D=dict(table['D'], alpha=table['D']['alpha'] * mult))
                        jobs.append(('diagnostic_decay_rate', str(mult), changed, {}))
                for control, value, changed, kw in jobs:
                    row, _ = fit(t, z, changed, model, start, **kw)
                    row.update(moment=name, control=control, control_value=value)
                    controls.append(row)

    # Проверка самого метода: реальный сигнал пересекает ноль, есть
    # затухающая неосциллирующая часть и ненулевая известная поправка фона.
    a, _ = design(t - 20, 'F3D', table)
    known = np.array([.3, -.2, .04, .07, .012, -.018, .08])
    synthetic = a @ known + .015
    test, _ = fit(t, synthetic, table, 'F3D', 20, shift=.015)
    if test['validation_relative_L2'] > 1e-10:
        raise ValueError('Не прошла проверка вещественной модели с фоном')

    def select(name, model, start):
        return next(r for r in rows if r['moment'] == name and r['model'] == model and r['start_M'] == start)

    cases = [select(name, model, start) for start in [0, 20, 70]
             for name in NATIVE for model in [NATIVE[name], 'F3O3', NATIVE[name] + 'D']]
    pairs = [('F1', 'F2'), ('F2', 'F3'), ('F1', 'F1O1'), ('F1O1', 'F1O2'),
             ('F1O2', 'F1O3'), ('F3', 'F3O1'), ('F3O1', 'F3O2'), ('F3O2', 'F3O3'),
             ('F1', 'F1D'), ('F2', 'F2D'), ('F3', 'F3D')]
    negative = []
    for name in NATIVE:
        for start in STARTS:
            for left, right in pairs:
                l, r = select(name, left, start), select(name, right, start)
                if r['validation_relative_L2'] > l['validation_relative_L2']:
                    negative.append({'moment': name, 'start_M': start, 'base': left, 'expanded': right,
                                     'base_validation_relative_L2': l['validation_relative_L2'],
                                     'expanded_validation_relative_L2': r['validation_relative_L2']})
    summary = {'version': VERSION, 'status': 'BH_YITIAN_A3_LATE_M0_DESCRIPTION_EARLY_PREDICTION_NOT_CLOSED',
               'native_models': NATIVE, 'models': MODELS, 'start_times_M': STARTS,
               'main_fit_count': len(rows), 'control_count': len(controls),
               'baseline': baseline_checks, 'representative_cases': cases,
               'nested_comparisons': len(NATIVE) * len(STARTS) * len(pairs),
               'expanded_model_worse_count': len(negative),
               'maximum_scaled_condition_number': max(r['scaled_condition_number'] for r in rows),
               'synthetic_prediction_relative_L2': test['validation_relative_L2'],
               'metric': '||delta-predicted_delta||_2/||delta||_2; canonical delta=Z-mean(Z[t>=400]); одинаковая нормировка во всех контролях фона',
               'validation_note': 'Отложенный интервал не используется для подгонки коэффициентов. Конечные Mf, chi и поздний фон известны из всей симуляции; это не полностью слепое предсказание.',
               'D_note': 'Экспонента с lambda=2*alpha220 — диагностическая гипотеза; её успех не идентифицирует физическую нелинейную моду.',
               'limits': ['Один уровень численного разрешения; нет оценки ошибки ранних моментов по двум сеткам.',
                          'Ошибка относится к малому возмущению, а не к полному моменту или всей геометрии.',
                          'Узлы пересечения нуля не исключаются; поточечная относительная ошибка не используется.',
                          'Проверки прореживания не заменяют сходимость пространственной сетки.',
                          'Окна перекрываются и не являются независимыми опытами.',
                          'Хорошее описание короткого участка не доказывает физическое присутствие каждой подогнанной моды.',
                          'Неполнота выбранных моделей не доказывает отклонение от ОТО или гиперформу.']}
    return summary, rows, controls, negative, delta


def make_plots(out, t, delta, spec, rows):
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10,
                         'axes.spines.top': False, 'axes.spines.right': False,
                         'axes.grid': True, 'grid.alpha': .18})
    fig, axs = plt.subplots(1, 3, figsize=(14.5, 4.8), layout='constrained')
    for ax, name in zip(axs, NATIVE):
        for model, label, color in [(NATIVE[name], 'Основные моды из статьи', '#286798'),
                                    ('F3O3', 'Три основные + три обертона', '#be6235'),
                                    (NATIVE[name] + 'D', 'Основные + плавное затухание', '#287e6c')]:
            rr = [r for r in rows if r['moment'] == name and r['model'] == model]
            ax.semilogy([r['start_M'] for r in rr], [100 * r['validation_relative_L2'] for r in rr],
                        'o-', ms=4, label=label, color=color)
        ax.set(title=rf'${name[0]}_{{{name[1:]}}}$', xlabel=r'Начало подгонки $t_0/M$',
               ylabel='Ошибка продолжения, %', ylim=(.5, 1500))
        ax.legend(fontsize=7.2, loc='upper right')
    fig.suptitle('A3: ранние осесимметричные моменты труднее предсказывать', fontsize=15, fontweight='bold')
    fig.supxlabel('Подгонка на 30M, проверка на следующих 20M. Нормировка — по возмущению над постоянным фоном.', fontsize=9)
    fig.savefig(out / 'BH_YITIAN_A3_prediction.png', dpi=155, bbox_inches='tight')
    plt.close(fig)

    fig, axs = plt.subplots(2, 3, figsize=(14.5, 7.3), layout='constrained')
    for col, name in enumerate(NATIVE):
        for row_index, start in enumerate([0, 70]):
            ax = axs[row_index, col]
            model = 'F3O3' if start == 0 else NATIVE[name]
            r, _ = fit(t, delta[name], spec['central'], model, start)
            mask = (t >= start) & (t <= start + 50)
            a, labels = design(t[mask] - start, model, spec['central'])
            c = np.array([r['coefficients_at_start'][label] for label in labels])
            scale = 1 if start == 0 else 1e6
            ax.plot(t[mask], delta[name][mask] * scale, color='#283746', lw=1.6, label='Данные')
            ax.plot(t[mask], a @ c * scale, '--', color='#b56332', lw=1.4, label='Модель')
            ax.axvspan(start + 30, start + 50, color='#dfe8f0', alpha=.8)
            ax.set(title=rf'${name[0]}_{{{name[1:]}}}$; ошибка проверки {100*r["validation_relative_L2"]:.1f}%',
                   xlabel=r'Время $t/M$', ylabel='Возмущение' if start == 0 else 'Возмущение × 10⁶')
            ax.legend(fontsize=8, loc='lower right' if start == 0 else 'upper right')
            if start == 0:
                inset = ax.inset_axes([.47, .47, .50, .43])
                check = t[mask] > 30
                inset.plot(t[mask][check], delta[name][mask][check] * 1e3,
                           color='#283746', lw=1.2)
                inset.plot(t[mask][check], (a @ c)[check] * 1e3,
                           '--', color='#b56332', lw=1.2)
                inset.set_title('Проверка 30–50M, ×10³', fontsize=7)
                inset.tick_params(labelsize=7)
    fig.suptitle('A3: сходство на подгонке и проверка за её пределами', fontsize=15, fontweight='bold')
    fig.supxlabel('Сверху: три основные моды и три обертона. Снизу: основные моды из статьи.\nЗатенение — отложенные данные; коэффициенты на них не подбирались.', fontsize=9)
    fig.savefig(out / 'BH_YITIAN_A3_examples.png', dpi=155, bbox_inches='tight')
    plt.close(fig)


def report(out, summary, rows, controls, spec):
    def get(name, model, start):
        return next(r for r in rows if r['moment'] == name and r['model'] == model and r['start_M'] == start)
    lines = ['# BH-YITIAN-A3: осесимметричная часть общего горизонта', '',
             'Дата: 15 сентября 2026 года. A3 продолжает A2. Эволюция QC0/SpEC и прежние A1/A2 не пересчитывались.', '',
             '**Результат:** поздние I20, I40, L30 описываются небольшими наборами основных керровских колебаний с ошибкой продолжения около 2.6–3.1% на отрезке 100–120M. На ранних интервалах успешная подгонка несколькими модами не обеспечивает хорошего продолжения. Добавочная неосциллирующая экспонента помогает на некоторых промежуточных отрезках, но не установлена как отдельная физическая мода.', '',
             '## Постоянный фон — часть геометрии', '',
             '$$Z_{\\ell0}(t)=Z_{\\ell0}^{\\infty}+\\delta Z_{\\ell0}(t),\\qquad Z_{\\ell0}^{\\infty}=\\operatorname{mean}_{400\\le t/M\\le500}Z_{\\ell0}(t).$$', '',
             'Исходные моменты сохранены; среднее вычитается только для анализа возмущений. При m=0 поворот exp(−imΩt) равен единице, а ряды вещественны. Поэтому особенности A3 нельзя объяснить забытым частотным сдвигом mΩ.', '',
             '| Момент | Поздний фон |', '| --- | ---: |']
    for name, b in summary['baseline'].items():
        lines.append(f"| {name} | {b['late_mean']:.12f} |")
    lines += ['', 'Согласие этих фонов с керровскими моментами уже проверено в A1; здесь соответствующие значения переиспользованы. Это безразмерные геометрические моменты, а не непосредственно масса и спин.', '',
              '## Модели и критерий сравнения', '',
              '$$\\delta Z_{\\mathrm{lin}}(t)=\\sum_k e^{-\\alpha_k(t-t_0)}[A_k\\cos\\omega_k(t-t_0)+B_k\\sin\\omega_k(t-t_0)].$$', '',
              'Частоты фиксированы численным спектром Керра при Mf=0.95162M и χ=0.68644. Для каждого колебания подбираются две вещественные амплитуды. Синусы и косинусы эквивалентны согласованной паре положительной и отрицательной частот; подгонка одной комплексной экспоненты к вещественному ряду здесь не использовалась. Точки пересечения нуля сохранены.', '',
              '| Модель | Компоненты | Число вещественных амплитуд |', '| --- | --- | ---: |']
    for model, keys in MODELS.items():
        lines.append(f"| {model} | {' + '.join(keys)} | {sum(1 if k=='D' else 2 for k in keys)} |")
    lines += ['', 'Основные модели для сопоставления со статьёй: I20 — F1 (200); L30 — F2 (200,300); I40 — F3 (200,300,400). Более богатые модели служат контролем. Наборы и окна исследовательские, без заявления о заранее зарегистрированном испытании.', '',
              '| Мода | ω, M⁻¹ | α, M⁻¹ |', '| --- | ---: | ---: |']
    for k, v in spec['central'].items():
        lines.append(f"| {k} | {v['omega']:.9f} | {v['alpha']:.9f} |")
    lines += ['', 'Обозначение D означает отдельную диагностическую экспоненту exp(−λ(t−t0)), λ=2α220. Это не линейная квазинормальная мода и не заявленное обнаружение нелинейной моды. Мотивация только алгебраическая: квадрат модуля затухающего колебания даёт exp(−2αt); произведение азимутальных множителей с m=2 и m=−2 допускает m=0. Из этой алгебры не следует, что конкретный момент обязан содержать такой член или иметь установленную амплитуду. Нужна отдельная проверка физических уравнений и устойчивости связи.', '',
              'Подгонка [t0,t0+30]M, проверка (t0+30,t0+50]M, t0=0,5,10,15,20,30,40,50,60,70M. Критерий:', '',
              '$$E=\\frac{\\|\\delta Z-\\delta Z_{\\mathrm{model}}\\|_2}{\\|\\delta Z\\|_2}.$$', '',
              'Все проценты ниже — 100E. Это ошибка возмущения на всём отрезке, не ошибка полной геометрии, не максимум поточечной ошибки и не статистическая вероятность. Коэффициенты не подбирались на отложенном отрезке. Поздний фон и конечные параметры известны по всей симуляции, поэтому проверка не является полностью слепым предсказанием будущего. Перекрывающиеся окна не считаются независимыми опытами.', '',
              '## Поздние моменты', '',
              '| Момент | Основная модель | Подгонка 70–100M: ошибка, % | Проверка 100–120M: ошибка, % | Абсолютная среднеквадратичная ошибка проверки |',
              '| --- | --- | ---: | ---: | ---: |']
    for name, model in NATIVE.items():
        r = get(name, model, 70)
        lines.append(f"| {name} | {model} | {100*r['fit_relative_L2']:.4f} | {100*r['validation_relative_L2']:.4f} | {r['validation_absolute_rms']:.4g} |")
    lines += ['', 'Для позднего I20 одной основной моды достаточно для описания в несколько процентов по этому критерию; I40 и L30 требуют смешения. Это согласуется с качественным выводом авторов, но наша ошибка продолжения не тождественна их мере несходства на подогнанном интервале.', '',
              '## Ранняя динамика: проверка против ложного успеха', '',
              '| Момент | Модель | Подгонка 0–30M: ошибка, % | Проверка 30–50M: ошибка, % |',
              '| --- | --- | ---: | ---: |']
    for name in NATIVE:
        for model in [NATIVE[name], 'F3O3']:
            r = get(name, model, 0)
            lines.append(f"| {name} | {model} | {100*r['fit_relative_L2']:.4f} | {100*r['validation_relative_L2']:.4f} |")
    lines += ['', 'Даже малая ошибка подгонки не разрешает объявить ранний режим объяснённым. Большая относительная ошибка проверки отчасти отражает малый размер уже затухшего сигнала; абсолютная ошибка и амплитуда проверочного сигнала сохранены в JSON. Результат относится к выбранным моделям и процедуре, а не является доказательством невозможности любого керровского или нелинейного описания.', '',
              '| Момент | Основная модель: проверка 50–70M, % | Та же модель + D: проверка 50–70M, % |', '| --- | ---: | ---: |']
    for name, model in NATIVE.items():
        a, b = get(name, model, 20), get(name, model + 'D', 20)
        lines.append(f"| {name} | {100*a['validation_relative_L2']:.4f} | {100*b['validation_relative_L2']:.4f} |")
    lines += ['', 'Добавление D помогает на этих промежуточных отрезках, но не закрывает самую раннюю динамику. Возможные причины остатка — нелинейная релаксация, неполный набор колебаний, изменение эффективных амплитуд, выбор временной координаты и численные эффекты. Одной успешной экспоненты недостаточно для выбора между ними.', '',
              '![Ошибка продолжения](BH_YITIAN_A3_prediction.png)', '',
              '![Примеры рядов](BH_YITIAN_A3_examples.png)', '',
              '## Проверки чувствительности', '',
              'Проверены окна позднего среднего 350–500, 400–500, 450–500M; шаг обучения 0.1,0.2,0.4M; округление Mf и χ на ±0.000005; длины обучения 20,30,40M на выбранных началах. Дополнительно проверены веса остатков exp(α200(t−t0)), уменьшающие преобладание больших ранних амплитуд. Эти веса — методический контроль, а не модель статистической погрешности.', '',
              '| Момент | Обычные веса: ранняя проверка F3O3, % | Веса с компенсацией затухания: ранняя проверка F3O3, % |', '| --- | ---: | ---: |']
    for name in NATIVE:
        r = get(name, 'F3O3', 0)
        c = next(c for c in controls if c['moment'] == name and c['model'] == 'F3O3' and c['start_M'] == 0 and c['control'] == 'training_weights')
        lines.append(f"| {name} | {100*r['validation_relative_L2']:.4f} | {100*c['validation_relative_L2']:.4f} |")
    lines += ['', 'Изменение процедуры может менять числа, поэтому отдельный остаток нельзя называть неизменной физической структурой. Полный набор результатов сохраняется, включая случаи ухудшения при расширении модели.', '',
              '| Поздняя основная модель | Фон 400–500M: ошибка проверки, % | Керровский фон при округлённом χ: ошибка проверки, % | Керровский фон при χ из L10: ошибка проверки, % |',
              '| --- | ---: | ---: | ---: |']
    for name, model in NATIVE.items():
        a = get(name, model, 70)
        rr = [c for c in controls if c['moment'] == name and c['model'] == model and c['start_M'] == 70 and c['control'] == 'theoretical_baseline']
        rounded = next(c for c in rr if c['control_value'] == 'kerr_reported_chi')
        inferred = next(c for c in rr if c['control_value'] == 'kerr_L10_inferred_chi')
        lines.append(f"| {name}, {model} | {100*a['validation_relative_L2']:.4f} | {100*rounded['validation_relative_L2']:.4f} | {100*inferred['validation_relative_L2']:.4f} |")
    lines += ['', 'Важное различие: округлённый χ может прекрасно описывать полный ненулевой момент, но давать заметную ошибку после вычитания почти равных величин. Для честного сравнения выше во всех случаях использован один знаменатель — норма канонического возмущения относительно среднего 400–500M; полные моменты восстанавливаются перед вычислением остатка. Согласие фона, восстановленного по L10, является внутренней проверкой той же симуляции.', '',
              '| Момент, модель + D | λ=0.8×2α220: проверка 50–70M, % | λ=2α220 | λ=1.2×2α220 |', '| --- | ---: | ---: | ---: |']
    for name, model in NATIVE.items():
        vals = []
        for mult in [.8, 1., 1.2]:
            c = next(c for c in controls if c['moment'] == name and c['model'] == model + 'D' and c['start_M'] == 20 and c['control'] == 'diagnostic_decay_rate' and c['control_value'] == str(mult))
            vals.append(f"{100*c['validation_relative_L2']:.4f}")
        lines.append(f"| {name} | {' | '.join(vals)} |")
    lines += ['', 'Эта проверка не является измерением λ или доверительным интервалом. Для вывода о квадратичной физической связи нужны устойчивость амплитуд, явная связь с родительскими компонентами и контроль численного разрешения.', '',
              f"Всего основных подгонок: {summary['main_fit_count']}; контрольных случаев: {summary['control_count']}. В {summary['expanded_model_worse_count']} из {summary['nested_comparisons']} вложенных сравнений расширение модели ухудшило продолжение. Максимальное число обусловленности нормированной матрицы: {summary['maximum_scaled_condition_number']:.4g}. Проверка на известном вещественном сигнале с постоянной поправкой фона дала ошибку продолжения {summary['synthetic_prediction_relative_L2']:.3g}.", '',
              'Прореживание не проверяет сходимость пространственной сетки. Позднее рассеяние не задаёт ошибку раннего сигнала. Из одного разрешения и целевого значения AMR нельзя получить достоверную погрешность каждого момента. Порог статистической значимости и вероятность новой структуры здесь не вычислялись.', '',
              '## Точка продолжения', '',
              '**A3 завершён как проверка этих трёх моментов; ранняя динамика не объявляется полностью объяснённой.** Следующий плановый шаг — A4: компоненты с m=4, начиная с I44, затем L54 и I64. Там необходимо вернуть авторский поворот exp(−4iΩt), проверить обычное смешение мод и лишь отдельно рассматривать дополнительные нелинейные шаблоны.', '',
              'Контроль QC0 95→105M остаётся в маршруте. Его сетка, контрольная точка и выполненные результаты сохраняются; параметры Yitian в QC0 не подставляются. После завершения HF3 — HF4 с проверкой того, что можно восстановить из угловых коэффициентов и что требует метрики поверхности. Научный фильм остаётся этапом после такой проверки.', '',
              '## Повторение', '',
              'Из папки распакованного пакета:', '', '```bash',
              'python BH_YITIAN_A3.py --input /путь/common_horizon_moments_data_and_scripts.zip --output repeated', '```', '',
              'Исходный архив не дублируется. Для обработки нужны numpy, scipy, h5py, matplotlib. Для пересчёта только спектра m=0 добавить --recompute-spectrum; тогда нужен qnm==0.4.4. Контекст A1/A2, спектр и точка продолжения включены в пакет.', '',
              'Источники методики: [Chen и соавторы, разделы IV C–D](https://arxiv.org/pdf/2208.02965); [документация qnm](https://qnm.readthedocs.io/en/latest/README.html). Диагностический член D и проверки продолжения — наш исследовательский вариант, а не вывод статьи о данном члене.', '']
    (out / 'BH_YITIAN_A3_report.md').write_text('\n'.join(lines), encoding='utf-8')


def main():
    root = Path(__file__).resolve().parent
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--input', type=Path, required=True)
    ap.add_argument('--output', type=Path, default=Path('BH_YITIAN_A3_results'))
    ap.add_argument('--context', type=Path, default=root / 'BH_YITIAN_A3_context.json')
    ap.add_argument('--spectrum', type=Path, default=root / 'BH_YITIAN_A3_spectrum.json')
    ap.add_argument('--recompute-spectrum', action='store_true')
    args = ap.parse_args()
    out = args.output.resolve()
    out.mkdir(parents=True, exist_ok=True)
    context = json.loads(args.context.read_text())
    t, raw, provenance = read_data(args.input)
    spec = compute_spectrum(context) if args.recompute_spectrum else json.loads(args.spectrum.read_text())
    if spec['mass_M'] != MF or spec['chi'] != CHI:
        raise ValueError('Не тот спектральный эталон')
    dump(out / 'BH_YITIAN_A3_spectrum.json', spec)
    summary, rows, controls, negative, delta = analyze(t, raw, spec, context)
    state = dict(context['resume_state_A2'])
    state.update(last_completed_Yitian_stage='BH-YITIAN-A3',
                 A3_status=summary['status'],
                 next_step='BH-YITIAN-A4: m=4, I44, затем L54/I64; поворот exp(-4i*Omega*t), обычное смешение мод, отдельный контроль нелинейных шаблонов.',
                 A3_early_closure_reached=False)
    provenance.update(version=VERSION, program_sha256=sha(__file__), context_sha256=sha(args.context),
                      python=platform.python_version(), numpy=np.__version__, scipy=scipy.__version__,
                      h5py=h5py.__version__, matplotlib=matplotlib.__version__,
                      m0_spectrum_recomputed=args.recompute_spectrum)
    for name, obj in [('BH_YITIAN_A3_summary.json', summary), ('BH_YITIAN_A3_all_fits.json', rows),
                      ('BH_YITIAN_A3_controls.json', controls), ('BH_YITIAN_A3_negative_controls.json', negative),
                      ('BH_YITIAN_A3_spectrum.json', spec), ('BH_YITIAN_A3_context.json', context),
                      ('BH_HF3_resume_state_A3.json', state), ('BH_YITIAN_A3_provenance.json', provenance)]:
        dump(out / name, obj)
    make_plots(out, t, delta, spec, rows)
    report(out, summary, rows, controls, spec)
    if Path(__file__).resolve() != out / Path(__file__).name:
        shutil.copy2(__file__, out / Path(__file__).name)
    (out / 'requirements.txt').write_text('numpy\nscipy\nh5py\nmatplotlib\n# Только для --recompute-spectrum:\n# qnm==0.4.4\n')
    print(json.dumps({'status': summary['status'], 'fits': len(rows), 'controls': len(controls),
                      'synthetic_check': summary['synthetic_prediction_relative_L2'],
                      'max_condition': summary['maximum_scaled_condition_number'],
                      'output': str(out)}, ensure_ascii=False, indent=2))


if __name__ == '__main__':
    main()
