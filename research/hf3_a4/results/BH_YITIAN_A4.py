#!/usr/bin/env python3
"""A4: компоненты m=4 горизонта, смешение и диагностический шаблон 220×220.

Использует уже проверенный вход A1–A3. Эволюцию пространства-времени не запускает.
Сохранённый спектр достаточен для повторения; qnm нужен только для его пересчёта.
"""
from __future__ import annotations
import argparse
import copy
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

VERSION = 'BH-YITIAN-A4-2026-09-15-v1'
H5_SHA = '9cdd268e85b783d8a832b72e46945c6b43bb2b3ea0a9e6665e57a5a1b7d56f0f'
MF, CHI, OMEGA = .95162, .68644, .20878380853617037
STARTS = [0, 5, 10, 15, 20, 30, 40, 50, 60, 70]
MOMENTS = {'I44': ('geo_mass/Adaptive', 4),
           'L54': ('geo_spin/Adaptive', 5),
           'I64': ('geo_mass/Adaptive', 6)}
MODES = [(4, 0), (5, 0), (6, 0), (7, 0), (8, 0), (4, 1), (4, 2), (4, 3)]
MODELS = {
    'F1': ['440'], 'F2': ['440', '540'],
    'F3': ['440', '540', '640'], 'F4': ['440', '540', '640', '740'],
    'F5': ['440', '540', '640', '740', '840'],
    'F1O3': ['440', '441', '442', '443'],
    'F2O1': ['440', '540', '441'],
    'F2O2': ['440', '540', '441', '442'],
    'F2O3': ['440', '540', '441', '442', '443'],
    'F4O3': ['440', '540', '640', '740', '441', '442', '443'],
    'Q': ['Q'], 'F1Q': ['440', 'Q'], 'F2Q': ['440', '540', 'Q'],
    'F3Q': ['440', '540', '640', 'Q'],
    'F4Q': ['440', '540', '640', '740', 'Q'],
    'F2O3Q': ['440', '540', '441', '442', '443', 'Q'],
    'F4O3Q': ['440', '540', '640', '740', '441', '442', '443', 'Q'],
}
REFERENCE = {'I44': 'F2', 'L54': 'F4', 'I64': 'F4'}
SOURCES = ['https://arxiv.org/pdf/2208.02965',
           'https://qnm.readthedocs.io/en/latest/README.html']


def dump(path, obj):
    Path(path).write_text(json.dumps(obj, ensure_ascii=False, indent=2,
                                   allow_nan=False) + '\n', encoding='utf-8')


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def cj(z):
    return {'real': float(np.real(z)), 'imag': float(np.imag(z))}


def un_cj(z):
    return z['real'] + 1j * z['imag']


def read_data(path):
    blob = Path(path).read_bytes()
    source_sha = hashlib.sha256(blob).hexdigest()
    notebook = None
    if zipfile.is_zipfile(io.BytesIO(blob)):
        with zipfile.ZipFile(io.BytesIO(blob)) as archive:
            names = [n for n in archive.namelist()
                     if n.endswith('.h5') and not n.startswith('__MACOSX/')]
            if len(names) != 1:
                raise ValueError('Ожидался один исходный HDF5')
            nbs = [n for n in archive.namelist()
                   if n.endswith('/I44.ipynb') and not n.startswith('__MACOSX/')]
            if nbs:
                nb = archive.read(nbs[0])
                notebook = {'name': nbs[0], 'sha256': hashlib.sha256(nb).hexdigest(),
                            'code': [''.join(c.get('source', [])) for c in json.loads(nb)['cells']
                                     if c.get('cell_type') == 'code']}
            blob = archive.read(names[0])
    if hashlib.sha256(blob).hexdigest() != H5_SHA:
        raise ValueError('Вход не совпал с подтверждённым HDF5 A1–A3')
    fields, parity = {}, {}
    with h5py.File(io.BytesIO(blob), 'r') as f:
        t = f['time_adaptive'][()]
        for name, (dataset, ell) in MOMENTS.items():
            fields[name] = f[dataset][:, ell*(ell+1)+4].conj()
            minus = f[dataset][:, ell*(ell+1)-4].conj()
            parity[name] = float(np.linalg.norm(minus-fields[name].conj()) /
                                 np.linalg.norm(fields[name]))
        fields['I22'] = f['geo_mass/Adaptive'][:, 8].conj()
    if t.shape != (5001,) or not np.allclose(t, np.arange(5001)*.1, atol=1e-10):
        raise ValueError('Неожиданная временная сетка')
    if not all(np.all(np.isfinite(z)) for z in fields.values()):
        raise ValueError('Некорректные числа во входе')
    return t, fields, {'input_name': Path(path).name, 'input_sha256': source_sha,
                       'hdf5_sha256': H5_SHA, 'm4_reality_relation_relative_L2': parity,
                       'author_I44_notebook': notebook}


def compute_spectrum(context, cache):
    import qnm
    cache.mkdir(parents=True, exist_ok=True)
    variants = [('central', MF, CHI),
                ('mass_minus_half_last_digit', MF-5e-6, CHI),
                ('mass_plus_half_last_digit', MF+5e-6, CHI),
                ('chi_minus_half_last_digit', MF, CHI-5e-6),
                ('chi_plus_half_last_digit', MF, CHI+5e-6)]
    result = {key: {} for key, _, _ in variants}
    for ell, n in MODES:
        key = f'{ell}4{n}'
        p = cache / f'mode_{key}.json'
        if p.exists():
            item = json.loads(p.read_text())
            if item['version'] != VERSION or item['qnm_version'] != qnm.__version__:
                raise ValueError('Несовместимое сохранение спектра')
        else:
            seq = qnm.spinsequence.KerrSpinSeq(s=-2, l=ell, m=4, n=n, a_max=.70)
            seq.do_find_sequence()
            item = {'version': VERSION, 'qnm_version': qnm.__version__, 'values': {}}
            for variant, mass, chi in variants:
                w, _, _ = seq(a=chi)
                cf, terms = seq.solver.get_cf_err()
                item['values'][variant] = {'omega': float(w.real/mass),
                                           'alpha': float(-w.imag/mass),
                                           'cf_indicator': float(cf), 'cf_terms': int(terms)}
            dump(p, item)
        for variant in result:
            result[variant][key] = item['values'][variant]
        print('Спектр сохранён:', key, flush=True)
    a2 = context['A2_spectrum']
    for variant, table in result.items():
        parent = a2['central'] if variant == 'central' else a2['rounding_variations'][variant]
        table['Q'] = {'omega': 2*parent['220']['omega'],
                      'alpha': 2*parent['220']['alpha'],
                      'meaning': 'Диагностическая сумма частот 220+220, не линейная QNM'}
    return {'version': VERSION, 'qnm_version': qnm.__version__, 's': -2, 'm': 4,
            'mass_M': MF, 'chi': CHI, 'units': 'M^-1, исходная масса ADM',
            'a_max': .70, 'l_max': 20, 'root_tol': float(np.sqrt(np.finfo(float).eps)),
            'cf_tol': 1e-10, 'central': result.pop('central'),
            'rounding_variations': result, 'sources': SOURCES,
            'note': 'Частоты QNM спинового веса -2 используются в модели статьи; '
                    'геометрические моменты разлагаются по обычным сферическим гармоникам.'}


def frequencies(table, keys):
    return np.array([table[k]['omega'] - 1j*table[k]['alpha'] for k in keys])


def design(tau, keys, table):
    return np.exp(-1j*np.outer(tau, frequencies(table, keys)))


def fit(t, canonical, table, model, start, duration=30, stride=1,
        weighting='uniform', baseline_delta=None):
    keys = MODELS[model]
    train = np.flatnonzero((t >= start-1e-9) & (t <= start+duration+1e-9))[::stride]
    valid = np.flatnonzero((t > start+duration+1e-9) & (t <= start+duration+20+1e-9))
    if not len(train) or not len(valid):
        raise ValueError('Пустой интервал')
    shift = np.zeros_like(canonical) if baseline_delta is None else baseline_delta
    target = canonical + shift
    a = design(t[train]-start, keys, table)
    b = design(t[valid]-start, keys, table)
    weights = (np.exp(table['440']['alpha']*(t[train]-start))
               if weighting == 'decay_compensation' else np.ones(len(train)))
    aw = a*weights[:, None]
    norms = np.linalg.norm(aw, axis=0)
    cs, _, rank, singular = np.linalg.lstsq(aw/norms, target[train]*weights, rcond=None)
    if rank != len(keys):
        raise ValueError('Матрица потеряла ранг')
    c = cs/norms
    p, v = a@c-shift[train], b@c-shift[valid]
    residual = canonical[valid]-v
    # Общая точка отсчёта нужна для сопоставления амплитуд между окнами.
    at20 = c*np.exp(-1j*frequencies(table, keys)*(20-start))
    row = {'model': model, 'modes': keys, 'start_M': start,
           'train_end_M': start+duration, 'validation_start_M': float(t[valid[0]]),
           'validation_end_M': start+duration+20, 'duration_M': duration,
           'training_samples': len(train), 'validation_samples': len(valid),
           'stride': stride, 'weighting': weighting, 'real_parameter_count': 2*len(keys),
           'fit_relative_L2': float(np.linalg.norm(canonical[train]-p)/np.linalg.norm(canonical[train])),
           'validation_relative_L2': float(np.linalg.norm(residual)/np.linalg.norm(canonical[valid])),
           'validation_absolute_rms': float(np.sqrt(np.mean(abs(residual)**2))),
           'validation_signal_rms': float(np.sqrt(np.mean(abs(canonical[valid])**2))),
           'scaled_condition_number': float(singular[0]/singular[-1]),
           'coefficients_at_start': {k: cj(x) for k,x in zip(keys,c)},
           'coefficients_at_reference_20M': {k: cj(x) for k,x in zip(keys,at20)}}
    return row, (train, valid, p, v)


def parent_square_control(t, z, table, base, start, parent_row, parent_table):
    """Родительская модель сохранена в A2; значения I22 на проверке не используются."""
    train = np.flatnonzero((t >= start-1e-9) & (t <= start+30+1e-9))
    valid = np.flatnonzero((t > start+30+1e-9) & (t <= start+50+1e-9))
    pk = parent_row['modes']
    pc = np.array([un_cj(parent_row['coefficients_at_start'][k]) for k in pk])
    indices = np.r_[train, valid]
    parent = design(t[indices]-start, pk, parent_table)@pc
    parent_square = parent**2
    keys = MODELS[base]
    a = np.column_stack((design(t[train]-start, keys, table), parent_square[:len(train)]))
    b = np.column_stack((design(t[valid]-start, keys, table), parent_square[len(train):]))
    norms = np.linalg.norm(a, axis=0)
    cs, _, rank, singular = np.linalg.lstsq(a/norms, z[train], rcond=None)
    if rank != len(keys)+1:
        raise ValueError('Родительский шаблон вырожден')
    c = cs/norms
    p, v = a@c, b@c
    return {'moment': None, 'model': base+'P', 'parent_model_A2': parent_row['model'],
            'start_M': start, 'train_end_M': start+30, 'validation_end_M': start+50,
            'real_parameter_count': 2*(len(keys)+1),
            'parent_parameters_fitted_in_A2': 2*len(pk),
            'fit_relative_L2': float(np.linalg.norm(z[train]-p)/np.linalg.norm(z[train])),
            'validation_relative_L2': float(np.linalg.norm(z[valid]-v)/np.linalg.norm(z[valid])),
            'validation_absolute_rms': float(np.sqrt(np.mean(abs(z[valid]-v)**2))),
            'scaled_condition_number': float(singular[0]/singular[-1]),
            'coupling_K': cj(c[-1]), 'K_absolute': float(abs(c[-1])),
            'K_phase_rad': float(np.angle(c[-1])),
            'coefficients_at_start': {k: cj(x) for k,x in zip(keys+['P'],c)},
            'validation_parent_data_used': False,
            'note': 'Квадрат продолжения модели I22 из A2; один новый комплексный коэффициент. '
                    'Это эмпирический шаблон, не выведенное уравнение связи моментов.'}


def analyze(t, raw, spec, context):
    table = spec['central']
    rot = np.exp(-4j*OMEGA*t)
    offsets = {name: raw[name][t >= 400].mean() for name in MOMENTS}
    data = {name: (raw[name]-offsets[name])*rot for name in MOMENTS}
    rows = []
    for name, z in data.items():
        for start in STARTS:
            for model in MODELS:
                row, _ = fit(t, z, table, model, start)
                row['moment'] = name
                rows.append(row)
    controls, negative, parent_rows = [], [], []
    for name, z in data.items():
        base = REFERENCE[name]
        for start in [0, 20, 50, 70]:
            for model in dict.fromkeys([base, base+'Q', 'F4O3', 'F4O3Q']):
                for floor in [350, 450, None]:
                    alternate = 0j if floor is None else raw[name][t >= floor].mean()
                    shift = (offsets[name]-alternate)*rot
                    row, _ = fit(t, z, table, model, start, baseline_delta=shift)
                    controls.append(dict(row, moment=name, control='late_mean', value=floor))
                for stride in [2, 4]:
                    row, _ = fit(t, z, table, model, start, stride=stride)
                    controls.append(dict(row, moment=name, control='stride', value=stride))
                row, _ = fit(t, z, table, model, start, weighting='decay_compensation')
                controls.append(dict(row, moment=name, control='weighting', value='decay_compensation'))
                for variation, tab in spec['rounding_variations'].items():
                    row, _ = fit(t, z, tab, model, start)
                    controls.append(dict(row, moment=name, control='parameter_rounding', value=variation))
                for duration in [20, 40]:
                    row, _ = fit(t, z, table, model, start, duration=duration)
                    controls.append(dict(row, moment=name, control='training_duration', value=duration))
            for model in [base+'Q', 'F4O3Q']:
                for parameter, values in [('omega_shift', [-.03, .03]),
                                          ('alpha_factor', [.8, 1.2])]:
                    for value in values:
                        tab = copy.deepcopy(table)
                        if parameter == 'omega_shift':
                            tab['Q']['omega'] += value
                        else:
                            tab['Q']['alpha'] *= value
                        row, _ = fit(t, z, tab, model, start)
                        controls.append(dict(row, moment=name, control='Q_'+parameter, value=value))
            for factor in [-1, 0, 1]:
                altered = (raw[name]-offsets[name])*np.exp(-4j*factor*OMEGA*t)
                row, _ = fit(t, altered, table, base, start)
                negative.append(dict(row, moment=name, control='rotation_factor', value=factor))
            # Два согласованных описания одной фазы обязаны давать ту же ошибку.
            shifted_table = copy.deepcopy(table)
            for value in shifted_table.values():
                value['omega'] -= 4*OMEGA
            row, _ = fit(t, raw[name]-offsets[name], shifted_table, base, start)
            canonical_row, _ = fit(t, z, table, base, start)
            negative.append({'moment': name, 'start_M': start, 'control': 'consistent_frame_change',
                             'validation_error_difference': abs(row['validation_relative_L2']-
                                                                canonical_row['validation_relative_L2'])})
        for start in STARTS:
            parent = next(r for r in context['A2_I22_fits']
                          if r['model'] == 'F2_O3' and r['start_M'] == start)
            row = parent_square_control(t, z, table, base, start, parent, context['A2_spectrum']['central'])
            row['moment'] = name
            parent_rows.append(row)

    lookup = {(r['moment'],r['start_M'],r['model']): r for r in rows}
    comparisons = []
    for name in data:
        for start in STARTS:
            for simple, expanded in [('F1','F2'),('F2','F3'),('F3','F4'),('F4','F5'),
                                     ('F2','F2O3'),('F4','F4O3'),('F2','F2Q'),
                                     ('F4','F4Q'),('F4O3','F4O3Q')]:
                a, b = lookup[name,start,simple], lookup[name,start,expanded]
                comparisons.append({'moment': name, 'start_M': start,
                                    'simple': simple, 'expanded': expanded,
                                    'training_improved': b['fit_relative_L2'] <= a['fit_relative_L2']+1e-12,
                                    'validation_improved': b['validation_relative_L2'] < a['validation_relative_L2']})
    coupling = []
    for name in data:
        for start in STARTS:
            row = lookup[name,start,REFERENCE[name]+'Q']
            parent = next(r for r in context['A2_I22_fits']
                          if r['model'] == 'F2_O3' and r['start_M'] == start)
            bq = un_cj(row['coefficients_at_reference_20M']['Q'])
            b220 = un_cj(parent['coefficients_at_reference_20M']['220'])
            k = bq/b220**2
            coupling.append({'moment': name, 'start_M': start, 'model': row['model'],
                             'K_Q_over_parent220_squared': cj(k), 'K_absolute': float(abs(k)),
                             'K_phase_rad': float(np.angle(k)),
                             'interpretation': 'Постоянство K проверяется, но не постулируется.'})
    # Проверка алгебры, комплексного знака и отделения постоянной поправки.
    ts = np.arange(0,50.001,.1)
    known = design(ts, MODELS['F2Q'], table)@np.array([.3+.2j, .04-.01j, .08+.03j])
    shift = (.005-.003j)*np.exp(-4j*OMEGA*ts)
    syn, _ = fit(ts, known-shift, table, 'F2Q', 0, baseline_delta=shift)
    if syn['validation_relative_L2'] > 1e-10:
        raise ValueError('Не пройдена проверка известного сигнала')
    eqerr = max(r['validation_error_difference'] for r in negative
                if r['control'] == 'consistent_frame_change')
    if eqerr > 1e-10:
        raise ValueError('Согласованная смена угловой системы изменила ошибку')
    geom = design(np.arange(0,30.001,.1), ['440','Q'], table)
    coherence = abs(np.vdot(geom[:,0],geom[:,1]))/(np.linalg.norm(geom[:,0])*np.linalg.norm(geom[:,1]))
    summary = {'version': VERSION, 'status': 'BH_YITIAN_A4_COMPLETED_INTERPRETATION_IN_REPORT',
               'moments': list(MOMENTS), 'models': MODELS, 'reference_models': REFERENCE,
               'reference_note': 'Для I44 пара 440+540 указана в статье; F4 для L54/I64 — '
                                 'выбранный нами общий контроль, а не дословная спецификация статьи.',
               'starts_M': STARTS, 'training_M': 30, 'validation_M': 20,
               'main_fit_count': len(rows), 'control_count': len(controls),
               'negative_control_count': len(negative), 'parent_square_count': len(parent_rows),
               'floor': {name: {'late_mean': cj(offsets[name]),
                                'late_fluctuation_rms': float(np.sqrt(np.mean(abs(raw[name][t>=400]-offsets[name])**2)))}
                         for name in data},
               'synthetic_validation_relative_L2': syn['validation_relative_L2'],
               'consistent_frame_change_max_error_difference': eqerr,
               'column_coherence_440_Q_30M': float(coherence),
               'max_scaled_condition_number': max(r['scaled_condition_number'] for r in rows),
               'nested_comparison_count': len(comparisons),
               'expanded_model_validation_worse_count': sum(not r['validation_improved'] for r in comparisons),
               'metric': '||Z-model||_2 / ||Z||_2 на отрезке; Z=(conj(HDF)-late_mean)*exp(-4i*Omega*t)',
               'validation_note': 'Амплитуды подбираются только на обучающем отрезке. Поздний фон '
                                  'и Mf,chi известны по полной симуляции. Окна перекрываются.',
               'Q_note': 'Частота и затухание равны удвоенным 220 из A2; это диагностический шаблон.',
               'P_note': 'Квадрат сохранённой модели I22 из A2 на том же обучающем отрезке. '
                         'Проверочные значения I22 и m=4 не используются при подгонке; '
                         'параметры родительской модели тоже учитываются в интерпретации.',
               'limits': ['Только одно численное разрешение; прореживание не заменяет сходимость сетки.',
                          'Остаток модели не доказывает новую физику или гиперформу.',
                          'Цель AMR 5e-8 не равна ошибке каждого момента.',
                          'Ни доверительные интервалы, ни статистическая значимость здесь не рассчитаны.',
                          'Скалярные моменты не задают сами по себе поверхность или её физическое вложение.']}
    return data, rows, controls, negative, parent_rows, coupling, comparisons, summary


def analyze_and_save(out, t, raw, spec, context, provenance, args):
    data, rows, controls, negative, parents, coupling, comparisons, summary = analyze(t, raw, spec, context)
    for name, value in [('all_fits', rows), ('controls', controls), ('negative_controls', negative),
                        ('parent_square', parents), ('coupling', coupling), ('nested_comparisons', comparisons),
                        ('summary', summary), ('context', context)]:
        dump(out/f'BH_YITIAN_A4_{name}.json', value)
    provenance.update(version=VERSION, program_sha256=sha(__file__),
                      context_sha256=sha(args.context), spectrum_sha256=sha(out/'BH_YITIAN_A4_spectrum.json'),
                      python=platform.python_version(), numpy=np.__version__, scipy=scipy.__version__,
                      h5py=h5py.__version__, no_new_evolution=True)
    dump(out/'BH_YITIAN_A4_provenance.json', provenance)
    print('Основных подгонок:',len(rows),'Контролей:',len(controls),flush=True)
    for name in MOMENTS:
        for start in [0,20,50,70]:
            print(name, start, [(r['model'],round(100*r['fit_relative_L2'],4),
                                round(100*r['validation_relative_L2'],4)) for r in rows
                               if r['moment']==name and r['start_M']==start and
                               r['model'] in ['F1','F2','F3','F4','F5','F2O3','F4O3','F2Q','F4Q','F4O3Q']],flush=True)
    if not args.analysis_only:
        make_deliverables(out, t, data, spec, context, provenance, rows, controls,
                          negative, parents, coupling, summary)


def make_plots(out, t, data, spec, rows, coupling):
    lookup = {(r['moment'], r['start_M'], r['model']): r for r in rows}
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10,
                         'axes.grid': True, 'grid.alpha': .2, 'figure.dpi': 130})
    fig, axes = plt.subplots(1, 3, figsize=(13.2, 4.5), layout='constrained')
    for ax, name in zip(axes, MOMENTS):
        base = REFERENCE[name]
        for model, label, color, style in [(base, base+': основные моды', '#2463aa', '-'),
                                           ('F4O3', 'F4O3: с обертонами', '#8453a3', '-.'),
                                           (base+'Q', base+'Q: с шаблоном Q', '#d58118', '--'),
                                           ('F4O3Q', 'F4O3Q: обертоны и Q', '#218268', '-')]:
            ax.semilogy(STARTS, [100*lookup[name,s,model]['validation_relative_L2'] for s in STARTS],
                        style, color=color, marker='o', markersize=3, label=label, lw=1.6)
        ax.set_title(name)
        ax.set_xlabel('Начало подгонки, t₀/M')
        ax.set_xlim(-1, 71)
        ax.legend(fontsize=7.1, loc='upper right')
    axes[0].set_ylabel('Ошибка продолжения E, %')
    fig.suptitle('A4: проверка следующего участка\nПодгонка 30M; затем проверка 20M без изменения амплитуд', fontsize=13)
    fig.savefig(out/'BH_YITIAN_A4_prediction.png', dpi=170)
    plt.close(fig)

    fig, axes = plt.subplots(2, 3, figsize=(13.2, 7.4), layout='constrained')
    for col, name in enumerate(MOMENTS):
        for ridx, start in enumerate([0, 70]):
            ax = axes[ridx,col]
            ids = np.flatnonzero((t >= start-1e-9) & (t <= start+50+1e-9))
            y = data[name][ids]
            exponent = int(np.floor(np.log10(max(abs(y)))))
            scale = 10.**(-exponent)
            ax.plot(t[ids], y.real*scale, color='#172538', label='Данные', lw=1.6)
            for model, label, color in [(REFERENCE[name], 'Основные моды', '#d58118'),
                                        ('F4O3Q', 'Моды, обертоны и Q', '#218268')]:
                row = lookup[name,start,model]
                co = np.array([un_cj(row['coefficients_at_start'][k]) for k in MODELS[model]])
                pred = design(t[ids]-start, MODELS[model], spec['central'])@co
                ax.plot(t[ids], pred.real*scale, '--', color=color, lw=1.1, label=label)
            ax.axvspan(start+30,start+50,color='#b9cee3',alpha=.26)
            ax.axvline(start+30,color='#677c90',ls=':',lw=1)
            ax.set_title(f'{name}: {start}–{start+50}M',fontsize=11)
            ax.set_ylabel(f'Re Z / 10^{exponent}')
            ax.set_xlabel('t/M')
            if ridx == 0:
                inset = ax.inset_axes([.50,.52,.47,.37])
                j = ids[t[ids] > start+30+1e-9]
                inset.plot(t[j], data[name][j].real*scale, color='#172538', lw=1.1)
                for model, color in [(REFERENCE[name],'#d58118'),('F4O3Q','#218268')]:
                    row = lookup[name,start,model]
                    co = np.array([un_cj(row['coefficients_at_start'][k]) for k in MODELS[model]])
                    p = design(t[j]-start, MODELS[model], spec['central'])@co
                    inset.plot(t[j], p.real*scale, '--',color=color,lw=.9)
                inset.tick_params(labelsize=6)
                inset.set_title('Проверочный участок',fontsize=7)
    handles, labels = axes[1,0].get_legend_handles_labels()
    fig.legend(handles,labels,loc='outside lower center',ncol=3,fontsize=9)
    fig.suptitle('A4: примеры рядов и их продолжения\nВыделенная область не используется для подгонки; ошибки вычислены по комплексному ряду',fontsize=12)
    fig.savefig(out/'BH_YITIAN_A4_examples.png',dpi=170)
    plt.close(fig)

    fig, axes = plt.subplots(1,2,figsize=(11.7,4.5),layout='constrained')
    colors = {'I44':'#2463aa','L54':'#d58118','I64':'#218268'}
    for name in MOMENTS:
        items = [r for r in coupling if r['moment']==name]
        axes[0].plot([r['start_M'] for r in items],[r['K_absolute'] for r in items],
                     'o-',color=colors[name],markersize=4,label=name)
        axes[1].plot([r['start_M'] for r in items],
                     np.unwrap([r['K_phase_rad'] for r in items]),
                     'o-',color=colors[name],markersize=4,label=name)
    axes[0].set_ylabel('|K|')
    axes[1].set_ylabel('Фаза K, рад; развёрнута для непрерывности')
    for ax in axes:
        ax.set_xlabel('Начало подгонки, t₀/M')
        ax.legend()
    fig.suptitle('A4: постоянная связь с квадратом моды 220 не установлена\nK = амплитуда Q / (амплитуда 220)²; единая точка отсчёта 20M',fontsize=12)
    fig.savefig(out/'BH_YITIAN_A4_coupling.png',dpi=170)
    plt.close(fig)


def make_report(out, spec, rows, controls, parents, coupling, summary, provenance):
    lookup = {(r['moment'],r['start_M'],r['model']):r for r in rows}
    def case(name,start,model):
        return lookup[name,start,model]
    def pc(r,field='validation_relative_L2'):
        return f"{100*r[field]:.4f}"
    lines = ['# BH-YITIAN-A4: компоненты m=4 общего горизонта', '',
             'Дата: 15 сентября 2026 года. Продолжение A3. Использован тот же подтверждённый архив Yitian. Эволюция QC0/SpEC и расчёты A1–A3 не повторялись.', '',
             '**Результат:** поздняя динамика I44, L54 и I64 согласуется с суммой керровских колебаний. Смешение особенно существенно для L54 и I64. Добавочные затухающие шаблоны улучшают некоторые промежуточные и поздние участки, однако конкретная квадратичная связь 220×220 не идентифицирована. Самая ранняя динамика не замкнута рассмотренными моделями.', '',
             '## Исходные определения', '',
             'Авторский блокнот I44.ipynb прочитан как текст; его команды не исполнялись. Обработка повторяет его соглашения:', '',
             r'$$Z_{\ell4}(t)=\left[H_{\ell4}^{*}(t)-c_{\ell4}\right]e^{-4i\Omega t},\quad c_{\ell4}=\operatorname{mean}_{400\leq t/M\leq500}H_{\ell4}^{*}(t),$$', '',
             f'Ω={OMEGA:.17f} M⁻¹. Здесь H — сохранённый коэффициент HDF5, а Z — момент после сопряжения, вычитания позднего смещения и поворота. Для сравнения используется масса Mf={MF}M и безразмерный спин χ={CHI}. Параметры QC0 сюда не подставляются.', '',
             '| Ряд | Набор HDF5 | Индекс k=ℓ(ℓ+1)+m |', '| --- | --- | ---: |',
             '| I44 | geo_mass/Adaptive | 24 |', '| L54 | geo_spin/Adaptive | 34 |',
             '| I64 | geo_mass/Adaptive | 46 |', '',
             'Сопоставление m=±4 проверяет вещественность исходного скалярного поля, а не даёт независимого численного разрешения. В отличие от A3, постоянные значения этих m≠0 моментов должны исчезать в осесимметричном пределе; здесь поздние небольшие смещения вычитаются как в авторском блокноте.', '',
             '| Момент | Модуль позднего среднего | Среднеквадратичное рассеяние 400–500M |',
             '| --- | ---: | ---: |']
    for name in MOMENTS:
        f = summary['floor'][name]
        lines.append(f"| {name} | {abs(un_cj(f['late_mean'])):.4e} | {f['late_fluctuation_rms']:.4e} |")
    lines += ['', 'Позднее рассеяние не используется как известная ошибка ранних данных.', '',
              '## Модели и способ проверки', '',
              r'$$Z_{\rm model}(t)=\sum_k C_k e^{-i(\omega_k-i\alpha_k)(t-t_0)}.$$', '',
              'Комплексные амплитуды определяются методом наименьших квадратов с нормированием столбцов. Подбираются одновременно действительная и мнимая части; частоты, затухания, Mf и χ под каждый ряд не подгоняются. QNM спинового веса −2 взяты для спектральной модели; сами моменты горизонта определены через обычные сферические гармоники.', '',
              '| Компонента | ω, M⁻¹ | α, M⁻¹ |', '| --- | ---: | ---: |']
    for key, value in spec['central'].items():
        lines.append(f"| {key} | {value['omega']:.9f} | {value['alpha']:.9f} |")
    lines += ['', 'F1–F5 — последовательные наборы основных мод: 440; 440+540; затем добавляются 640, 740 и 840. O1–O3 добавляют обертоны 441–443. Полная таблица моделей:', '',
              '| Модель | Состав | Число вещественных амплитуд |', '| --- | --- | ---: |']
    for key, modes in MODELS.items():
        lines.append(f"| {key} | {' + '.join(modes)} | {2*len(modes)} |")
    lines += ['',
              'В разделе IV D статьи авторы указывают пару 440+540 для позднего I44 и существенное смешение для L54/I64. Поэтому в качестве опорных сравнений здесь взяты F2 для I44 и общий расширенный набор F4 для L54/I64. F4 в этих двух случаях — наш контрольный выбор; статья не перечисляет для них точный набор в этом абзаце. [Chen и соавторы](https://arxiv.org/pdf/2208.02965).', '',
              'Начала t0=0,5,10,15,20,30,40,50,60,70M. Амплитуды определяются по [t0,t0+30]M, затем без изменения продолжаются на (t0+30,t0+50]M. Включены все отсчёты, в том числе вблизи нулей действительной или мнимой части.', '',
              r'$$E=\frac{\|Z-Z_{\rm model}\|_2}{\|Z\|_2}.$$', '',
              'Проценты означают 100E на всём указанном отрезке, а не поточечную погрешность геометрии и не вероятность гипотезы. Поздний фон и конечные параметры известны из всей симуляции, поэтому проверка не является полностью слепым прогнозом. Проверочные данные не участвуют в подборе амплитуд, но наборы моделей и выбор представительных окон исследовательские; независимого окончательного испытания для выбора моделей здесь нет. Перекрывающиеся окна не независимы.', '',
              '## Поздняя динамика и смешение', '',
              'Подгонка 70–100M, проверка 100–120M. Во всех столбцах показана ошибка продолжения, %:', '',
              '| Момент | F1: 440 | F2: 440+540 | F3: +640 | F4: +740 | F5: +840 |',
              '| --- | ---: | ---: | ---: | ---: | ---: |']
    for name in MOMENTS:
        lines.append('| '+name+' | '+' | '.join(pc(case(name,70,model)) for model in ['F1','F2','F3','F4','F5'])+' |')
    lines += ['', 'Сильное улучшение L54 при добавлении 540 и I64 при добавлении 540/640 подтверждает практическую важность смешения. Дальнейшее расширение основного набора не даёт сопоставимого выигрыша на этом проверочном интервале. Это не означает, что остальные компоненты строго отсутствуют.', '',
              '| Момент | Опорная модель | Подгонка, % | Проверка, % | Абсолютная среднеквадратичная ошибка проверки |',
              '| --- | --- | ---: | ---: | ---: |']
    for name in MOMENTS:
        row = case(name,70,REFERENCE[name])
        lines.append(f"| {name} | {REFERENCE[name]} | {pc(row,'fit_relative_L2')} | {pc(row)} | {row['validation_absolute_rms']:.4e} |")
    lines += ['', 'Расширенная модель F4O3Q даёт на той же проверке:', '',
              '| Момент | F4O3Q: ошибка, % |', '| --- | ---: |']
    for name in MOMENTS:
        lines.append(f'| {name} | {pc(case(name,70,"F4O3Q"))} |')
    lines += ['', 'Эта модель имеет 16 вещественных амплитуд; её успех не доказывает физическое присутствие каждой компоненты. Сначала следует проверять, устойчиво ли выделяются их индивидуальные амплитуды.', '',
              '## Самое раннее окно', '',
              'Подгонка 0–30M, проверка 30–50M:', '',
              '| Момент | Модель | Подгонка, % | Проверка, % |', '| --- | --- | ---: | ---: |']
    for name in MOMENTS:
        for model in [REFERENCE[name], 'F4O3', 'F4O3Q']:
            row = case(name,0,model)
            lines.append(f"| {name} | {model} | {pc(row,'fit_relative_L2')} | {pc(row)} |")
    lines += ['', 'Дополнительные моды уменьшают ошибку подгонки, но не гарантируют лучшего продолжения. В частности, Q ухудшает раннее продолжение F4O3 для L54 и I64. Большие относительные ошибки относятся к уже затухшему сигналу на проверке; абсолютные ошибки и его амплитуда сохранены для каждого случая. Неудача выбранных моделей не доказывает невозможности более полного описания в ОТО.', '',
              '![Ошибка продолжения](BH_YITIAN_A4_prediction.png)', '',
              '![Примеры рядов](BH_YITIAN_A4_examples.png)', '',
              '## Что именно проверяет квадратичный шаблон', '',
              r'$$Q(t)=\exp[-2i(\omega_{220}-i\alpha_{220})(t-t_0)].$$', '',
              'Если одно затухающее колебание возвести в квадрат, частота и затухание удваиваются, а азимутальный индекс складывается: 2+2=4. Это алгебраическая мотивация Q. Уравнение эволюции горизонтовых моментов и коэффициент связи из неё не следуют. Здесь Q не называется отдельной обнаруженной модой.', '',
              'Для промежуточного окна подгонки 20–50M и проверки 50–70M:', '',
              '| Момент | Основные моды, % | Основные + Q, % | Основные + квадрат модели I22 из A2, % |',
              '| --- | ---: | ---: | ---: |']
    for name in MOMENTS:
        p = next(r for r in parents if r['moment']==name and r['start_M']==20)
        lines.append(f"| {name} | {pc(case(name,20,REFERENCE[name]))} | {pc(case(name,20,REFERENCE[name]+'Q'))} | {pc(p)} |")
    lines += ['', 'Последний столбец использует уже сохранённую модель I22 из A2 (220+320+221+222+223), подогнанную только на том же отрезке. Её продолжение возводится в квадрат; добавляется один комплексный коэффициент K к основной модели m=4. Фактические значения I22 на проверочном отрезке не используются. Родительская модель имеет свои десять вещественных амплитуд, что также учитывается в интерпретации; это не сравнение моделей одинаковой общей сложности.', '',
              'Чтобы проверить специфичность Q, его частота сдвигалась на ±0.03M⁻¹ при том же затухании. Это диагностические смещения, не физическая погрешность частоты. На том же окне:', '',
              '| Момент | ωQ−0.03: ошибка, % | ωQ=2ω220 | ωQ+0.03: ошибка, % |',
              '| --- | ---: | ---: | ---: |']
    for name in MOMENTS:
        vals=[]
        for value in [-.03,.03]:
            vals.append(next(r for r in controls if r['moment']==name and r['model']==REFERENCE[name]+'Q'
                             and r['start_M']==20 and r['control']=='Q_omega_shift' and r['value']==value))
        lines.append(f"| {name} | {pc(vals[0])} | {pc(case(name,20,REFERENCE[name]+'Q'))} | {pc(vals[1])} |")
    lines += ['', 'Именно удвоенная частота не выделяется как единственный подходящий вариант. Кроме того, коэффициент K=амплитуда Q/(амплитуда 220)², приведённый к общей точке отсчёта 20M, не остаётся постоянным при сдвиге окна:', '',
              '| Момент | |K| при t0=20M | |K| при t0=40M | |K| при t0=70M |',
              '| --- | ---: | ---: | ---: |']
    for name in MOMENTS:
        values=[next(r for r in coupling if r['moment']==name and r['start_M']==s)['K_absolute'] for s in [20,40,70]]
        lines.append('| '+name+' | '+' | '.join(f'{x:.5f}' for x in values)+' |')
    lines += ['', 'Меняется и комплексная фаза. При позднем сдвиге окна физически ожидаемый квадратичный вклад быстро ослабевает, поэтому его коэффициент особенно чувствителен к другим остаткам модели. Изменчивость K означает, что данная процедура не устанавливает постоянную связь; это не доказательство отсутствия нелинейного взаимодействия.', '',
              f"Нормированное скалярное произведение столбцов 440 и Q на 30M равно {summary['column_coherence_440_Q_30M']:.6f}. Шаблоны заметно похожи, что ограничивает их раздельную интерпретацию в коротком окне.", '',
              '![Коэффициент диагностической связи](BH_YITIAN_A4_coupling.png)', '',
              '## Контроль процедуры', '',
              'Проверены средние за 350–500, 400–500, 450–500M и отсутствие вычитания; шаг подгонки 0.1,0.2,0.4M; округление Mf и χ на ±0.000005; длины обучения 20,30,40M. При изменении длины меняется и проверочный интервал: это проверка чувствительности, а не сравнение на неизменном отрезке. Дополнительно проверены веса остатков exp(α440(t−t0)) и изменение затухания Q на множители 0.8 и 1.2.', '',
              'Средние в этих рядах очень малы, и смена позднего окна практически не меняет показанные результаты. Прореживание и веса могут менять ранние ошибки; поэтому конкретный остаток не следует автоматически считать неизменной физической структурой. Все варианты сохранены, включая ухудшения.', '',
              'Проверена согласованная смена угловой системы: если убрать поворот данных и одновременно заменить частоты на ω−4Ω, ошибка должна остаться прежней. Также сохранены намеренно несогласованные варианты без поворота и с обратным знаком.', '',
              f"Максимальное изменение ошибки при согласованном преобразовании: {summary['consistent_frame_change_max_error_difference']:.3e}. Проверка известного комплексного сигнала с поправкой фона: {summary['synthetic_validation_relative_L2']:.3e}. Это проверки вычислительной процедуры, не независимые физические эксперименты.", '',
              f"Выполнено {summary['main_fit_count']} основных подгонок, {summary['control_count']} вариантов чувствительности, {summary['negative_control_count']} проверок фазы и {summary['parent_square_count']} проверок квадрата родительской модели. В {summary['expanded_model_validation_worse_count']} из {summary['nested_comparison_count']} вложенных сравнений расширение модели не улучшило продолжение. Максимальное число обусловленности нормированной основной матрицы: {summary['max_scaled_condition_number']:.1f}.", '',
              'Имеется один численный уровень разрешения; целевая ошибка AMR не служит оценкой ошибки каждого момента. Доверительные интервалы, статистическая значимость и вероятность гиперформы здесь не вычислялись.', '',
              '## Обновлённая точка продолжения', '',
              '**A4 завершён для I44, L54 и I64.** Позднее смешение мод подтверждено в рамках выполненной проверки; раннее описание и физическая идентификация нелинейного вклада остаются открытыми. Отсутствие полного раннего спектрального описания не запрещает изучать измеренную геометрию, но её визуализацию нельзя выдавать за доказанную новую физику.', '',
              'Следующий шаг основного маршрута — подготовить оставшийся контроль QC0 95→105M из подтверждённой точки it=48640, затем проверить шов 95M и третий период на r=40M. Исходный prepare_qc0_t95_recovery.py уже найден; заново искать его у пользователя не требуется. Восстановление и эволюция в A4 не запускались. Ориентир QC0 остаётся его собственным; параметры Yitian его не заменяют.', '',
              'После технического завершения HF3 — HF4: установить, какие скалярные поля можно восстановить из моментов и какого геометрического материала не хватает для физической поверхности. Научный 3D-фильм следует после этой проверки. Имеющиеся моменты относятся к общему кажущемуся горизонту; они сами по себе не задают глобальный горизонт событий или внутреннюю геометрию чёрной дыры.', '',
              'Для отдельного доказательства квадратичной связи понадобятся физически обоснованный закон связи, устойчивые амплитуды и данные другого разрешения. Эта открытая задача зафиксирована, но бесконечное добавление подгоночных мод не становится обязательным условием перехода к HF4.', '',
              '## Воспроизведение и происхождение', '',
              'Программа, вычисленный спектр, контекст A2/A3, все результаты, графики и точка продолжения входят в пакет. Исходный архив не дублируется.', '',
              '```bash', 'python BH_YITIAN_A4.py --input /путь/common_horizon_moments_data_and_scripts.zip --output repeated', '```', '',
              'Необходимы numpy, scipy, h5py и matplotlib. qnm==0.4.4 нужен только для --recompute-spectrum; новые частоты сохраняются по мере расчёта. --render-only восстанавливает отчёт и графики из готовых численных JSON в папке --output.', '',
              f"SHA256 исходного HDF5: `{H5_SHA}`.", '',
              'Методический источник: [Chen и соавторы, раздел IV D](https://arxiv.org/pdf/2208.02965). Вычисление спектра: [документация qnm](https://qnm.readthedocs.io/en/latest/README.html). Процедура продолжения, Q, квадрат родительской модели и проверки их специфичности являются нашим анализом предоставленных рядов.', '']
    # Вертикальные черты внутри обозначения модуля экранируются для Markdown-таблицы.
    lines = [line.replace('| |K| при', '| Модуль K при') if line.startswith('| Момент | |K|') else line for line in lines]
    (out/'BH_YITIAN_A4_report.md').write_text('\n'.join(lines),encoding='utf-8')


def make_deliverables(out, t, data, spec, context, provenance, rows, controls,
                      negative, parents, coupling, summary):
    summary['status'] = 'BH_YITIAN_A4_LATE_M4_MIXING_SUPPORTED_EARLY_CLOSURE_NOT_REACHED_QUADRATIC_IDENTIFICATION_NOT_ESTABLISHED'
    summary['late_reference_cases'] = [r for r in rows if r['start_M']==70 and r['model']==REFERENCE[r['moment']]]
    dump(out/'BH_YITIAN_A4_summary.json',summary)
    make_plots(out,t,data,spec,rows,coupling)
    make_report(out,spec,rows,controls,parents,coupling,summary,provenance)
    state = copy.deepcopy(context['resume_state_A3'])
    state.update(last_completed_Yitian_stage='BH-YITIAN-A4', A4_status=summary['status'],
                 A4_early_closure_reached=False, A4_quadratic_identification_established=False,
                 next_step='QC0: подготовить 95→105M из it=48640, аудит шва 95M, третий период r40 с численным спектром.',
                 existing_preparation_script='hf3_continuation/prior/Пространства/prepare_qc0_t95_recovery.py',
                 next_after_QC0='HF4: аудит геометрической реконструкции из моментов, затем научная визуализация.',
                 open_research='Ранняя динамика m=0/m=4; нелинейные связи не идентифицированы. '
                               'Это не основание приписывать остаткам новую физику.',
                 A4_source_HDF5_sha256=H5_SHA)
    dump(out/'BH_HF3_resume_state_A4.json',state)
    (out/'requirements.txt').write_text('numpy\nscipy\nh5py\nmatplotlib\n# Только для --recompute-spectrum:\n# qnm==0.4.4\n',encoding='utf-8')
    (out/'README_RU.txt').write_text(
        'BH-YITIAN-A4\n\nНачните с BH_YITIAN_A4_report.md.\n'
        'Численные JSON содержат все варианты, а не только успешные.\n'
        'BH_HF3_resume_state_A4.json — точка продолжения; A1–A3 и эволюция не повторяются.\n\n'
        'Повторение с исходным архивом Yitian:\n'
        'python BH_YITIAN_A4.py --input /путь/common_horizon_moments_data_and_scripts.zip --output repeated\n\n'
        'Частоты уже включены. Их пересчёт: добавить --recompute-spectrum (нужен qnm==0.4.4).\n'
        'Только оформление готовых JSON: --render-only с прежней папкой --output.\n'
        'Эволюция Einstein Toolkit/SpEC не запускается.\n',encoding='utf-8')
    dest = out/'BH_YITIAN_A4.py'
    if Path(__file__).resolve() != dest.resolve():
        shutil.copy2(__file__, dest)
    manifest = {p.name: {'bytes': p.stat().st_size,'sha256': sha(p)}
                for p in sorted(out.iterdir()) if p.is_file() and p.name != 'BH_YITIAN_A4_manifest.json'}
    dump(out/'BH_YITIAN_A4_manifest.json',manifest)
    bundle = out.parent/'BH_YITIAN_A4_results.zip'
    with zipfile.ZipFile(bundle,'w',zipfile.ZIP_DEFLATED) as archive:
        for p in sorted(out.iterdir()):
            if p.is_file():
                archive.write(p,p.name)
    print(summary['status'],flush=True)
    print('Пакет:',bundle,flush=True)


def main():
    here = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', required=True, type=Path)
    parser.add_argument('--output', required=True, type=Path)
    parser.add_argument('--context', type=Path, default=here/'BH_YITIAN_A4_context.json')
    parser.add_argument('--spectrum', type=Path, default=here/'BH_YITIAN_A4_spectrum.json')
    parser.add_argument('--recompute-spectrum', action='store_true')
    parser.add_argument('--spectrum-only', action='store_true')
    parser.add_argument('--analysis-only', action='store_true')
    parser.add_argument('--render-only', action='store_true')
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    context = json.loads(args.context.read_text())
    t, raw, provenance = read_data(args.input)
    if args.recompute_spectrum:
        spec = compute_spectrum(context, args.output/'spectrum_checkpoints')
    else:
        spec = json.loads(args.spectrum.read_text())
    dump(args.output/'BH_YITIAN_A4_spectrum.json', spec)
    if args.spectrum_only:
        print('Спектр A4 готов.', flush=True)
        return
    if args.render_only:
        def saved(name):
            return json.loads((args.output/f'BH_YITIAN_A4_{name}.json').read_text())
        provenance.update(version=VERSION,program_sha256=sha(__file__),
                          context_sha256=sha(args.context), spectrum_sha256=sha(args.output/'BH_YITIAN_A4_spectrum.json'),
                          python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__,
                          h5py=h5py.__version__,no_new_evolution=True,rendered_from_saved_results=True)
        dump(args.output/'BH_YITIAN_A4_provenance.json',provenance)
        data={k:(raw[k]-raw[k][t>=400].mean())*np.exp(-4j*OMEGA*t) for k in MOMENTS}
        make_deliverables(args.output,t,data,spec,context,provenance,saved('all_fits'),saved('controls'),
                          saved('negative_controls'),saved('parent_square'),saved('coupling'),saved('summary'))
        return
    analyze_and_save(args.output, t, raw, spec, context, provenance, args)


if __name__ == '__main__':
    main()
