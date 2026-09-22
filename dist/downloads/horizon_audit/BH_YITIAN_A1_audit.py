#!/usr/bin/env python3
"""Воспроизводимый первичный разбор архива мультиполей Yitian.

Пример: python BH_YITIAN_A1_audit.py --input common_horizon_moments_data_and_scripts.zip
Зависимости: numpy, scipy, h5py, matplotlib. Исходные блокноты не исполняются.
Результаты — проверки входа и описательные оценки, без статистического
утверждения об открытии новой физики или независимой проверке сходимости сетки.
"""
from __future__ import annotations

import argparse
import hashlib
import io
import json
import platform
import zipfile
from pathlib import Path

import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import scipy
from scipy.special import eval_legendre

OMEGA = 0.20878380853617037  # Из I44.ipynb и L32.ipynb, без подгонки.
CHI_REPORTED = 0.68644      # Округлённое значение из письма Yitian.
MASS_REPORTED = 0.95162     # Единица M: начальная полная масса ADM, по статье.
QNM_OMEGA = 0.5535         # Округлённая таблица II, arXiv:2208.02965.
QNM_ALPHA = 1.0 / 11.707   # Обратное табличное время затухания.
PAPER = "https://arxiv.org/pdf/2208.02965"
KERR_PAPER = "https://arxiv.org/pdf/2602.05823"
ISOLATED_PAPER = "https://arxiv.org/pdf/gr-qc/0401114"


def idx(ell, m):
    return ell * (ell + 1) + m


def complex_json(z):
    return {"real": float(np.real(z)), "imag": float(np.imag(z))}


def dump(path, obj):
    path.write_text(json.dumps(obj, ensure_ascii=False, indent=2, allow_nan=False) + "\n", encoding="utf-8")


def read_input(path):
    raw = path.read_bytes()
    meta = {"input_name": path.name, "input_sha256": hashlib.sha256(raw).hexdigest(),
            "input_bytes": len(raw), "notebook_sources": []}
    if zipfile.is_zipfile(io.BytesIO(raw)):
        with zipfile.ZipFile(io.BytesIO(raw)) as z:
            bad = z.testzip()
            if bad:
                raise ValueError(f"Ошибка целостности ZIP: {bad}")
            names = [n for n in z.namelist() if not n.startswith("__MACOSX/") and not n.endswith("/")]
            h5names = [n for n in names if n.endswith(".h5")]
            if len(h5names) != 1:
                raise ValueError("Ожидался один HDF5-файл")
            data = z.read(h5names[0])
            meta.update({"zip_crc_ok": True, "payload_files": names, "hdf5_name": h5names[0]})
            for name in names:
                if name.endswith(".ipynb"):
                    nbraw = z.read(name)
                    nb = json.loads(nbraw)
                    meta["notebook_sources"].append({"name": name, "sha256": hashlib.sha256(nbraw).hexdigest(),
                        "cells": ["".join(c.get("source", [])) for c in nb["cells"]]})
    else:
        data = raw
        meta["hdf5_name"] = path.name
    meta["hdf5_sha256"] = hashlib.sha256(data).hexdigest()
    with h5py.File(io.BytesIO(data), "r") as f:
        inventory = []
        def visit(name, obj):
            inventory.append({"name": name, "type": "dataset" if isinstance(obj, h5py.Dataset) else "group",
                "shape": list(obj.shape) if isinstance(obj, h5py.Dataset) else None,
                "attributes": {k: str(v) for k, v in obj.attrs.items()}})
        f.visititems(visit)
        meta["inventory"] = inventory
        meta["root_attributes"] = {k: str(v) for k, v in f.attrs.items()}
        t = f["time_adaptive"][()]
        fields = {"I": f["geo_mass/Adaptive"][()], "L": f["geo_spin/Adaptive"][()]}
    # Авторские блокноты используют conj() для m != 0. Сохраняем также
    # малую мнимую часть m=0 для проверки, вместо её молчаливого удаления.
    return t, {name: a.conj() for name, a in fields.items()}, meta


def kerr_moments(chi, order=256):
    """Только стационарный осесимметричный эталон; не формула раннего горизонта.

    b=a/r_+=chi/(1+sqrt(1-chi^2));
    I_l0+i L_l0 = pi*(1+b^2)^2 integral Y_l0(x)/(1-i*b*x)^3 dx.
    Знак L соответствует положительному L10 в данных.
    """
    b = chi / (1.0 + np.sqrt(1.0 - chi * chi))
    x, w = np.polynomial.legendre.leggauss(order)
    return np.array([np.pi * (1+b*b)**2 * np.sum(w * (1-1j*b*x)**-3 *
        np.sqrt((2*ell+1)/(4*np.pi)) * eval_legendre(ell, x)) for ell in range(9)])


def phase_envelope_fit(t, z, start, end):
    use = (t >= start) & (t <= end)
    tt = t[use]
    zz = z[use]
    if len(tt) < 3 or np.any(np.abs(zz) == 0):
        raise ValueError("Недостаточно точек для оценки фазы и огибающей")
    x = tt - tt.mean()
    phase = np.unwrap(np.angle(zz))
    logamp = np.log(np.abs(zz))
    pp = np.polyfit(x, phase, 1)
    pa = np.polyfit(x, logamp, 1)
    model = np.exp(np.polyval(pa, x) + 1j*np.polyval(pp, x))
    return {"start_M": start, "end_M": end, "n": len(tt),
        "omega_signed_M_inverse": float(-pp[0]), "alpha_M_inverse": float(-pa[0]),
        "phase_residual_rms_rad": float(np.sqrt(np.mean((phase-np.polyval(pp, x))**2))),
        "log_amplitude_residual_rms": float(np.sqrt(np.mean((logamp-np.polyval(pa, x))**2))),
        "complex_relative_L2_residual": float(np.linalg.norm(zz-model)/np.linalg.norm(zz))}


def analyze(t, fields):
    expected_t = np.arange(5001) * 0.1
    checks = {"expected_shapes": bool(t.shape == (5001,) and all(a.shape == (5001,81) for a in fields.values())),
        "finite_values": bool(np.isfinite(t).all() and all(np.isfinite(a).all() for a in fields.values())),
        "strictly_increasing_time": bool(np.all(np.diff(t) > 0)),
        "expected_time_grid": bool(t.shape == expected_t.shape and np.allclose(t, expected_t, rtol=0, atol=1e-10))}
    if not all(checks.values()):
        raise ValueError(f"Вход отличается от проверенного набора: {checks}")
    late = t >= 400
    asym = {name: a[late].mean(axis=0) for name,a in fields.items()}
    modal = []
    parity = {}
    for name, a in fields.items():
        allow, forbid, reality = [], [], []
        for ell in range(9):
            for m in range(-ell, ell+1):
                k = idx(ell,m)
                permitted = ell % 2 == (0 if name == "I" else 1) and m % 2 == 0
                (allow if permitted else forbid).append(k)
                z = a[:,k]
                c = asym[name][k]
                modal.append({"kind": name,"ell":ell,"m":m,"column":k,"symmetry_allowed":permitted,
                    "initial":complex_json(z[0]),"final":complex_json(z[-1]),"peak_abs":float(np.max(abs(z))),
                    "late_mean":complex_json(c),"late_rms_about_mean":float(np.sqrt(np.mean(abs(z[late]-c)**2)))})
                if m > 0:
                    reality.append(np.max(abs(a[:,idx(ell,-m)]-(-1)**m*z.conj())))
        peak = np.max(abs(a[:,allow]))
        fpeak = np.max(abs(a[:,forbid]))
        parity[name] = {"forbidden_peak_abs":float(fpeak),"allowed_peak_abs":float(peak),
            "global_peak_ratio":float(fpeak/peak), "reality_max_abs_residual":float(max(reality)),
            "normalization_note":"Отношение глобальных максимумов, не поточечная относительная ошибка и не метрика HF3-R2."}
    v = asym["L"][idx(1,0)].real / np.sqrt(3*np.pi)
    chi_inferred = float(2*v/(1+v*v))
    theory = kerr_moments(CHI_REPORTED)
    inferred_theory = kerr_moments(chi_inferred)
    # Повтор квадратуры нужен для оценки погрешности интегрирования высоких l.
    theory_512 = kerr_moments(CHI_REPORTED,512)
    kr = []
    for ell in range(9):
        kind = "I" if ell % 2 == 0 else "L"
        obs = float(asym[kind][idx(ell,0)].real)
        ref = float(theory[ell].real if kind=="I" else theory[ell].imag)
        fitted = float(inferred_theory[ell].real if kind=="I" else inferred_theory[ell].imag)
        kr.append({"kind":kind,"ell":ell,"m":0,"late_mean":obs,"kerr_reported_chi":ref,
            "relative_difference_reported_chi":abs(obs-ref)/abs(ref),
            "kerr_L10_inferred_chi":fitted,"absolute_difference_inferred_chi":abs(obs-fitted),
            "relative_difference_inferred_chi":abs(obs-fitted)/abs(fitted),
            "used_to_infer_chi":ell==1})
    z = fields["I"][:,idx(2,2)]
    c = asym["I"][idx(2,2)]
    bar = z-c
    tilted = bar*np.exp(-2j*OMEGA*t)
    windows = [(10,30),(20,50),(30,60),(40,80),(50,100),(70,120)]
    fits = []
    for start,end in windows:
        row = phase_envelope_fit(t, tilted,start,end)
        before = phase_envelope_fit(t,bar,start,end)
        row["omega_before_rotation"] = before["omega_signed_M_inverse"]
        row["rotation_shift_identity_error"] = row["omega_signed_M_inverse"]-before["omega_signed_M_inverse"]-2*OMEGA
        row["relative_omega_difference_from_rounded_table"] = abs(row["omega_signed_M_inverse"]-QNM_OMEGA)/QNM_OMEGA
        row["relative_alpha_difference_from_rounded_table"] = abs(row["alpha_M_inverse"]-QNM_ALPHA)/QNM_ALPHA
        use = (t>=start)&(t<=end)
        row["minimum_corrected_amplitude_over_abs_offset"] = float(np.min(abs(bar[use]))/abs(c))
        fits.append(row)
    sensitivity = []
    for start in [350,400,450]:
        cc = z[t>=start].mean()
        rr = phase_envelope_fit(t,(z-cc)*np.exp(-2j*OMEGA*t),40,80)
        sensitivity.append({"baseline_start_M":start,"baseline_end_M":500,
            "offset":complex_json(cc),"offset_delta_from_400":float(abs(cc-c)),"fit":rr})
    audit = {"status":"BH_YITIAN_A1_INPUT_VERIFIED_LATE_MOMENTS_KERR_CONSISTENT_DESCRIPTIVE_RINGDOWN",
        "scope":"Первичный разбор предоставленного набора; не воспроизведение эволюции SpEC, не новая проверка QC0 и не обнаружение GW150914.",
        "checks":checks,"samples":len(t),"time_min_M":float(t[0]),"time_max_M":float(t[-1]),
        "time_step_M":float(np.median(np.diff(t))),"moments_per_kind":81,"total_complex_samples":int(sum(a.size for a in fields.values())),
        "late_samples":int(late.sum()),"phase_I22_at_zero":complex_json(z[0]),
        "I00_max_abs_error_sqrt_pi":float(np.max(abs(fields['I'][:,0]-np.sqrt(np.pi)))),
        "parity":parity,"rotation_Omega_M_inverse":OMEGA,"I22_offset":complex_json(c),"I22_offset_abs":float(abs(c)),
        "rotation_modulus_relative_max_error":float(np.max(abs(abs(tilted)-abs(bar)))/np.max(abs(bar))),
        "reported_mass_M":MASS_REPORTED,"reported_chi":CHI_REPORTED,"chi_inferred_from_L10":chi_inferred,
        "chi_absolute_difference_from_reported":abs(chi_inferred-CHI_REPORTED),
        "mass_independently_measured":False,"kerr_axisymmetric_comparison":kr,
        "kerr_quadrature_256_512_max_abs_difference":float(np.max(abs(theory-theory_512))),
        "QNM_reference":{"omega":QNM_OMEGA,"alpha":QNM_ALPHA,"tau":11.707,"source":PAPER,"precision":"Округлённая таблица II; не точный независимый расчёт частот."},
        "phase_envelope_fits_I22":fits,"baseline_sensitivity":sensitivity,
        "limitations":["Один файл мультиполей: независимая проверка сходимости сетки недоступна.",
            "Целевая ошибка AMR 5e-8 не является погрешностью каждого момента.",
            "HDF5 не содержит атрибутов с идентификатором симуляции, площадью, массой и спином.",
            "SXS:BBH:0389 указан в статье; принадлежность архива связывается с письмом, параметрами и скриптами, а не с метаданными HDF5.",
            "В архиве отсутствуют Psi4, h, поверхностная сетка, q_AB, X^a и A(t).",
            "Одноэкспоненциальные оценки описательны; остаток может включать смешение мод, калибровочные и численные эффекты.",
            "Данные t>150M после вычитания уровня не следует считать доказательством надёжного физического сигнала.",
            "Позднее согласие не доказывает точный Керр в ранней динамике."]}
    return audit,modal,asym


def figures(out,t,fields,asym,audit):
    plt.rcParams.update({"font.family":"DejaVu Sans","font.size":10,"axes.spines.top":False,
        "axes.spines.right":False,"axes.grid":True,"grid.alpha":0.18,"figure.facecolor":"white"})
    blue, orange, green = "#19588b", "#c36a1b", "#237460"
    fig,axs=plt.subplots(2,2,figsize=(13.4,8.5),layout="constrained")
    fig.suptitle("Общий горизонт: проверка предоставленных мультиполей",fontsize=17,fontweight="bold")
    z=fields["I"][:,idx(2,2)]; bar=z-asym["I"][idx(2,2)]; tilted=bar*np.exp(-2j*OMEGA*t)
    ax=axs[0,0]
    ax.semilogy(t,abs(z),color=blue,label=r"Исходный $|I_{22}|$")
    ax.semilogy(t,abs(bar),color=orange,label=r"После вычитания $|\bar I_{22}|$")
    ax.axhline(abs(asym["I"][idx(2,2)]),color="#69717a",ls=":",label="Модуль вычитаемого уровня")
    ax.axvspan(150,180,color="#e5e8eb",alpha=.7)
    ax.set(xlim=(0,180),ylim=(1e-8,2),xlabel=r"Время $t/M$",ylabel="Безразмерный модуль",
           title="Затухание и постоянный остаток")
    ax.legend(fontsize=9,loc="upper right")
    ax=axs[0,1]; mask=(t>=20)&(t<=100)
    for vals,label,col in [(bar,r"До поворота: $\bar I_{22}$",blue),(tilted,r"После поворота: $\tilde I_{22}$",orange)]:
        ph=np.unwrap(np.angle(vals[mask])); ax.plot(t[mask],ph-ph[0],label=label,color=col)
    ax.set(xlabel=r"Время $t/M$",ylabel="Изменение фазы, рад",title="Выбор угловой системы меняет частоту")
    ax.legend(fontsize=9)
    ax=axs[1,0]
    for kind,ell,col in [("I",2,blue),("I",4,orange),("L",3,green)]:
        k=idx(ell,0); ax.plot(t,fields[kind][:,k].real,color=col,label=rf"${kind}_{{{ell}0}}$")
        ax.axhline(asym[kind][k].real,color=col,ls=":",alpha=.8)
    ax.set(xlim=(0,80),xlabel=r"Время $t/M$",ylabel="Безразмерный момент",title="Осесимметричные моменты остаются ненулевыми")
    ax.legend(fontsize=9)
    ax=axs[1,1]; rows=audit['kerr_axisymmetric_comparison'][1:]
    ax.bar(np.arange(8),[r['relative_difference_reported_chi']*100 for r in rows],color=[blue if r['kind']=='I' else green for r in rows],width=.65)
    ax.set_xticks(np.arange(8),[rf"${r['kind']}_{{{r['ell']}0}}$" for r in rows])
    ax.set(ylabel="Различие с Керром, %",title=r"Поздние значения: сравнение при $\chi=0{,}68644$")
    ax.ticklabel_format(axis="y",style="sci",scilimits=(-3,-3),useMathText=True)
    fig.text(.5,-.014,"Расчёт по архиву Yitian. Поворот и позднее вычитание взяты из его скриптов.\nСогласие поздних моментов не заменяет проверку численной сходимости.",ha="center",fontsize=9,color="#4a535b")
    fig.savefig(out/'BH_YITIAN_A1_overview.png',dpi=155,bbox_inches="tight")
    plt.close(fig)
    rows=audit['phase_envelope_fits_I22']
    fig,axs=plt.subplots(1,2,figsize=(12,4.3),layout="constrained")
    fig.suptitle(r"Описательные оценки $\tilde I_{22}$ в нескольких временных окнах",fontsize=14)
    labels=[f"{r['start_M']}–{r['end_M']}" for r in rows]
    for ax,key,ref,title in [(axs[0],'omega_signed_M_inverse',QNM_OMEGA,"Частота"),(axs[1],'alpha_M_inverse',QNM_ALPHA,"Коэффициент затухания")]:
        ax.plot(np.arange(len(rows)),[r[key] for r in rows],"o-",color=blue,label="Оценка из данных")
        ax.axhline(ref,color=orange,ls="--",label="Керровская мода (2, 2, 0), таблица II")
        ax.set_xticks(np.arange(len(rows)),labels,rotation=25)
        ax.set(xlabel="Границы окна, M",ylabel=r"$M^{-1}$",title=title)
        ax.legend(fontsize=8)
    fig.text(.5,-.065,"Наклоны развёрнутой фазы и логарифма модуля; равные веса точек.\nЭто не доверительные интервалы и не доказательство отсутствия других мод.",ha="center",fontsize=9)
    fig.savefig(out/'BH_YITIAN_A1_windows.png',dpi=155,bbox_inches="tight")
    plt.close(fig)


def report(out,a,meta):
    krows="\n".join(f"| {r['kind']}{r['ell']}0 | {r['late_mean']:.11g} | {r['kerr_reported_chi']:.11g} | {100*r['relative_difference_reported_chi']:.6g} |" for r in a['kerr_axisymmetric_comparison'])
    frows="\n".join(f"| {r['start_M']}–{r['end_M']} | {r['omega_before_rotation']:.6f} | {r['omega_signed_M_inverse']:.6f} | {r['alpha_M_inverse']:.6f} | {100*r['relative_omega_difference_from_rounded_table']:.3f} | {100*r['relative_alpha_difference_from_rounded_table']:.3f} |" for r in a['phase_envelope_fits_I22'])
    p=a['parity']; c=a['I22_offset']
    text=rf"""# Разбор данных общего горизонта Yitian — BH-YITIAN-A1

Дата разбора: 15 сентября 2026 года. Исходный архив: `{meta['input_name']}`.

**Результат:** файл целостен и пригоден для анализа. Поздние осесимметричные моменты согласуются с керровским эталоном; корректно повёрнутый I22 имеет частоту и затухание, близкие к ожидаемым. Это воспроизводимая проверка предоставленных рядов. Расчёт эволюции SpEC здесь заново не выполнялся. HF3/QC0 и наблюдения GW150914 этим разбором не перепроверены.

## 1. Что прочитано и проверено

Архив содержит `Moments_500M_L30.h5` и четыре блокнота: `I40.ipynb`, `I44.ipynb`, `L30.ipynb`, `L32.ipynb`. Пример I22 в блокнотах отсутствует, но соответствующий столбец данных есть. Имя `L30` в имени HDF5 не означает ограничение данных моментом (3,0).

| Величина | Результат |
| --- | --- |
| Время | 0–500 M, шаг 0,1 M, 5001 отсчёт |
| Массовые моменты | `geo_mass/Adaptive`, 5001 × 81, комплексные |
| Спиновые моменты | `geo_spin/Adaptive`, 5001 × 81, комплексные |
| Диапазон индексов | 0 ≤ ℓ ≤ 8, −ℓ ≤ m ≤ ℓ |
| Индекс столбца с нуля | k = ℓ(ℓ+1)+m |
| Пропуски, бесконечности, повтор времени | Не обнаружены |
| Всего комплексных значений | {a['total_complex_samples']} |
| Позднее окно среднего | 400–500 M, 1001 отсчёт |

81 столбец на семейство не означает 81 независимую физическую степень свободы: есть сопряжённые пары и симметрийные ограничения. Содержимое блокнотов прочитано как текст; их программный код не исполнялся.

Технические имена, контрольные суммы и полный перечень объектов HDF5 сохранены в `BH_YITIAN_A1_provenance.json`. В HDF5 нет атрибутов, удостоверяющих параметры симуляции.

## 2. Происхождение и объект

В статье Yitian указана конфигурация SXS:BBH:0389: равные массы, исходно невращающиеся чёрные дыры, почти круговая орбита. Общий горизонт здесь — внешняя маргинально захваченная поверхность; это не восстановленный глобальный горизонт событий. Гармоники переносятся назад с позднего горизонта вдоль вектора сшивки Xᵃ. [Статья, разделы II–III]({PAPER}).

Связь архива с этой работой поддерживают письмо, совпадающие параметры и скрипты. Встроенного идентификатора SXS в HDF5 нет. Письмо задаёт t=0 при первом появлении общего горизонта, Mf=0,95162 M, χf=0,68644 и целевую ошибку адаптивного уточнения около 5×10⁻⁸. Это не означает, что каждый мультиполь измерен с абсолютной или относительной ошибкой 5×10⁻⁸.

Архив содержит только послеслияние общего горизонта. В нём нет двух отдельных горизонтов до слияния, волновых рядов Ψ₄ или h, площади A(t), метрики поверхности qAB, поверхностной сетки и вектора Xᵃ. Полную связь «горизонт–волна» на одном физическом случае по этому архиву отдельно проверить нельзя.

## 3. Три представления: точное воспроизведение авторских скриптов

Для обозначения любого массового или спинового ряда используем Zℓm. Последовательность для m≠0:

\[
Z_{{\ell m}}(t)=\operatorname{{conj}}\bigl(H[:,\ell(\ell+1)+m]\bigr),
\qquad c_{{\ell m}}=\frac1{{1001}}\sum_{{t_j\ge400M}}Z_{{\ell m}}(t_j),
\]
\[
\bar Z_{{\ell m}}(t)=Z_{{\ell m}}(t)-c_{{\ell m}},
\qquad \tilde Z_{{\ell m}}(t)=\bar Z_{{\ell m}}(t)e^{{-im\Omega t}},
\qquad \Omega={OMEGA}\,M^{{-1}}.
\]

Число Ω взято непосредственно из I44.ipynb и L32.ipynb, не подбиралось по полученному результату. Операция сопряжения существенна: её пропуск меняет знак фазы. При соглашении Z∝exp[−αt−iωt] поворот даёт ω̃=ω̄+mΩ и сохраняет |Z̄|=|Z̃|. Следовательно, α огибающей сохраняется. Вычитание постоянного уровня, в отличие от чистого поворота, меняет модуль исходного ряда.

Для m=0 авторские примеры берут вещественную часть; поворот тождественен единице. В названиях переменных и подписях I40/L30 встречаются оба обозначения bar/tilde для одного и того же вычитания. Численного противоречия здесь нет.

**Нельзя удалять поздние константы при проверке полной геометрии.** Для осесимметричных I20, I40, L10, L30 и других разрешённых моментов поздняя константа — часть конечного состояния. Вычитание полезно для изучения отклонений от него. Все исходные значения сохранены в анализе.

Для I22 получено c={c['real']:.12g}{c['imag']:+.12g}i, |c|={a['I22_offset_abs']:.6g}. Это наблюдаемый постоянный остаток, а не автоматически новая физическая структура. Простое вычитание остатка не доказывает достоверность слабого хвоста. В статье ограничивают физический анализ приблизительно t≤150 M. [Раздел IV A]({PAPER}).

## 4. Симметрия и нормировка

Проверены разрешённые секторы: I — чётные ℓ и чётные m; L — нечётные ℓ и чётные m. Результаты по всем временам:

| Проверка | Массовые I | Спиновые L |
| --- | --- | --- |
| Максимум модуля запрещённого коэффициента | {p['I']['forbidden_peak_abs']:.6g} | {p['L']['forbidden_peak_abs']:.6g} |
| Отношение к глобальному максимуму разрешённого сектора | {p['I']['global_peak_ratio']:.6g} | {p['L']['global_peak_ratio']:.6g} |
| Максимальная ошибка Zℓ,−m=(−1)^m Zℓm* | {p['I']['reality_max_abs_residual']:.6g} | {p['L']['reality_max_abs_residual']:.6g} |

Отношение глобальных максимумов не является поточечной относительной погрешностью и не совпадает с метрикой прежнего HF3-R2. Числа разных проверок нельзя непосредственно ранжировать как качество двух симуляций.

Максимальное абсолютное отклонение I00 от √π равно {a['I00_max_abs_error_sqrt_pi']:.6g}. Это проверка нормировки и интегральной кривизны. Геометрические Iℓm, Lℓm безразмерны; I00 нельзя читать как массу, а L10 — как χ без преобразования. [Определения изолированного горизонта]({ISOLATED_PAPER}).

## 5. Поздний Керр: проверка без вычитания констант

В осесимметричном стационарном эталоне положим b=χ/(1+√(1−χ²)). При ориентации, в которой L10 положителен, используем

\[
I_{{\ell0}}^{{\rm K}}+iL_{{\ell0}}^{{\rm K}}
=\pi(1+b^2)^2\int_{{-1}}^1
\frac{{\sqrt{{(2\ell+1)/(4\pi)}}P_\ell(x)}}{{(1-ibx)^3}}\,dx.
\]

Это осесимметричное определение горизонтовых моментов, а не формула мультиполей поля на бесконечности. Оно следует из стационарного керровского Ψ₂ и поверхностной меры; к раннему динамическому горизонту эта стационарная формула не применяется. [Керровские горизонтовые моменты, осесимметричное определение, раздел VI A]({KERR_PAPER}).

Квадратура Гаусса–Лежандра: 256 узлов. Проверка 512 узлами изменила теоретические значения не более чем на {a['kerr_quadrature_256_512_max_abs_difference']:.3g} по абсолютной величине.

Сравнение средних за 400–500 M с эталоном при сообщённом χ=0,68644:

| Момент | Среднее из данных | Керр | Относительное различие, % |
| --- | --- | --- | --- |
{krows}

Кроме того, L10=√(3π)b позволяет восстановить χ=2b/(1+b²). Из ряда получено **χ={a['chi_inferred_from_L10']:.11f}**, отличие от округлённого значения письма — {a['chi_absolute_difference_from_reported']:.3g}. Оно укладывается в точность представления χ пятью знаками после запятой.

Если использовать этот χ, извлечённый только из L10, то остальные разрешённые осесимметричные моменты ℓ=2…8 согласуются с эталоном с абсолютными отклонениями порядка 10⁻¹¹. L10 в этой второй проверке не является независимым подтверждением: он использован для определения χ. Это проверка взаимного согласования поздних моментов, а не утверждение о физической точности 10⁻¹¹; систематическая ошибка одной симуляции отдельно не оценена.

Масса Mf=0,95162 M в этом разборе остаётся величиной из письма. По одним поздним безразмерным моментам без площади абсолютный масштаб массы не определён. Динамические частоты позволяют поставить отдельную модельную задачу оценки массы, но здесь она не решалась.

## 6. Частота после поворота

Равновесная проверка дополнена описательными оценками I22. На каждом окне методом наименьших квадратов с равными весами оценены наклоны развёрнутой фазы и ln|Z|. Условность модели: один ведущий затухающий комплексный экспоненциальный сигнал. Никакая частота при оценке наклонов не фиксировалась.

Округлённый керровский ориентир из таблицы II статьи: ω220≈0,5535 M⁻¹, τ220≈11,707 M, α220=1/τ≈{QNM_ALPHA:.8f} M⁻¹. [Таблица II]({PAPER}).

| Окно, M | ω до поворота | ω после поворота | α | Различие ω, % | Различие α, % |
| --- | --- | --- | --- | --- | --- |
{frows}

Это показывает, почему пропуск множителя exp(−imΩt) дал бы большую ложную частотную несогласованность. Смена позднего окна для среднего на 350–500 или 450–500 M проверена; результаты и остатки регрессии доступны в JSON. Зависимость оценок от времени остаётся. Полный многомодовый анализ, оценки неопределённости и поиск дополнительных устойчивых компонент здесь не выполнены. Малое различие с одной модой не доказывает, что весь сигнал одномодов.

![Проверки данных](BH_YITIAN_A1_overview.png)

![Зависимость оценок от окна](BH_YITIAN_A1_windows.png)

## 7. Что это меняет для HF3 и HF4

Прежний вывод «масса и спин отличаются от SpEC, значит Kerr-согласие R6B сомнительно» был необоснован. Равенство физических начальных условий QC0 и данного случая не установлено. Одинаковое название режима «слияние двух чёрных дыр» не делает два расчёта одной задачей. В каждом расчёте частоты следует проверять по его собственным массе, спину, единицам и времени.

Значения QC0 Mf=0,977641203 и χ=0,648087073 известны из предыдущего диалога. Самих данных и программы R6B в этом архиве нет, поэтому сейчас их вычисление не проверено. Подставлять вместо них параметры Yitian автоматически нельзя. Моменты QuasiLocalMeasures qlm_mp_mN / qlm_mp_jN также нельзя отождествлять с Iℓm/Lℓm без формального соответствия определений.

Разумный следующий анализ этого архива: устойчивость ранних и промежуточных компонент к смене окна, числа мод и угловой системы, с отдельным контролем численного уровня. Сначала следует воспроизвести стандартные керровские смешения и затухание; остаток после конкретной подгонки сам по себе не определяет новую структуру.

Для будущего фильма уже доступны временные коэффициенты геометрических полей. Однако рисунок поля на вспомогательной сфере, радиальная деформация сферы и изометрическое вложение физического горизонта — разные операции. Полный научный фильм слияния, внутреннюю геометрию чёрной дыры или гиперформу этот архив однозначно не задаёт. При ℓ≤8 неизмеренные старшие гармоники также остаются неизвестными.

Продолжение локального QC0 с 95 до 105 M сохраняет смысл отдельного контроля поздней волны на r=40 M. Из-за нового архива оно не стало ошибочным и не заменяется данными SpEC. В данном разборе новый расчёт эволюции не запускался.

## 8. Повторение

Установить зависимости из `requirements.txt`, поместить исходный архив рядом с программой и выполнить:

```bash
python BH_YITIAN_A1_audit.py --input common_horizon_moments_data_and_scripts.zip --output results
```

Программа заново создаст JSON, отчёт и рисунки. Необходимо сохранить исходный архив: он не включён повторно в пакет результатов. Код принимает также сам `Moments_500M_L30.h5`; в таком случае контрольная сумма входа относится к HDF5, а перечень блокнотов будет отсутствовать.

Сведения о версиях среды приведены в `BH_YITIAN_A1_provenance.json`. Контрольная сумма HDF5 SHA-256: `{meta['hdf5_sha256']}`.
"""
    (out/'BH_YITIAN_A1_report.md').write_text(text,encoding="utf-8")


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input",type=Path,default=Path("common_horizon_moments_data_and_scripts.zip"))
    parser.add_argument("--output",type=Path,default=Path("BH_YITIAN_A1_output"))
    args=parser.parse_args()
    args.output.mkdir(parents=True,exist_ok=True)
    t,fields,meta=read_input(args.input)
    audit,modal,asym=analyze(t,fields)
    meta['runtime']={"python":platform.python_version(),"numpy":np.__version__,"scipy":scipy.__version__,
        "h5py":h5py.__version__,"matplotlib":matplotlib.__version__}
    meta['references']=[PAPER,ISOLATED_PAPER,KERR_PAPER]
    meta['source_code_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    dump(args.output/'BH_YITIAN_A1_summary.json',audit)
    dump(args.output/'BH_YITIAN_A1_mode_statistics.json',modal)
    dump(args.output/'BH_YITIAN_A1_provenance.json',meta)
    figures(args.output,t,fields,asym,audit)
    report(args.output,audit,meta)
    print(audit['status'])
    print(f"Данные: {len(t)} времён, 2 семейства по 81 моменту")
    print(f"Спин из позднего L10: {audit['chi_inferred_from_L10']:.11f}")
    print(f"Результаты: {args.output.resolve()}")


if __name__=="__main__":
    main()
