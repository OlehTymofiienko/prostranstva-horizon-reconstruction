#!/usr/bin/env python3
"""Summarize saved B5 fits without refitting or rerunning previous stages."""
from pathlib import Path
import csv, json, hashlib, shutil, zipfile
from datetime import datetime, timezone
import numpy as np
from BH_HF4_B5 import FIELDS, LEVELS, sha, write_csv, write_json
ROOT=Path(__file__).resolve().parent
OUT=ROOT/'results'
source=ROOT/'input/BH_HF4_B4_fields_geometry_masks.npz'
z=np.load(source); fits=np.load(OUT/'BH_HF4_B5_fits_residuals.npz')
s=json.loads((OUT/'BH_HF4_B5_summary.json').read_text())
checks=json.loads((OUT/'BH_HF4_B5_checks.json').read_text())
rows=[]
for it in fits['iterations']:
    w=z[f'physical_weights_{it}']
    for f in FIELDS:
        r=fits[f'residual_{f}_L8_{it}']; norm=w*r*r; order=np.argsort(norm)[::-1]
        row={'iteration':int(it),'t_M':float(it/512),'field':f,'L':8,
             'top2_nodes_squared_norm_fraction':float(norm[order[:2]].sum()/norm.sum()),
             'top2_nodes_physical_area_fraction':float(w[order[:2]].sum()),
             'top10_nodes_squared_norm_fraction':float(norm[order[:10]].sum()/norm.sum())}
        for j,idx in enumerate(order[:2],1):
            row.update({f'node{j}_index':int(idx),f'node{j}_theta_deg':float((idx//72+.5)*180/35),
                        f'node{j}_phi_deg':float(idx%72*5),f'node{j}_residual':float(r[idx])})
        rows.append(row)
write_csv(OUT/'BH_HF4_B5_residual_node_concentration.csv',rows)
first=next(r for r in rows if r['iteration']==9536 and r['field']=='repsi2')
names={'repsi2':'Re Ψ₂','impsi2':'Im Ψ₂','abs_npsigma':'|σ|'}
table='\n'.join('| '+names[f]+' | '+' | '.join(f"{100*s['fields'][f][str(L)]['R2']['min']:.3f}–{100*s['fields'][f][str(L)]['R2']['max']:.3f}%" for L in LEVELS)+' |' for f in FIELDS)
last='\n'.join(f"| {names[f]} | {100*s['fields'][f]['8']['R2']['last']:.5f}% | {100*s['fields'][f]['8']['residual_RMS_over_original_SD']['last']:.3f}% |" for f in FIELDS)
overlap='\n'.join(f"| {names[f]} | {s['fields'][f]['8']['dice']['min']:.4f}–{s['fields'][f]['8']['dice']['max']:.4f} | {100*s['fields'][f]['8']['residual_squared_norm_fraction_in_original_mask']['min']:.3f}–{100*s['fields'][f]['8']['residual_squared_norm_fraction_in_original_mask']['max']:.3f}% |" for f in FIELDS)
report=f'''# HF4-B5 — угловое описание полей и локализация остатка

22 сентября 2026 года. **Этап B5 завершён по сохранённому протоколу.**

Основная угловая картина трёх полей QC0 хорошо описывается сферическими гармониками в исходных координатах. При этом остаток неоднороден: в самом раннем кадре Re Ψ₂ существенную часть невязки дают два узла. Это описательный результат, а не обнаружение дополнительного физического объекта и не доказательство отсутствия какой-либо дополнительной структуры.

## Что было сделано

Использован неизменённый результат B4: 55 кадров общего **кажущегося** горизонта QC0 от 18,625M до 22M, шаг 0,0625M, по 2520 угловых узлов (35×72). Физические веса площади, поля и исходные маски прочитаны из готового файла. Эволюция Einstein Toolkit, площадь, маски B4, кольцевание HF3 и ветви HF1/Yitian не пересчитывались.

SHA-256 входа: `{s['input_sha256']}`.

SHA-256 ранее сохранённого протокола: `{checks['protocol_sha256']}`.

Для каждого времени и каждого поля независимо найдена взвешенная аппроксимация в пространствах 0≤ℓ≤L, L=2,4,6,8. Они содержат соответственно 9,25,49,81 вещественный коэффициент. Всего 660 подгонок полей (220 общих матриц для трёх правых частей).

Вещественный базис: Yℓ0; √2 Re Yℓm и √2 Im Yℓm при m>0, стандартная нормировка на единичной сфере и фаза Кондона–Шортли. Используется `scipy.special.sph_harm_y(ell,m,theta,phi)`. θ — полярный угол, φ — азимут. Исходные оси и угловые индексы QC0 сохранены.

Минимизируется Σᵢwᵢ(Fᵢ−Fᴸᵢ)². Основная характеристика:

R² = 1 − Σᵢwᵢ(Fᵢ−Fᴸᵢ)² / Σᵢwᵢ(Fᵢ−⟨F⟩w)².

Это доля **взвешенной угловой вариации**, описанная моделью на этих же узлах. Это не доля физической энергии, не вероятность верности модели, не оценка точности симуляции и не проверка предсказания на независимых данных. Более широкое вложенное пространство по построению не ухудшает такую подгонку. Вырожденных знаменателей нет.

Сферические гармоники относятся к координатной сфере. Они не являются найденными собственными функциями внутреннего лапласиана деформированного горизонта. Коэффициенты B5 не отождествляются с геометрическими мультиполями Yitian. Для |σ| раскладывается только вещественный модуль, без восстановления комплексной фазы и без спин-взвешенного разложения.

## Доля описанной вариации по всем 55 кадрам

| Поле | L≤2 | L≤4 | L≤6 | L≤8 |
|---|---:|---:|---:|---:|
{table}

Диапазоны охватывают все времена; это не доверительные интервалы. До L=8 объясняется не менее 96,388% вариации Re Ψ₂, 99,074% Im Ψ₂ и 99,815% |σ|. Уже L=2 описывает большую часть крупной картины, но качество отличается по полям и времени.

На последнем кадре t=22M:

| Поле | R² при L≤8 | RMS остатка / стандартное отклонение исходного поля |
|---|---:|---:|
{last}

Относительная невязка |σ| к концу интервала немного растёт, хотя абсолютный RMS остатка снижается с {s['fields']['abs_npsigma']['8']['residual_RMS']['first']:.8f} до {s['fields']['abs_npsigma']['8']['residual_RMS']['last']:.8f}. Поэтому нельзя заменять результат утверждением о монотонном улучшении всех полей.

## Где остаётся невязка

Исходная маска B4 выделяет примерно 5% физической площади с наибольшим отклонением поля от взвешенной медианы. Она не изменена. Новая маска остатка использует такое же правило для r=F−Fᴸ: 95-й взвешенный квантиль |r−medianw(r)|, включая все равные порогу узлы. Точное правило зафиксировано в `BH_HF4_B5_definition.json` перед подгонкой.

Измерены взвешенное перекрытие Dice и доля Σwr² внутри прежней маски. Здесь Σwr² — квадрат нормы численного остатка, **не физическая энергия**.

| Поле, L≤8 | Dice старой и остаточной масок | Доля квадрата нормы остатка внутри старой маски |
|---|---:|---:|
{overlap}

У Im Ψ₂ и |σ| максимум остатка обычно располагается вне прежней выделенной области. У Re Ψ₂ картина меняется со временем; нельзя объявлять все остатки одной устойчивой локализованной сущностью. Нулевое перекрытие масок не означает нулевой остаток в старой области: сравниваются только верхние хвосты двух распределений.

**Особенность первого кадра.** При t=18,625M для Re Ψ₂ R²=96,38844%, а RMS остатка составляет 19,004% исходного стандартного отклонения. В старой маске находится 83,476% квадрата нормы остатка. Дополнительная проверка уже сохранённых остатков показывает: два узла при θ=90°, φ=40° и 220° дают **{100*first['top2_nodes_squared_norm_fraction']:.3f}%** всей Σwr², занимая **{100*first['top2_nodes_physical_area_fraction']:.4f}%** дискретной физической меры поверхности.

Это важное ограничение интерпретации первого кадра. Высокая общая доля описанной вариации не устраняет отдельные большие локальные отклонения. Такая зависимость от двух узлов требует проверки исходной численной диагностики и разрешения прежде, чем ей можно приписать самостоятельный физический смысл. В B5 их не удаляли, не сглаживали и не подгоняли повторно после исключения. Эта дополнительная характеристика локализации приведена описательно, без статистического порога обнаружения.

На рисунках `BH_HF4_B5_maps_it9536.png` и `BH_HF4_B5_maps_it11264.png` показаны исходное поле минус среднее, аппроксимация минус то же среднее и остаток. Во всех трёх столбцах строки единая цветовая шкала. Оранжевый контур — исходная маска B4. Это карты углов, а не новое восстановление формы. Все остатки и маски для остальных времён сохранены в NPZ.

## Численная проверка

- Все веса положительны; максимальное отклонение их суммы от единицы {checks['weight_sum_max_error']:.3g}.
- Все 220 матриц имеют полный столбцовый ранг. Числа обусловленности: {checks['condition_min']:.6f}–{checks['condition_max']:.6f}.
- Аналитические Y00, Y10 и теорема сложения проверены; максимальная ошибка {checks['addition_theorem_max_error']:.3g}.
- Максимальная относительная невязка нормальных уравнений {checks['max_relative_normal_equation_error']:.3g}; взвешенное среднее остатка согласуется с нулём.
- Независимые SVD- и QR-решения для крайних времён, L=8, отличаются в предсказаниях не более {checks['svd_qr_prediction_max_relative_difference']:.3g} относительно масштаба исходных значений.
- Вложенные подгонки не ухудшают квадрат невязки. Контрольная сумма входа после анализа не изменилась.

Эти проверки относятся к корректности линейной обработки, а не к сходимости симуляции по разрешению.

## Научный итог и завершённость

**B5 закрывает запланированный вопрос:** количественно установлено, сколько вариации описывают выбранные угловые пространства и как остаток соотносится с прежними областями. Основная картина допускает компактное угловое описание. Ненулевая невязка остаётся, особенно в первых кадрах Re Ψ₂, но сама по себе не устанавливает новую физику или отдельный «латентный объект».

Теперь можно завершить выпуск 1.0 как проект **реконструкции и описательной проверки доступных данных**, включающий HF1, QC0/HF3/HF4 и анализ мультиполей Yitian. Это три отдельных численных конфигурации: SXS:BBH:0305, QC0 и SXS:BBH:0389. Их нельзя склеивать в одну физическую историю.

Открытые границы выпуска:

1. Нет единого проверенного ряда геометрии QC0 от 0 до 105M: поздняя эволюция/кольцевание и ранняя геометрия имеют разное покрытие.
2. Нет проверки этих локальных остатков на втором разрешении, другой координатной системе и независимой тетраде. Источник узловых особенностей первого кадра не установлен.
3. HF1 показывает области по своему критерию градиента, QC0 — по отклонению полей. Их физическое тождество не доказано.
4. Конечный набор мультиполей Yitian не задаёт единственную полную 3D-геометрию. Ранняя динамика и уникальная идентификация нелинейных компонент остаются незамкнутыми выводами предыдущих этапов, а не молча объявляются доказанными.
5. Визуализируются кажущиеся горизонты или явно обозначенные скалярные карты. Это не найденный горизонт событий и не оптический фильм с трассировкой света.

Эти границы входят в результаты проекта. Они не требуют скрывать завершённую работу и не превращаются в утверждение о полностью решённой физической задаче.

## Воспроизведение и состав

`BH_HF4_B5.py` — только новый анализ B5; `finish_B5.py` — сводка и упаковка сохранённых результатов. В архив включён входной NPZ B4, его геометрическая таблица и исходный протокол. `README.md` содержит команды. При обычном продолжении повторно запускать программы не нужно.

Основные таблицы: `BH_HF4_B5_metrics.csv`, `BH_HF4_B5_coefficients.csv`, `BH_HF4_B5_design_checks.csv`, `BH_HF4_B5_singular_values.csv`, `BH_HF4_B5_residual_node_concentration.csv`. Массивы аппроксимаций, остатков и их масок — `BH_HF4_B5_fits_residuals.npz`. Времена и единицы сохранены из QC0.

Математическое соглашение реализации: [официальная документация SciPy sph_harm_y](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.sph_harm_y.html). Нормировка дополнительно проверена аналитически в этом запуске; численные выводы отчёта получены из указанных файлов проекта.
'''
(OUT/'BH_HF4_B5_report.md').write_text(report)
state={'stage':'BH-HF4-B5','status':'COMPLETE_DESCRIPTIVE_COORDINATE_HARMONIC_DECOMPOSITION',
       'saved_UTC':datetime.now(timezone.utc).isoformat(),'completed':['HF4-B4 saved geometry reused','HF4-B5 all 660 field fits','residual localization and conditioning checks'],
       'input_sha256':s['input_sha256'],'previous_results_recomputed':False,
       'resume_action':'Package and publish reconstruction research release 1.0; do not rerun B4 or B5.',
       'scientific_scope':'Available-data reconstruction and descriptive analysis, not complete all-time geometry or evidence of new physics',
       'publication_status':'IN_PREPARATION_NOT_YET_PUBLISHED',
       'unresolved':['QC0 continuous geometry 0..105M absent','resolution/gauge/tetrad robustness of localized residuals untested','initial RePsi2 residual dominated by two nodes; origin not determined','Yitian early nonlinear closure not established'],
       'report':'results/BH_HF4_B5_report.md','input':'input/BH_HF4_B4_fields_geometry_masks.npz',
       'results':'results/BH_HF4_B5_fits_residuals.npz'}
write_json(OUT/'BH_HF4_resume_state_B5.json',state)
manifest=[{'path':str(p.relative_to(ROOT)),'bytes':p.stat().st_size,'sha256':sha(p)} for p in sorted(ROOT.rglob('*')) if p.is_file() and '__pycache__' not in p.parts and p.name not in ('BH_HF4_B5_manifest.json','BH_HF4_B5_results.zip','run.log')]
write_json(ROOT/'BH_HF4_B5_manifest.json',manifest)
with zipfile.ZipFile(OUT/'BH_HF4_B5_results.zip','w',zipfile.ZIP_DEFLATED,compresslevel=6) as archive:
    for record in manifest:archive.write(ROOT/record['path'],'hf4_b5/'+record['path'])
    archive.write(ROOT/'BH_HF4_B5_manifest.json','hf4_b5/BH_HF4_B5_manifest.json')
print(json.dumps({'report':str(OUT/'BH_HF4_B5_report.md'),'archive_bytes':(OUT/'BH_HF4_B5_results.zip').stat().st_size,'first_frame_top2':first},ensure_ascii=False))
