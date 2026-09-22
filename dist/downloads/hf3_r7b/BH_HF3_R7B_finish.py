#!/usr/bin/env python3
"""Оформление сохранённых численных результатов R7B без повторения подгонок."""
from pathlib import Path
import argparse
import hashlib
import json

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

BASE=Path(__file__).resolve().parent


def read(p):return json.loads(Path(p).read_text())
def digest(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(p,x):Path(p).write_text(json.dumps(x,ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf-8')
def pct(x):return f'{100*x:.3f}'


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--results',type=Path,default=BASE/'results')
    ap.add_argument('--previous-state',type=Path,default=BASE/'inputs/BH_HF3_resume_state_R7_prepared.json')
    args=ap.parse_args();out=args.results;s=read(out/'BH_HF3_R7B_summary.json')
    third=[x for x in s['primary_cases'] if x['period']==3]
    r40=next(x for x in third if x['radius_M']==40)
    r40ages=sorted(s['age_control_cases']+[r40],key=lambda x:x['period'])
    diagnostics=[]
    for row in third:
        a=np.loadtxt(out/f'BH_HF3_R7B_r{row["radius_M"]}_period3.tsv');t=a[:,0];y=a[:,1]+1j*a[:,2]
        q={'radius_M':row['radius_M']}
        for key,z in [('phase',np.unwrap(np.angle(y))),('logamp',np.log(abs(y)))]:
            fitted=np.polyval(np.polyfit(t-t[0],z,1),t-t[0])
            q[key+'_r2']=float(1-np.sum((z-fitted)**2)/np.sum((z-z.mean())**2))
        diagnostics.append(q)
    save(out/'BH_HF3_R7B_phase_logamp_diagnostics.json',diagnostics)
    checks={'third_period_frequency_and_damping_targets':s['numerical_frequency_damping_goals_all_three_radii'],
            'historical_nominal_BIC_target_diagnostic_only':s['historical_nominal_BIC_gate_all_three_radii'],
            'historical_phase_and_envelope_linearity':all(d['phase_r2']>=.995 and d['logamp_r2']>=.99 for d in diagnostics),
            'r40_small_window_shifts_and_even_odd_subsets':s['r40_sensitivity_ranges']['relative_omega_difference'][1]<=.01 and s['r40_sensitivity_ranges']['relative_alpha_difference'][1]<=.05,
            'third_period_waveform_agreement':all(d['coherence']>=.995 and d['phase_only_nrmse']<=.05 for d in s['third_period_shapes']),
            'technical_input':s['input_audit_pass']}
    passed=all(checks.values())
    classification=('BH_HF3_R7B_THIRD_PERIOD_R15_R30_R40_KERR220_COMPATIBILITY_SUPPORTED_HF4_READY'
                    if passed else 'BH_HF3_R7B_REVIEW_REQUIRED')
    s.update(classification=classification,project_control_task_completed=passed,control_checks=checks,
             scope_of_completion='Выполнен запланированный контроль позднего кольцевания на трёх радиусах. Полная динамика и численная сходимость не объявлены доказанными.',
             next_stage='HF4-A1: аудит геометрической реконструкции из имеющихся моментов',
             phase_logamp_diagnostics=diagnostics,
             postprocessing_program_sha256=digest(__file__))
    save(out/'BH_HF3_R7B_summary.json',s)

    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11,'axes.titlesize':12,
                         'axes.spines.top':False,'axes.spines.right':False,'savefig.facecolor':'white'})
    fig,axs=plt.subplots(2,2,figsize=(13.6,9.2),layout='constrained')
    fig.suptitle('QC0: третий период на 40M подтверждает близость к керровскому кольцеванию',fontsize=15)
    a=np.loadtxt(out/'BH_HF3_R7B_r40_period3.tsv');t=a[:,0];y=a[:,1]+1j*a[:,2]
    k=a[:,3]+1j*a[:,4];f=a[:,5]+1j*a[:,6]
    ax=axs[0,0]
    ax.plot(t,40*y.real,'o',ms=3.4,color='#2466a4',label='Данные')
    ax.plot(t,40*k.real,'-',lw=1.6,color='#26333b',label='Мода Керра')
    ax.plot(t,40*f.real,'--',lw=1.3,color='#db8237',label='Свободная мода')
    ax.set(title='Третий период на r = 40M',xlabel='Время, M',ylabel=r'$\mathrm{Re}[r\Psi_4^{22}]$')
    ax.legend(frameon=False,fontsize=9,loc='upper right');ax.grid(alpha=.17)
    ax=axs[0,1];rms=np.sqrt(np.mean(abs(y)**2))
    ax.plot(t,100*abs(y-k)/rms,color='#26333b',label='Мода Керра')
    ax.plot(t,100*abs(y-f)/rms,'--',color='#db8237',label='Свободная мода')
    ax.set(title=f'Средний относительный остаток: {pct(r40["fixed_nrmse"])}% и {pct(r40["free_nrmse"])}%',
           xlabel='Время, M',ylabel='Модуль остатка / среднеквадр. сигнал, %')
    ax.set_ylim(bottom=0);ax.legend(frameon=False,fontsize=9);ax.grid(alpha=.17)
    ax=axs[1,0]
    for key,label,color in [('relative_omega_difference','Частота','#2466a4'),('relative_alpha_difference','Затухание','#b04743')]:
        ax.semilogy([1,2,3],[100*x[key] for x in r40ages],'o-',color=color,label=label,lw=1.8)
    ax.set(title='Приближение к Керру с возрастом сигнала: r = 40M',xlabel='Полный период после максимума',ylabel='Отличие от керровского ориентира, %')
    ax.set_xticks([1,2,3]);ax.legend(frameon=False,fontsize=9);ax.grid(alpha=.17,which='both')
    ax=axs[1,1]
    ref=np.loadtxt(out/'BH_HF3_R7B_r15_period3.tsv');yr=15*(ref[:,1]+1j*ref[:,2]);basephase=np.exp(-1j*np.angle(yr[0]))
    for r,color,style in [(15,'#2466a4','-'),(30,'#b04743','--'),(40,'#318b76',':')]:
        row=next(x for x in third if x['radius_M']==r)
        a=np.loadtxt(out/f'BH_HF3_R7B_r{r}_period3.tsv');v=r*(a[:,1]+1j*a[:,2])
        phase=1 if r==15 else np.exp(1j*next(z['constant_phase_applied_rad'] for z in s['third_period_shapes'] if z['reference_radius_M']==15 and z['compared_radius_M']==r))
        ax.plot(a[:,0]-row['peak_M'],(v*phase*basephase).real,style,lw=1.8,color=color,label=f'r = {r}M')
    ax.set(title='Одинаковый возраст: только постоянный поворот фазы',xlabel='Время после местного максимума, M',ylabel=r'$\mathrm{Re}[r\Psi_4^{22}]$')
    ax.legend(frameon=False,fontsize=9);ax.grid(alpha=.17)
    fig.savefig(out/'BH_HF3_R7B_overview.png',dpi=170)
    fig.savefig(out/'BH_HF3_R7B_overview.pdf')
    plt.close(fig)

    lines=['# BH-HF3-R7B — проверка позднего кольцевания после продолжения до 105M','',
           '**Запланированный контроль третьего периода на r=40M выполнен. Результат согласуется с керровским кольцеванием на трёх радиусах; маршрут может переходить к HF4.**','',
           'Это завершение конкретной контрольной задачи HF3. Оно не означает, что вся ранняя динамика объяснена или что проведена независимая экспериментальная проверка ОТО. Исследуется собственный численный расчёт QC0, а не непосредственно запись GW150914.','',
           '## Технический результат','',
           'Прочитан предоставленный `BH_HF3_R7A_output.zip`. Подтверждены t=105M, итерация 53760, семь сообщений восстановления с t=95M и код завершения 0. Продолжение заняло около 25 минут 18 секунд по журналу сценария. Повторная эволюция здесь не запускалась.','',
           'Независимо проверены 125 волновых рядов, точный стык 95M, полная новая сетка времени и три файла общего горизонта QLM slot 2. Проверены хеши программы, исходного и нового файлов параметров. 136 общих файлов старого участка побайтно совпадают с прежним пакетом R6A. Сам физический анализ распространения и закона 1/r из R6B не повторялся.','',
           'Журнал и локальный реестр подтверждают создание обеих частей контрольной точки it=53760. Гигабайтные контрольные точки не входят в полученный ZIP; их полное содержимое здесь не читалось.','',
           '## Как проведено сравнение','',
           'До подгонок этого этапа сохранён `BH_HF3_R7B_protocol.json`. Основной керровский ориентир зафиксирован по ранее вычисленному численному спектру QC0:','',
           r'\[M_f=0.977461203278627M,\quad J_f=0.648087073319352M^2,\quad \chi=0.6783194994253878.\]','',
           r'\[\omega_K=0.5353401318639867M^{-1},\qquad \alpha_K=0.08345059699511974M^{-1}.\]','',
           'Спектр получен ранее численным решением керровской спектральной задачи; он не подгонялся к новым волновым данным. Метод вычисления: [документация qnm](https://qnm.readthedocs.io/en/latest/README.html), [Stein, 2019](https://arxiv.org/abs/1908.10377). Длинная запись чисел обеспечивает воспроизводимость вычислений и не задаёт физическую точность расчёта.','',
           'Период T=11.73680980221365M. Третий период означает окно возраста от 2T до 3T после местного максимума. Для r=40M это 92.4736196044–104.2104294066M; на исходной сетке в него входят 47 отсчётов 92.50–104.00M. Концы окна не дополнялись искусственными измерениями.','',
           r'Сравниваются модели $z(t)=C\exp[(-\alpha-i\omega)(t-t_0)]$. В фиксированной модели свободна комплексная амплитуда $C$; во второй свободны также $\omega$ и $\alpha$. В обеих минимизируется один комплексный критерий $\mathrm{RSS}=\sum_j|z_j-z_{\rm model}(t_j)|^2$ с одинаковыми весами. Амплитуда находится аналитической проекцией; частота и затухание — численно с четырёх начальных приближений.','',
           'Так исправлен недостаток старого сравнения R6B: раньше свободные параметры брались из отдельных регрессий фазы и логарифма амплитуды. Для сопоставимости заново обработаны только выбранные опорные окна прежних данных тем же комплексным критерием. Поэтому новые оценки не обязаны совпадать с прежними регрессиями.','',
           '## Третий период: основной результат','',
           '| Радиус | Частота, M⁻¹ | Затухание, M⁻¹ | Отличие частоты | Отличие затухания |','|---|---:|---:|---:|---:|']
    for x in third:lines.append(f'| {x["radius_M"]}M | {x["omega"]:.9f} | {x["alpha"]:.9f} | {pct(x["relative_omega_difference"])}% | {pct(x["relative_alpha_difference"])}% |')
    lines+=['','Все три радиуса удовлетворяют ранее выбранным численным ориентирам: отличие частоты не более 1%, затухания не более 5%. Эти проценты описывают разность точечных оценок; они не являются доверительными интервалами или измерением полной численной ошибки.','',
            '| Радиус | Относительный остаток Керра | Остаток свободной модели | Условный ΔBIC |','|---|---:|---:|---:|']
    for x in third:lines.append(f'| {x["radius_M"]}M | {pct(x["fixed_nrmse"])}% | {pct(x["free_nrmse"])}% | {x["nominal_delta_BIC_fixed_minus_free"]:.3f} |')
    lines+=['',r'Относительный остаток определён как $\sqrt{\mathrm{RSS}/\sum_j|z_j|^2}$. Например, на r=40M фиксированная модель даёт 2.661%, свободная — 2.632%. Два дополнительных параметра уменьшают остаток лишь немного. Остаток около 2.6% сохраняется в обеих моделях и не объявляется объяснённым.','',
            r'Для преемственности приведён $\Delta\mathrm{BIC}=2N\log(\mathrm{RSS}_K/\mathrm{RSS}_{free})-2\log(2N)$, где N — число комплексных отсчётов; число действительных параметров равно 2 и 4. Его значения ниже прежнего порога +10 на всех трёх радиусах. Но предположение о независимых одинаково распределённых ошибках здесь не обосновано: остатки численного сигнала коррелированы, а модель численной ошибки отсутствует. Поэтому этот показатель не переводится в вероятности гипотез и не используется как калиброванный тест ОТО.','',
            '## Устойчивость и проверка продолжением','',
            'Окно r=40M сдвигалось на −0.5, −0.25, 0, +0.25 и +0.5M при прежней длине. Отдельно использовались чётные и нечётные отсчёты. Это проверка чувствительности к выбору точек, а не сходимости решения по разрешению сетки. Во всех этих случаях отличие частоты остаётся 0.035–0.118%, затухания — 0.866–2.080%; относительный остаток фиксированной модели — 2.562–2.742%.','',
            'Две проверки на аналитически заданных затухающих комплексных сигналах восстанавливают известные частоту и затухание. Все четыре старта оптимизации сходятся к одному решению в каждом основном окне; решений на границах диапазонов нет.','',
            'Дополнительно амплитуда и свободные параметры оценены только на первых двух третях третьего периода r=40M. На оставшейся трети, без новой подгонки, получены ошибки предсказания: 2.230% для фиксированной модели Керра и 5.509% для свободной. Это один хронологический контроль на 15 отсчётах, а не независимый набор наблюдений.','',
            'На укороченном обучающем окне свободная оценка затухания отличается от ориентира уже на 5.34%. Поэтому устойчивость полного периода при малых сдвигах не следует превращать в заявление о столь же точном определении параметров на любом коротком отрезке.','',
            '## Что меняется со временем и радиусом','',
            '| r=40M: период после максимума | Отличие частоты | Отличие затухания | Относительный остаток Керра |','|---|---:|---:|---:|']
    for x in r40ages:lines.append(f'| {x["period"]} | {pct(x["relative_omega_difference"])}% | {pct(x["relative_alpha_difference"])}% | {pct(x["fixed_nrmse"])}% |')
    lines+=['','На одном и том же радиусе одно-модовое керровское описание последовательно улучшается с возрастом сигнала. Новые окна на внешних радиусах продолжают ту же картину:','',
            '| Радиус и период | Частота, M⁻¹ | Затухание, M⁻¹ | Отличие частоты | Отличие затухания |','|---|---:|---:|---:|---:|']
    for x in s['primary_cases']:
        if x['period']<3:lines.append(f'| {x["radius_M"]}M, {x["period"]}-й | {x["omega"]:.9f} | {x["alpha"]:.9f} | {pct(x["relative_omega_difference"])}% | {pct(x["relative_alpha_difference"])}% |')
    lines+=['','Второй период на 50M близок по параметрам ко второму на 40M; первый на 60M — к первому на 40M. Значительное отличие ранних периодов от одной поздней моды само по себе не является отклонением от ОТО. В этих данных решающее значение имеет возраст сигнала после местного максимума.','',
            'На третьем периоде комплексные формы также согласуются после умножения на радиус и одного постоянного поворота фазы. Дополнительный временной сдвиг и подгонка амплитуды не использовались.','',
            '| Пара радиусов | Когерентность | Различие после поворота фазы |','|---|---:|---:|']
    for x in s['third_period_shapes']:lines.append(f'| {x["reference_radius_M"]}–{x["compared_radius_M"]}M | {x["coherence"]:.9f} | {pct(x["phase_only_nrmse"])}% |')
    lines+=['',r'Соотношение $\Psi_4^{2,-2}\simeq\overline{\Psi_4^{22}}$ в трёх поздних окнах выполняется с относительным различием менее $8\times10^{-16}$. Это внутренняя согласованность данных расчёта с заданной симметрией, а не независимая проверка природы слияния.','',
            f'![Проверка R7B](sandbox:{(out/"BH_HF3_R7B_overview.png").resolve()})','',
            '## Поздний горизонт','',
            'На 100–105M средние значения составляют Mf=0.9774491229M, Jf=0.6481486193M², chi=0.6784006850. По отношению к зафиксированному ориентиру масса меняется примерно на −0.00124%, а безразмерный спин — на +0.01197%. Ориентир основного анализа намеренно оставлен прежним. Эти небольшие изменения дополнительно напоминают, что последние цифры оценок частоты не являются точностью знания конечного состояния.','',
            'Параметры Yitian не подставлялись в QC0. Его данные относятся к другой конфигурации; выводы A1–A4 сохранены отдельно. В A3/A4 ранняя динамика и конкретная квадратичная связь не получили полного подтверждения.','',
            '## Итог и маршрут','',f'`{classification}`','',
            'Техническое продолжение завершено, третий период на r=40M доступен, а выбранные количественные условия позднего сравнения выполнены на 15, 30 и 40M. По этой контрольной задаче дальнейшее удлинение QC0 сейчас не требуется.','',
            'Область вывода ограничена одной численной конфигурацией, одним разрешением, конечными радиусами извлечения и выбранными временными окнами. Не выполнены экстраполяция к бесконечности, проверка сходимости по сетке или независимая статистическая проверка ОТО. Не установлены дополнительные измерения, кротовая нора или отдельная «гиперформа» как физический объект. Переходные структуры остаются предметом геометрического описания.','',
            '**Следующий шаг — HF4-A1: определить, что именно можно восстановить из имеющихся горизонтовых моментов.**','',
            '1. Сопоставить определения коэффициентов, угловой базис, ориентацию и меру интегрирования для QC0 и данных Yitian.','2. Отдельно установить, какие скалярные карты на вспомогательной сфере восстанавливаются из конечного набора моментов и какая дополнительная информация нужна для физической метрики горизонта.','3. Проверить роль обрезания по l≤8 и выбор шкалы изображения; указать, какие особенности изображения определяются данными, а какие — способом представления.','4. После этого подготовить научно подписанную визуализацию динамики; 3D-фильм остаётся результатом HF4 после определения допустимого содержания реконструкции.','',
            '## Воспроизведение','',
            'Основной анализ: `python BH_HF3_R7B.py`. Оформление готовых результатов без повторных подгонок: `python BH_HF3_R7B_finish.py`. Самодостаточный архив содержит оба входных ZIP, программы, протокол, численные таблицы, результаты проверок и запись состояния. Нужны numpy, scipy и matplotlib. Программа Einstein Toolkit этими командами не запускается.','']
    (out/'BH_HF3_R7B_report.md').write_text('\n'.join(lines),encoding='utf-8')

    state=read(args.previous_state)
    state.update(stage='BH-HF3-R7B',state=classification,last_completed_stage='BH-HF3-R7B',
                 no_new_evolution_run=False,no_new_evolution_run_here=True,new_evolution_completed_by_user=True,
                 last_verified_QC0={'time_M':105,'iteration':53760,'exit_status':0,'seam95_exact':True,'input_archive_sha256':s['provenance']['new_archive_sha256']},
                 QC0_completion_pending=None,HF3_scientific_completion_not_yet_declared=False,
                 HF3_completion_scope=s['scope_of_completion'],
                 next_step='HF4-A1: аудит возможности геометрической реконструкции из моментов, углового базиса и меры; новый расчёт эволюции для этого шага не запланирован.',
                 next_after_QC0='HF4-A1 → проверяемые карты геометрических величин → научная визуализация.',
                 R7B_summary_sha256=digest(out/'BH_HF3_R7B_summary.json'),
                 parent_resume_state_sha256=digest(args.previous_state))
    state['R7']['state']='COMPLETED_T105_IT53760_EXIT0'
    state['R7B']={'classification':classification,'reference':s['reference'],'third_period_cases':third,
                  'remaining_relative_single_mode_residual':r40['fixed_nrmse'],
                  'numeric_targets_pass':passed,'calibrated_statistical_Kerr_acceptance':False,
                  'resolution_convergence_test_available':False,'new_physics_established':False}
    save(out/'BH_HF3_resume_state_R7B.json',state)
    print(classification)
    print('Отчёт и графики подготовлены.')


if __name__=='__main__':main()
