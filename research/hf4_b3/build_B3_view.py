#!/usr/bin/env python3
"""Build B3 presentation from saved results, without repeating the analysis."""
from pathlib import Path
import csv,json,re
import numpy as np

HERE=Path(__file__).resolve().parent
OUT=HERE/'results'
B2=HERE.parent/'hf4_b2'


def build():
    frozen=np.load(OUT/'BH_HF4_B3_intrinsic_geometry_and_masks.npz')
    data=json.loads((B2/'results/BH_HF4_B2_view_data.json').read_text())
    for frame in data['frames']:
        if 'fields' not in frame:continue
        it=frame['iteration']
        for f in data['fields']:
            frame['masks'][f]=np.flatnonzero(frozen[f'mask_{f}_{it}']).tolist()
    text=json.dumps(data,separators=(',',':'),allow_nan=False)
    (OUT/'BH_HF4_B3_view_data.json').write_text(text)
    template=(B2/'viewer_template.html').read_text().replace('HF4-B2','HF4-B3').replace('HF4B2','HF4B3')
    template=template.replace('примерно 5% координатной площади сетки с наибольшим отклонением поля от медианы данного кадра. Это новое выделение B2 в QC0; маска HF1 из SXS:BBH:0305 сюда не переносилась.',
       'примерно 5% восстановленной физической площади с наибольшим отклонением поля от взвешенной медианы данного кадра. Использована внутренняя метрика из тетрады QC0. Это маска B3; тождество структуре HF1 не установлено.')
    template=template.replace('Для замыкания изображения добавлены две маленькие нейтральные крышки, не участвующие в анализе.',
       'Для замыкания изображения добавлены две маленькие нейтральные крышки. Они не участвуют в анализе; физическая мера интегрируется на исходной угловой сетке.')
    template=template.replace('<details><summary>Как понимать эту сцену</summary>',
       '<details><summary>Как понимать эту сцену</summary><p>Форма и координаты поверхности сохранены из расчёта. Физическая площадь используется для выделения областей. Метод прошёл контроль на пяти поздних состояниях; отдельная сверка ранних площадей с выводом того же запуска ещё предстоит.</p>')
    vendor=B2/'vendor/plotly-3.3.1.min.js'
    html=template.replace('__DATA__',text.replace('</','<\\/')).replace('__PLOTLY__',vendor.read_text())
    assert '__DATA__' not in html and '__PLOTLY__' not in html
    assert 'примерно 5% координатной площади' not in html
    (OUT/'BH_HF4_B3_fields_3D.html').write_text(html)

    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10})
    def rows(name):
        with (OUT/name).open() as f:return list(csv.DictReader(f))
    geo=rows('BH_HF4_B3_geometry.csv');loc=rows('BH_HF4_B3_localization.csv');comp=rows('BH_HF4_B3_mask_comparison.csv')
    late=json.loads((OUT/'BH_HF4_B3_late_area_controls.json').read_text())
    fig,axs=plt.subplots(2,2,figsize=(12,8),facecolor='white')
    colors={'repsi2':'#ad4c2b','impsi2':'#246999','abs_npsigma':'#7b4ca3'}
    labels={'repsi2':'Re Ψ₂ · продольная область','impsi2':'Im Ψ₂ · полярная область','abs_npsigma':'|σ| · полярная область'}
    def blocks(ax,rs,key,**kwargs):
        for ids in [[9536,9696],[10368,10528]]:
            r=[x for x in rs if ids[0]<=int(x['iteration'])<=ids[1]]
            ax.plot([float(x['tau_M']) for x in r],[float(x[key]) for x in r],marker='o',ms=3,**kwargs)
            kwargs.pop('label',None)
    a=axs[0,0]
    blocks(a,geo,'area_physical_fejer',color='#246999',label='Восстановленная физическая площадь')
    blocks(a,geo,'area_coordinate_fd4_fejer',color='#7d8894',label='Евклидова координатная площадь')
    a.set(title='Разные меры площади одной поверхности',ylabel='Площадь / M²',ylim=(0,44));a.legend(fontsize=8,loc='center right')
    a=axs[0,1]
    for field in colors:
        r=[x for x in loc if x['field']==field and x['measure']=='fejer' and float(x['quantile'])==.95]
        blocks(a,r,'long_axis_fraction' if field=='repsi2' else 'polar_fraction',color=colors[field],label=labels[field])
    a.set(title='Где расположены области B3',ylabel='Доля выделенной физической площади',ylim=(.65,1.035));a.legend(fontsize=8,loc='lower right')
    a=axs[1,0]
    for field in colors:
        blocks(a,[x for x in comp if x['field']==field],'B2_B3_Dice_physical_measure',color=colors[field],label=labels[field].split(' · ')[0])
    a.set(title='Совпадение масок B2 и B3',ylabel='Взвешенное перекрытие Dice',ylim=(.88,1.01));a.legend(fontsize=8)
    a=axs[1,1];a.semilogy([x['time_M'] for x in late],[abs(x['relative_difference_fejer']) for x in late],'o-',color='#246999')
    a.set(title='Контроль по сохранённой площади QLM',xlabel='Время контрольной точки / M',ylabel='Относительное расхождение')
    for a in [axs[0,0],axs[0,1],axs[1,0]]:
        a.axvspan(.3125,1.625,color='#d5dce4',alpha=.5,zorder=-2)
        a.set_xlabel('Возраст общего горизонта τ / M')
        a.text(.49,.43 if a is axs[0,1] else .05,'Нет кадров в пакете',transform=a.transAxes,ha='center',fontsize=8,color='#556171')
    for a in axs.ravel():a.grid(alpha=.18)
    fig.suptitle('HF4-B3 · общая локализация сохраняется при физическом взвешивании',fontsize=15,y=.99)
    fig.text(.07,.025,'Проверка восстановления метрики не является проверкой сходимости симуляции. Прямая сверка ранней площади ещё нужна.',fontsize=9,color='#536171')
    fig.tight_layout(rect=[.015,.055,1,.95]);fig.savefig(OUT/'BH_HF4_B3_overview.png',dpi=170);plt.close(fig)
    verification={'display_mask_indices_match_saved_B3':True,'surface_coordinates_and_field_values_reused_from_B2':True,
                  'new_evolution_or_missing_frame_interpolation':False,'browser_WebGL_render_tested':False}
    (OUT/'BH_HF4_B3_view_checks.json').write_text(json.dumps(verification,indent=2)+'\n')
    print(json.dumps({'viewer_bytes':(OUT/'BH_HF4_B3_fields_3D.html').stat().st_size,'saved_figure':True}))


if __name__=='__main__':build()
