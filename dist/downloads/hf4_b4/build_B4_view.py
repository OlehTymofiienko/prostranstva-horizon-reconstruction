#!/usr/bin/env python3
"""Present already saved B4 results; no numerical area/mask reconstruction."""
from pathlib import Path
import json,re
import numpy as np
from BH_HF4_B4 import shape_gp,csvrows

HERE=Path(__file__).resolve().parent;OUT=HERE/'results';BASE=HERE.parent
FIELDS=['repsi2','impsi2','abs_npsigma']


def build():
    data=np.load(OUT/'BH_HF4_B4_fields_geometry_masks.npz')
    previous=json.loads((BASE/'hf4_b2/results/BH_HF4_B2_view_data.json').read_text())
    geo={int(x['iteration']):x for x in csvrows(OUT/'BH_HF4_B4_geometry.csv')}
    frames=[];gp_tris=None;limit=0;range_abs={f:0 for f in FIELDS};mask_checks=0
    for it in range(0,11265,32):
        frame={'iteration':it,'t_M':it/512,'tau_M':(it-9536)/512,'objects':[],'missing_individuals':[]}
        for h in [1,2]:
            p=HERE/f'input/shapes/h.t{it}.ah{h}.gp'
            if not p.is_file():frame['missing_individuals'].append(h);continue
            P,T,_=shape_gp(p,stride=2)
            if gp_tris is None:gp_tris=T
            else:assert np.array_equal(gp_tris,T)
            limit=max(limit,float(np.abs(P).max()))
            frame['objects'].append({'horizon':h,'mesh':'gp','xyz':np.round(P,8).ravel().tolist()})
        if it>=9536:
            P=data[f'points_{it}'];limit=max(limit,float(np.abs(P).max()))
            frame['objects'].append({'horizon':3,'mesh':'qlm','xyz':np.round(P,8).ravel().tolist()})
            frame['area']=float(geo[it]['area_physical_fejer']);frame['span_ratio']=float(geo[it]['coordinate_span_max_over_min'])
            frame['fields']={};frame['masks']={}
            for f in FIELDS:
                F=data[f'field_{f}_{it}'];range_abs[f]=max(range_abs[f],float(np.abs(F).max()))
                frame['fields'][f]=np.round(F,8).tolist()
                frame['masks'][f]=np.flatnonzero(data[f'mask_{f}_{it}']).tolist();mask_checks+=1
        frames.append(frame)
    D={'frames':frames,'fields':FIELDS,'meshes':{'qlm':previous['triangles'],'gp':gp_tris.tolist()},
       'ranges':{f:[0 if f=='abs_npsigma' else -range_abs[f],range_abs[f]] for f in FIELDS},
       'nt':35,'np':72,'birth_index':9536//32,'birth_iteration':9536,'full_bound':float(np.ceil(limit*11)/10),
       'common_frames':55,'individual_gp_stride_for_display':2,'individual_original_patch_nodes':19}
    serialized=json.dumps(D,separators=(',',':'),allow_nan=False)
    (OUT/'BH_HF4_B4_view_data.json').write_text(serialized)
    vendor=(BASE/'hf4_b2/vendor/plotly-3.3.1.min.js').read_text()
    template=(HERE/'viewer_template.html').read_text()
    controller=re.search(r'<script id="controller">(.*?)</script>',template,re.S)[1]
    (OUT/'BH_HF4_B4_controller.js').write_text(controller)
    html=template.replace('__DATA__',serialized.replace('</','<\\/')).replace('__PLOTLY__',vendor)
    assert '__DATA__' not in html and '__PLOTLY__' not in html
    (OUT/'BH_HF4_B4_fields_3D.html').write_text(html)
    # Plot scientific diagnostics from saved tables only.
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10})
    ar=csvrows(OUT/'BH_HF4_B4_area_comparison.csv');lo=csvrows(OUT/'BH_HF4_B4_localization.csv');pairs=csvrows(OUT/'BH_HF4_B4_temporal_pairs.csv')
    tau=np.array([float(r['tau_M']) for r in ar]);fig,ax=plt.subplots(2,2,figsize=(12,8))
    colors={'repsi2':'#b14e31','impsi2':'#246f9b','abs_npsigma':'#774ba0'}
    labels={'repsi2':'Re Ψ₂ · продольные концы','impsi2':'Im Ψ₂ · полярные области','abs_npsigma':'|σ| · полярные области'}
    a=ax[0,0];a.plot(tau,[float(r['area_recovered']) for r in ar],color='#246f9b',label='Из тетрады')
    a.plot(tau,[float(r['area_QLM']) for r in ar],ls='',marker='o',markevery=5,ms=4,mfc='none',color='#c37637',label='Сохранённая QLM')
    a.set(title='Физическая площадь общего горизонта',ylabel='Площадь / M²');a.legend(fontsize=9)
    a=ax[0,1];a.semilogy(tau,[max(abs(float(r['relative_recovered_minus_QLM'])),1e-16) for r in ar],color='#246f9b',label='Восстановление и QLM')
    a.semilogy(tau,[abs(float(r['relative_AHFinderDirect_minus_QLM'])) for r in ar],color='#c37637',label='AHFinderDirect и QLM')
    a.set(title='Согласование способов определения площади',ylabel='Относительное расхождение');a.legend(fontsize=9)
    a=ax[1,0]
    for f in FIELDS:
        rows=[r for r in lo if r['field']==f and float(r['quantile'])==.95]
        key='long_axis_fraction' if f=='repsi2' else 'polar_fraction'
        values=[float(r[key]) if r[key]!='' else np.nan for r in rows]
        a.plot([float(r['tau_M']) for r in rows],values,color=colors[f],label=labels[f],lw=1.6)
    a.set(title='Локализация на полном полученном участке',ylabel='Доля выделенной площади',ylim=(.65,1.04));a.legend(fontsize=8,loc='lower right')
    a=ax[1,1]
    for f in FIELDS:
        rows=[r for r in pairs if r['field']==f]
        a.plot([(int(r['iteration_to'])-9536)/512 for r in rows],[float(r['Dice_physical_measure']) for r in rows],color=colors[f],label=labels[f].split(' · ')[0])
    a.set(title='Совпадение масок соседних кадров',ylabel='Перекрытие Dice',ylim=(.77,1.025));a.legend(fontsize=9)
    for a in ax.ravel():a.set_xlabel('Возраст общего горизонта τ / M');a.grid(alpha=.2)
    fig.suptitle('HF4-B4 · 55 кадров общего горизонта · ранняя площадь проверена',fontsize=15,y=.99)
    fig.text(.055,.025,'Прежние 12 кадров переиспользованы. В двух кадрах продольная ось недоступна; маски и полярные оценки сохранены.',fontsize=9,color='#526171')
    fig.tight_layout(rect=[.015,.055,1,.96]);fig.savefig(OUT/'BH_HF4_B4_overview.png',dpi=170);plt.close(fig)
    checks={'viewer_common_frames':55,'viewer_time_slots':len(frames),'empty_slots':[x['iteration'] for x in frames if not x['objects']],
        'display_mask_count':mask_checks,'common_geometry_precision_decimals':8,'individual_angular_display_subsampling_only':True,
        'all_time_samples_included':True,'intermediate_time_or_shape_interpolation':False,'browser_WebGL_render_tested':False}
    (OUT/'BH_HF4_B4_view_checks.json').write_text(json.dumps(checks,indent=2)+'\n')
    print(json.dumps({'viewer_bytes':len(html.encode()),'data_bytes':len(serialized),'time_slots':len(frames),'common_frames':55,'bound':D['full_bound']}))


if __name__=='__main__':build()
