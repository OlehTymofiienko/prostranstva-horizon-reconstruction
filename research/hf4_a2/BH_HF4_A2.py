#!/usr/bin/env python3
"""HF4-A2: reproducible time atlas of horizon scalar fields. No evolution.

Python >=3.10, numpy, scipy>=1.15, matplotlib, h5py; ffmpeg for the video.
The accompanying atlas_template.html is needed to produce the offline viewer.
"""
from __future__ import annotations
import argparse
import base64
import hashlib
import json
import platform
from pathlib import Path
import shutil

import h5py
import numpy as np
import scipy
from scipy.special import sph_harm_y
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import AsinhNorm, Normalize
from matplotlib.animation import FFMpegWriter

SHA = '9cdd268e85b783d8a832b72e46945c6b43bb2b3ea0a9e6665e57a5a1b7d56f0f'
OMEGA = 0.20878380853617037
ELL = np.array([l for l in range(9) for m in range(-l,l+1)])
EM = np.array([m for l in range(9) for m in range(-l,l+1)])
TIMES = [0, 5, 10, 20, 40, 80, 120, 150]


def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def dump(p, a): p.write_text(json.dumps(a, ensure_ascii=False, indent=2, allow_nan=False)+'\n')
def b64(a): return base64.b64encode(np.asarray(a,dtype='<f8').tobytes()).decode('ascii')


def pack(c):
    """Real harmonic representation; 81 reals, l block starts at l*l.
    m0: Re(c_l0); each positive m: 2Re(c_lm), -2Im(c_lm).
    Corresponding basis: Y_l0, Re(Y_lm), Im(Y_lm).
    """
    out = np.empty(c.shape[:-1]+(81,))
    for l in range(9):
        out[...,l*l] = c[...,l*(l+1)].real
        for m in range(1,l+1):
            out[...,l*l+2*m-1] = 2*c[...,l*(l+1)+m].real
            out[...,l*l+2*m] = -2*c[...,l*(l+1)+m].imag
    return out


def complex_basis(theta, phi):
    th, ph = np.broadcast_arrays(theta, phi)
    return np.array([sph_harm_y(l,m,th.ravel(),ph.ravel()) for l,m in zip(ELL,EM)])


def real_basis(theta, phi):
    th, ph = np.broadcast_arrays(theta, phi)
    ys=[]
    for l in range(9):
        ys.append(sph_harm_y(l,0,th,ph).real.ravel())
        for m in range(1,l+1):
            y=sph_harm_y(l,m,th,ph)
            ys.extend([y.real.ravel(),y.imag.ravel()])
    return np.array(ys)


def frame_audit(t, fields, means, delta):
    u = np.exp(-1j*t[:,None]*OMEGA*EM)
    results={}
    th=np.array([.17,.52,1.23,1.57,2.33,2.91])
    ph=np.array([-2.31,-.44,.25,1.32,2.71,-1.12])
    y=complex_basis(th,ph)
    for name,c,mu,d in zip(['K','B'],fields,means,delta):
        norms=np.linalg.norm(d,axis=1)
        rot=d*u
        s={'max_relative_norm_change':float(np.max(np.abs(np.linalg.norm(rot,axis=1)/norms-1))),
           'max_abs_covariant_background_subtraction_error':float(np.max(np.abs((c*u-mu*u)-rot)))}
        residual=[]
        for tm in TIMES:
            i=int(round(tm*10))
            residual.append(np.max(np.abs(rot[i]@y-d[i]@complex_basis(th,ph-OMEGA*t[i]))))
        s['coordinate_shift_identity_max_abs_error']=float(max(residual))
        # Deliberately wrong order as a diagnostic, never used for output maps.
        all_u=np.exp(-1j*TIME_ALL[:,None]*OMEGA*EM)
        bad_mu=(ALL_FIELDS[['K','B'].index(name)]*all_u)[TIME_ALL>=400].mean(axis=0)
        bad_delta=c*u-bad_mu
        s['wrong_order_example']={str(tm):float(np.linalg.norm(bad_delta[int(tm*10)]-rot[int(tm*10)])/norms[int(tm*10)]) for tm in [80,120,150]}
        s['wrong_order_is_not_used']=True
        results[name]=s
    assert max(x['max_relative_norm_change'] for x in results.values())<1e-12
    assert max(x['coordinate_shift_identity_max_abs_error'] for x in results.values())<1e-11
    return results


def make_plot(out,t,norms,differences,baseaudit):
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11,'axes.spines.top':False,
                         'axes.spines.right':False,'axes.grid':True,'grid.alpha':.18})
    fig,ax=plt.subplots(2,2,figsize=(12.8,8.8),layout='constrained')
    fig.suptitle('HF4-A2 · Что меняется после вычитания позднего состояния',fontsize=17)
    for i,(name,label,color) in enumerate([('K','Кривизна','#225fad'),('B','Поле вращения','#9c3977')]):
        ax[0,0].semilogy(t,norms[i]/max(norms[i]),label=label,color=color)
        for cut,ls in [(4,'--'),(6,'-')]:
            ax[0,1].plot(t,100*differences[f'{name}_delta_L{cut}_vs_L8'],ls,color=color,label=f'{label}, ℓ ≤ {cut}')
        ax[1,i].semilogy(t,100*differences[f'{name}_full_L6_vs_L8'],color='#a6b8ca',label='Полное поле')
        ax[1,i].semilogy(t,100*differences[f'{name}_delta_L6_vs_L8'],color=color,label='После вычитания позднего среднего')
        ax[1,i].set(title=label+': ℓ ≤ 6 против ℓ ≤ 8',xlabel='Время SpEC, M',ylabel='Относительная разность по норме L², %')
        ax[1,i].legend(fontsize=9)
    ax[0,0].set(title='Величина отклонений от позднего среднего',xlabel='Время SpEC, M',ylabel='Норма / её максимум на 0–150 M')
    ax[0,0].legend(fontsize=9)
    ax[0,1].set(title='Чувствительность возмущений к отсечению',xlabel='Время SpEC, M',ylabel='Относительная разность по норме L², %')
    ax[0,1].legend(fontsize=9)
    fig.savefig(out/'BH_HF4_A2_checks.png',dpi=165)
    plt.close(fig)


def make_movie(out, times, maps, limits, t, norms, fps=15):
    if shutil.which('ffmpeg') is None: raise RuntimeError('ffmpeg is required for --video')
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11,'axes.grid':False})
    fig=plt.figure(figsize=(12.8,9.6),dpi=100)
    gs=fig.add_gridspec(3,2,height_ratios=[1,1,.60],left=.065,right=.94,top=.865,bottom=.15,hspace=.53,wspace=.28)
    fig.suptitle('HF4-A2 · Динамика скалярных полей общего горизонта',fontsize=17,y=.985)
    title=fig.text(.5,.938,'',ha='center',fontsize=12)
    images=[]
    panel_names=['Полная кривизна  K','Полное поле вращения  B','Отклонение кривизны  δK','Отклонение поля вращения  δB']
    for i in range(4):
        ax=fig.add_subplot(gs[i//2,i%2])
        lim=limits[i]
        norm=Normalize(*lim) if i<2 else AsinhNorm(linear_width=1e-5,vmin=lim[0],vmax=lim[1])
        im=ax.imshow(maps[i][0],origin='upper',extent=(-180,180,-90,90),aspect='auto',
                     cmap='viridis' if i==0 else 'RdBu_r',norm=norm,interpolation='bilinear')
        ax.set_title(panel_names[i],pad=9,fontsize=12)
        ax.set_xticks([-180,0,180]);ax.set_yticks([-90,0,90]);ax.tick_params(labelsize=9)
        ax.set_xlabel('Азимут, °',fontsize=9);ax.set_ylabel('Широта, °',fontsize=9)
        cb=fig.colorbar(im,ax=ax,pad=.025,fraction=.045)
        cb.ax.tick_params(labelsize=8)
        if i>=2: cb.set_ticks([lim[0],-.01,-.00001,0,.00001,.01,lim[1]])
        images.append(im)
    ax=fig.add_subplot(gs[2,:])
    for i,label in enumerate(['δK','δB']):
        ax.semilogy(t,norms[i]/max(norms[i]),label=label,color=['#225fad','#9c3977'][i])
    line=ax.axvline(0,color='#222',lw=1)
    ax.set(xlim=(0,150),ylim=(1e-6,2),xlabel='Время SpEC после появления общего горизонта, M',ylabel='Норма / максимум')
    ax.grid(alpha=.2);ax.legend(loc='upper right',fontsize=9)
    fig.text(.5,.069,'Внизу вычтено среднее 400–500 M; сжатые шкалы цвета asinh неизменны во всех кадрах.',ha='center',fontsize=10)
    fig.text(.5,.035,'Карты в авторских углах · ℓ ≤ 8 · цвет показывает поле; физическая форма поверхности здесь не задана.',ha='center',fontsize=10,color='#475569')
    writer=FFMpegWriter(fps=fps,codec='libx264',bitrate=2800,extra_args=['-pix_fmt','yuv420p','-movflags','+faststart','-threads','2'])
    with writer.saving(fig,str(out/'BH_HF4_A2_dynamics.mp4'),dpi=100):
        for j,tm in enumerate(times):
            title.set_text(f't = {tm:5.1f} M   ·   исходные переносимые углы   ·   постоянные шкалы цвета')
            for i,im in enumerate(images): im.set_data(maps[i][j])
            line.set_xdata([tm,tm]);writer.grab_frame()
            if j in [0,60,120,180,240,300]: print(f'video frame {j+1}/{len(times)}',flush=True)
            if j==20: fig.savefig(out/'BH_HF4_A2_preview.png',dpi=130)
    plt.close(fig)


def main():
    global ALL_FIELDS,TIME_ALL
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--moments',type=Path,required=True)
    ap.add_argument('--parent-state',type=Path,required=True)
    ap.add_argument('--out',type=Path,default=Path('results'))
    ap.add_argument('--video',action='store_true')
    a=ap.parse_args();out=a.out;out.mkdir(parents=True,exist_ok=True)
    assert sha(a.moments)==SHA
    parent=json.loads(a.parent_state.read_text())
    assert parent['last_completed_stage']=='BH-HF4-A1'
    with h5py.File(a.moments) as f:
        TIME_ALL=f['time_adaptive'][:]
        ALL_FIELDS=np.array([2*f['geo_mass/Adaptive'][:].conj(),f['geo_spin/Adaptive'][:].conj()])
    allc=ALL_FIELDS;ta=TIME_ALL
    assert allc.shape==(2,5001,81) and np.isfinite(allc).all()
    select=ta<=150;t=ta[select];c=allc[:,select]
    mu=allc[:,ta>=400].mean(axis=1)
    d=c-mu[:,None,:]
    norms=np.linalg.norm(d,axis=2)
    fullnorms=np.linalg.norm(c,axis=2)
    differences={}
    baseline={}
    snapshots=[]
    for i,name in enumerate(['K','B']):
        for cut in [4,6]:
            differences[f'{name}_delta_L{cut}_vs_L8']=np.linalg.norm(d[i][:,ELL>cut],axis=1)/norms[i]
            differences[f'{name}_full_L{cut}_vs_L8']=np.linalg.norm(c[i][:,ELL>cut],axis=1)/fullnorms[i]
        baseline[name]={}
        for lo,hi in [(300,400),(400,450),(450,500)]:
            alt=allc[i,(ta>=lo)&(ta<=hi)].mean(axis=0)
            dist=float(np.linalg.norm(alt-mu[i]))
            baseline[name][f'{lo}_{hi}']={'L2_shift':dist,'fraction_of_delta_at_150M':dist/norms[i,-1]}
        baseline[name]['RMS_time_variation_400_500']=float(np.sqrt(np.mean(np.linalg.norm(allc[i,ta>=400]-mu[i],axis=1)**2)))
        baseline[name]['is_not_total_error_bound']=True
    for tm in TIMES:
        j=int(round(tm*10));s={'t_M':float(t[j])}
        for i,name in enumerate(['K','B']):
            s[name+'_delta_L2']=float(norms[i,j]);s[name+'_fraction_of_max']=float(norms[i,j]/max(norms[i]))
            s.update({key:float(v[j]) for key,v in differences.items() if key.startswith(name)})
        snapshots.append(s)
    audit=frame_audit(t,c,mu,d)
    indices=np.arange(0,len(t),5);tf=t[indices]
    # Store pre-subtracted perturbations in double precision: never subtract rounded display values.
    packed=np.array([pack(c[0,indices]),pack(c[1,indices]),pack(d[0,indices]),pack(d[1,indices])])
    lat=np.linspace(90,-90,49);lon=np.linspace(-180,180,97)
    rb=real_basis(np.deg2rad(90-lat)[:,None],np.deg2rad(lon)[None,:])
    maps=(packed@rb).reshape(4,len(tf),len(lat),len(lon))
    # Same scale for all times, frames and cutoffs, with a 10% display margin.
    ext=np.zeros((4,2));ext[:,0]=np.inf;ext[:,1]=-np.inf
    for cut in [4,6,8]:
        arr=packed[:,:,:((cut+1)**2)]@rb[:((cut+1)**2)]
        ext[:,0]=np.minimum(ext[:,0],arr.min(axis=(1,2)))
        ext[:,1]=np.maximum(ext[:,1],arr.max(axis=(1,2)))
    limits=[]
    limits.append([float(np.floor(ext[0,0]*1.1*10)/10),float(np.ceil(ext[0,1]*1.1*10)/10)])
    for i in range(1,4):
        amp=float(np.ceil(max(abs(ext[i,0]),abs(ext[i,1]))*1.1*100)/100);limits.append([-amp,amp])
    # Independent complex-basis check of the real representation used by the viewer.
    sample=[0,20,160,300]
    yy=complex_basis(np.deg2rad(90-lat)[:,None],np.deg2rad(lon)[None,:])
    packing_errors=[]
    for i,v in enumerate([c[0],c[1],d[0],d[1]]):
        packing_errors.append(float(np.max(np.abs((v[indices[sample]]@yy).real-maps[i,sample].reshape(len(sample),-1)))))
    assert max(packing_errors)<1e-11
    cutsum={key:{'max_0_150':float(max(v)),'max_5_120':float(max(v[(t>=5)&(t<=120)]))} for key,v in differences.items()}
    summary={'stage':'BH-HF4-A2','date_UTC':'2026-09-16',
       'classification':'BH_HF4_A2_TIME_ATLAS_VALIDATED_FRAME_COVARIANCE_PASSED_EARLY_SPIN_PERTURBATIONS_REQUIRE_HIGHER_ANGULAR_MODES',
       'source_HDF5_sha256':sha(a.moments),'parent_state_sha256':sha(a.parent_state),
       'definitions':parent['reconstruction_conventions'],
       'background_window_M':[400,500], 'analysis_window_M':[0,150],
       'sample_step_M':.1,'atlas_frame_step_M':.5,'atlas_frames':len(tf),'atlas_grid':[49,97],
       'phase_rotation_omega':OMEGA,'rotation_rule':'c_lm -> c_lm exp(-im Omega t); both signal and original background transform at time t',
       'frame_audit':audit,'baseline_sensitivity':baseline,'snapshots':snapshots,'cutoff_summary':cutsum,
       'real_basis_packing_max_abs_errors':packing_errors,'fixed_color_limits':limits,
       'color_scale_order':['K','B','delta_K','delta_B'],'asinh_linear_width':1e-5,
       'delta_norm_maxima':{name:{'value':float(max(norms[i])),'time_M':float(t[np.argmax(norms[i])])} for i,name in enumerate(['K','B'])},
       'limits':['No numerical resolution comparison supplied','Differences between l<=4,6,8 do not bound l>8','Late baseline variation is not the total numerical error','Maps do not specify a radial shape or coordinate embedding','SpEC and QC0 are not combined'],
       'new_evolution_run':False,'new_physics_established':False,
       'environment':{'python':platform.python_version(),'numpy':np.__version__,'scipy':scipy.__version__,'h5py':h5py.__version__,'matplotlib':matplotlib.__version__}}
    dump(out/'BH_HF4_A2_summary.json',summary)
    np.savez_compressed(out/'BH_HF4_A2_analysis.npz',time_M=t,delta_K_L2=norms[0],delta_B_L2=norms[1],background_K=mu[0],background_B=mu[1],**differences)
    np.savez_compressed(out/'BH_HF4_A2_atlas_coefficients.npz',time_M=tf,packed_coefficients=packed,longitude_deg=lon,latitude_deg=lat)
    validation_points=[]
    for frame,field,cut,rot,row,col in [(0,0,8,False,24,48),(20,3,6,False,15,31),(160,2,8,True,26,61),(300,3,4,True,9,13)]:
        coeff=[c[0],c[1],d[0],d[1]][field][indices[frame]].copy()
        coeff[ELL>cut]=0
        phi=np.deg2rad(lon[col])-(OMEGA*tf[frame] if rot else 0)
        val=float((coeff@complex_basis(np.deg2rad(90-lat[row]),phi)).real[0])
        validation_points.append({'frame':frame,'field':field,'cutoff':cut,'rotate':rot,'row':row,'col':col,'expected':val})
    palette={name:(plt.get_cmap(name)(np.linspace(0,1,256))[:,:3]*255).astype(int).tolist() for name in ['viridis','RdBu_r']}
    payload={'times':tf.tolist(),'packed_b64':b64(packed),'omega':OMEGA,'limits':limits,'width':len(lon),'height':len(lat),
        'palettes':palette,'norms':norms[:,indices].tolist(),'cutoffs':{k:v[indices].tolist() for k,v in differences.items()},
        'validation_points':validation_points,'source_hash':SHA}
    template=Path(__file__).with_name('atlas_template.html').read_text()
    (out/'BH_HF4_A2_atlas.html').write_text(template.replace('__ATLAS_DATA__',json.dumps(payload,ensure_ascii=False,separators=(',',':'))))
    make_plot(out,t,norms,differences,baseline)
    print(json.dumps({'summary':str(out/'BH_HF4_A2_summary.json'),'snapshots':snapshots,'frame_audit':audit,'limits':limits,'cutoff_summary':cutsum},ensure_ascii=False),flush=True)
    if a.video: make_movie(out,tf,maps,limits,t,norms)


if __name__=='__main__': main()
