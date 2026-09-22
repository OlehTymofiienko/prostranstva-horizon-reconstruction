#!/usr/bin/env python3
"""R7B: проверка полученных данных QC0. Новая эволюция не запускается."""
from __future__ import annotations
import argparse
import hashlib
import importlib.util
import io
import json
import platform
import re
import zipfile
from pathlib import Path

import numpy as np
import scipy
from scipy.optimize import least_squares

BASE = Path(__file__).resolve().parent
VERSION = 'BH-HF3-R7B-2026-09-15-V1'


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def dump(p, obj):
    Path(p).write_text(json.dumps(obj, ensure_ascii=False, indent=2, allow_nan=False)+'\n', encoding='utf-8')


def array(blob):
    return np.loadtxt(io.BytesIO(blob), comments='#', ndmin=2)


def columns(blob):
    for line in blob.decode().splitlines():
        if line.startswith('# data columns:'):
            return {name:int(j)-1 for j,name in re.findall(r'(\d+):([^\s]+)', line)}
    raise ValueError('Не найден заголовок колонок')


def merge(arrays, time_column=0):
    rows = {}
    for a in arrays:
        for row in a:
            t=float(row[time_column])
            if t in rows and not np.array_equal(rows[t],row,equal_nan=True):
                raise ValueError(f'Несовпадение повторной строки при t={t}')
            rows[t]=row
    return np.array([rows[t] for t in sorted(rows)])


def signal(r, m, oldzip, newzip):
    name=f'mp_Psi4_l2_m{m}_r{r:.2f}.asc'
    chunks=[array(oldzip.read(n)) for n in oldzip.namelist() if n.endswith('/'+name) and '/provenance/' not in n]
    chunks.append(array(newzip.read('t95_t105/'+name)))
    a=merge(chunks)
    return a[:,0], a[:,1]+1j*a[:,2]


def profile(t, y, omega, alpha):
    u=t-t[0]
    v=np.exp((-alpha-1j*omega)*u)
    c=np.vdot(v,y)/np.vdot(v,v)
    model=c*v
    rss=float(np.vdot(y-model,y-model).real)
    return model,c,rss


def fit(t, y, reference):
    if len(t)<8 or not np.all(np.diff(t)>0) or not np.all(np.isfinite(y)):
        raise ValueError('Непригодное окно подгонки')
    norm=float(np.linalg.norm(y))
    if norm<=0:raise ValueError('Нулевой сигнал')
    yn=y/norm
    omega_log=-np.polyfit(t-t[0],np.unwrap(np.angle(y)),1)[0]
    alpha_log=-np.polyfit(t-t[0],np.log(abs(y)),1)[0]
    seeds=[[reference['omega'],reference['alpha']],
           [np.clip(omega_log,.201,.899),np.clip(alpha_log,.0011,.249)],
           [.4,.05],[.7,.15]]
    def residual(p):
        model,_,_=profile(t,yn,*p)
        d=model-yn
        return np.concatenate([d.real,d.imag])
    trials=[least_squares(residual,x,bounds=([.2,.001],[.9,.25]),
                          xtol=1e-12,ftol=1e-12,gtol=1e-12,max_nfev=2000)
            for x in seeds]
    best=min(trials,key=lambda x:np.dot(x.fun,x.fun))
    omega,alpha=map(float,best.x)
    free,c_free,rssf=profile(t,y,omega,alpha)
    fixed,c_fixed,rssk=profile(t,y,reference['omega'],reference['alpha'])
    if rssf>rssk+norm**2*1e-12:raise ValueError('Свободный минимум хуже вложенной фиксированной модели')
    n=2*len(t)
    bic=float(n*np.log(max(rssk,1e-300)/max(rssf,1e-300))-2*np.log(n))
    out={'samples_complex':len(t),'first_sample_M':float(t[0]),'last_sample_M':float(t[-1]),
         'omega':omega,'alpha':alpha,'signed_phase_slope':float(-omega_log),
         'omega_logphase':float(omega_log),'alpha_logamplitude':float(alpha_log),
         'relative_omega_difference':abs(omega/reference['omega']-1),
         'relative_alpha_difference':abs(alpha/reference['alpha']-1),
         'fixed_nrmse':float(np.sqrt(rssk)/norm),'free_nrmse':float(np.sqrt(rssf)/norm),
         'fixed_rss':rssk,'free_rss':rssf,'signal_energy':norm**2,
         'fixed_amplitude_re_im':[float(c_fixed.real),float(c_fixed.imag)],
         'free_amplitude_re_im':[float(c_free.real),float(c_free.imag)],
         'nominal_delta_BIC_fixed_minus_free':bic,'nominal_BIC_n_real':n,
         'nominal_BIC_interpretation':'Некалиброванный описательный показатель при фиктивном предположении независимых одинаковых ошибок.',
         'all_optimizers_successful':bool(all(x.success for x in trials)),
         'optimizer_solution_spread':np.ptp(np.stack([x.x for x in trials]),axis=0).tolist(),
         'at_parameter_boundary':bool(omega<.20001 or omega>.89999 or alpha<.00101 or alpha>.24999)}
    return out, fixed, free


def self_check(reference):
    t=np.arange(47)*.25+92.5
    records=[]
    for omega,alpha,c in [(reference['omega'],reference['alpha'],.7-.9j),(.51,.074,-.23+.18j)]:
        y=c*np.exp((-alpha-1j*omega)*(t-t[0]))
        result,_,_=fit(t,y,reference)
        passed=bool(abs(result['omega']-omega)<1e-8 and abs(result['alpha']-alpha)<1e-8 and result['free_nrmse']<1e-9)
        records.append({'true_omega':omega,'true_alpha':alpha,'recovered_omega':result['omega'],
                        'recovered_alpha':result['alpha'],'relative_residual':result['free_nrmse'],'pass':passed})
    if not all(r['pass'] for r in records):raise ValueError('Не пройдена проверка на известных аналитических сигналах')
    return records


def technical_audit(oldzip,newzip,extracted,helper):
    summary=json.loads(newzip.read('provenance/BH_HF3_R7A_summary.json'))
    status=json.loads(newzip.read('provenance/BH_HF3_R7_run_status.json'))
    manifest=json.loads(newzip.read('provenance/qc0-horizon-multipoles-lite-t105-recovery_manifest.json'))
    parname='qc0-horizon-multipoles-lite-t105-recover.par'
    source=newzip.read('provenance/'+helper.OLD_PAR_NAME).decode()
    expected,_=helper.build_parameter(source,Path('/home/alex/simulations/qc0_lite_t95_recover'),Path('/home/alex/simulations/qc0_lite_t105_recover'))
    checks={'uploaded_classification_pass':summary['technical_pass'],
            'source_parameter_hash':hashlib.sha256(source.encode()).hexdigest()==helper.SOURCE_SHA,
            'new_parameter_exactly_as_prepared':newzip.read('provenance/'+parname).decode()==expected,
            'parameter_hash_matches_manifest':sha(extracted/'provenance'/parname)==manifest['generated_sha256'],
            'executed_program_matches_prepared':sha(extracted/'provenance/prepare_qc0_t105_recovery.py')==sha(BASE/'prepare_qc0_t105_recovery.py'),
            'wrapper_exit_0':status['process_exit_code']==0 and status['guard_error'] is None,
            'wrapper_reached_105':status['last_time_M']==105 and status['last_iteration']==53760}
    log=newzip.read('provenance/ET_QC0_lite_t105_recover.log').decode(errors='replace')
    restart=re.findall(r'restarting simulation .* at iteration (\d+) \(simulation time ([0-9.]+)\)',log)
    checks.update(restart_95_verified=bool(restart and all(int(it)==48640 and float(t)==95 for it,t in restart)),
                  native_termination_105='Terminating due to cctk_final_time at t = 105.000000' in log,
                  native_final_checkpoint='Dumping termination checkpoint at iteration 53760, simulation time 105' in log,
                  native_exit_0=bool(re.search(r'^\s*Exit status:\s*0\s*$',log,re.MULTILINE)))
    waves=[]
    for r in helper.RADII:
        for ell in range(5):
            for m in range(-ell,ell+1):
                name=f'mp_Psi4_l{ell}_m{m}_r{r:.2f}.asc'
                waves.append(helper.audit_mode(extracted/'t70_t95'/name,extracted/'t95_t105'/name))
    qlm=[helper.audit_qlm(extracted/'t70_t95'/n,extracted/'t95_t105'/n) for n in helper.QLM_FILES]
    # Сравнение байтов старых данных, а не повторный физический анализ R6A.
    old_shared=[]
    for n in newzip.namelist():
        if n.startswith('t70_t95/'):
            candidates=[p for p in oldzip.namelist() if p.endswith('/'+n)]
            if candidates:old_shared.append({'file':n,'same':newzip.read(n)==oldzip.read(candidates[0])})
    checks['all_shared_old_bytes_unchanged']=bool(old_shared and all(p['same'] for p in old_shared))
    checks['all_125_wave_seams_pass']=len(waves)==125 and all(w['pass'] for w in waves)
    checks['all_3_qlm_files_pass']=all(q['pass'] for q in qlm)
    checks['checkpoint_files_reported_by_local_runner']=len(summary['new_checkpoint_files'])==2 and all(x['header_is_hdf5'] for x in summary['new_checkpoint_files'])
    checks['source_checkpoint_preservation_reported_by_local_runner']=summary['source_checkpoints_unchanged']
    return {'checks':checks,'pass':all(checks.values()),'wave_seams':waves,'qlm':qlm,
            'old_files_compared':old_shared,'checkpoint_scope':'Журнал и локальный реестр; гигабайтных контрольных точек в архиве нет.',
            'elapsed_seconds':status['elapsed_seconds'],'new_checkpoint_inventory':summary['new_checkpoint_files']}


def horizon(newzip,reference,out):
    name='quasilocalmeasures-qlm_scalars..asc'
    blobs=[newzip.read(s+'/'+name) for s in ['t70_t95','t95_t105']]
    col=columns(blobs[0]); a=merge([array(b) for b in blobs],8)
    t=a[:,8]; mass=a[:,col['qlm_mass[2]']];spin=a[:,col['qlm_spin[2]']];chi=spin/mass**2
    records=[]
    for left,right in [(90,95),(95,100),(100,105)]:
        mask=(t>=left)&(t<=right)
        r={'window_M':[left,right],'samples':int(mask.sum())}
        for key,y in [('Mf',mass),('Jf',spin),('chi',chi)]:
            r[key]={'mean':float(np.mean(y[mask])),'median':float(np.median(y[mask])),
                    'relative_range':float(np.ptp(y[mask])/np.mean(y[mask])),
                    'relative_mean_change_from_frozen':float(np.mean(y[mask])/reference[key]-1)}
        records.append(r)
    np.savetxt(out/'BH_HF3_R7B_horizon.tsv',np.column_stack([t,mass,spin,chi]),delimiter='\t',header='t_M\tMf\tJf\tchi',fmt='%.17g')
    return records


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input',type=Path,default=BASE/'inputs/BH_HF3_R7A_output.zip')
    parser.add_argument('--previous',type=Path,default=BASE/'inputs/BH_HF3_R6A_audit_bundle.zip')
    parser.add_argument('--protocol',type=Path,default=BASE/'results/BH_HF3_R7B_protocol.json')
    parser.add_argument('--output',type=Path,default=BASE/'results')
    args=parser.parse_args();out=args.output;out.mkdir(parents=True,exist_ok=True)
    protocol=json.loads(args.protocol.read_text());ref=protocol['reference'];T=protocol['period_M']
    spec=importlib.util.spec_from_file_location('r7helper',BASE/'prepare_qc0_t105_recovery.py')
    helper=importlib.util.module_from_spec(spec);spec.loader.exec_module(helper)
    checks=self_check(ref);dump(out/'BH_HF3_R7B_fit_checks.json',checks)
    with zipfile.ZipFile(args.previous) as oldzip, zipfile.ZipFile(args.input) as newzip:
        if newzip.testzip() is not None:raise ValueError('Ошибка CRC нового архива')
        extracted=out.parent/'input';extracted.mkdir(exist_ok=True)
        for n in newzip.namelist():
            if not (extracted/n).resolve().is_relative_to(extracted.resolve()):raise ValueError('Некорректный путь ZIP')
        newzip.extractall(extracted)
        audit=technical_audit(oldzip,newzip,extracted,helper)
        dump(out/'BH_HF3_R7B_input_audit.json',audit)
        if not audit['pass']:raise ValueError('Технические условия R7B не пройдены')
        print('Технический аудит подтверждён.',flush=True)
        h=horizon(newzip,ref,out);dump(out/'BH_HF3_R7B_horizon.json',h)
        cache={}
        def get(r,m=2):
            if (r,m) not in cache:cache[r,m]=signal(r,m,oldzip,newzip)
            return cache[r,m]
        cases=[]
        for r,period in protocol['primary_cases']+protocol['age_control_cases']:
            t,y=get(r);peak=protocol['peaks_M'][str(r)]
            left,right=peak+(period-1)*T,peak+period*T
            mask=(t>=left)&(t<=right)
            result,k,f=fit(t[mask],y[mask],ref)
            result.update(radius_M=r,period=period,requested_window_M=[left,right],peak_M=peak,
                          case_kind='primary' if [r,period] in protocol['primary_cases'] else 'age_control')
            tm,ym=get(r,-2)
            if not np.array_equal(tm,t):raise ValueError('Несовпадение времён m=±2')
            result['minus2_conjugacy_relative_difference']=float(np.linalg.norm(ym[mask]-y[mask].conj())/np.linalg.norm(y[mask]))
            cases.append(result)
            np.savetxt(out/f'BH_HF3_R7B_r{r}_period{period}.tsv',
                       np.column_stack([t[mask],y[mask].real,y[mask].imag,k.real,k.imag,f.real,f.imag]),
                       delimiter='\t',header='t_M\tPsi4_re\tPsi4_im\tKerr_re\tKerr_im\tfree_re\tfree_im',fmt='%.17g')
            print(f'r={r}, период={period}: omega={result["omega"]:.9f}, alpha={result["alpha"]:.9f}, '
                  f'ошибки={100*result["relative_omega_difference"]:.3f}% / {100*result["relative_alpha_difference"]:.3f}%, '
                  f'остатки={100*result["fixed_nrmse"]:.3f}% / {100*result["free_nrmse"]:.3f}%, '
                  f'условный dBIC={result["nominal_delta_BIC_fixed_minus_free"]:.3f}',flush=True)
        dump(out/'BH_HF3_R7B_cases.json',cases)
        sensitivity=[]
        t,y=get(40);left=protocol['peaks_M']['40']+2*T;right=left+T
        for shift in protocol['sensitivity_r40_third']['window_shift_M']:
            ids=np.flatnonzero((t>=left+shift)&(t<=right+shift))
            result,_,_=fit(t[ids],y[ids],ref)
            sensitivity.append(dict(result,window_shift_M=shift,subset='all'))
        ids=np.flatnonzero((t>=left)&(t<=right))
        for offset,name in [(0,'even'),(1,'odd')]:
            ii=ids[offset::2];result,_,_=fit(t[ii],y[ii],ref)
            sensitivity.append(dict(result,window_shift_M=0,subset=name))
        dump(out/'BH_HF3_R7B_sensitivity.json',sensitivity)
        split=left+2*T/3
        train=(t>=left)&(t<=split);test=(t>split)&(t<=right)
        result,_,_=fit(t[train],y[train],ref)
        u=t[test]-t[train][0]
        holdout={'training_requested_M':[left,split],'validation_requested_M':[split,right],
                 'training_fit':result,'validation_samples':int(test.sum())}
        predicted=[]
        for model,omega,alpha in [('fixed',ref['omega'],ref['alpha']),('free',result['omega'],result['alpha'])]:
            c=complex(*result[model+'_amplitude_re_im'])
            yp=c*np.exp((-alpha-1j*omega)*u)
            holdout[model+'_prediction_nrmse']=float(np.linalg.norm(yp-y[test])/np.linalg.norm(y[test]))
            predicted.extend([yp.real,yp.imag])
        np.savetxt(out/'BH_HF3_R7B_r40_holdout.tsv',np.column_stack([t[test],y[test].real,y[test].imag,*predicted]),
                   delimiter='\t',header='t_M\tPsi4_re\tPsi4_im\tKerr_prediction_re\tKerr_prediction_im\tfree_prediction_re\tfree_prediction_im',fmt='%.17g')
        dump(out/'BH_HF3_R7B_holdout.json',holdout)
        age_series={}
        for r in [15,30,40]:
            tt,yy=get(r);age=tt-protocol['peaks_M'][str(r)];sel=(age>=2*T)&(age<=3*T)
            age_series[r]=(age[sel],r*yy[sel])
        shape=[]
        for a,b in [(15,30),(15,40),(30,40)]:
            x,v=age_series[a];xt,w=age_series[b]
            if not np.array_equal(x,xt):raise ValueError('Не совпадает сетка возраста кольцевания')
            phase=np.exp(-1j*np.angle(np.vdot(v,w)))
            corr=abs(np.vdot(v,w))/(np.linalg.norm(v)*np.linalg.norm(w))
            shape.append({'reference_radius_M':a,'compared_radius_M':b,'samples':len(x),
                          'coherence':float(corr),'norm_amplitude_ratio':float(np.linalg.norm(w)/np.linalg.norm(v)),
                          'constant_phase_applied_rad':float(np.angle(phase)),
                          'phase_only_nrmse':float(np.linalg.norm(phase*w-v)/np.linalg.norm(v))})
        dump(out/'BH_HF3_R7B_third_period_shape.json',shape)
        summary={'stage':'BH-HF3-R7B','version':VERSION,'input_audit_pass':audit['pass'],
                 'reference':ref,'period_M':T,'primary_cases':[r for r in cases if r['case_kind']=='primary'],
                 'age_control_cases':[r for r in cases if r['case_kind']=='age_control'],
                 'r40_sensitivity_ranges':{k:[min(s[k] for s in sensitivity),max(s[k] for s in sensitivity)]
                                          for k in ['omega','alpha','relative_omega_difference','relative_alpha_difference','fixed_nrmse','free_nrmse','nominal_delta_BIC_fixed_minus_free']},
                 'holdout':holdout,'third_period_shapes':shape,'horizon_windows':h,
                 'provenance':{'new_archive_sha256':sha(args.input),'previous_archive_sha256':sha(args.previous),
                               'protocol_sha256':sha(args.protocol),'program_sha256':sha(__file__),
                               'python':platform.python_version(),'numpy':np.__version__,'scipy':scipy.__version__},
                 'numerical_frequency_damping_goals_all_three_radii':all(r['relative_omega_difference']<=.01 and r['relative_alpha_difference']<=.05 for r in cases if r['period']==3),
                 'historical_nominal_BIC_gate_all_three_radii':all(r['nominal_delta_BIC_fixed_minus_free']<=10 for r in cases if r['period']==3),
                 'no_new_evolution_run_here':True,'resolution_convergence_test_available':False,
                 'calibrated_statistical_Kerr_acceptance':False,'new_physics_established':False}
        dump(out/'BH_HF3_R7B_summary.json',summary)
        print('Численный анализ R7B сохранён.',flush=True)


if __name__=='__main__':
    main()
