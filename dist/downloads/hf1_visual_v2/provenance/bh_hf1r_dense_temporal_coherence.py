#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""BH-HF1R — Dense Temporal Coherence Audit.

Follow-up to positive BH-HF1 on SXS:BBH:0305/Lev6 common horizon C.
Tests whether sparse localized high-gradient components form temporally coherent
tracks. Still coordinate/gauge dependent; not domain wall, not physical phase,
not QNM-subtracted, not resolution/frame/slicing robust yet.
"""
from __future__ import annotations
import argparse, importlib.util, json
from pathlib import Path
import h5py, numpy as np, pandas as pd, matplotlib.pyplot as plt
from scipy.spatial import cKDTree

STAGE='BH-HF1R'; SID='SXS:BBH:0305'; LEV=6; HORIZON='C'
WINDOW_END_M=12.0
KNN=8; ACTIVE_Q=.95; N_SHUFFLE=199; RNG_SEED=150914
P_GATE=.01; FRAC_GATE=.35
GRID_N=4096; GRID_KNN=8
PAIR_STEP_MAX=25.0; PAIR_DICE_MIN=.15
GAP_FLOOR=.75; GAP_FACTOR=2.5
TRACK_N_MIN=4; TRACK_DUR_MIN=1.0; TRACK_STEP_MED_MAX=15.0; TRACK_DICE_MED_MIN=.25
PRIMARY_FIELDS_MIN=2
RE_MAX=10.0; RB_ANTI_DEV_MAX=15.0; SECONDARY_N_MIN=3; SECONDARY_FRAC_MIN=.60
FIELDS=['Rstar','Estar','Bstar']

FROZEN_PROTOCOL={
 'stage':STAGE,'title':'Dense Temporal Coherence Audit','status':'FROZEN_BEFORE_LOCAL_EXECUTION',
 'simulation':SID,'lev':LEV,'horizon':HORIZON,'dense_window_M':[0.0,WINDOW_END_M],
 'snapshot_rule':{'shuffle_p_max':P_GATE,'largest_component_fraction_min':FRAC_GATE,
                  'active_threshold_quantile':ACTIVE_Q,'shuffle_trials':N_SHUFFLE},
 'tracking_rule':{'fixed_grid_points':GRID_N,'component_grid_dilation':'one 8-NN hop',
                  'pair_centroid_step_max_deg':PAIR_STEP_MAX,'pair_dilated_dice_min':PAIR_DICE_MIN,
                  'max_gap_M':'max(0.75, 2.5 * median_native_cadence)',
                  'track_min_snapshots':TRACK_N_MIN,'track_min_duration_M':TRACK_DUR_MIN,
                  'track_median_centroid_step_max_deg':TRACK_STEP_MED_MAX,
                  'track_median_dilated_dice_min':TRACK_DICE_MED_MIN,
                  'primary_min_fields_pass':PRIMARY_FIELDS_MIN},
 'secondary_cross_field_rule':{'status':'secondary; motivated by sparse HF1; not primary',
                  'RE_angle_max_deg':RE_MAX,'RB_antipodal_deviation_max_deg':RB_ANTI_DEV_MAX,
                  'min_simultaneous_positive_snapshots':SECONDARY_N_MIN,
                  'min_fraction_relation_pass':SECONDARY_FRAC_MIN},
 'forbidden_claims':['not a physical phase measurement','not a domain wall detection',
                     'not an intrinsic/covariant gradient measurement','not gauge invariant',
                     'not slicing invariant','not numerical-resolution robust yet','not QNM-subtracted',
                     'does not establish new physics','does not establish AARFI or APK']}

def imp(path,name):
    s=importlib.util.spec_from_file_location(name,str(path)); m=importlib.util.module_from_spec(s); s.loader.exec_module(m); return m

def robust_z(x):
    x=np.asarray(x,float); med=np.median(x); mad=np.median(np.abs(x-med)); s=1.4826*mad
    if not np.isfinite(s) or s<=1e-14: s=np.std(x)
    if not np.isfinite(s) or s<=1e-14: s=1.0
    return (x-med)/s

def unit_dirs(coords,center):
    x=np.asarray(coords,float)-np.asarray(center,float)[None,:]; r=np.linalg.norm(x,axis=1)
    good=np.isfinite(r)&(r>0)
    if good.sum()<max(20,int(.9*len(r))): raise RuntimeError('Too many invalid radii')
    return x[good]/r[good,None],good

def knn_graph(n,k):
    tr=cKDTree(n); _,idx=tr.query(n,k=min(k+1,len(n))); nei=idx[:,1:]
    dot=np.einsum('ij,ikj->ik',n,n[nei]); ang=np.arccos(np.clip(dot,-1,1))
    return nei,ang

def grad_proxy(z,nei,ang):
    return np.median(np.abs(z[nei]-z[:,None])/np.maximum(ang,1e-8),axis=1)

def largest_component(active,nei):
    ids=np.flatnonzero(active)
    if not len(ids): return 0,np.array([],int)
    aset=set(map(int,ids)); seen=set(); best=[]
    for s in ids:
        s=int(s)
        if s in seen: continue
        stack=[s]; seen.add(s); comp=[]
        while stack:
            u=stack.pop(); comp.append(u)
            for v in nei[u]:
                v=int(v)
                if v in aset and v not in seen: seen.add(v); stack.append(v)
        if len(comp)>len(best): best=comp
    return len(best),np.asarray(best,int)

def morphology(z,n,nei,ang,rng):
    g=grad_proxy(z,nei,ang); thr=np.quantile(g,ACTIVE_Q); active=g>=thr
    ln,li=largest_component(active,nei); an=int(active.sum()); frac=ln/max(an,1)
    null=[]
    for _ in range(N_SHUFFLE):
        zp=rng.permutation(z); gp=grad_proxy(zp,nei,ang); ap=gp>=np.quantile(gp,ACTIVE_Q)
        lnp,_=largest_component(ap,nei); null.append(lnp/max(int(ap.sum()),1))
    null=np.asarray(null); p=(1+np.sum(null>=frac-1e-15))/(N_SHUFFLE+1)
    if len(li):
        v=np.sum(n[li],axis=0); vn=np.linalg.norm(v); c=v/vn if vn>0 else np.full(3,np.nan)
    else: c=np.full(3,np.nan)
    return {'gradient_q95':float(thr),'active_points':an,'largest_component_points':int(ln),
            'largest_component_fraction_of_active':float(frac),'shuffle_p':float(p),
            'centroid_x':float(c[0]),'centroid_y':float(c[1]),'centroid_z':float(c[2]),
            'candidate_snapshot':bool(p<=P_GATE and frac>=FRAC_GATE),'_li':li}

def fibonacci_sphere(N):
    i=np.arange(N,dtype=float); phi=(1+np.sqrt(5))/2; z=1-2*(i+.5)/N
    th=2*np.pi*i/phi; r=np.sqrt(np.maximum(0,1-z*z))
    return np.c_[r*np.cos(th),r*np.sin(th),z]

def mask_on_grid(comp_dirs,grid,tree):
    m=np.zeros(len(grid),bool)
    if len(comp_dirs): _,ii=tree.query(comp_dirs,k=1); m[np.asarray(ii,int)]=True
    return m

def dilate(m,nei):
    out=m.copy(); ids=np.flatnonzero(m)
    if len(ids): out[np.unique(nei[ids].ravel())]=True
    return out

def dice(a,b):
    den=a.sum()+b.sum(); return 0.0 if den==0 else float(2*np.logical_and(a,b).sum()/den)

def angle_deg(a,b):
    if not (np.all(np.isfinite(a)) and np.all(np.isfinite(b))): return np.nan
    return float(np.degrees(np.arccos(np.clip(np.dot(a,b),-1,1))))

def center_at(h2,hh,t):
    tx,y=h2.h5_series(hh,HORIZON,'CoordCenterInertial'); return h2.interp_cols(tx,y[:,:3],[t])[0]

def mch_at(h2,hh,t):
    tx,y=h2.h5_series(hh,HORIZON,'ChristodoulouMass'); vals=y if y.ndim==1 else y[:,0]
    return float(np.interp(t,tx,vals))

def fvals(field,R,E,B,m):
    return R if field=='Rstar' else (m*m*E if field=='Estar' else m*m*B)

def plot_outputs(out,snap):
    fig=plt.figure(figsize=(11,6))
    for f in FIELDS:
        g=snap[snap.field==f].sort_values('offset_M')
        plt.plot(g.offset_M,g.largest_component_fraction_of_active,marker='o',label=f)
    plt.axhline(FRAC_GATE,ls='--',label='component gate=0.35')
    plt.xlabel('t - tC [M]'); plt.ylabel('largest component / active points')
    plt.title('BH-HF1R dense localization timeline'); plt.legend(); plt.tight_layout()
    fig.savefig(out/'BH_HF1R_localization_timeline.png',dpi=180); plt.close(fig)

    fig=plt.figure(figsize=(10,5.5)); ax=fig.add_subplot(111,projection='mollweide')
    for f in FIELDS:
        g=snap[(snap.field==f)&snap.candidate_snapshot].sort_values('offset_M')
        if len(g):
            x,y,z=g.centroid_x.to_numpy(),g.centroid_y.to_numpy(),g.centroid_z.to_numpy()
            ax.plot(np.arctan2(y,x),np.arcsin(np.clip(z,-1,1)),marker='o',label=f)
    ax.grid(True,alpha=.3); ax.set_title('BH-HF1R positive-component centroid tracks'); ax.legend(loc='lower left')
    fig.tight_layout(); fig.savefig(out/'BH_HF1R_centroid_tracks.png',dpi=180); plt.close(fig)

def run(args):
    out=Path(args.out); out.mkdir(parents=True,exist_ok=True)
    h2=imp(args.h2_script,'h2r2_frozen'); rng=np.random.default_rng(RNG_SEED)
    grid=fibonacci_sphere(GRID_N); gtree=cKDTree(grid); gnei,_=knn_graph(grid,GRID_KNN)
    rows=[]; comps={}
    with h5py.File(args.horizons,'r') as hh, h5py.File(args.dump,'r') as hd:
        idx=np.asarray(h2.tdm_index(hd,HORIZON,'DimlessRicciScalar.tdm')); times=idx[:,0].astype(float)
        tC=float(times[0]); sel=np.flatnonzero((times>=tC)&(times<=tC+WINDOW_END_M))
        if len(sel)<TRACK_N_MIN: raise RuntimeError('Too few dense snapshots')
        dts=np.diff(times[sel]); cadence=float(np.median(dts[dts>0])); maxgap=max(GAP_FLOOR,GAP_FACTOR*cadence)
        print(f'tC={tC:.12f} M; snapshots={len(sel)}; cadence={cadence:.6f} M; maxgap={maxgap:.6f} M')
        for k,i in enumerate(sel,1):
            t=float(times[i]); off=t-tC
            R,xR=h2.read_scalar_snapshot(hd,HORIZON,'DimlessRicciScalar.tdm',int(i))
            E,xE=h2.read_scalar_snapshot(hd,HORIZON,'WeylE_NN.tdm',int(i))
            B,xB=h2.read_scalar_snapshot(hd,HORIZON,'WeylB_NN.tdm',int(i))
            C,xC=h2.read_coord_snapshot(hd,HORIZON,int(i))
            if not (xR==xE==xB==xC): raise RuntimeError(f'Extent mismatch at snapshot {i}')
            c=center_at(h2,hh,t); mch=mch_at(h2,hh,t); n,good=unit_dirs(C,c); nei,ang=knn_graph(n,KNN)
            R,E,B=np.asarray(R)[good],np.asarray(E)[good],np.asarray(B)[good]
            for f in FIELDS:
                z=robust_z(fvals(f,R,E,B,mch)); m=morphology(z,n,nei,ang,rng)
                raw=mask_on_grid(n[m['_li']],grid,gtree); dm=dilate(raw,gnei)
                comps[(f,int(i))]={'dil':dm,'c':np.array([m['centroid_x'],m['centroid_y'],m['centroid_z']])}
                rows.append({'stage':STAGE,'sxs_id':SID,'lev':LEV,'horizon':HORIZON,
                             'snapshot_index':int(i),'field':f,'actual_time_M':t,'offset_M':off,
                             'n_surface_points':int(len(z)),'Mch':mch,
                             **{q:v for q,v in m.items() if not q.startswith('_')}})
            if k%10==0 or k==len(sel): print(f'processed {k}/{len(sel)}',flush=True)
    snap=pd.DataFrame(rows)
    pairs=[]
    for f in FIELDS:
        g=snap[(snap.field==f)&snap.candidate_snapshot].sort_values('actual_time_M')
        rr=list(g.to_dict('records'))
        for a,b in zip(rr[:-1],rr[1:]):
            ca,cb=comps[(f,int(a['snapshot_index']))],comps[(f,int(b['snapshot_index']))]
            dt=b['actual_time_M']-a['actual_time_M']; step=angle_deg(ca['c'],cb['c']); dd=dice(ca['dil'],cb['dil'])
            linked=bool(dt<=maxgap and step<=PAIR_STEP_MAX and dd>=PAIR_DICE_MIN)
            pairs.append({'field':f,'snapshot_i':int(a['snapshot_index']),'snapshot_j':int(b['snapshot_index']),
                          'offset_i_M':a['offset_M'],'offset_j_M':b['offset_M'],'dt_M':dt,
                          'centroid_step_deg':step,'dilated_dice':dd,'linked':linked})
    pair=pd.DataFrame(pairs)

    tracks=[]; track_json={}
    for f in FIELDS:
        g=snap[(snap.field==f)&snap.candidate_snapshot].sort_values('actual_time_M')
        rec=list(g.to_dict('records')); best=[]; cur=[]
        lookup={(int(r.snapshot_i),int(r.snapshot_j)):r for _,r in pair[pair.field==f].iterrows()}
        for r in rec:
            if not cur: cur=[r]; continue
            a=cur[-1]; e=lookup.get((int(a['snapshot_index']),int(r['snapshot_index'])))
            if e is not None and bool(e.linked): cur.append(r)
            else:
                if len(cur)>len(best): best=cur
                cur=[r]
        if len(cur)>len(best): best=cur
        edges=[]
        for a,b in zip(best[:-1],best[1:]):
            e=lookup.get((int(a['snapshot_index']),int(b['snapshot_index'])))
            if e is not None: edges.append(e)
        N=len(best); dur=0.0 if N<2 else best[-1]['actual_time_M']-best[0]['actual_time_M']
        medstep=np.nan if not edges else float(np.median([e.centroid_step_deg for e in edges]))
        meddice=np.nan if not edges else float(np.median([e.dilated_dice for e in edges]))
        passed=bool(N>=TRACK_N_MIN and dur>=TRACK_DUR_MIN and np.isfinite(medstep) and medstep<=TRACK_STEP_MED_MAX and np.isfinite(meddice) and meddice>=TRACK_DICE_MED_MIN)
        tracks.append({'field':f,'n_snapshots':N,'duration_M':dur,'median_centroid_step_deg':medstep,
                       'median_dilated_dice':meddice,'temporal_coherence_pass':passed})
        track_json[f]={'n_snapshots':N,'duration_M':dur,'median_centroid_step_deg':medstep,
                       'median_dilated_dice':meddice,'temporal_coherence_pass':passed,
                       'snapshot_indices':[int(r['snapshot_index']) for r in best],
                       'offsets_M':[float(r['offset_M']) for r in best]}
    track=pd.DataFrame(tracks)

    sets=[set(snap[(snap.field==f)&snap.candidate_snapshot].snapshot_index.astype(int)) for f in FIELDS]
    common=sorted(set.intersection(*sets)) if sets else []; cf=[]
    for i in common:
        def cen(f):
            r=snap[(snap.field==f)&(snap.snapshot_index==i)].iloc[0]
            return np.array([r.centroid_x,r.centroid_y,r.centroid_z]),float(r.offset_M)
        cr,off=cen('Rstar'); ce,_=cen('Estar'); cb,_=cen('Bstar')
        re=angle_deg(cr,ce); rb=angle_deg(cr,cb); dev=abs(180-rb)
        cf.append({'snapshot_index':i,'offset_M':off,'RE_angle_deg':re,'RB_angle_deg':rb,
                   'RB_antipodal_deviation_deg':dev,
                   'secondary_relation_pass':bool(re<=RE_MAX and dev<=RB_ANTI_DEV_MAX)})
    cf=pd.DataFrame(cf); secfrac=0.0 if not len(cf) else float(cf.secondary_relation_pass.mean())
    secpass=bool(len(cf)>=SECONDARY_N_MIN and secfrac>=SECONDARY_FRAC_MIN)

    nfields=int(track.temporal_coherence_pass.sum())
    cls='BH_HF1R_TEMPORALLY_COHERENT_TRANSIENT_MORPHOLOGY' if nfields>=PRIMARY_FIELDS_MIN else 'BH_HF1R_NO_TEMPORAL_COHERENCE_BY_FROZEN_TRACKING'

    snap.to_csv(out/'BH_HF1R_snapshot_metrics.csv',index=False)
    pair.to_csv(out/'BH_HF1R_pairwise_tracking.csv',index=False)
    track.to_csv(out/'BH_HF1R_track_summary.csv',index=False)
    cf.to_csv(out/'BH_HF1R_cross_field_geometry.csv',index=False)
    plot_outputs(out,snap)

    summary={'stage':STAGE,'classification':cls,'simulation':SID,'lev':LEV,'horizon':HORIZON,
             'common_horizon_scalar_birth_time_M':float(snap.actual_time_M.min()),
             'dense_window_M':[0.0,WINDOW_END_M],'native_snapshot_count':int(snap.snapshot_index.nunique()),
             'median_native_cadence_M':cadence,'max_link_gap_M':maxgap,
             'fields_passing_temporal_coherence':nfields,'required_fields_for_primary_pass':PRIMARY_FIELDS_MIN,
             'field_tracks':track_json,
             'secondary_cross_field_test':{'status':'SUPPORTED' if secpass else 'NOT_SUPPORTED',
                 'simultaneous_positive_snapshots':int(len(cf)),'fraction_relation_pass':secfrac,
                 'RE_angle_gate_deg':RE_MAX,'RB_antipodal_deviation_gate_deg':RB_ANTI_DEV_MAX,
                 'note':'Secondary pattern was motivated by sparse HF1 and does not enter primary classification.'},
             'methodological_status':'coordinate/gauge-dependent temporal morphology audit; not intrinsic/covariant, not QNM-subtracted',
             'forbidden_claims':FROZEN_PROTOCOL['forbidden_claims'],
             'next_if_pass':'BH-HF2 QNM/multipole closure; only then resolution/frame/slicing robustness.',
             'next_if_fail':'Stop this morphology route.'}
    (out/'BH_HF1R_summary.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2),encoding='utf-8')
    (out/'BH_HF1R_preregistered_protocol.json').write_text(json.dumps(FROZEN_PROTOCOL,ensure_ascii=False,indent=2),encoding='utf-8')
    print(json.dumps(summary,ensure_ascii=False,indent=2)); print(f'\nOutputs written to: {out.resolve()}')

def self_test():
    g=fibonacci_sphere(512); assert g.shape==(512,3); assert np.allclose(np.linalg.norm(g,axis=1),1,atol=1e-12)
    nei,_=knn_graph(g,8); tree=cKDTree(g); m=mask_on_grid(g[:20],g,tree)
    assert dice(dilate(m,nei),dilate(m,nei))>.99
    assert abs(angle_deg(np.array([1.,0,0]),np.array([0.,1,0]))-90)<1e-10
    print('BH-HF1R self-test: PASS')

def parser():
    p=argparse.ArgumentParser(); p.add_argument('--self-test',action='store_true')
    p.add_argument('--h2-script',default='/mnt/c/Users/Admin/gw150914_GW_XA4_4E2H2R2_surface_field_recurrence_nonreturn.py')
    p.add_argument('--horizons',default='/mnt/c/Users/Admin/Horizons.h5')
    p.add_argument('--dump',default='/mnt/c/Users/Admin/HorizonsDump.h5')
    p.add_argument('--out',default='/mnt/c/Users/Admin/bh_hf1r_output')
    return p

if __name__=='__main__':
    args=parser().parse_args()
    if args.self_test: self_test()
    else:
        for q in [args.h2_script,args.horizons,args.dump]:
            if not Path(q).exists(): raise SystemExit(f'Missing required file: {q}')
        run(args)
