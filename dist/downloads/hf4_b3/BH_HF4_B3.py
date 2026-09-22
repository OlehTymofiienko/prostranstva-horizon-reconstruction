#!/usr/bin/env python3
"""HF4-B3: recover intrinsic surface measure from the stored NP dyad.

Read-only analysis. No evolution, ringdown fitting, or B1/B2 result rerun.
"""
from pathlib import Path
import argparse
import csv
import hashlib
import json
import re
import numpy as np

HERE = Path(__file__).resolve().parent
BASE = HERE.parent


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def vtk(path):
    s=path.read_text()
    h=re.search(r'POINTS\s+(\d+)\s+\w+\s*\n',s)
    P=np.fromstring(s[h.end():s.index('POLYGONS',h.end())],sep=' ').reshape(int(h[1]),3)
    hh=list(re.finditer(r'SCALARS\s+(\S+)\s+\w+\s+1\s*\nLOOKUP_TABLE\s+default\s*\n',s))
    F={h[1]:np.fromstring(s[h.end():hh[k+1].start() if k+1<len(hh) else len(s)],sep=' ') for k,h in enumerate(hh)}
    assert P.shape==(2520,3) and all(v.shape==(2520,) for v in F.values())
    F['abs_npsigma']=np.hypot(F['renpsigma'],F['imnpsigma'])
    m=np.stack([F['rem'+str(k)]+1j*F['imm'+str(k)] for k in range(4)],axis=-1)
    l=np.stack([F['l'+str(k)] for k in range(4)],axis=-1)
    n=np.stack([F['n'+str(k)] for k in range(4)],axis=-1)
    return P.reshape(35,72,3),m.reshape(35,72,4),l.reshape(35,72,4),n.reshape(35,72,4),F


def derivatives(P,order=4):
    """Same centered FD stencil as QLM, Cartesian components have scalar parity.

    Across a pole: theta -> -theta, phi -> phi+pi; Cartesian x/y/z
    themselves have positive scalar parity on this coordinate covering.
    """
    nt,np_=P.shape[:2]; g=2 if order==4 else 1
    ext=np.concatenate([np.roll(P[:g][::-1],np_//2,axis=1),P,
                        np.roll(P[-g:][::-1],np_//2,axis=1)],axis=0)
    dt=np.pi/nt;dp=2*np.pi/np_
    if order==4:
        Et=(-ext[4:]+8*ext[3:-1]-8*ext[1:-3]+ext[:-4])/(12*dt)
        Ep=(-np.roll(P,-2,axis=1)+8*np.roll(P,-1,axis=1)-8*np.roll(P,1,axis=1)+np.roll(P,2,axis=1))/(12*dp)
    else:
        Et=(ext[2:]-ext[:-2])/(2*dt)
        Ep=(np.roll(P,-1,axis=1)-np.roll(P,1,axis=1))/(2*dp)
    return np.stack([Et,Ep],axis=-1)


def metric(E,m,l=None,n=None):
    """E[..., i,A]; m Cartesian contravariant. q^AB=2 Re(m^A conj(m^B))."""
    eu=np.einsum('...iA,...iB->...AB',E,E)
    rhs=np.einsum('...iA,...i->...A',E,m[...,1:])
    mc=np.linalg.solve(eu,rhs[...,None])[...,0]
    mproject=np.einsum('...iA,...A->...i',E,mc)
    residual=np.linalg.norm(mproject-m[...,1:],axis=-1)/np.linalg.norm(m[...,1:],axis=-1)
    qi=2*np.real(mc[..., :,None]*mc.conj()[...,None,:])
    q=np.linalg.inv(qi)
    rho=np.sqrt(np.linalg.det(q));rho_e=np.sqrt(np.linalg.det(eu))
    out={'q':q,'rho':rho,'rho_e':rho_e,'tangent_residual':residual,
         'q_eigenvalues':np.linalg.eigvalsh(q),'m0_abs_max':float(np.max(np.abs(m[...,0])))}
    if l is not None:
        gi=-(l[..., :,None]*n[...,None,:]+n[..., :,None]*l[...,None,:])+2*np.real(m[..., :,None]*m.conj()[...,None,:])
        g=np.linalg.inv(gi)
        qp=np.einsum('...iA,...ij,...jB->...AB',E,g[...,1:,1:],E)
        out['full_tetrad_q_relative_difference']=np.linalg.norm(q-qp,axis=(-2,-1))/np.linalg.norm(qp,axis=(-2,-1))
        out['spacetime_metric_eigenvalues']=np.linalg.eigvalsh(g)
    return out


def quadrature(rho,kind='midpoint'):
    nt,np_=rho.shape;dt=np.pi/nt;dp=2*np.pi/np_
    if kind=='midpoint':return rho*dt*dp
    # Fejer-I weights: integrate in x=cos(theta), i.e. rho/sin(theta).
    th=(np.arange(nt)+.5)*dt
    fw=np.ones(nt)
    for k in range(1,nt//2+1):fw-=2*np.cos(2*k*th)/(4*k*k-1)
    fw*=2/nt
    assert np.all(fw>0) and abs(fw.sum()-2)<1e-13
    return rho*(fw/np.sin(th))[:,None]*dp


def checkpoint(root,index):
    meta=json.loads((root/f'hdf_metadata/hdf_{index:03}.json').read_text())
    arrays=np.load(root/f'arrays/hdf_{index:03}.npz')
    keys={a['dataset'].split()[0].lower():a['key'] for a in meta['arrays'] if ' tl=0' in a['dataset']}
    def get(name):
        a=arrays[keys['quasilocalmeasures::qlm_'+name+'[2]']]
        if a.dtype.names:a=a['real']+1j*a['imag']
        return a
    def scalar(name):return get(name).item()
    gt=int(scalar('nghoststheta'));gp=int(scalar('nghostsphi'))
    nt=int(scalar('ntheta'));np_=int(scalar('nphi'))
    def interior(name):return get(name)[:np_,:nt].T[gt:nt-gt,gp:np_-gp]
    P=np.stack([interior(c) for c in ['x','y','z']],axis=-1)
    m=np.stack([interior('m'+str(k)) for k in range(4)],axis=-1)
    l=np.stack([interior('l'+str(k)) for k in range(4)],axis=-1)
    n=np.stack([interior('n'+str(k)) for k in range(4)],axis=-1)
    assert scalar('have_valid_data')==1
    assert abs(scalar('delta_theta')-np.pi/P.shape[0])<1e-14
    assert abs(scalar('delta_phi')-2*np.pi/P.shape[1])<1e-14
    # Exact coordinate ghosts are available in checkpoints: compare our pole
    # reconstruction with them, independently of the scalar area comparison.
    full=np.stack([get(c).T for c in ['x','y','z']],axis=-1)
    Et=(-full[gt+2:nt-gt+2,gp:np_-gp]+8*full[gt+1:nt-gt+1,gp:np_-gp]-8*full[gt-1:nt-gt-1,gp:np_-gp]+full[gt-2:nt-gt-2,gp:np_-gp])/(12*scalar('delta_theta'))
    Ep=(-full[gt:nt-gt,gp+2:np_-gp+2]+8*full[gt:nt-gt,gp+1:np_-gp+1]-8*full[gt:nt-gt,gp-1:np_-gp-1]+full[gt:nt-gt,gp-2:np_-gp-2])/(12*scalar('delta_phi'))
    Eref=np.stack([Et,Ep],axis=-1)
    return P,m,l,n,{'time_M':scalar('time'),'iteration':int(scalar('iteration')),
                   'area_QLM':scalar('area'),'ghost_derivative_difference':float(np.max(np.abs(Eref-derivatives(P)))),
                   'source':meta['path'],'npz_sha256':sha(root/f'arrays/hdf_{index:03}.npz')}


def controls(root,out):
    rows=[]
    for idx in [0,2,4,6,8]:
        P,m,l,n,row=checkpoint(root,idx)
        r=metric(derivatives(P),m,l,n)
        for k in ['midpoint','fejer']:
            row['area_'+k]=float(quadrature(r['rho'],k).sum())
            row['relative_difference_'+k]=row['area_'+k]/row['area_QLM']-1
        row['max_tangent_residual']=float(r['tangent_residual'].max())
        row['min_q_eigenvalue']=float(r['q_eigenvalues'].min())
        row['full_tetrad_q_difference_max']=float(r['full_tetrad_q_relative_difference'].max())
        row['m0_abs_max']=r['m0_abs_max'];rows.append(row)
    out.mkdir(parents=True,exist_ok=True)
    (out/'BH_HF4_B3_late_area_controls.json').write_text(json.dumps(rows,indent=2,allow_nan=False)+'\n')
    print(json.dumps(rows,indent=2))
    return rows


def write_csv(path,rows):
    with path.open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)


def wquantile(v,w,q):
    ix=np.argsort(v,kind='stable');c=np.cumsum(w[ix])
    return float(v[ix[min(np.searchsorted(c,q*c[-1]),len(ix)-1)]])


def synthetic_checks():
    nt,np_=35,72;th=(np.arange(nt)+.5)*np.pi/nt;ph=np.arange(np_)*2*np.pi/np_
    T,H=np.meshgrid(th,ph,indexing='ij')
    R=1.7;P=R*np.stack([np.sin(T)*np.cos(H),np.sin(T)*np.sin(H),np.cos(T)],axis=-1)
    Et=R*np.stack([np.cos(T)*np.cos(H),np.cos(T)*np.sin(H),-np.sin(T)],axis=-1)
    Ep=R*np.stack([-np.sin(T)*np.sin(H),np.sin(T)*np.cos(H),np.zeros_like(T)],axis=-1)
    E=np.stack([Et,Ep],axis=-1);G=np.diag([1.3,2.1,3.2])
    ip=lambda a,b:np.einsum('...i,ij,...j->...',a,G,b)
    u=Ep/np.sqrt(ip(Ep,Ep))[...,None]
    v=Et-ip(Et,u)[...,None]*u;v/=np.sqrt(ip(v,v))[...,None]
    m=np.zeros(P.shape[:-1]+(4,),complex);m[...,1:]=(u+1j*v)/np.sqrt(2)
    qdirect=np.einsum('...iA,ij,...jB->...AB',E,G,E)
    q=metric(E,m)['q']
    phase=.37+.25*np.cos(T)+.17*np.sin(H)
    qrot=metric(E,m*np.exp(1j*phase)[...,None])['q']
    qerr=float(np.max(np.abs(q-qdirect))/np.max(np.abs(qdirect)))
    perr=float(np.max(np.abs(qrot-q))/np.max(np.abs(q)))
    weights=quadrature(np.sin(T),'fejer')[:,0]/(2*np.pi/np_)
    errs=[]
    for k in range(nt):
        exact=0 if k%2 else 2/(1-k*k)
        errs.append(abs(float(weights@np.cos(k*th))-exact))
    serr=abs(float(quadrature(R*R*np.sin(T),'fejer').sum())/(4*np.pi*R*R)-1)
    assert qerr<1e-13 and perr<1e-13 and max(errs)<1e-13 and serr<1e-13
    return {'anisotropic_known_metric_relative_error':qerr,'local_dyad_phase_metric_relative_error':perr,
            'fejer_polynomial_moments_up_to_degree':nt-1,'fejer_moment_error_max':max(errs),
            'analytic_round_sphere_area_relative_error':serr,'tests_do_not_establish_simulation_convergence':True}


def analyze(input_dir,b2_dir,out,late):
    out.mkdir(parents=True,exist_ok=True)
    definition={
        'stage':'BH-HF4-B3','primary_measure':'dyad metric, FD4 tangents matching QLM, Fejer-I in cos(theta) and periodic phi',
        'grid':'theta=(i+0.5)*pi/35; phi=j*2*pi/72; no measured polar vertices',
        'primary_quantile':.95,'sensitivity_quantiles':[.90,.98],
        'deviation':'absolute deviation from median under each stated measure',
        'frozen_reference':'B2 saved fields, primary masks, coordinate-area weights, coordinate centers and axes',
        'auxiliary_measure_1':'FD4 physical density with midpoint integration',
        'auxiliary_measure_2':'saved B2 triangle weights times pointwise physical/Euclidean FD4 density ratio',
        'spatial_order_2_role':'deliberate mismatch diagnostic, not another simulation resolution',
        'early_independent_QLM_area_in_received_VTK_package':False,
        'late_control_role':'validate reconstruction against separately stored QLM scalar; no ringdown recalculation',
        'statistical_significance_test':False,'new_evolution':False,
    }
    (out/'BH_HF4_B3_definition.json').write_text(json.dumps(definition,indent=2)+'\n')
    checks=synthetic_checks()
    db=json.loads((b2_dir/'BH_HF4_B2_view_data.json').read_text())
    frozen=np.load(b2_dir/'BH_HF4_B2_fields_and_masks.npz')
    b2summary=json.loads((b2_dir/'BH_HF4_B2_summary.json').read_text())
    hashes={a['file']:a['sha256'] for a in b2summary['input_files']}
    geo=[];rows=[];comparisons=[];pairs=[];saved={};previous={};common=[]
    for frame in db['frames']:
        if 'fields' not in frame:continue
        it=frame['iteration'];common.append(it);path=input_dir/f'surface03_{it:08}.vtk'
        assert sha(path)==hashes[path.name]
        P,m,l,n,F=vtk(path);E=derivatives(P);r=metric(E,m,l,n)
        assert np.array_equal(P.reshape(-1,3),frozen[f'points_{it}'])
        assert r['m0_abs_max']==0 and r['tangent_residual'].max()<1e-10
        assert np.isfinite(r['q']).all() and r['q_eigenvalues'].min()>0
        assert np.all(r['spacetime_metric_eigenvalues'][...,0]<0) and np.all(r['spacetime_metric_eigenvalues'][...,1:]>0)
        nt,np_=P.shape[:2];theta=(np.arange(nt)+.5)*np.pi/nt;phi=np.arange(np_)*2*np.pi/np_
        T,H=np.meshgrid(theta,phi,indexing='ij')
        nhat=np.stack([np.sin(T)*np.cos(H),np.sin(T)*np.sin(H),np.cos(T)],axis=-1)
        origins=P-F['shape'].reshape(nt,np_,1)*nhat
        orspread=float(np.max(np.abs(origins-origins.mean(axis=(0,1)))))
        assert orspread<1e-10
        wc=frozen[f'weights_{it}']
        center=np.average(P.reshape(-1,3),axis=0,weights=wc)
        # B2 axes are rounded in its display JSON (absolute error <=5e-9).
        # Keep those axes fixed and check that no vertex lies near region edges.
        local=(P.reshape(-1,3)-center)@np.asarray(frame['axes'])
        local/=np.max(np.abs(local),axis=0)
        region_margin=float(np.min(np.abs(np.abs(local)-.65)))
        assert region_margin>1e-7
        weights={k:quadrature(r['rho'],k).ravel() for k in ['fejer','midpoint']}
        ratio=(r['rho']/r['rho_e']).ravel()
        weights['matched_triangle_support']=wc*ratio
        area=float(weights['fejer'].sum());area_e=float(quadrature(r['rho_e'],'fejer').sum())
        r2=metric(derivatives(P,2),m)
        rotated=metric(E,m*np.exp(1j*(.37+.25*np.cos(T)+.17*np.sin(H)))[...,None])
        A=np.exp(.2*np.cos(T))
        boosted=metric(E,m,l*A[...,None],n/A[...,None])
        g=dict(iteration=it,t_M=frame['t_M'],tau_M=frame['tau_M'],area_physical_fejer=area,
               area_physical_midpoint=float(weights['midpoint'].sum()),area_coordinate_fd4_fejer=area_e,
               physical_to_coordinate_area=area/area_e,physical_to_coordinate_density_min=float(ratio.min()),
               physical_to_coordinate_density_max=float(ratio.max()),max_tangent_residual=float(r['tangent_residual'].max()),
               min_q_eigenvalue=float(r['q_eigenvalues'].min()),full_tetrad_q_difference_max=float(r['full_tetrad_q_relative_difference'].max()),
               local_phase_density_relative_difference=float(np.max(np.abs(rotated['rho']/r['rho']-1))),
               local_boost_full_tetrad_q_difference_max=float(boosted['full_tetrad_q_relative_difference'].max()),
               angular_grid_radial_fit_max_error=orspread,region_boundary_margin=region_margin,
               second_order_mismatch_area_relative=float(quadrature(r2['rho'],'fejer').sum()/area-1),
               second_order_mismatch_tangent_residual=float(r2['tangent_residual'].max()))
        geo.append(g)
        wp=weights['fejer']/area
        saved[f'q_{it}']=r['q'];saved[f'physical_weights_{it}']=wp
        saved[f'physical_density_{it}']=r['rho'];saved[f'physical_coordinate_density_ratio_{it}']=ratio
        for field in definition.get('fields',['repsi2','impsi2','abs_npsigma']):
            values=F[field];assert np.array_equal(values,frozen[f'field_{field}_{it}'])
            masks={}
            for kind,unscaled in weights.items():
                W=unscaled/unscaled.sum();med=wquantile(values,W,.5);D=np.abs(values-med)
                for quant in [.90,.95,.98]:
                    cutoff=wquantile(D,W,quant);mask=D>=cutoff;mass=W[mask].sum()
                    row={'iteration':it,'t_M':frame['t_M'],'tau_M':frame['tau_M'],'field':field,'measure':kind,
                         'quantile':quant,'median':med,'deviation_cutoff':cutoff,'vertices':int(mask.sum()),
                         'selected_measure_fraction':float(mass),
                         'long_axis_fraction':float(W[mask&(np.abs(local[:,0])>=.65)].sum()/mass),
                         'polar_fraction':float(W[mask&(np.abs(local[:,2])>=.65)].sum()/mass)}
                    rows.append(row)
                    if quant==.95:masks[kind]=mask
            mask=masks['fejer'];old=frozen[f'mask_{field}_{it}']
            def dice(a,b,W):return float(2*W[a&b].sum()/(W[a].sum()+W[b].sum()))
            comparisons.append({'iteration':it,'t_M':frame['t_M'],'tau_M':frame['tau_M'],'field':field,
              'B2_B3_Dice_physical_measure':dice(old,mask,wp),'B2_B3_Dice_coordinate_measure':dice(old,mask,wc),
              'B2_mask_physical_area_fraction':float(wp[old].sum()),
              'B3_fejer_midpoint_mask_Dice_physical':dice(mask,masks['midpoint'],wp),
              'B3_fejer_matched_support_mask_Dice_physical':dice(mask,masks['matched_triangle_support'],wp)})
            saved[f'mask_{field}_{it}']=mask
            if field in previous:
                pit,pm,pw=previous[field];avg=(pw+wp)/2
                pairs.append({'field':field,'iteration_from':pit,'iteration_to':it,'gap_M':(it-pit)/512,
                              'contiguous':it-pit==32,'Dice_physical_measure':dice(pm,mask,avg)})
            previous[field]=(it,mask,wp)
    write_csv(out/'BH_HF4_B3_geometry.csv',geo);write_csv(out/'BH_HF4_B3_localization.csv',rows)
    write_csv(out/'BH_HF4_B3_mask_comparison.csv',comparisons);write_csv(out/'BH_HF4_B3_temporal_pairs.csv',pairs)
    np.savez_compressed(out/'BH_HF4_B3_intrinsic_geometry_and_masks.npz',**saved)
    summaries={}
    for field in ['repsi2','impsi2','abs_npsigma']:
        rr=[x for x in rows if x['field']==field and x['measure']=='fejer' and x['quantile']==.95]
        cc=[x for x in comparisons if x['field']==field];tt=[x for x in pairs if x['field']==field and x['contiguous']]
        ran=lambda rs,key:[min(x[key] for x in rs),max(x[key] for x in rs)]
        summaries[field]={'long_axis_fraction_range':ran(rr,'long_axis_fraction'),'polar_fraction_range':ran(rr,'polar_fraction'),
                         'B2_B3_Dice_physical_range':ran(cc,'B2_B3_Dice_physical_measure'),
                         'fejer_midpoint_mask_Dice_range':ran(cc,'B3_fejer_midpoint_mask_Dice_physical'),
                         'fejer_matched_support_mask_Dice_range':ran(cc,'B3_fejer_matched_support_mask_Dice_physical'),
                         'adjacent_mask_Dice_range':ran(tt,'Dice_physical_measure')}
    summary={'stage':'BH-HF4-B3','classification':'BH_HF4_B3_INTRINSIC_MEASURE_RECOVERED_LATE_AREA_CONTROL_PASSED_LOCALIZATION_SUPPORTED_EARLY_AREA_CONTROL_PENDING',
             'common_iterations':common,'definition':definition,'fields':summaries,'checks':checks,
             'late_area_control_max_relative_difference':max(abs(x['relative_difference_fejer']) for x in late),
             'max_early_tangent_residual':max(x['max_tangent_residual'] for x in geo),
             'min_early_q_eigenvalue':min(x['min_q_eigenvalue'] for x in geo),
             'max_local_phase_density_relative_difference':max(x['local_phase_density_relative_difference'] for x in geo),
             'early_independent_area_control_complete':False,'simulation_convergence_established':False,
             'HF1_identification_established':False,'B2_results_recomputed':False,
             'frozen_B2_npz_sha256':sha(b2_dir/'BH_HF4_B2_fields_and_masks.npz'),
             'source_hashes':[{'file':f'surface03_{it:08}.vtk','sha256':hashes[f'surface03_{it:08}.vtk']} for it in common]}
    (out/'BH_HF4_B3_summary.json').write_text(json.dumps(summary,indent=2,allow_nan=False)+'\n')
    print(json.dumps(summaries,indent=2));return summary


if __name__=='__main__':
    ap=argparse.ArgumentParser()
    ap.add_argument('--checkpoint-input',type=Path,default=BASE/'hf4_a3_review/input')
    ap.add_argument('--output',type=Path,default=HERE/'results')
    ap.add_argument('--input',type=Path,default=BASE/'hf4_latest/full')
    ap.add_argument('--b2',type=Path,default=BASE/'hf4_b2/results')
    args=ap.parse_args()
    cp=args.output/'BH_HF4_B3_late_area_controls.json'
    late=json.loads(cp.read_text()) if cp.exists() else controls(args.checkpoint_input,args.output)
    analyze(args.input,args.b2,args.output,late)
