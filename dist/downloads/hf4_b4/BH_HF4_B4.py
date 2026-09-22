#!/usr/bin/env python3
"""HF4-B4: audit same-run scalar areas and extend only the 43 new frames.

The previous 12 metrics/masks and 30 contiguous pair results are reused.
No evolution, fitting, synthetic re-tests, or B1/B2/B3 recalculation.
"""
from pathlib import Path
import csv, hashlib, importlib.util, json, re
import numpy as np

HERE=Path(__file__).resolve().parent
BASE=HERE.parent
INP=HERE/'input'
OUT=HERE/'results'
FIELDS=['repsi2','impsi2','abs_npsigma']
spec=importlib.util.spec_from_file_location('frozen_B3',BASE/'hf4_b3/BH_HF4_B3.py')
B3=importlib.util.module_from_spec(spec);spec.loader.exec_module(B3)


def dump(path,obj):
    path.write_text(json.dumps(obj,ensure_ascii=False,indent=2,allow_nan=False)+'\n')


def csvrows(path):
    with path.open() as f:return list(csv.DictReader(f))


def writecsv(path,rows):
    keys=list(dict.fromkeys(k for r in rows for k in r))
    with path.open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=keys);writer.writeheader();writer.writerows(rows)


def read_asc(name):
    p=INP/'diagnostics'/name;lines=p.read_text().splitlines()
    header=next(x for x in lines if x.startswith('# data columns:'))
    cols={m[2]:int(m[1])-1 for m in re.finditer(r'(\d+):([^\s]+)',header)}
    a=np.loadtxt(p,ndmin=2)
    assert len(set(a[:,0]))==len(a) and np.array_equal(a[:,0],np.arange(0,11265,32))
    return {'rows':{int(r[0]):r for r in a},'cols':cols,
            'run_id':next(x.split(': ',1)[1] for x in lines if x.startswith('# Run ID:'))}


def val(table,it,name):return float(table['rows'][it][table['cols'][name]])


def shape_gp(path,stride=1):
    """Read native six cubed-sphere patches; no time or angular interpolation.

    Full-resolution meshes are used for centroids. A subset of nodes (stride=2)
    is sufficient for the separate plain-color individual-horizon viewer.
    """
    text=path.read_text();points=[];triangles=[];specs=[]
    for block in text.split('### ')[1:]:
        nr=int(re.search(r'N_rho = (\d+)',block)[1]);ns=int(re.search(r'N_sigma = (\d+)',block)[1])
        rows=[line for line in block.splitlines() if line and not line.startswith('#') and len(line.split())==6]
        a=np.fromstring(' '.join(rows),sep=' ').reshape(nr,ns,6)
        idx=np.unique(np.r_[np.arange(0,nr,stride),nr-1]);jdx=np.unique(np.r_[np.arange(0,ns,stride),ns-1])
        sub=a[idx][:,jdx,3:6];offset=len(points);points.extend(sub.reshape(-1,3))
        h,w=sub.shape[:2];specs.append([h,w])
        for i in range(h-1):
            for j in range(w-1):
                k=offset+i*w+j;triangles.extend([[k,k+1,k+w+1],[k,k+w+1,k+w]])
    P=np.asarray(points);T=np.asarray(triangles,dtype=int)
    assert len(specs)==6 and np.isfinite(P).all()
    A=np.linalg.norm(np.cross(P[T[:,1]]-P[T[:,0]],P[T[:,2]]-P[T[:,0]]),axis=1)/2
    assert np.all(A>0)
    C=np.average(P[T].mean(axis=1),axis=0,weights=A)
    return P,T,C


def coordinate_weights(P,tris):
    a,b,c=P[tris[:,0]],P[tris[:,1]],P[tris[:,2]]
    area=np.linalg.norm(np.cross(b-a,c-a),axis=1)/2
    w=np.zeros(len(P))
    for k in range(3):np.add.at(w,tris[:,k],area/3)
    return w/w.sum()


def frame_axes(it):
    paths=[INP/f'shapes/h.t{it}.ah{h}.gp' for h in [1,2]]
    if not all(p.is_file() for p in paths):return None
    c1,c2=[shape_gp(p)[2] for p in paths]
    u=c2-c1;u[2]=0;u/=np.linalg.norm(u);w=np.array([0.,0.,1.])
    return np.stack([u,np.cross(w,u),w],axis=1)


def main():
    OUT.mkdir(exist_ok=True)
    # This small JSON is committed before new numerical work.
    definition={'stage':'BH-HF4-B4','input_run':'qc0_early_geometry_t22',
      'common_iterations':list(range(9536,11265,32)),'frozen_method':'HF4-B3 FD4 dyad metric, Fejer-I physical measure',
      'primary_quantile':.95,'sensitivity_quantiles':[.90,.98],
      'old_frames_policy':'reuse saved B3 q, physical weights, masks and localization rows',
      'new_frames_policy':'compute only missing 43 iterations',
      'area_control':'QLM slot 2 same iteration/time; AHFinderDirect area secondary cross-method control',
      'axis_extension':'native AH1/AH2 gp triangle-area centroids; overlap control against frozen B2 axes',
      'missing_individuals_policy':'no invented long axis; leave long-axis fractions empty, retain z-polar fractions',
      'area_comparison_not_simulation_error_estimate':True,'no_new_evolution':True}
    dump(OUT/'BH_HF4_B4_definition.json',definition)
    tables={key:read_asc(f'quasilocalmeasures-qlm_{key}..asc') for key in ['state','scalars','grid_int','grid_real']}
    assert len(set(t['run_id'] for t in tables.values()))==1
    par=(INP/'parameters/qc0-horizon-early-geometry-t22.par').read_text()
    assert re.search(r'QuasiLocalMeasures::surface_index\[2\]\s*=\s*2',par)
    assert re.search(r'which_surface_to_store_info\s*\[3\]\s*=\s*2',par)
    assert re.search(r'QuasiLocalMeasures::spatial_order\s*=\s*4',par)
    state=tables['state'];scalars=tables['scalars'];gr=tables['grid_real'];gi=tables['grid_int']
    for it in range(0,9536,32):assert val(state,it,'qlm_have_valid_data[2]')==0
    olddir=BASE/'hf4_b3/results';oldsum=json.loads((olddir/'BH_HF4_B3_summary.json').read_text())
    oldset=set(oldsum['common_iterations']);oldnpz=np.load(olddir/'BH_HF4_B3_intrinsic_geometry_and_masks.npz')
    b2=np.load(BASE/'hf4_b2/results/BH_HF4_B2_fields_and_masks.npz')
    db=json.loads((BASE/'hf4_b2/results/BH_HF4_B2_view_data.json').read_text());tris=np.asarray(db['triangles'])
    oldframes={f['iteration']:f for f in db['frames'] if 'fields' in f}
    oldgeo={int(r['iteration']):r for r in csvrows(olddir/'BH_HF4_B3_geometry.csv')}
    oldloc=[r for r in csvrows(olddir/'BH_HF4_B3_localization.csv') if r['measure']=='fejer']
    oldpairs={(int(r['iteration_from']),int(r['iteration_to']),r['field']):r for r in csvrows(olddir/'BH_HF4_B3_temporal_pairs.csv') if r['contiguous']=='True'}
    bh=np.loadtxt(INP/'diagnostics/BH_diagnostics.ah3.gp');bh={int(r[0]):r for r in bh}
    geo=[];loc=[];areas=[];axischecks=[];pairs=[];saved={};done=[]
    previous={};newcount=0;pair_reused=0
    for it in range(9536,11265,32):
        time=it/512;tau=(it-9536)/512
        assert val(state,it,'qlm_have_valid_data[2]')==1 and val(state,it,'qlm_calc_error[2]')==0
        assert val(scalars,it,'qlm_time[2]')==time
        # First valid sample has iteration-history sentinel -1 in this export.
        counter=val(state,it,'qlm_iteration[2]');assert counter==it or (it==9536 and counter==-1)
        assert val(gi,it,'qlm_ntheta[2]')==39 and val(gi,it,'qlm_nphi[2]')==76
        assert val(gi,it,'qlm_nghoststheta[2]')==val(gi,it,'qlm_nghostsphi[2]')==2
        assert abs(val(gr,it,'qlm_delta_theta[2]')-np.pi/35)<1e-13
        assert abs(val(gr,it,'qlm_delta_phi[2]')-2*np.pi/72)<1e-13
        newaxis=frame_axes(it)
        if it in oldset:
            P=b2[f'points_{it}'];wc=b2[f'weights_{it}'];F={f:b2[f'field_{f}_{it}'] for f in FIELDS}
            q=oldnpz[f'q_{it}'];rho=oldnpz[f'physical_density_{it}'];wp=oldnpz[f'physical_weights_{it}']
            area=float(oldgeo[it]['area_physical_fejer']);masks={f:oldnpz[f'mask_{f}_{it}'] for f in FIELDS}
            oldaxis=np.asarray(oldframes[it]['axes']);center=np.average(P,axis=0,weights=wc)
            if newaxis is not None:
                u=oldaxis[:,0];ang=np.arccos(np.clip(newaxis[:,0]@u/np.linalg.norm(u),-1,1))*180/np.pi
                a=(P-center)@oldaxis;a/=np.max(np.abs(a),axis=0)
                b=(P-center)@newaxis;b/=np.max(np.abs(b),axis=0)
                differing=int(np.count_nonzero((np.abs(a)>=.65)!=(np.abs(b)>=.65)))
                axischecks.append({'iteration':it,'axis_angle_degrees':float(ang),'differing_region_memberships':differing})
            axis=oldaxis
            for r in oldloc:
                if int(r['iteration'])==it:loc.append({**r,'origin':'reused_B3','axis_source':'frozen_B2'})
            g={k:float(v) for k,v in oldgeo[it].items()};g['iteration']=it
        else:
            shaped,m,l,n,F=B3.vtk(INP/f'surfaces/surface03_{it:08}.vtk');P=shaped.reshape(-1,3)
            E=B3.derivatives(shaped);r=B3.metric(E,m,l,n)
            assert r['m0_abs_max']==0 and r['tangent_residual'].max()<1e-10
            assert np.isfinite(r['q']).all() and r['q_eigenvalues'].min()>0
            assert np.all(r['spacetime_metric_eigenvalues'][...,0]<0) and np.all(r['spacetime_metric_eigenvalues'][...,1:]>0)
            q=r['q'];rho=r['rho'];weights=B3.quadrature(rho,'fejer').ravel();area=float(weights.sum());wp=weights/area
            wc=coordinate_weights(P,tris);center=np.average(P,axis=0,weights=wc);axis=newaxis
            theta=(np.arange(35)+.5)*np.pi/35;phi=np.arange(72)*2*np.pi/72;T,H=np.meshgrid(theta,phi,indexing='ij')
            radial=np.stack([np.sin(T)*np.cos(H),np.sin(T)*np.sin(H),np.cos(T)],axis=-1)
            origin=np.array([val(gr,it,f'qlm_origin_{c}[2]') for c in ['x','y','z']])
            grid_error=float(np.max(np.abs(shaped-F['shape'].reshape(35,72,1)*radial-origin)))
            assert grid_error<1e-10
            zcoord=P[:,2]-center[2];polar=np.abs(zcoord/np.max(np.abs(zcoord)))>=.65
            if axis is not None:
                local=(P-center)@axis;local/=np.max(np.abs(local),axis=0);long=np.abs(local[:,0])>=.65
            masks={}
            for f in FIELDS:
                values=F[f];assert np.isfinite(values).all();median=B3.wquantile(values,wp,.5);deviation=np.abs(values-median)
                for quant in [.90,.95,.98]:
                    cutoff=B3.wquantile(deviation,wp,quant);mask=deviation>=cutoff;mass=float(wp[mask].sum())
                    loc.append({'iteration':it,'t_M':time,'tau_M':tau,'field':f,'measure':'fejer','quantile':quant,
                        'median':median,'deviation_cutoff':cutoff,'vertices':int(mask.sum()),'selected_measure_fraction':mass,
                        'long_axis_fraction':float(wp[mask&long].sum()/mass) if axis is not None else '',
                        'polar_fraction':float(wp[mask&polar].sum()/mass),'origin':'new_B4',
                        'axis_source':'native_AH1_AH2_gp' if axis is not None else 'missing_AH1_AH2_z_only'})
                    if quant==.95:masks[f]=mask
            g={'iteration':it,'t_M':time,'tau_M':tau,'area_physical_fejer':area,
               'area_physical_midpoint':float(B3.quadrature(rho,'midpoint').sum()),
               'area_coordinate_fd4_fejer':float(B3.quadrature(r['rho_e'],'fejer').sum()),
               'max_tangent_residual':float(r['tangent_residual'].max()),'min_q_eigenvalue':float(r['q_eigenvalues'].min()),
               'full_tetrad_q_difference_max':float(r['full_tetrad_q_relative_difference'].max()),
               'angular_grid_radial_fit_max_error':grid_error}
            newcount+=1
        g.update(origin='reused_B3' if it in oldset else 'new_B4',long_axis_available=axis is not None,
                 coordinate_span_max_over_min=float(np.ptp(P,axis=0).max()/np.ptp(P,axis=0).min()) if it not in oldset else oldframes[it]['B1_span_ratio'])
        geo.append(g)
        qal=val(scalars,it,'qlm_area[2]');assert np.isfinite(qal) and qal>0
        assert abs(bh[it][1]-time)<=.00050001
        areas.append({'iteration':it,'time_M':time,'tau_M':tau,'area_recovered':area,'area_QLM':qal,
          'relative_recovered_minus_QLM':(area-qal)/qal,'area_AHFinderDirect':float(bh[it][25]),
          'relative_AHFinderDirect_minus_QLM':float(bh[it][25]/qal-1),'QLM_history_iteration':int(counter),
          'QLM_valid':True,'QLM_calc_error':0,'origin':g['origin']})
        saved[f'points_{it}']=P;saved[f'coordinate_weights_{it}']=wc;saved[f'q_{it}']=q
        saved[f'physical_weights_{it}']=wp;saved[f'physical_density_{it}']=rho
        if axis is not None:saved[f'axes_{it}']=axis
        for f in FIELDS:
            saved[f'field_{f}_{it}']=F[f];saved[f'mask_{f}_{it}']=masks[f]
            if f in previous:
                pit,pm,pw=previous[f];key=(pit,it,f)
                if key in oldpairs:
                    row={**oldpairs[key],'origin':'reused_B3'};pair_reused+=1
                else:
                    W=(pw+wp)/2;row={'field':f,'iteration_from':pit,'iteration_to':it,'gap_M':(it-pit)/512,
                        'contiguous':True,'Dice_physical_measure':float(2*W[pm&masks[f]].sum()/(W[pm].sum()+W[masks[f]].sum())),
                        'origin':'new_B4'}
                pairs.append(row)
            previous[f]=(it,masks[f],wp)
        done.append(it)
        if len(done)%10==0:print('Verified common frames:',len(done),'new metrics:',newcount,flush=True)
    assert newcount==43 and pair_reused==30
    np.savez_compressed(OUT/'BH_HF4_B4_fields_geometry_masks.npz',**saved)
    writecsv(OUT/'BH_HF4_B4_geometry.csv',geo);writecsv(OUT/'BH_HF4_B4_area_comparison.csv',areas)
    writecsv(OUT/'BH_HF4_B4_localization.csv',loc);writecsv(OUT/'BH_HF4_B4_temporal_pairs.csv',pairs)
    dump(OUT/'BH_HF4_B4_axis_checks.json',axischecks)
    def span(rows,key):
        values=[float(r[key]) for r in rows if r.get(key) not in ['',None]]
        return [min(values),max(values)]
    fields={}
    for f in FIELDS:
        r=[x for x in loc if x['field']==f and float(x['quantile'])==.95]
        pp=[x for x in pairs if x['field']==f]
        fields[f]={'primary_mask_long_axis_fraction_range':span(r,'long_axis_fraction'),
            'primary_mask_polar_fraction_range':span(r,'polar_fraction'),
            'selected_area_fraction_range':span(r,'selected_measure_fraction'),
            'adjacent_mask_Dice_range':span(pp,'Dice_physical_measure'),
            'long_axis_frames':sum(x['long_axis_fraction']!='' for x in r),'polar_frames':len(r),
            'quantile_sensitivity':{str(q):{'long_axis_fraction_range':span([x for x in loc if x['field']==f and float(x['quantile'])==q],'long_axis_fraction'),
                'polar_fraction_range':span([x for x in loc if x['field']==f and float(x['quantile'])==q],'polar_fraction')} for q in [.9,.95,.98]}}
    summary={'stage':'BH-HF4-B4','classification':'BH_HF4_B4_COMPLETE_AH3_SERIES_EARLY_AREA_AGREEMENT_CONFIRMED_LOCALIZATION_EXTENDED',
      'common_frame_count':55,'old_B3_frames_reused':12,'new_frames_computed':43,'old_contiguous_pair_results_reused':30,
      'new_contiguous_pair_results':len(pairs)-pair_reused,'same_run_id':scalars['run_id'],
      'QLM_common_slot':2,'QLM_early_validity_checks_passed':True,
      'first_valid_history_counter_sentinel':{'iteration':9536,'qlm_iteration':-1,'qlm_time':18.625,'have_valid_data':1,'calc_error':0,'interpretation':'Recorded explicitly; same-time area corroborated. Not used as a dynamical time derivative.'},
      'first_common_time_M':18.625,'preceding_valid_sample_without_common_time_M':18.5625,
      'last_common_time_M':22,'common_step_M':.0625,'common_time_span_M':3.375,
      'physical_area_first_last':[areas[0]['area_recovered'],areas[-1]['area_recovered']],
      'relative_area_increase':areas[-1]['area_recovered']/areas[0]['area_recovered']-1,
      'area_monotonically_increasing_at_all_sampled_frames':bool(np.all(np.diff([a['area_recovered'] for a in areas])>0)),
      'max_abs_relative_area_difference_QLM':max(abs(a['relative_recovered_minus_QLM']) for a in areas),
      'AHFinderDirect_vs_QLM_relative_range':span(areas,'relative_AHFinderDirect_minus_QLM'),
      'max_abs_relative_area_difference_AHFinderDirect_QLM':max(abs(a['relative_AHFinderDirect_minus_QLM']) for a in areas),
      'max_tangent_residual_new_frames':max(g['max_tangent_residual'] for g in geo if g['origin']=='new_B4'),
      'max_axis_overlap_difference_degrees':max(a['axis_angle_degrees'] for a in axischecks),
      'axis_overlap_differing_region_memberships':sum(a['differing_region_memberships'] for a in axischecks),
      'long_axis_unavailable_iterations':[g['iteration'] for g in geo if not g['long_axis_available']],
      'fields':fields,'new_evolution_run':False,'previous_results_recomputed':False,
      'simulation_convergence_established':False,'all_tetrad_or_slicing_independence_established':False,
      'material_transport_established':False,'HF1_identification_established':False,'new_physics_established':False}
    dump(OUT/'BH_HF4_B4_summary.json',summary)
    print(json.dumps(summary,ensure_ascii=False,indent=2))


if __name__=='__main__':main()
