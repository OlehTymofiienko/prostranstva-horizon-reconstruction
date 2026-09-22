#!/usr/bin/env python3
"""Reconstruct missing deterministic masks; compare to saved HF1R/HF1 metrics.
No permutations, QNM fits, evolution, or new scientific classification.
"""
import hashlib,json,time
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.spatial import ConvexHull
from frozen_mask_functions import robust_z,unit_dirs,knn_graph,grad_proxy,largest_component

P=Path(__file__).resolve().parent;IN=P/'input';OUT=P/'results';OUT.mkdir(exist_ok=True)

def main():
    a={p.stem:np.load(p,mmap_mode='r',allow_pickle=False) for p in (IN/'arrays').glob('*.npy')}
    report=json.loads((IN/'BH_HF1_visual_collection.json').read_text())
    for key,v in report['arrays'].items():
        assert hashlib.sha256(a[key].tobytes()).hexdigest()==v['sha256_float64_C_order']
        assert np.isfinite(a[key]).all()
    ref=pd.read_csv(IN/'BH_HF1R_snapshot_metrics.csv').set_index(['snapshot_index','field'])
    sparse=pd.read_csv(IN/'BH_HF1_snapshot_metrics.csv')
    assert len(ref)==4692 and len(a['snapshot_index'])==1566
    assert np.all(np.diff(a['actual_time_M'])>0)
    assert np.array_equal(np.diff(a['offsets']),np.prod(a['extents'],axis=1))
    for file,key in [('BH_HF1R_snapshot_metrics.csv','reference_csv_sha256'),('BH_HF1_snapshot_metrics.csv','sparse_csv_sha256')]:
        assert hashlib.sha256((IN/file).read_bytes()).hexdigest()==report[key]
        original=P/'provenance'/file
        assert original.read_bytes()==(IN/file).read_bytes()
    fields=['Rstar','Estar','Bstar'];raw=['DimlessRicciScalar','WeylE_NN','WeylB_NN']
    masks=np.zeros((len(a['coords_measurement_frame']),3),dtype=np.uint8)
    decisions=np.zeros((len(a['snapshot_index']),3),dtype=bool)
    stored_p=np.zeros((len(a['snapshot_index']),3));fractions=np.zeros_like(stored_p)
    thresholds=np.zeros_like(stored_p);triangles=[];toff=[0];checks=[];geom=[]
    for j,sid in enumerate(a['snapshot_index']):
        lo,hi=a['offsets'][j:j+2];C=np.asarray(a['coords_measurement_frame'][lo:hi])
        n,good=unit_dirs(C,a['center_inertial'][j]);assert good.all()
        nei,ang=knn_graph(n,8)
        hull=ConvexHull(n);tri=hull.simplices.copy()
        cn=np.cross(n[tri[:,1]]-n[tri[:,0]],n[tri[:,2]]-n[tri[:,0]])
        flip=np.einsum('ij,ij->i',cn,n[tri[:,0]])<0
        tri[flip]=tri[flip][:,[0,2,1]]
        # Hull in direction space gives connectivity; positions remain original C.
        assert len(hull.vertices)==len(C)
        edges=np.sort(np.concatenate([tri[:,[0,1]],tri[:,[1,2]],tri[:,[2,0]]]),axis=1)
        unique,counts=np.unique(edges,axis=0,return_counts=True)
        assert np.all(counts==2) and len(C)-len(unique)+len(tri)==2
        cn=np.cross(C[tri[:,1]]-C[tri[:,0]],C[tri[:,2]]-C[tri[:,0]])
        doublearea=np.linalg.norm(cn,axis=1)
        assert doublearea.min()>1e-10
        angular_edges=np.degrees(np.arccos(np.clip(np.einsum('ij,ij->i',n[unique[:,0]],n[unique[:,1]]),-1,1)))
        geom.append({'snapshot_index':int(sid),'offset_M':float(a['offset_M'][j]),'vertices':len(C),'triangles':len(tri),
                     'max_angular_edge_deg':float(angular_edges.max()),'min_triangle_doublearea_coordinate':float(doublearea.min()),
                     'coordinate_center_inertial_x':float(a['center_inertial'][j,0]),
                     'coordinate_center_inertial_y':float(a['center_inertial'][j,1]),
                     'closed_manifold':True})
        triangles.append(tri.astype(np.int32));toff.append(toff[-1]+len(tri))
        for k,(f,rname) in enumerate(zip(fields,raw)):
            x=np.asarray(a[rname][lo:hi]);x=x if k==0 else float(a['Mch'][j])**2*x
            z=robust_z(x);g=grad_proxy(z,nei,ang);q=float(np.quantile(g,.95));active=g>=q
            ln,li=largest_component(active,nei);an=int(active.sum())
            cen=n[li].sum(0);cen/=np.linalg.norm(cen)
            masks[lo+li,k]=1
            if (int(sid),f) in ref.index:
                rr=ref.loc[(int(sid),f)];reference_kind='HF1R'
            else:
                ss=sparse[(sparse.field==f)&(np.abs(sparse.actual_time_M-float(a['actual_time_M'][j]))<1e-8)]
                assert len(ss)==1;rr=ss.iloc[0];reference_kind='HF1_sparse'
            dc=float(np.max(abs(cen-rr[['centroid_x','centroid_y','centroid_z']].to_numpy(float))))
            dq=abs(q-float(rr.gradient_q95));df=abs(ln/max(1,an)-float(rr.largest_component_fraction_of_active))
            matched=bool(ln==rr.largest_component_points and an==rr.active_points and dc<1e-9 and dq<1e-9 and df<1e-12)
            checks.append({'snapshot_index':int(sid),'field':f,'reference':reference_kind,'matched':matched,
                           'centroid_max_abs_diff':dc,'gradient_q95_abs_diff':dq,'fraction_abs_diff':df,
                           'active_points':an,'largest_component_points':ln})
            decisions[j,k]=bool(rr.candidate_snapshot);stored_p[j,k]=float(rr.shuffle_p)
            fractions[j,k]=float(rr.largest_component_fraction_of_active);thresholds[j,k]=q
        if j%250==0:print(f'Восстановление масок и поверхности: {j}/{len(a["snapshot_index"])}',flush=True)
    check=pd.DataFrame(checks);check.to_csv(OUT/'BH_HF1_V2_mask_verification.csv',index=False)
    pd.DataFrame(geom).to_csv(OUT/'BH_HF1_V2_mesh_verification.csv',index=False)
    bad=check[~check.matched]
    if len(bad):
        print(bad.head(12).to_string(index=False));raise RuntimeError(f'Несовпадений: {len(bad)}; визуализация исходной находки заблокирована.')
    np.savez_compressed(OUT/'BH_HF1_V2_masks_mesh.npz',masks=masks,decisions=decisions,stored_p=stored_p,
        fractions=fractions,gradient_thresholds=thresholds,triangles=np.concatenate(triangles),triangle_offsets=np.array(toff))
    validation={'stage':'BH-HF1-VISUAL-V2','status':'ALL_FROZEN_MASK_METRICS_MATCHED',
        'input_archive_sha256':hashlib.sha256((P.parent/'upload/BH_HF1_visual_data.zip').read_bytes()).hexdigest(),
        'snapshots':len(a['snapshot_index']),'dense_snapshots':1564,'mask_comparisons':len(check),'all_matched':True,
        'max_centroid_component_difference':float(check.centroid_max_abs_diff.max()),
        'max_gradient_q95_difference':float(check.gradient_q95_abs_diff.max()),
        'all_component_and_active_node_counts_match':True,'all_meshes_closed_two_manifold':True,
        'max_mesh_angular_edge_deg':max(x['max_angular_edge_deg'] for x in geom),
        'surface_positions':'unmodified CoordsMeasurementFrame native Cartesian nodes',
        'surface_connectivity':'convex hull of normalized original HF1R directions; radial triangulation, no convex hull of Cartesian positions',
        'frame_handling':'no frame transformation assumed; scalar/node association is by shared snapshot and point index; prior center used only to recover frozen directions and mesh connectivity',
        'statistical_classifications':'copied unchanged from saved tables','null_tests_repeated':False,
        'evolution_run':False,'source_code_sha256':hashlib.sha256((P/'provenance/bh_hf1r_dense_temporal_coherence.py').read_bytes()).hexdigest()}
    (OUT/'BH_HF1_V2_verification.json').write_text(json.dumps(validation,ensure_ascii=False,indent=2))
    print(json.dumps(validation,ensure_ascii=False,indent=2),flush=True)

if __name__=='__main__':main()
