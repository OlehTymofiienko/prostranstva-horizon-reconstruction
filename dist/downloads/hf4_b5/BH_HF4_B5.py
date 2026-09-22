#!/usr/bin/env python3
"""HF4-B5: descriptive coordinate-harmonic fits to frozen HF4-B4 fields.

This program never re-evolves QC0, rebuilds geometry, or redefines B4 masks.
Run from any directory; inputs default to the adjacent saved B4 results.
"""
import os
for _key in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_key, '1')
import argparse
import csv
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
import numpy as np
import scipy
from scipy.special import sph_harm_y
from scipy.linalg import lstsq
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent
FIELDS = ('repsi2', 'impsi2', 'abs_npsigma')
LABELS = {'repsi2': r'Re $\Psi_2$', 'impsi2': r'Im $\Psi_2$', 'abs_npsigma': r'$|\sigma|$'}
LEVELS = (2, 4, 6, 8)

def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def write_json(path, obj):
    Path(path).write_text(json.dumps(obj, ensure_ascii=False, indent=2, allow_nan=False)+'\n')

def write_csv(path, rows):
    with Path(path).open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)

def weighted_quantile(v, w, q):
    ix = np.argsort(v, kind='stable'); c = np.cumsum(w[ix])
    return float(v[ix[min(np.searchsorted(c, q*c[-1]), len(ix)-1)]])

def real_basis(theta, phi):
    columns, labels = [], []
    for ell in range(9):
        columns.append(sph_harm_y(ell, 0, theta, phi).real)
        labels.append([ell, 0, 'real'])
        for m in range(1, ell+1):
            y = sph_harm_y(ell, m, theta, phi)
            columns.extend((np.sqrt(2)*y.real, np.sqrt(2)*y.imag))
            labels.extend(([ell, m, 'cos'], [ell, m, 'sin']))
    return np.column_stack(columns), labels

def ranges(values):
    values = np.asarray(values, dtype=float)
    return {'min': float(values.min()), 'median': float(np.median(values)), 'max': float(values.max()),
            'first': float(values[0]), 'last': float(values[-1])}

def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--input', type=Path, default=ROOT.parent/'hf4_b4/results/BH_HF4_B4_fields_geometry_masks.npz')
    p.add_argument('--protocol', type=Path, default=ROOT.parent/'hf4_b4/results/BH_HF4_B5_protocol.json')
    p.add_argument('--geometry', type=Path, default=ROOT.parent/'hf4_b4/results/BH_HF4_B4_geometry.csv')
    p.add_argument('--output', type=Path, default=ROOT/'results')
    a = p.parse_args(); out = a.output; out.mkdir(parents=True, exist_ok=True)
    protocol = json.loads(a.protocol.read_text()); assert sha(a.input) == protocol['input_sha256']
    assert tuple(protocol['fields']) == FIELDS and tuple(protocol['angular_degree_limits']) == LEVELS
    iterations = protocol['frames']; assert len(iterations) == 55
    geometry = {int(r['iteration']): r for r in csv.DictReader(a.geometry.open())}
    times = np.array([float(geometry[it]['t_M']) for it in iterations])
    assert np.all(np.diff(times) > 0) and np.allclose(times, np.array(iterations)/512, rtol=0, atol=1e-14)
    theta, phi = np.meshgrid((np.arange(35)+.5)*np.pi/35, np.arange(72)*2*np.pi/72, indexing='ij')
    X, labels = real_basis(theta.ravel(), phi.ravel())
    assert X.shape == (2520,81)
    # Analytic conventions and the addition theorem check the angular basis.
    y00_error = float(np.max(np.abs(X[:,0]-1/np.sqrt(4*np.pi))))
    y10_error = float(np.max(np.abs(X[:,1]-np.sqrt(3/(4*np.pi))*np.cos(theta.ravel()))))
    addition_error = max(float(np.max(np.abs(np.sum(X[:,ell**2:(ell+1)**2]**2,axis=1)-(2*ell+1)/(4*np.pi)))) for ell in range(9))
    assert max(y00_error, y10_error, addition_error) < 1e-13
    conventions = {
        'stage': 'BH-HF4-B5', 'registered_protocol_sha256': sha(a.protocol),
        'basis': 'Y_l0; sqrt(2) Re Y_lm, sqrt(2) Im Y_lm, m>0; unit-sphere normalization with Condon-Shortley phase',
        'coordinates': 'Original QLM theta/phi, fixed simulation axes; no alignment or fitted rotation',
        'measure': 'Unmodified saved B4 physical_weights; normalized already',
        'fit': 'SVD least squares of sqrt(w) X and sqrt(w) F; three fields solved together per time and L',
        'residual_mask': 'Top tail at weighted q=0.95 of abs(residual-weighted_median(residual)); include all ties',
        'original_mask': 'Read saved B4 mask verbatim; never rebuilt',
        'residual_energy': 'sum(w*residual**2); a squared numerical field norm, NOT physical energy',
        'concentration': 'sum_inside_original_mask(w*residual**2)/sum_all(w*residual**2)',
        'dice': '2*weighted_area(intersection)/(weighted_area(original)+weighted_area(residual))',
        'degeneracy': 'R2 undefined if weighted variance <= 1e-24 * max(weighted raw mean square, float tiny)',
        'scope': 'Descriptive projection on coordinate harmonics; not intrinsic eigenfunctions, Kerr multipoles, QNMs or significance test',
        'shear': 'Ordinary scalar decomposition of |sigma|, not spin-weighted decomposition of complex sigma',
        'libraries': {'numpy': np.__version__, 'scipy': scipy.__version__, 'matplotlib': matplotlib.__version__},
    }
    # Persist the exact operational choices before fitting.
    write_json(out/'BH_HF4_B5_definition.json', conventions)
    z = np.load(a.input, allow_pickle=False)
    metric_rows, design_rows, singular_rows, coefficient_rows = [], [], [], []
    saved = {'iterations': np.array(iterations), 'times_M': times, 'theta': theta, 'phi': phi}
    audit = {'input_sha256': sha(a.input), 'protocol_sha256': sha(a.protocol), 'geometry_sha256': sha(a.geometry),
             'y00_error': y00_error, 'y10_error': y10_error, 'addition_theorem_max_error': addition_error}
    max_weight_sum_error = 0.; max_normal_error = 0.; max_qr_error = 0.; max_mean_error = 0.
    min_nested_gain = float('inf')
    for ti, (it, time) in enumerate(zip(iterations, times)):
        w = z[f'physical_weights_{it}']; F = np.column_stack([z[f'field_{f}_{it}'] for f in FIELDS])
        masks = np.column_stack([z[f'mask_{f}_{it}'] for f in FIELDS])
        assert w.shape == (2520,) and F.shape == masks.shape == (2520,3)
        assert np.all(np.isfinite(w)) and np.all(w > 0) and np.all(np.isfinite(F)) and masks.dtype == bool
        weight_error = abs(w.sum()-1); max_weight_sum_error = max(max_weight_sum_error, float(weight_error)); assert weight_error < 1e-12
        mu = w@F; var = w@((F-mu)**2); raw_sq = w@(F**2)
        assert np.all(var > 1e-24*np.maximum(raw_sq, np.finfo(float).tiny))
        previous_sse = var.copy()
        sw = np.sqrt(w); target = sw[:,None]*F
        for L in LEVELS:
            n = (L+1)**2; basis = X[:,:n]; design = sw[:,None]*basis
            coeff, _, rank, singular = lstsq(design, target, cond=None, lapack_driver='gelsd')
            assert rank == n
            fitted = basis@coeff; residual = F-fitted
            sse = w@(residual**2); rms = np.sqrt(sse); r2 = 1-sse/var
            min_nested_gain = min(min_nested_gain, float(np.min((previous_sse-sse)/var)))
            assert np.all(sse <= previous_sse+1e-12*var); previous_sse = sse
            normal = design.T@(sw[:,None]*residual)
            normal_error = float(np.max(np.abs(normal))/max(np.max(np.abs(design.T@target)),1e-300))
            max_normal_error = max(max_normal_error, normal_error)
            max_mean_error = max(max_mean_error,float(np.max(np.abs(w@residual))/max(np.sqrt(raw_sq).max(),1e-300)))
            assert normal_error < 1e-11
            if ti in (0,54) and L == 8:
                alternative, _, alt_rank, _ = lstsq(design, target, cond=None, lapack_driver='gelsy')
                error = float(np.max(np.abs(basis@alternative-fitted))/np.max(np.abs(F)))
                max_qr_error = max(max_qr_error,error); assert alt_rank==n and error < 1e-11
            design_rows.append({'iteration':it,'t_M':time,'L':L,'columns':n,'rank':rank,
                                'condition':float(singular[0]/singular[-1]),'s_min':float(singular[-1]),
                                's_max':float(singular[0]),'relative_normal_equation_error':normal_error})
            singular_rows.extend({'iteration':it,'t_M':time,'L':L,'index':j,'singular_value':s} for j,s in enumerate(singular))
            for fi, f in enumerate(FIELDS):
                r = residual[:,fi]; mask = masks[:,fi]; rmed=weighted_quantile(r,w,.5)
                rcut=weighted_quantile(np.abs(r-rmed),w,.95); rmask=np.abs(r-rmed)>=rcut
                original_area=float(w[mask].sum()); residual_area=float(w[rmask].sum()); intersection=float(w[mask & rmask].sum())
                energy_in=float(w[mask]@(r[mask]**2)); energy_fraction=energy_in/sse[fi]
                row={'iteration':it,'t_M':float(time),'tau_M':float(time-times[0]),'field':f,'L':L,
                     'weighted_mean':float(mu[fi]),'weighted_variance':float(var[fi]),'R2':float(r2[fi]),
                     'residual_RMS':float(rms[fi]),'residual_RMS_over_original_SD':float(rms[fi]/np.sqrt(var[fi])),
                     'squared_residual_norm':float(sse[fi]),'original_mask_area':original_area,
                     'residual_mask_area':residual_area,'intersection_area':intersection,
                     'dice':2*intersection/(original_area+residual_area),
                     'jaccard':intersection/(original_area+residual_area-intersection),
                     'original_mask_recall':intersection/original_area,
                     'residual_squared_norm_fraction_in_original_mask':float(energy_fraction),
                     'residual_concentration_over_area':float(energy_fraction/original_area),
                     'residual_RMS_inside_original_mask':float(np.sqrt(energy_in/original_area)),
                     'residual_RMS_outside_original_mask':float(np.sqrt(max(sse[fi]-energy_in,0)/(1-original_area))),
                     'residual_weighted_median':rmed,'residual_mask_cutoff':rcut}
                metric_rows.append(row)
                coefficient_rows.extend({'iteration':it,'t_M':time,'field':f,'L':L,'ell':lab[0],'m':lab[1],'part':lab[2],'coefficient':float(c)} for lab,c in zip(labels[:n],coeff[:,fi]))
                saved[f'fit_{f}_L{L}_{it}']=fitted[:,fi]
                saved[f'residual_{f}_L{L}_{it}']=r
                saved[f'residual_mask_{f}_L{L}_{it}']=rmask
        if ti%10 == 0 or ti==54: print(f'B5: {ti+1}/55 frames, t={time:g}M',flush=True)
    audit.update(frames=55, nodes_per_frame=2520, field_fits=len(metric_rows), design_matrices=len(design_rows),
                 weight_sum_max_error=max_weight_sum_error, max_relative_normal_equation_error=max_normal_error,
                 max_relative_residual_mean=max_mean_error, svd_qr_prediction_max_relative_difference=max_qr_error,
                 min_nested_variance_fraction_gain=min_nested_gain, full_rank_all=True,
                 condition_min=min(r['condition'] for r in design_rows),condition_max=max(r['condition'] for r in design_rows),
                 original_input_unchanged=(sha(a.input)==protocol['input_sha256']), previous_science_recomputed=False)
    assert audit['original_input_unchanged']
    write_json(out/'BH_HF4_B5_checks.json',audit)
    write_csv(out/'BH_HF4_B5_metrics.csv',metric_rows)
    write_csv(out/'BH_HF4_B5_design_checks.csv',design_rows)
    write_csv(out/'BH_HF4_B5_singular_values.csv',singular_rows)
    write_csv(out/'BH_HF4_B5_coefficients.csv',coefficient_rows)
    np.savez_compressed(out/'BH_HF4_B5_fits_residuals.npz',**saved)
    summary={'stage':'BH-HF4-B5','status':'COMPLETE_DESCRIPTIVE_COORDINATE_HARMONIC_DECOMPOSITION',
             'simulation':'QC0 early geometry t22','time_range_M':[float(times[0]),float(times[-1])],
             'frames':55,'fields':{},'input_sha256':sha(a.input),'condition_range':[audit['condition_min'],audit['condition_max']]}
    for f in FIELDS:
        summary['fields'][f]={}
        for L in LEVELS:
            rows=[r for r in metric_rows if r['field']==f and r['L']==L]
            summary['fields'][f][str(L)]={key:ranges([r[key] for r in rows]) for key in (
                'R2','residual_RMS','residual_RMS_over_original_SD','dice',
                'residual_squared_norm_fraction_in_original_mask','residual_concentration_over_area')}
    write_json(out/'BH_HF4_B5_summary.json',summary)
    make_plots(out,z,saved,metric_rows,times,iterations)
    write_json(out/'BH_HF4_B5_compute_state.json',{
        'stage':'BH-HF4-B5','computed_UTC':datetime.now(timezone.utc).isoformat(),
        'numerical_status':'COMPLETE','input_sha256':sha(a.input),'program_sha256':sha(__file__),
        'outputs':['BH_HF4_B5_metrics.csv','BH_HF4_B5_fits_residuals.npz','BH_HF4_B5_summary.json'],
        'previous_stages_recomputed':False})
    print(json.dumps({'status':summary['status'],'checks':audit,'R2_L8':{f:summary['fields'][f]['8']['R2'] for f in FIELDS}},indent=2))

def make_plots(out,z,saved,metrics,times,iterations):
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    colors=['#ac5423','#227c96','#386e42','#6d58a7']
    fig,ax=plt.subplots(3,3,figsize=(14,10),constrained_layout=True)
    for j,f in enumerate(FIELDS):
        for L,c in zip(LEVELS,colors):
            rows=[r for r in metrics if r['field']==f and r['L']==L]
            ax[0,j].plot(times,[100*r['R2'] for r in rows],label=f'L={L}',c=c)
            ax[1,j].semilogy(times,[r['residual_RMS_over_original_SD'] for r in rows],c=c)
            ax[2,j].plot(times,[100*r['residual_squared_norm_fraction_in_original_mask'] for r in rows],c=c)
        ax[0,j].set_title(LABELS[f]);ax[0,j].legend(ncol=2,fontsize=8)
        ax[2,j].axhline(5,c='gray',ls=':',lw=1)
        for row in ax[:,j]:row.grid(alpha=.16);row.set_xlabel('QC0 simulation time t/M')
    ax[0,0].set_ylabel('Explained angular variance (%)')
    ax[1,0].set_ylabel('Residual RMS / original standard deviation')
    ax[2,0].set_ylabel('Squared residual norm in original mask (%)')
    fig.suptitle('HF4-B5 | 55 saved common-horizon frames | physical-area weighted fits\nCoordinate spherical harmonics; squared residual norm is not physical energy',fontsize=14)
    fig.savefig(out/'BH_HF4_B5_overview.png',dpi=170);plt.close(fig)
    for it in (iterations[0],iterations[-1]):
        fig,axes=plt.subplots(3,3,figsize=(14,10),constrained_layout=True)
        w=z[f'physical_weights_{it}']
        for j,f in enumerate(FIELDS):
            F=z[f'field_{f}_{it}'];mu=w@F;fit=saved[f'fit_{f}_L8_{it}'];res=saved[f'residual_{f}_L8_{it}']
            limit=np.max(np.abs(F-mu));old=z[f'mask_{f}_{it}'].reshape(35,72)
            for k,v in enumerate((F-mu,fit-mu,res)):
                art=axes[j,k].imshow(v.reshape(35,72),origin='upper',extent=(0,360,180,0),aspect='auto',vmin=-limit,vmax=limit,cmap='RdBu_r')
                axes[j,k].contour(np.arange(72)*5,(np.arange(35)+.5)*180/35,old,levels=[.5],colors=['#ec9a1b'],linewidths=.8)
                axes[j,k].set_xlabel('phi (deg)');axes[j,k].set_ylabel('theta (deg)')
                axes[j,k].set_title(f'{LABELS[f]}: '+['original minus mean','L <= 8 minus same mean','original minus fit'][k])
                fig.colorbar(art,ax=axes[j,k],shrink=.7)
        fig.suptitle(f'QC0, t={it/512:g}M | same color scale within each row\nOrange contour: saved B4 selection; angular maps, not reconstructed shapes',fontsize=13)
        fig.savefig(out/f'BH_HF4_B5_maps_it{it}.png',dpi=155);plt.close(fig)

if __name__=='__main__':main()
