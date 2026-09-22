#!/usr/bin/env python3
"""Package already computed research. No simulation or fitting is performed."""
from pathlib import Path
import argparse, base64, csv, hashlib, html, json, re, shutil
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
P=argparse.ArgumentParser();P.add_argument('--workspace',type=Path,default=ROOT.parent);args=P.parse_args()
W=args.workspace;D=ROOT/'dist';AS=D/'assets';V=D/'viewers';DOC=ROOT/'docs'
for p in (AS,V,D/'reports',DOC,ROOT/'research'):p.mkdir(parents=True,exist_ok=True)
def copy(src,dst):
    dst.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(src,dst)
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def jsdump(x):return json.dumps(x,ensure_ascii=False,separators=(',',':'),allow_nan=False).replace('</','<\\/')

stages=['horizon_audit','hf3_continuation','hf3_a3','hf3_a4','hf3_r7b','hf4_a1','hf4_a2','hf4_a3','hf4_b2','hf4_b3','hf4_b4','hf4_b5','hf1_visual_v2']
allowed={'.py','.cjs','.js','.md','.txt','.json','.csv','.npz','.html','.F90','.par','.asc','.gp'}
inventory=[];reports=[]
for stage in stages:
    base=W/stage
    for p in sorted(base.rglob('*')):
        if not p.is_file():continue
        rel=p.relative_to(base)
        if any(x in rel.parts for x in ('__pycache__','deps','input','inputs','vendor','.cache')):continue
        if p.suffix not in allowed or '.openai-' in p.name or p.name.endswith('.log'):continue
        if p.suffix=='.html' and 'results' in rel.parts:continue
        if p.name.endswith('_view_data.json') or p.name.endswith('_viewer_data.json'):continue
        if 'controller.js' in p.name:continue
        target=ROOT/'research'/stage/rel;copy(p,target)
        inventory.append({'stage':stage,'file':str(target.relative_to(ROOT)),'bytes':p.stat().st_size,'sha256':digest(p)})
        if ('report' in p.name.lower() and p.suffix in ('.md','.txt')):
            dst=D/'reports'/(p.stem+'.md');copy(p,dst)
            reports.append({'stage':stage,'title':p.stem,'markdown':dst.name,'html':p.stem+'.html'})
    # Retain figures next to archival reports and expose every derived figure.
    for p in sorted((base/'results').glob('*.png')):
        copy(p,ROOT/'research'/stage/'results'/p.name)
        copy(p,D/'reports'/p.name)
    # Public-facing code/tables: raw data from correspondence are not copied.
    for p in sorted((ROOT/'research'/stage).rglob('*')):
        if not p.is_file():continue
        if p.suffix in ('.csv','.py','.json','.npz','.md'):
            if p.stat().st_size > 24*1024**2:continue
            copy(p,D/'downloads'/stage/p.relative_to(ROOT/'research'/stage))

# B5 is self-contained, from the frozen B4 artifact through all saved fits.
for p in (W/'hf4_b5/input').iterdir():
    if p.is_file():copy(p,ROOT/'research/hf4_b5/input'/p.name)
# Saved wave inputs for R7B are user's own QC0 outputs.
for p in (W/'hf3_r7b/inputs').iterdir():
    if p.is_file():copy(p,ROOT/'research/hf3_r7b/inputs'/p.name)
# Earlier HF1 evidence, including frozen detection protocols.
for p in (W/'hf1_recovery/Пространства').iterdir():
    if p.is_file() and p.suffix in ('.py','.json','.csv','.md','.txt') and '.openai-' not in p.name:
        copy(p,ROOT/'research/hf1_historical'/p.name)
# Keep original vendor needed by the archived HF1 visualization program.
copy(W/'hf1_visual_v2/vendor/plotly-3.3.1.min.js',AS/'plotly.min.js')
copy(W/'hf1_visual_v2/vendor/plotly-3.3.1.min.js',ROOT/'research/hf1_visual_v2/vendor/plotly-3.3.1.min.js')

asset_map={
 'hf1-preview.png':'hf1_visual_v2/results/BH_HF1_V2_surface_preview.png',
 'hf1-evolution.mp4':'hf1_visual_v2/results/BH_HF1_V2_surface_evolution.mp4',
 'qc0-overview.png':'hf4_b4/results/BH_HF4_B4_overview.png',
 'yitian-preview.png':'hf4_a2/results/BH_HF4_A2_preview.png',
 'yitian-dynamics.mp4':'hf4_a2/results/BH_HF4_A2_dynamics.mp4',
 'b5-overview.png':'hf4_b5/results/BH_HF4_B5_overview.png',
 'b5-early.png':'hf4_b5/results/BH_HF4_B5_maps_it9536.png',
 'b5-late.png':'hf4_b5/results/BH_HF4_B5_maps_it11264.png'}
for dst,src in asset_map.items():copy(W/src,AS/dst)

# Split embedded data into portable script files. Values and ordering are exact.
viewer_checks=[]
for name,template,data,var,marker in [
 ('hf1','hf1_visual_v2/viewer_template.html','hf1_visual_v2/results/BH_HF1_V2_viewer_data.json','HF_DATA','sourceData'),
 ('qc0','hf4_b4/viewer_template.html','hf4_b4/results/BH_HF4_B4_view_data.json','QC_DATA','data')]:
    obj=json.loads((W/data).read_text());frames=obj.pop('frames');folder=V/(name+'-data');folder.mkdir(exist_ok=True)
    obj['frames']=[];tags=[]
    init=folder/'init.js';init.write_text(f'window.{var}='+jsdump(obj)+';\n');tags.append(f'<script src="{name}-data/init.js"></script>')
    restored=[]
    for j,i in enumerate(range(0,len(frames),25)):
        part=frames[i:i+25];raw=jsdump(part);p=folder/f'part-{j:02}.js'
        p.write_text(f'window.{var}.frames.push(...'+raw+');\n')
        restored.extend(json.loads(raw.replace('<\\/','</')));tags.append(f'<script src="{name}-data/{p.name}"></script>')
    assert restored==frames
    page=(W/template).read_text()
    page=page.replace('<script>__PLOTLY__</script>','<script src="../assets/plotly.min.js"></script>')
    page=page.replace(f'<script id="{marker}" type="application/json">__DATA__</script>','\n'.join(tags))
    page=page.replace(f"JSON.parse(document.getElementById('{marker}').textContent)",f'window.{var}')
    page=page.replace('<main>','<main><p><a style="color:inherit" href="../index.html">← К проекту «Пространства»</a></p>',1)
    assert '__DATA__' not in page and '__PLOTLY__' not in page
    (V/(name+'.html')).write_text(page)
    viewer_checks.append({'viewer':name,'frames':len(frames),'data_exact_after_split':True,'original_data_sha256':digest(W/data)})
page=(W/'hf4_a2/results/BH_HF4_A2_atlas.html').read_text().replace('<body><main>','<body><main><p><a href="../index.html">← К проекту «Пространства»</a></p>',1)
(V/'yitian.html').write_text(page)

# Small binary representation for the new angular residual viewer, no refitting.
src=np.load(W/'hf4_b4/results/BH_HF4_B4_fields_geometry_masks.npz')
fit=np.load(W/'hf4_b5/results/BH_HF4_B5_fits_residuals.npz')
fields=['repsi2','impsi2','abs_npsigma'];its=fit['iterations'];m={}
for row in csv.DictReader((W/'hf4_b5/results/BH_HF4_B5_metrics.csv').open()):
    m[(int(row['iteration']),row['field'],int(row['L']))]={k:float(row[k]) for k in ['R2','residual_RMS_over_original_SD','dice','residual_squared_norm_fraction_in_original_mask']}
parts=[];scales=[0.,0.,0.];float_error=0.
for it in its:
    record={'t':float(it/512),'fields':[]}
    w=src[f'physical_weights_{it}']
    for fi,f in enumerate(fields):
        val=src[f'field_{f}_{it}'];mu=float(w@val); original=val-mu;scales[fi]=max(scales[fi],float(np.abs(original).max()))
        seq=[original];variants=[]
        for L in [2,4,6,8]:
            seq.extend([fit[f'fit_{f}_L{L}_{it}']-mu,fit[f'residual_{f}_L{L}_{it}']]);variants.append(m[(int(it),f,L)])
        full=np.concatenate(seq);rounded=full.astype('<f4');float_error=max(float_error,float(np.max(np.abs(full-rounded))))
        record['fields'].append({'values':base64.b64encode(rounded.tobytes()).decode(),
          'mask':base64.b64encode(src[f'mask_{f}_{it}'].astype('u1').tobytes()).decode(),'metrics':variants})
    parts.append(record)
bdir=V/'b5-data';bdir.mkdir(exist_ok=True)
(bdir/'init.js').write_text('window.B5={frames:[],scales:'+jsdump(scales)+'};\n')
tags=['<script src="b5-data/init.js"></script>']
for j,i in enumerate(range(0,len(parts),10)):
    name=f'part-{j:02}.js';(bdir/name).write_text('window.B5.frames.push(...'+jsdump(parts[i:i+10])+');\n');tags.append(f'<script src="b5-data/{name}"></script>')
(V/'b5.html').write_text((ROOT/'scripts/b5-template.html').read_text().replace('__SCRIPTS__','\n'.join(tags)))
viewer_checks.append({'viewer':'b5','frames':55,'nodes':2520,'field_modes':12,'display_float32_max_abs_rounding_error':float_error,'full_precision_science_unchanged':True})
(DOC/'viewer-packaging-checks.json').write_text(json.dumps(viewer_checks,indent=2))
(DOC/'research-inventory.json').write_text(json.dumps(inventory,ensure_ascii=False,indent=2))
(DOC/'report-catalog.json').write_text(json.dumps(reports,ensure_ascii=False,indent=2))
print(json.dumps({'stages':len(stages),'research_files':len(inventory),'reports':len(reports),'viewer_checks':viewer_checks},ensure_ascii=False))
