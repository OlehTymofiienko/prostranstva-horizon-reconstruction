#!/usr/bin/env python3
import json
from pathlib import Path
import numpy as np
from scene_data import ROOT,Data

def main():
    d=Data();out=ROOT/'results'
    selected={d.previous(float(t)) for t in np.linspace(0,float(d.times[1563]),241)}
    for t in [0,2,5,10,20,50]:selected.add(d.nearest(t))
    transitions=np.flatnonzero(np.any(np.diff(d.m['decisions'][:1564].astype(int),axis=0),axis=1))
    selected.update(map(int,transitions));selected.update(map(int,transitions+1));selected.add(1563)
    records=[]
    for j in sorted(selected):
        lo,hi=d.a['offsets'][j:j+2];s,e=d.m['triangle_offsets'][j:j+2]
        record={'row':j,'id':int(d.ids[j]),'t':float(d.times[j]),'dense':j<1564,
                'xyz':np.round(d.a['coords_measurement_frame'][lo:hi],8).ravel().tolist(),
                'tri':d.m['triangles'][s:e].ravel().tolist(),
                'masks':[np.flatnonzero(d.m['masks'][lo:hi,k]).tolist() for k in range(3)],
                'ok':d.m['decisions'][j].tolist(),'p':d.m['stored_p'][j].tolist()}
        records.append(record)
    payload={'simulation':'SXS:BBH:0305','lev':6,'horizon':'C','denseEnd':float(d.times[1563]),
             'all_snapshots':1566,'coordinate_rounding_decimals_for_view_only':8,'frames':records}
    jsdata=json.dumps(payload,separators=(',',':'),ensure_ascii=False)
    (out/'BH_HF1_V2_viewer_data.json').write_text(jsdata)
    vendor=(ROOT/'vendor/plotly-3.3.1.min.js').read_text()
    page=(ROOT/'viewer_template.html').read_text().replace('__PLOTLY__',vendor).replace('__DATA__',jsdata)
    (out/'BH_HF1_V2_surface_3D.html').write_text(page)
    print('Снимков в просмотрщике:',len(records),'размер, МБ:',len(page.encode())/1e6,flush=True)

if __name__=='__main__':main()
