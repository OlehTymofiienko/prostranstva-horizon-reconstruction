import json
from pathlib import Path
import numpy as np
import pandas as pd
from scene_data import ROOT,Data

def area(x,t):
    return np.linalg.norm(np.cross(x[t[:,1]]-x[t[:,0]],x[t[:,2]]-x[t[:,0]]),axis=1).sum()/2

def main():
    out=ROOT/'results';d=Data();v=json.loads((out/'BH_HF1_V2_viewer_data.json').read_text())
    maxerr=0
    for f in v['frames']:
        j=f['row'];lo,hi=d.a['offsets'][j:j+2];s,e=d.m['triangle_offsets'][j:j+2]
        x=np.asarray(f['xyz']).reshape(-1,3)
        maxerr=max(maxerr,float(abs(x-d.a['coords_measurement_frame'][lo:hi]).max()))
        assert maxerr<=5.01e-9
        assert np.array_equal(np.asarray(f['tri']).reshape(-1,3),d.m['triangles'][s:e])
        for k in range(3):
            assert np.array_equal(f['masks'][k],np.flatnonzero(d.m['masks'][lo:hi,k]))
            assert f['ok'][k]==d.m['decisions'][j,k] and f['p'][k]==d.m['stored_p'][j,k]
    # Cross-check JavaScript subdivision against independently invoked Python.
    relerr=0
    for r in json.loads((out/'viewer_geometry_check.json').read_text()):
        p=d.get(r['row'],r['k']);x=np.asarray(r['xyz']);tri=np.asarray(r['triangles'])
        assert np.array_equal(tri,p['triangles']) and np.array_equal(r['selected'],p['selected_faces'])
        assert np.array_equal(r['edges'],p['edges'])
        assert np.max(abs(x-p['xyz']))<=5.01e-9
        before=area(p['native_xyz'],p['native_triangles']);after=area(p['xyz'],p['triangles'])
        relerr=max(relerr,abs(after-before)/before)
    assert relerr<1e-12
    frames=pd.read_csv(out/'BH_HF1_V2_video_frames.csv')
    assert len(frames)==1200
    assert np.allclose(frames.offset_M,d.times[frames.source_row],atol=1e-12,rtol=0)
    early=frames[frames.chapter=='evolution'];orbit=frames[frames.chapter=='orbit']
    assert early.offset_M.is_monotonic_increasing and early.camera_azimuth_A_deg.nunique()==1
    assert orbit.offset_M.nunique()==1 and orbit.camera_azimuth_A_deg.nunique()>1
    result={'all_viewer_masks_match_source':True,'all_viewer_meshes_match_saved_mesh':True,
        'max_viewer_coordinate_rounding_error_M':maxerr,
        'javascript_and_python_subdivision_match':True,'subdivision_coordinate_area_relative_error_max':relerr,
        'video_frame_count':len(frames),'video_frame_times_match_source':True,
        'fixed_camera_during_evolution':True,'fixed_time_during_orbit':True,
        'native_times_interpolated':False,'native_geometry_smoothed_or_exaggerated':False}
    (out/'BH_HF1_V2_display_verification.json').write_text(json.dumps(result,indent=2))
    print(json.dumps(result,indent=2))

if __name__=='__main__':main()
