"""Coordinate surface and a piecewise-linear display of the recovered node mask."""
from pathlib import Path
from functools import lru_cache
import numpy as np

ROOT=Path(__file__).resolve().parent
DATA=ROOT/'input/arrays'
FIELDS=['Rstar','Estar','Bstar']
COLORS=['#f6a52b','#4ac3ec','#e082ce']

def project(x,az=-58.,el=24.):
    az,el=np.deg2rad([az,el])
    right=np.array([-np.sin(az),np.cos(az),0.])
    up=np.array([-np.sin(el)*np.cos(az),-np.sin(el)*np.sin(az),np.cos(el)])
    eye=np.array([np.cos(el)*np.cos(az),np.cos(el)*np.sin(az),np.sin(el)])
    return np.asarray(x)@np.array([right,up,eye]).T

def split_surface(x,faces,mask):
    """Split each mixed triangle at mask=1/2. All positions lie on input facets.

    No vertex smoothing, radial deformation, geometric exaggeration, or mask
    dilation. Original selected vertices are retained as a separate trace.
    """
    vertices=list(x);out=[];inside=[];bound=[];bound_parent=[];mids={}
    def mid(a,b):
        key=tuple(sorted([int(a),int(b)]))
        if key not in mids:
            mids[key]=len(vertices);vertices.append((x[a]+x[b])/2)
        return mids[key]
    for ii,tri in enumerate(faces):
        b=mask[tri].astype(bool);count=int(b.sum())
        if count in [0,3]:out.append(tri);inside.append(count==3);continue
        # Unique selected (k=1) or unique unselected (k=2) corner.
        s=int(np.flatnonzero(b if count==1 else ~b)[0])
        a,bb,c=map(int,np.roll(tri,-s));ab=mid(a,bb);ac=mid(a,c)
        out.extend([(a,ab,ac),(bb,c,ac),(bb,ac,ab)])
        inside.extend([count==1,count==2,count==2])
        bound.append([ab,ac]);bound_parent.append(ii)
    return np.asarray(vertices),np.asarray(out,dtype=np.int32),np.asarray(inside,bool),np.asarray(bound,dtype=np.int32).reshape(-1,2),np.asarray(bound_parent,dtype=np.int32)

class Data:
    def __init__(self):
        self.a={p.stem:np.load(p,mmap_mode='r',allow_pickle=False) for p in DATA.glob('*.npy')}
        with np.load(ROOT/'results/BH_HF1_V2_masks_mesh.npz',allow_pickle=False) as z:
            self.m={key:z[key] for key in z.files}
        self.times=self.a['offset_M'];self.ids=self.a['snapshot_index']
    @lru_cache(maxsize=4)
    def get(self,j,k=0):
        lo,hi=self.a['offsets'][j:j+2];s,e=self.m['triangle_offsets'][j:j+2]
        x=np.asarray(self.a['coords_measurement_frame'][lo:hi]);faces=self.m['triangles'][s:e]
        mask=self.m['masks'][lo:hi,k].astype(bool)
        xyz,tri,selected,edges,parents=split_surface(x,faces,mask)
        normals=np.cross(xyz[tri[:,1]]-xyz[tri[:,0]],xyz[tri[:,2]]-xyz[tri[:,0]])
        normals/=np.linalg.norm(normals,axis=1)[:,None]
        oldnormals=np.cross(x[faces[:,1]]-x[faces[:,0]],x[faces[:,2]]-x[faces[:,0]])
        vn=np.zeros_like(x)
        for q in range(3):np.add.at(vn,faces[:,q],oldnormals)
        vn/=np.linalg.norm(vn,axis=1)[:,None]
        return {'xyz':xyz,'triangles':tri,'selected_faces':selected,'edges':edges,
                'edge_normals':oldnormals[parents],'normals':normals,
                'nodes':x[mask],'node_normals':vn[mask],'native_xyz':x,'native_triangles':faces,
                'native_mask':mask,'valid':bool(self.m['decisions'][j,k]),
                'p':float(self.m['stored_p'][j,k]),'fraction':float(self.m['fractions'][j,k])}
    def nearest(self,t):return int(np.argmin(abs(self.times-t)))
    def previous(self,t):return int(np.clip(np.searchsorted(self.times[:1564],t,side='right')-1,0,1563))
