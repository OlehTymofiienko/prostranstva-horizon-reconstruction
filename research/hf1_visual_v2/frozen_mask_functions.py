# Functions copied verbatim from frozen HF1R source. No statistical tests.

import numpy as np

from scipy.spatial import cKDTree

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
