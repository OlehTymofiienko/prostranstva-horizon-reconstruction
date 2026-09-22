#!/usr/bin/env python3
import argparse,csv,json,subprocess,time
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection,LineCollection
from matplotlib.colors import to_rgb
from scene_data import ROOT,Data,project,COLORS

OUT=ROOT/'results';BG='#071421';FG='#e8eff6'
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':12,'text.color':FG,
                     'axes.labelcolor':FG,'xtick.color':'#a5b7c8','ytick.color':'#a5b7c8'})

class Renderer:
    def __init__(self,data):
        self.d=data;self.fig=plt.figure(figsize=(16,9),dpi=100,facecolor=BG)
        self.title=self.fig.text(.04,.945,'Выделенная область на общем горизонте',fontsize=27,weight='bold')
        self.fig.text(.04,.902,'SXS:BBH:0305 · Lev6 · исходная находка HF1/HF1R · поле кривизны R*',fontsize=14,color='#a5b7c8')
        self.clock=self.fig.text(.04,.846,'',fontsize=22,weight='bold')
        self.mode=self.fig.text(.48,.853,'',fontsize=14,color='#a5b7c8')
        self.ax=[];self.poly=[];self.lines=[];self.dots=[];self.triads=[]
        for j,x in enumerate([.027,.514]):
            ax=self.fig.add_axes([x,.226,.46,.60],facecolor=BG)
            ax.set(xlim=(-2.65,2.65),ylim=(-2.65,2.65),aspect='equal');ax.axis('off')
            p=PolyCollection([],edgecolors='none',antialiased=True,zorder=1)
            l=LineCollection([],colors='#fff4d8',linewidths=1.1,zorder=2)
            ax.add_collection(p);ax.add_collection(l)
            dots=ax.scatter([],[],s=4,color='#fff2d3',linewidths=0,zorder=3)
            self.ax.append(ax);self.poly.append(p);self.lines.append(l);self.dots.append(dots)
            self.fig.text(x+.23,.782,'Ракурс A' if j==0 else 'Противоположный ракурс',ha='center',fontsize=14,color='#c1d0df')
            ax.plot([-2.25,-1.25],[-1.95,-1.95],color='#9baabd',lw=2)
            ax.text(-1.75,-2.20,'1 M (координаты)',ha='center',fontsize=10,color='#9baabd')
            ll=[ax.plot([],[],color=c,lw=1.4,zorder=5)[0] for c in ['#f08c7f','#8ddaab','#87aff7']]
            tt=[ax.text(0,0,n,color=c,fontsize=10,zorder=5) for n,c in zip('xyz',['#f08c7f','#8ddaab','#87aff7'])]
            self.triads.append((ll,tt))
        self.fig.text(.50,.255,'Один и тот же горизонт с двух сторон; масштаб поверхности постоянный.',ha='center',fontsize=15,weight='bold')
        self.status=self.fig.text(.04,.208,'',fontsize=15)
        self.fig.text(.04,.166,'Светлые точки — выделенные узлы. Граница области проведена между ними по сетке.',fontsize=13,color='#a5b7c8')
        ax=self.fig.add_axes([.10,.089,.835,.040],facecolor=BG)
        passed=data.m['decisions'][:1564,0]
        ax.scatter(data.times[:1564][passed],np.zeros(passed.sum()),s=5,marker='s',c=COLORS[0],linewidths=0)
        ax.set(xlim=(-.1,12.1),ylim=(-1,1),yticks=[],xticks=np.arange(0,13,2))
        ax.spines[['top','left','right']].set_visible(False);ax.spines['bottom'].set_color('#5b7188')
        self.cursor=ax.axvline(0,color='white',lw=1.6)
        self.fig.text(.10,.136,'Отсчёты R*, прошедшие сохранённый критерий HF1R',fontsize=10,color='#a5b7c8')
        self.fig.text(.93,.040,'t − t₀, M',ha='right',fontsize=11,color='#a5b7c8')
        self.fig.text(.04,.013,'Координатная реконструкция общего кажущегося горизонта. Серая область: критерий не выполнен; это не доказательство исчезновения поля.',fontsize=10,color='#91a6bb')

    def update(self,j,az=-58.,mode='evolution'):
        s=self.d.get(j,0);t=float(self.d.times[j]);ok=s['valid']
        self.clock.set_text(f't − t₀ = {t:7.4f} M')
        self.mode.set_text({'evolution':'Время меняется; камера неподвижна',
            'orbit':'Время остановлено; движется только камера',
            'late':'Отдельный поздний контрольный снимок'}[mode])
        status='критерий выполнен' if ok else 'критерий не выполнен'
        self.status.set_text(f'Область R*: {status}   ·   выделено узлов: {len(s["nodes"])}   ·   сохранённое p = {s["p"]:.3f}')
        self.status.set_color(COLORS[0] if ok else '#b0bcc8')
        self.cursor.set_visible(t<=12)
        self.cursor.set_xdata([t,t])
        material=np.array(to_rgb(COLORS[0] if ok else '#758598'))
        for view in range(2):
            aa=az+180*view;el=24 if view==0 else -24
            xyz=project(s['xyz'],aa,el);normals=project(s['normals'],aa,el)
            visible=np.flatnonzero(normals[:,2]>0)
            depth=xyz[s['triangles'][visible],2].mean(1)
            visible=visible[np.argsort(depth)]
            polys=xyz[s['triangles'][visible],:2]
            light=np.array([-.35,.65,.85]);light/=np.linalg.norm(light)
            diffuse=np.clip(normals[visible]@light,0,1)
            shade=.30+.70*diffuse
            colors=np.tile(np.array(to_rgb('#78889b')),(len(visible),1))
            colors[s['selected_faces'][visible]]=material
            colors*=shade[:,None]
            self.poly[view].set_verts(polys);self.poly[view].set_facecolor(colors)
            eok=project(s['edge_normals'],aa,el)[:,2]>0
            self.lines[view].set_segments(xyz[s['edges'][eok],:2])
            self.lines[view].set_color('#fff0d1' if ok else '#b2c0d0')
            nok=project(s['node_normals'],aa,el)[:,2]>0
            self.dots[view].set_offsets(project(s['nodes'][nok],aa,el)[:,:2])
            self.dots[view].set_color('#fff0d1' if ok else '#b2c0d0')
            ll,tt=self.triads[view];q=project(np.eye(3)*.38,aa,el)[:,:2];base=np.array([1.95,-1.92])
            for k in range(3):
                ll[k].set_data([base[0],base[0]+q[k,0]],[base[1],base[1]+q[k,1]])
                tt[k].set_position(base+q[k]*1.15)

    def buffer(self):
        self.fig.canvas.draw();return self.fig.canvas.buffer_rgba()

def main():
    p=argparse.ArgumentParser();p.add_argument('--preview',action='store_true');args=p.parse_args()
    d=Data();r=Renderer(d);r.update(d.nearest(2));r.fig.savefig(OUT/'BH_HF1_V2_surface_preview.png',dpi=100)
    if args.preview:return
    fps=20;frames=[]
    for t in np.linspace(0,float(d.times[1563]),36*fps):frames.append((d.previous(float(t)),-58.,'evolution'))
    for az in np.linspace(-58,302,12*fps):frames.append((d.nearest(2),float(az),'orbit'))
    for t in [10,20,50]:
        frames.extend([(d.nearest(t),-58.,'late')]*4*fps)
    cmd=['ffmpeg','-y','-loglevel','warning','-f','rawvideo','-vcodec','rawvideo','-pix_fmt','rgba',
         '-s','1600x900','-r',str(fps),'-i','-','-an','-c:v','libx264','-preset','fast','-crf','19',
         '-pix_fmt','yuv420p','-movflags','+faststart','-threads','2',str(OUT/'BH_HF1_V2_surface_evolution.mp4')]
    log=(OUT/'render.log').open('w');ff=subprocess.Popen(cmd,stdin=subprocess.PIPE,stderr=log)
    prev=None;mapping=[];start=time.monotonic()
    try:
        for k,key in enumerate(frames):
            j,az,mode=key
            if key!=prev:r.update(j,az,mode);buf=r.buffer();prev=key
            ff.stdin.write(buf)
            mapping.append(dict(video_frame=k,video_seconds=k/fps,source_row=j,
                source_snapshot_index=int(d.ids[j]),offset_M=float(d.times[j]),camera_azimuth_A_deg=az,chapter=mode))
            if k%200==0:print(f'Фильм: {k}/{len(frames)} кадров; {time.monotonic()-start:.0f} с',flush=True)
        ff.stdin.close();rc=ff.wait()
        if rc:raise RuntimeError(f'ffmpeg: {rc}')
    finally:
        if ff.poll() is None:ff.kill()
        log.close()
    with (OUT/'BH_HF1_V2_video_frames.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=mapping[0]);w.writeheader();w.writerows(mapping)
    print('Фильм сохранён',flush=True)

if __name__=='__main__':main()
