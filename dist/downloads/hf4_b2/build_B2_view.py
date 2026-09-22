#!/usr/bin/env python3
"""Build the standalone B2 viewer and a scientific preview from saved B2 data.

Does not rerun any analysis. The Plotly source is reused from the existing B1
viewer, or from the vendor copy included with this result package.
"""
import argparse
import json
import re
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent


def build(out, source_html=None):
    data_text = (out / 'BH_HF4_B2_view_data.json').read_text()
    data = json.loads(data_text)
    vendor = HERE / 'vendor' / 'plotly-3.3.1.min.js'
    if not vendor.exists():
        if source_html is None:
            raise FileNotFoundError('Provide --source-html with the saved B1 viewer')
        scripts = re.findall(r'<script\b[^>]*>(.*?)</script>', source_html.read_text(), re.S | re.I)
        candidates = [s for s in scripts if 'plotly.js v3.3.1' in s[:1000]]
        if len(candidates) != 1:
            raise ValueError('Expected exactly one embedded Plotly.js v3.3.1 library')
        vendor.parent.mkdir(parents=True, exist_ok=True)
        vendor.write_text(candidates[0])
    template = (HERE / 'viewer_template.html').read_text()
    assert template.count('__DATA__') == template.count('__PLOTLY__') == 1
    html = template.replace('__DATA__', data_text.replace('</', '<\\/')).replace('__PLOTLY__', vendor.read_text())
    (out / 'BH_HF4_B2_fields_3D.html').write_text(html)
    controller = re.search(r'<script id="controller">(.*?)</script>', html, re.S)[1]
    (out / 'BH_HF4_B2_controller.js').write_text(controller)

    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.colors import Normalize
    from matplotlib.cm import ScalarMappable
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection

    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10})
    fields = [('repsi2', 'Re Ψ₂'), ('impsi2', 'Im Ψ₂'), ('abs_npsigma', '|σ|')]
    selected = [next(f for f in data['frames'] if f['iteration'] == it) for it in [9536, 10528]]
    fig = plt.figure(figsize=(14, 8.2), facecolor='white')
    fig.suptitle('HF4-B2 · расположение сильных отклонений на общем горизонте QC0', fontsize=17, y=.98)
    azim, elev = 33, 23
    camera = np.array([np.cos(np.radians(elev))*np.cos(np.radians(azim)),
                       np.cos(np.radians(elev))*np.sin(np.radians(azim)), np.sin(np.radians(elev))])
    nt, nphi = data['nt'], data['np']
    for row, frame in enumerate(selected):
        obj = next(o for o in frame['objects'] if o['horizon'] == 3)
        P = np.asarray(obj['xyz']).reshape(nt, nphi, 3)
        # Front-facing subset is used only to avoid drawing rear markers through
        # the opaque surface in this static preview. HTML stores the full mask.
        dt = np.gradient(P, axis=0)
        dp = (np.roll(P, -1, axis=1)-np.roll(P, 1, axis=1))/2
        normals = np.cross(dt, dp)
        relative = P-P.mean(axis=(0,1))
        normals *= np.where(np.sum(normals*relative, axis=2)>=0, 1, -1)[...,None]
        visible = np.sum(normals*camera, axis=2).ravel()>0
        for col, (field, label) in enumerate(fields):
            ax = fig.add_subplot(2, 3, 1+row*3+col, projection='3d', computed_zorder=False)
            V = np.asarray(frame['fields'][field]).reshape(nt, nphi)
            # Repeat longitude zero for the periodic closing strip.
            closed = np.concatenate([P, P[:, :1]], axis=1)
            vc = np.concatenate([V, V[:, :1]], axis=1)
            norm = Normalize(*data['ranges'][field])
            cmap = plt.get_cmap('cividis' if field == 'abs_npsigma' else 'RdBu')
            ax.plot_surface(*np.moveaxis(closed, 2, 0), facecolors=cmap(norm(vc)),
                            rstride=1, cstride=1, shade=False, linewidth=0, antialiased=True, zorder=1)
            cap_faces = []
            for ring in [P[0], P[-1]]:
                center = ring.mean(axis=0)
                cap_faces.extend([[center, ring[j], ring[(j+1)%nphi]] for j in range(nphi)])
            ax.add_collection3d(Poly3DCollection(cap_faces, facecolors='#aeb9c5', linewidth=0, zorder=2))
            ids = np.asarray(frame['masks'][field]); pts = P.reshape(-1,3)[ids[visible[ids]]]
            ax.scatter(*pts.T, s=7, c='#ffcf33', edgecolors='#715513', linewidths=.25,
                       depthshade=False, zorder=3)
            ax.set(xlim=(-1.08,1.08), ylim=(-1.08,1.08), zlim=(-1.08,1.08),
                   xlabel='x / M', ylabel='y / M', zlabel='z / M')
            ax.set_box_aspect((1,1,1)); ax.set_proj_type('ortho'); ax.view_init(elev, azim)
            ax.set_xticks([-1,0,1]); ax.set_yticks([-1,0,1]); ax.set_zticks([-1,0,1])
            ax.tick_params(labelsize=8, pad=1)
            ax.set_title(f'{label} · τ = {frame["tau_M"]:.4f} M', pad=4, fontsize=12)
            cbar = fig.colorbar(ScalarMappable(norm=norm, cmap=cmap), ax=ax, shrink=.5, pad=.02, aspect=17)
            cbar.ax.tick_params(labelsize=8)
    fig.text(.035, .047, 'Жёлтый: верхние ≈5% отклонений от медианы по координатной площади; один ракурс и одинаковый масштаб.', fontsize=10)
    fig.text(.035, .019, 'Показаны два сохранённых момента, а не непрерывная эволюция. Координатная поверхность; полярные крышки — только для изображения.', fontsize=9, color='#536171')
    fig.subplots_adjust(left=.015, right=.985, bottom=.095, top=.90, hspace=.12, wspace=.06)
    fig.savefig(out / 'BH_HF4_B2_overview.png', dpi=160)
    plt.close(fig)
    print(json.dumps({'html_bytes': (out/'BH_HF4_B2_fields_3D.html').stat().st_size,
                      'frames': len(data['frames']), 'preview': 'BH_HF4_B2_overview.png'}))


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--output', type=Path, default=HERE/'results')
    ap.add_argument('--source-html', type=Path)
    args = ap.parse_args()
    build(args.output, args.source_html)
