"""Base map for the figures (land, provincial boundaries, rivers; Natural Earth 1:10m), in the km projection of the
corpora (see normalize.xy). The layers are prepared by src/11_basemap.py."""
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon, Rectangle
from matplotlib.collections import PatchCollection
from . import paths

C40 = np.cos(np.radians(40.0))
INK = '#2b2a33'; LAND = '#f1efe8'; LINE = '#b9b4a6'; WATER = '#3f6e8c'
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 8.5, 'text.color': INK, 'axes.labelcolor': INK,
                     'xtick.color': INK, 'ytick.color': INK, 'axes.edgecolor': '#9a968a', 'figure.facecolor': 'white',
                     'axes.facecolor': 'white', 'savefig.facecolor': 'white', 'savefig.dpi': 200})
_B = None


def layers():
    global _B
    if _B is None:
        _B = json.load(open(paths.derived('basemap.json'), encoding='utf-8'))
    return _B


def km(lon, lat):
    return np.asarray(lon) * 111.32 * C40, np.asarray(lat) * 110.57


def extent(ax, ext):
    x0, y0 = km(ext[0], ext[2]); x1, y1 = km(ext[1], ext[3])
    ax.set_xlim(x0, x1); ax.set_ylim(y0, y1)


def base(ax, land=True, prov=True, rivers=('Duero', 'Tajo', 'Ebro', 'Guadiana', 'Guadalquivir', 'Miño', 'Garona'),
         ext=(-9.6, 3.6, 35.9, 45.4), zorder=0):
    B = layers()
    if land:
        for poly in B['land']:
            p = np.array(poly); x, y = km(p[:, 0], p[:, 1])
            ax.add_patch(Polygon(np.c_[x, y], closed=True, fc=LAND, ec='#8f8a7d', lw=0.4, zorder=zorder))
    if prov:
        for seg in B['prov']:
            p = np.array(seg); x, y = km(p[:, 0], p[:, 1])
            ax.plot(x, y, c=LINE, lw=0.25, zorder=zorder + 3)
        for seg in B['border']:
            p = np.array(seg); x, y = km(p[:, 0], p[:, 1])
            ax.plot(x, y, c='#7d786c', lw=0.5, zorder=zorder + 3)
    for r in B['rivers']:
        if any(r['n'].startswith(n) for n in rivers):
            p = np.array(r['c']); x, y = km(p[:, 0], p[:, 1])
            ax.plot(x, y, c=WATER, lw=0.6, alpha=0.8, zorder=zorder + 4)
    extent(ax, ext)
    ax.set_aspect('equal'); ax.axis('off')


def cells(ax, cx_km, cy_km, values, cell=20.0, cmap='magma_r', vmin=None, vmax=None, zorder=2):
    rects = [Rectangle((x - cell / 2, y - cell / 2), cell, cell) for x, y in zip(cx_km, cy_km)]
    pc = PatchCollection(rects, cmap=cmap, lw=0, zorder=zorder)
    pc.set_array(np.asarray(values, float)); pc.set_clim(vmin, vmax)
    ax.add_collection(pc)
    return pc


def label(ax, lon, lat, text, dx=6, dy=4, fs=7, dot=True):
    x, y = km(lon, lat)
    if dot:
        ax.plot(x, y, 'o', c=INK, ms=2.2, zorder=10)
    ax.text(x + dx, y + dy, text, fontsize=fs, color=INK, zorder=10)


def pretty(m):
    """'ena$' -> '-ena'; '^beni' -> 'beni-'."""
    return ('-' + m[:-1]) if m.endswith('$') else (m[1:] + '-') if m.startswith('^') else m


def save(fig, name):
    fig.savefig(paths.figure(name), bbox_inches='tight', facecolor='white')
    plt.close(fig)
