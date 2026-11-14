"""Figures of the article and of the supplementary material, drawn from the files in results/.
Output: figures/fig1_data.png … fig8_scan.png and figures/figS1_level2.png … figS4_ena_provinces.png (PNG, 200 dpi).

With the environment variable FIG_LAYOUT=print, figures 1–8 are drawn in a narrow layout (7.2 in wide, no text
below 7 pt) for a journal page; the default layout is unchanged.
"""
import os
import numpy as np, pandas as pd
from scipy.ndimage import gaussian_filter
from lib import paths
from lib import mapbase as mb
from lib.mapbase import pretty
from lib.reference import AREAS, FORTS, RIVER_LL, CONQ
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, ListedColormap, LinearSegmentedColormap
from matplotlib.patches import Circle

R = paths.RESULTS
LD = os.path.join(R, 'layers')
EXT = (-9.8, 3.6, 35.9, 45.4)
AREA_COLORS = ['#4c78a8', '#9ecae9', '#f58518', '#ffbf79', '#54a24b', '#88d27a', '#b79a20', '#f2cf5b', '#439894',
               '#83bcb6', '#e45756', '#ff9d98', '#79706e', '#bab0ac', '#d67195', '#fcbfd2', '#b279a2', '#d6a5c9',
               '#9e765f', '#d8b5a5', '#6c6c9e']
VLAB = {'sig': 'selection by significance', 'eff': 'selection by effect size', 'two': 'two-level analysis'}
VCOL = {'sig': '#9a968a', 'eff': '#3f6e8c', 'two': '#a8322a'}
CORPUS_COL = dict(zip(['S', 'M', 'G', 'I'], ['#3f6e8c', '#d9762b', '#5a8f4d', '#a8322a']))
CORPUS_MK = dict(zip(['S', 'M', 'G', 'I'], ['o', 's', '^', 'D']))
PRINT = os.environ.get('FIG_LAYOUT') == 'print'
PW = 7.2  # width of the print layout, inches


def fs(x):
    """Font size: unchanged in the default layout, at least 7 pt in the print layout."""
    return max(7.0, x) if PRINT else x


def print_title(top, n, recovers='', head=''):
    """Panel title of the print layout: elements on two lines, so that they fit a narrow panel."""
    el = [pretty(x) for x in str(top).split()[:n]]; h = (len(el) + 1) // 2
    lines = ([head] if head else []) + [' '.join(el[:h]), ' '.join(el[h:])]
    if recovers:
        lines.append('→ ' + recovers.replace('Central-Southern ', 'Central-Southern\n    '))
    return '\n'.join(lines)


def pretty_list(s, n=6):
    return ' '.join(pretty(x) for x in str(s).split()[:n])


def grid_of(key):
    g = np.load(os.path.join(LD, f'{key}_grid.npz'), allow_pickle=True)
    return (g['cx'] + 0.5) * 20 + float(g['x0']), (g['cy'] + 0.5) * 20 + float(g['y0']), g['lab']


def layers_of(key, v, k):
    t = pd.read_csv(os.path.join(LD, f'{key}_{v}_k{k}_layers.csv'), keep_default_na=False)
    z = np.load(os.path.join(LD, f'{key}_{v}_k{k}_H.npz'))
    t['row_all'] = np.arange(len(t))
    acc = t[t.accepted.astype(str) == 'True'].copy(); acc['row'] = np.arange(len(acc))
    return t, acc, z['H'], z['Hall']


def clean(ax):
    for sp in ['top', 'right']:
        ax.spines[sp].set_visible(False)


def layer_panel(ax, cx, cy, h, title, ext=EXT, fs=6.6):
    h = h / h.max(); keep = h > 0.03
    mb.base(ax, prov=False, rivers=(), ext=ext)
    mb.cells(ax, cx[keep], cy[keep], h[keep], cmap='magma_r', vmin=0, vmax=1)
    ax.set_title(title, fontsize=fs, loc='left')


# ------------------------------------------------------------------ 1. data and reference areas

def fig1():
    fig, axs = plt.subplots(1, 3, figsize=(PW, 3.1) if PRINT else (12.5, 4.4))
    for ax, key, title in [(axs[0], 'S', '(a) Corpus S' if PRINT else '(a) S: OpenStreetMap and GeoNames settlements'),
                           (axs[1], 'SI', '(b) Corpus S + I' if PRINT else '(b) S + I: with the official IGN/INE gazetteer')]:
        d = paths.load_corpus(key, cols=['x', 'y'])
        x0, y0 = d.x.min(), d.y.min()
        gx = ((d.x - x0) // 20).astype(int); gy = ((d.y - y0) // 20).astype(int)
        c = pd.Series(1, index=pd.MultiIndex.from_arrays([gy, gx])).groupby(level=[0, 1]).sum()
        mb.base(ax, prov=False, ext=EXT)
        pc = mb.cells(ax, (c.index.get_level_values(1) + 0.5) * 20 + x0, (c.index.get_level_values(0) + 0.5) * 20 + y0,
                      c.values / 400 * 1000, cmap='viridis_r')
        pc.set_norm(LogNorm(vmin=10, vmax=10000))
        ax.set_title(f'{title}\n({len(d):,} records)', fontsize=8.8, loc='left')
    cb = fig.colorbar(pc, ax=axs[:2], shrink=0.5, pad=0.02, location='bottom', aspect=45)
    cb.set_label('records per 1,000 km² (20 km cells, log scale)', fontsize=8)
    ax = axs[2]
    cx, cy, lab = grid_of('SI'); names = list(AREAS)
    mb.base(ax, prov=False, ext=EXT)
    val = np.array([names.index(a) if a in names else -1 for a in lab]); ok = val >= 0
    mb.cells(ax, cx[ok], cy[ok], val[ok], cmap=ListedColormap(AREA_COLORS), vmin=-0.5, vmax=20.5)
    ex = np.array([a == 'EXCL' for a in lab])
    mb.cells(ax, cx[ex], cy[ex], np.zeros(ex.sum()), cmap=ListedColormap(['#d9d6cc']), vmin=0, vmax=1)
    for i in range(len(names)):
        m = val == i
        if m.sum():
            ax.text(np.median(cx[m]), np.median(cy[m]), str(i + 1), fontsize=fs(6.5), ha='center', va='center', zorder=12)
    ax.set_title('(c) The 21 reference areas\n(grey: not scored)' if PRINT else
                 '(c) The 21 reference dialect areas\n(grey: transition provinces, not scored)', fontsize=8.8, loc='left')
    mb.save(fig, 'fig1_data.png')


# ------------------------------------------------------------------ 2. first-level layers

def fig2():
    t, acc, H, _ = layers_of('SI', 'two', 18)
    sel = acc[acc.level.astype(int) == 1]
    cx, cy, _ = grid_of('SI')
    ncol = 5 if PRINT else 6; nrow = int(np.ceil(len(sel) / ncol))
    fig, axs = plt.subplots(nrow, ncol, figsize=(PW, 1.95 * nrow) if PRINT else (2.15 * ncol, 2.3 * nrow))
    axs = np.atleast_1d(axs).ravel()
    for ax, (_, r) in zip(axs, sel.iterrows()):
        title = print_title(r.top, 6, r.recovers) if PRINT else pretty_list(r.top, 6) + (f'\n→ {r.recovers}' if r.recovers else '')
        layer_panel(ax, cx, cy, H[r.row], title, fs=fs(6.6))
    for ax in axs[len(sel):]:
        ax.axis('off')
    plt.tight_layout(); mb.save(fig, 'fig2_layers.png')


# ------------------------------------------------------------------ 3. evaluation

def fig3():
    if PRINT:
        fig = plt.figure(figsize=(PW, 7.6)); g3 = fig.add_gridspec(2, 2, height_ratios=[1, 1.4])
        axs = [fig.add_subplot(g3[0, 0]), fig.add_subplot(g3[0, 1]), fig.add_subplot(g3[1, :])]
    else:
        fig, axs = plt.subplots(1, 3, figsize=(13, 4.6), gridspec_kw=dict(width_ratios=[1, 1, 1.25]))
    for ax, key, title in [(axs[0], 'S', '(a) Corpus S'), (axs[1], 'SI', '(b) Corpus S + I')]:
        s = pd.read_csv(os.path.join(LD, f'summary_{key}.csv'))
        for v in ['sig', 'eff', 'two']:
            sv = s[s.variant == v].sort_values('k')
            ax.plot(sv.k, sv.recovered, '-o', c=VCOL[v], ms=3.5, label=f'{VLAB[v]} (mean {sv.recovered.mean():.1f})')
            ax.plot(sv.k, sv.null_mean, ':', c=VCOL[v], lw=1)
        ax.set_ylim(0, 17.5 if PRINT else 14); ax.set_xlabel('number of first-level layers, k'); ax.set_ylabel('reference areas recovered (of 21)')
        ax.legend(fontsize=7.2, frameon=False, loc='upper left'); ax.set_title(title + '; dotted lines: chance', fontsize=9, loc='left')
        clean(ax)
    ax = axs[2]; names = list(AREAS); cols, M = [], []
    for key in ['S', 'SI', 'M', 'I', 'G']:
        for v in ['sig', 'two']:
            a = pd.read_csv(os.path.join(LD, f'{key}_{v}_k18_areas.csv')).set_index('area').reindex(names)
            M.append(np.where(a.recovered.isna(), np.nan, a.recovered.astype(float))); cols.append(f'{key}\n{v}')
    M = np.array(M).T
    ax.imshow(np.where(np.isnan(M), 0.5, M), cmap=ListedColormap(['#f1efe8', '#d9d6cc', '#1f6f5c']), vmin=0, vmax=1, aspect='auto')
    for i in range(M.shape[0]):
        for j in range(M.shape[1]):
            if np.isnan(M[i, j]):
                ax.text(j, i, '–', ha='center', va='center', fontsize=fs(6), color='#9a968a')
    ax.set_yticks(range(len(names))); ax.set_yticklabels([f'{i + 1}. {n}' for i, n in enumerate(names)], fontsize=fs(6.8))
    ax.set_xticks(range(len(cols))); ax.set_xticklabels(cols, fontsize=fs(6.8))
    for j in range(1, len(cols), 2):
        ax.axvline(j + 0.5, c='white', lw=2)
    ax.set_title('(c) Recovered at k = 18 (green); – = not covered', fontsize=9, loc='left')
    plt.tight_layout(); mb.save(fig, 'fig3_evaluation.png')


# ------------------------------------------------------------------ 4. the south

def fig4():
    panels = []
    tS, accS, HS, _ = layers_of('S', 'eff', 18)
    r = accS[accS.top.str.contains(r'\^guad', regex=True)].iloc[0]
    cxS, cyS, _ = grid_of('S')
    panels.append((cxS, cyS, HS[r.row], print_title(r.top, 6, '', '(a) S, first level') if PRINT else
                   '(a) S, first level: ' + pretty_list(r.top, 6)))
    t, acc, H, _ = layers_of('SI', 'two', 18)
    south = acc[(acc.level.astype(int) == 2) & (acc.lat.astype(float) < 40.5)]
    cx, cy, _ = grid_of('SI')
    for i, (_, r) in enumerate(south.iterrows()):
        panels.append((cx, cy, H[r.row], print_title(r.top, 6, r.recovers, f'({chr(98 + i)}) S + I, second level') if PRINT else
                       f'({chr(98 + i)}) S+I, second level: ' + pretty_list(r.top, 5) + (f'\n→ {r.recovers}' if r.recovers else '')))
    if PRINT:
        nrow = int(np.ceil(len(panels) / 3)); fig, axs = plt.subplots(nrow, 3, figsize=(PW, 2.75 * nrow))
        for ax in np.atleast_1d(axs).ravel()[len(panels):]:
            ax.axis('off')
    else:
        fig, axs = plt.subplots(1, len(panels), figsize=(2.5 * len(panels), 2.8))
    for ax, (x, y, h, title) in zip(np.atleast_1d(axs).ravel(), panels):
        layer_panel(ax, x, y, h, title, ext=(-9.8, 3.6, 35.9, 44.2), fs=fs(6.6))
    plt.tight_layout(); mb.save(fig, 'fig4_south.png')


# ------------------------------------------------------------------ 5. the Castilian wedge

def fig5():
    C = os.path.join(R, 'wedge')
    fig = plt.figure(figsize=(PW, 6.2) if PRINT else (11.5, 7.4))
    gsp = fig.add_gridspec(2, 2, height_ratios=[1.15, 1], hspace=0.3 if PRINT else 0.16, wspace=0.05)
    ax = fig.add_subplot(gsp[0, 0])
    t, _, _, Hall = layers_of('S', 'sig', 18)
    r = t[t.top.str.split().map(lambda xs: 'iello$' in xs)].iloc[0]
    cx, cy, _ = grid_of('S')
    h = Hall[r.row_all] / Hall[r.row_all].max(); keep = h > 0.02
    mb.base(ax, ext=(-9.6, 3.6, 36.0, 44.6))
    mb.cells(ax, cx[keep], cy[keep], h[keep], cmap='magma_r', vmin=0, vmax=1)
    if PRINT:
        mb.label(ax, -5.84, 43.36, 'Oviedo', dx=-70, dy=62, fs=7); mb.label(ax, -0.41, 42.14, 'Huesca', dx=45, dy=-75, fs=7)
        mb.label(ax, -3.70, 42.34, 'Burgos', fs=7)
    else:
        mb.label(ax, -5.84, 43.36, 'Oviedo', fs=6.5); mb.label(ax, -0.41, 42.14, 'Huesca', fs=6.5); mb.label(ax, -3.70, 42.34, 'Burgos', fs=6.5)
    ax.set_title('(a) Unsupervised layer of -iello\n(corpus S, selection by significance)', fontsize=8 if PRINT else 9, loc='left')
    ax = fig.add_subplot(gsp[0, 1])
    t = pd.concat([pd.read_csv(os.path.join(C, f'fh_types_{k}.csv')) for k in ['M', 'S', 'I']]).drop_duplicates(['w', 'cx', 'cy'])
    gx = t.cx.values - t.cx.min(); gy = t.cy.values - t.cy.min(); Wd, Hd = gx.max() + 1, gy.max() + 1
    nf = np.zeros((Hd, Wd)); nn = np.zeros((Hd, Wd))
    np.add.at(nn, (gy, gx), 1); np.add.at(nf, (gy, gx), (t.fh == 'f').values.astype(float))
    sf, sn = gaussian_filter(nf, 1.5), gaussian_filter(nn, 1.5); has = gaussian_filter((nn > 0).astype(float), 0.7) > 0.05
    ok = (sn > 1.5) & has
    yy, xx = np.nonzero(ok)
    mb.base(ax, ext=(-9.6, 3.6, 36.0, 44.6))
    pc = mb.cells(ax, (xx + t.cx.min() + 0.5) * 20, (yy + t.cy.min() + 0.5) * 20, (sf / np.maximum(sn, 1e-9))[ok], cmap='RdBu', vmin=0, vmax=1)
    cb = fig.colorbar(pc, ax=ax, shrink=0.55, pad=0.01); cb.set_label('share of F- forms among F-/h- pairs', fontsize=7.5)
    ax.set_title('(b) Latin F- kept (blue) or h- (red)\n(corpora M + S + I)' if PRINT else
                 '(b) Latin F- kept (blue) or turned into h- (red)\n(corpora M + S + I)', fontsize=8 if PRINT else 9, loc='left')
    ax = fig.add_subplot(gsp[1, :])
    tr = pd.read_csv(os.path.join(C, 'iello_transect.csv'))
    for k in ['S', 'M', 'G', 'I']:
        g = tr[tr.corpus == k]
        ax.plot(g.lon, g.share, '-' + CORPUS_MK[k], c=CORPUS_COL[k], ms=3.5, lw=1.2, label=f'-iello share, {k}' if PRINT else f'-iello / (-iello + -illo), {k}')
    for k in ['S', 'M', 'G', 'I']:
        g = pd.read_csv(os.path.join(C, f'fh_transect_{k}.csv')); g = g[g.n >= 8]
        ax.plot(g.lon, g.share, ':' + CORPUS_MK[k], c=CORPUS_COL[k], ms=3, lw=1, alpha=0.8, label=f'F- share, {k}' if PRINT else f'F- / (F- + h-), {k}')
    for x, l in [(-6.3, 'Asturias–León'), (-3.25, 'Castile'), (-0.3, 'Aragon')]:
        ax.text(x, 1.05, l, ha='center', fontsize=8.5)
    ax.axvspan(-4.5, -2.0, color='#d9a441', alpha=0.12)
    ax.set_xlim(-9.5, 3.5); ax.set_ylim(-0.02, 1.12)
    ax.set_xlabel('longitude (band 41.5–43.6° N; Basque provinces excluded from the F- counts)' if PRINT else
                  'longitude (band 41.5–43.6° N; Basque provinces excluded from the F-/h- counts)')
    ax.set_ylabel('share of the conservative form')
    ax.legend(ncol=4, fontsize=7, frameon=False, loc='lower center', bbox_to_anchor=(0.5, -0.45 if PRINT else -0.36))
    clean(ax)
    ax.set_title('(c) Two isoglosses across the north: -iello among -iello + -illo (solid), F- among F- + h- (dotted)'
                 if PRINT else '(c) Two isoglosses across the north, four corpora', fontsize=8 if PRINT else 9, loc='left')
    mb.save(fig, 'fig5_wedge.png')


# ------------------------------------------------------------------ 6. -ena

def fig6():
    E = os.path.join(R, 'ena')
    if PRINT:
        fig = plt.figure(figsize=(PW, 6.6))
        g6 = fig.add_gridspec(2, 2, width_ratios=[1.35, 1], height_ratios=[1.45, 1], wspace=0.25, hspace=0.55)
        gsp = {0: g6[0, 0], 1: g6[0, 1], 3: g6[1, :].subgridspec(1, 2, width_ratios=[0.5, 1])[0, 1]}
    else:
        fig = plt.figure(figsize=(13, 4.3)); gsp = fig.add_gridspec(1, 4, width_ratios=[1.3, 0.95, 0.62, 0.8], wspace=0.3)
    ax = fig.add_subplot(gsp[0])
    d = paths.load_corpus('SIo', cols=['x', 'y', 'region', 'tok'])
    d = d[d.region.str.startswith('ES-') & ~d.region.isin({'ES-SS', 'ES-BI', 'ES-VI', 'ES-NA'})]
    rows = [(t, r, x, y) for tk, r, x, y in zip(d.tok, d.region, d.x, d.y) for t in tk.split() if len(t) >= 5 and t[-3:] in ('ena', 'ana')]
    T = pd.DataFrame(rows, columns=['t', 'r', 'x', 'y']).groupby(['t', 'r'], as_index=False).agg(x=('x', 'mean'), y=('y', 'mean'))
    T['e'] = T.t.str.endswith('ena')
    mb.base(ax, ext=(-9.6, 3.6, 35.9, 44.0))
    a = T[~T.e]; e = T[T.e]
    ax.scatter(a.x, a.y, s=4, c='#9a968a', lw=0, alpha=0.7, label=f'-ana ({len(a)} names)', zorder=5)
    ax.scatter(e.x, e.y, s=9, c='#7a3b69', lw=0, label=f'-ena ({len(e)} names)', zorder=6)
    ax.legend(loc='upper right', fontsize=7.5, frameon=True, framealpha=0.9, edgecolor='none')
    ax.set_title('(a) Opaque settlement names, S + I\n(Basque provinces omitted)', fontsize=8.6, loc='left')
    ax = fig.add_subplot(gsp[1])
    r1 = pd.read_csv(os.path.join(E, 'lift_all_provinces.csv')); r2 = pd.read_csv(os.path.join(E, 'without_basque.csv'))
    lab = {'So': 'S', 'Go': 'G', 'Mo': 'M', 'Io': 'I', 'SIo': 'S + I'}
    y = np.arange(len(r1))[::-1]
    ax.errorbar(r1.ena_lift, y + 0.2, xerr=[r1.ena_lift - r1.ena_lo, r1.ena_hi - r1.ena_lift], fmt='o', c='#7a3b69', ms=4, capsize=2, label='-ena, all provinces')
    ax.errorbar(r2.ena_lift, y, xerr=[r2.ena_lift - r2.ena_lo, r2.ena_hi - r2.ena_lift], fmt='o', mfc='white', c='#7a3b69', ms=4, capsize=2, label='-ena, without Basque provinces')
    ax.errorbar(r1.ana_lift, y - 0.2, xerr=[r1.ana_lift - r1.ana_lo, r1.ana_hi - r1.ana_lift], fmt='s', c='#9a968a', ms=3.5, capsize=2, label='-ana, all provinces')
    ax.axvline(1, c='#9a968a', ls=':', lw=1)
    ax.set_yticks(y); ax.set_yticklabels([lab[k] for k in r1.corpus], fontsize=8)
    ax.set_xlabel('concentration in core Baetica (lift)')
    ax.legend(fontsize=fs(6.6), frameon=False, loc='upper center', bbox_to_anchor=(0.5, -0.17), ncol=1)
    clean(ax)
    ax.set_title('(b) Five corpora', fontsize=9, loc='left')
    ax = fig.add_subplot(gsp[3])
    g = r2[r2.corpus == 'SIo'].iloc[0]
    groups = [c[:-len(' share')] for c in r2.columns if c.endswith(' share')]
    vals = [float(g[f'{c} share']) for c in groups]
    col = ['#7a3b69' if ('Baetica' in n or 'Jaén' in n or 'Nasrid' in n) else '#d9762b' if ('Ebro' in n or 'Murcia' in n) else '#9a968a' for n in groups]
    ax.barh(range(len(vals))[::-1], vals, color=col)
    ax.set_yticks(range(len(vals))[::-1]); ax.set_yticklabels(groups, fontsize=7)
    ax.set_xlabel('-ena share of -ena + -ana names (S + I)')
    clean(ax)
    ax.set_title('(c) By region (Basque provinces omitted)', fontsize=9, loc='left')
    mb.save(fig, 'fig6_ena.png')


# ------------------------------------------------------------------ 7. Arabic names and the Middle March

def fig7():
    A = os.path.join(R, 'arabic')
    s = np.load(os.path.join(A, 'surface_SI.npz'))
    summ = pd.read_csv(os.path.join(A, 'frontier_summary.csv')).set_index('corpus').loc['SI']
    if PRINT:
        fig = plt.figure(figsize=(PW, 7.9)); gsp = fig.add_gridspec(2, 1, height_ratios=[1.9, 1], hspace=0.16)
    else:
        fig = plt.figure(figsize=(11.5, 5.0)); gsp = fig.add_gridspec(1, 2, width_ratios=[1.45, 1], wspace=0.12)
    ax = fig.add_subplot(gsp[0])
    mb.base(ax, rivers=('Duero', 'Tajo', 'Ebro', 'Guadiana', 'Guadalquivir'), ext=(-9.4, 3.4, 36.0, 43.9))
    rate = s['rate']; has = s['tot'] > 0; cell = float(s['cell'])
    ys, xs = np.nonzero(np.isfinite(rate) & has)
    cx = (xs + 0.5) * cell + float(s['x0']); cy = (ys + 0.5) * cell + float(s['y0'])
    cm = LinearSegmentedColormap.from_list('f', ['#f4f1e8', '#e3c27a', '#c96a3a', '#7d2a2a', '#2e1520'])
    pc = mb.cells(ax, cx, cy, np.clip(rate[ys, xs], 0, 60), cell=cell, cmap=cm, vmin=0, vmax=60)
    cb = fig.colorbar(pc, ax=ax, shrink=0.6, pad=0.01); cb.set_label('‰ distinct names with an Arabic marker', fontsize=7.5)
    for frac, ls, lw, c in [(25, ':', 1.1, '#555555'), (75, '--', 1.1, '#555555'), (50, '-', 2.2, '#111111')]:
        fr = pd.read_csv(os.path.join(A, f'frontier_SI_{frac}.csv'))
        x, y = mb.km(fr.lon, fr.phi)
        ax.plot(x, y, c=c, lw=lw, ls=ls, zorder=9, label=f'frontier, threshold at {frac} %')
    F = sorted(FORTS, key=lambda t: t[2])
    ax.plot(*mb.km([f[2] for f in F], [f[1] for f in F]), c='#1f6f5c', lw=1.3, ls='--', label='Umayyad Middle March strongholds', zorder=9)
    off = {'Coria': (4, 5), 'Talavera': (4, -14), 'Madrid': (-34, -14), 'Talamanca': (-62, 6), 'Guadalajara': (6, -12),
           'Gormaz': (-34, 8), 'Atienza': (-40, -14), 'Medinaceli': (6, 2)}
    for n, la, lo in FORTS:
        xx, yy = mb.km(lo, la); ax.plot(xx, yy, 'o', c='#1f6f5c', ms=3, zorder=10)
        ax.text(xx + off[n][0], yy + off[n][1], n, fontsize=fs(6.3), zorder=11,
                bbox=dict(fc='white', ec='none', alpha=0.7, pad=0.4) if PRINT else None)
    ax.plot(*mb.km([q[2] for q in RIVER_LL], [q[1] for q in RIVER_LL]), c='#3f6e8c', lw=2.4, label='Duero', zorder=8)
    if PRINT:
        ax.legend(loc='lower left', fontsize=7, frameon=True, framealpha=0.85, edgecolor='none')
    else:
        ax.legend(loc='lower left', fontsize=7, frameon=False)
    ax.set_title('(a) Share of Arabic place names (S + I) and the frontier traced from it', fontsize=8.8, loc='left')
    ax = fig.add_subplot(gsp[1])
    prof = pd.read_csv(os.path.join(A, 'profiles.csv'))
    for band, c, ls in [('west (6.5–5.1° W)', '#3f6e8c', '-'), ('west and centre (6.5–3.5° W)', '#7fa7c0', '--'),
                        ('east (3.0–2.2° W)', '#a8322a', '-')]:
        p = prof[(prof.corpus == 'SI') & (prof.band == band)].sort_values('lat')
        ax.plot(p.lat, p.per_mille, ls, c=c, lw=1.6, label=band)
    ax.set_ylim(0, 36)
    for xv, n in [(41.55, 'Duero'), (40.4, 'Central\nSystem'), (39.9, 'Tagus')]:
        ax.axvline(xv, ls=':', c='#9a968a', lw=0.8); ax.text(xv + 0.04, 33.5, n, fontsize=7, va='top')
    L50 = float(summ.threshold)
    ax.axhline(L50, c='#111111', lw=0.8, ls='--'); ax.text(43.35, L50 + 0.6, 'threshold at 50 %', fontsize=7, ha='right')
    ax.set_xlabel('latitude'); ax.set_ylabel('‰ distinct names with an Arabic marker')
    ax.legend(fontsize=7.5, frameon=False, loc='upper left' if PRINT else 'best')
    clean(ax)
    ax.set_title('(b) North–south profiles across the Meseta', fontsize=9, loc='left')
    mb.save(fig, 'fig7_march.png')


# ------------------------------------------------------------------ 8. scan and -Vque

FAM_COL = {'Basque': '#54a24b', 'Galician': '#4c78a8', 'Portuguese': '#9ecae9', 'Occitan': '#b279a2', 'Catalan': '#e45756',
           'Asturian': '#f58518', 'Arabic': '#2e1520'}
VCOLS = {'a': '#a8322a', 'e': '#d9762b', 'i': '#3f6e8c', 'u': '#5a8f4d'}
CENTRE = {'ES-' + c for c in 'TO CR CU AB GU M'.split()}
SOUTH = {'ES-' + c for c in 'AL CA CO GR H J MA SE MU'.split()}


def fig8():
    B = os.path.join(R, 'scan')
    if PRINT:
        fig = plt.figure(figsize=(PW, 6.9)); g8 = fig.add_gridspec(2, 2, height_ratios=[1.9, 1], wspace=0.06, hspace=0.42)
        gsp = {0: g8[0, 0], 1: g8[0, 1], 2: g8[1, :]}
    else:
        fig = plt.figure(figsize=(14, 4.9)); g8 = fig.add_gridspec(1, 4, width_ratios=[1.2, 1.2, 0.3, 0.7], wspace=0.1)
        gsp = {0: g8[0], 1: g8[1], 2: g8[3]}
    ax = fig.add_subplot(gsp[0])
    m = pd.read_csv(os.path.join(B, 'top50_classified.csv'), keep_default_na=False)
    ext8 = (-9.8, 4.6, 35.9, 45.4)
    mb.base(ax, prov=False, ext=ext8)
    for _, r in m.iterrows():
        x, y = mb.km(float(r.lon), float(r.lat))
        col = FAM_COL.get(r.family.split(' ')[0], '#7d786c')
        ax.add_patch(Circle((x, y), float(r.radius), fc=col, ec=col, alpha=0.12, lw=0.6, zorder=6))
        ax.plot(x, y, 'o', ms=2.5, c=col, zorder=7)
    for f, c in list(FAM_COL.items()) + [('other documented strata', '#7d786c')]:
        ax.plot([], [], 'o', c=c, label=f, ms=4)
    if PRINT:
        ax.legend(loc='upper center', bbox_to_anchor=(0.5, 0.0), fontsize=7, frameon=False, ncol=3, columnspacing=0.8, handletextpad=0.2)
    else:
        ax.legend(loc='upper center', bbox_to_anchor=(0.5, 0.0), fontsize=7, frameon=False, ncol=4, columnspacing=1.0, handletextpad=0.2)
    mb.extent(ax, ext8)
    ax.set_title('(a) The 50 strongest clusters\n(opaque names, S + I)' if PRINT else '(a) The 50 strongest clusters (opaque names, S + I)',
                 fontsize=8.6, loc='left')
    ax = fig.add_subplot(gsp[1])
    V = pd.read_csv(os.path.join(R, 'aque_eque', 'vque_SIo.csv'))
    e8 = (-9.4, 3.4, 36.0, 43.9)
    mb.base(ax, rivers=('Duero', 'Tajo', 'Guadiana', 'Guadalquivir', 'Ebro'), ext=e8)
    for v in 'aeiu':
        s = V[V.v == v]; x, y = mb.km(s.lon, s.lat)
        ax.scatter(x, y, s=14, c=VCOLS[v], label=f'-{v}que ({len(s)})', zorder=6, edgecolor='white', lw=0.3)
    sec = pd.read_csv(os.path.join(B, 'secondary_clusters.csv'))
    for mm, cl, ls, txt, dy in [('eque$', 1, '-', '-eque cluster', 1.0), ('aque$', 2, '--', '-aque, second cluster', -1.15)]:
        r = sec[(sec.corpus == 'SIo') & (sec.element == mm) & (sec.cluster == cl)]
        if not len(r):
            continue
        r = r.iloc[0]; x, y = mb.km(r.lon, r.lat)
        ax.add_patch(Circle((x, y), r.radius, fill=False, ec='#111111', lw=1.1, ls=ls, zorder=8))
        ax.text(x + r.radius * 0.75, y + dy * r.radius * 0.75, txt, fontsize=fs(6.8), zorder=9,
                bbox=dict(fc='white', ec='none', alpha=0.8, pad=0.5))
    mb.label(ax, -4.024, 39.857, 'Toledo', fs=fs(6.5))
    if PRINT:
        ax.legend(loc='upper center', bbox_to_anchor=(0.5, 0.0), fontsize=7, frameon=False, ncol=2, columnspacing=0.8, handletextpad=0.2)
    else:
        ax.legend(loc='upper center', bbox_to_anchor=(0.5, 0.0), fontsize=7, frameon=False, ncol=4, columnspacing=1.0, handletextpad=0.2)
    mb.extent(ax, e8)
    ax.set_title('(b) Opaque names in -Vque\nand their scan clusters' if PRINT else '(b) Opaque names in -Vque and their scan clusters',
                 fontsize=8.6, loc='left')
    ax = fig.add_subplot(gsp[2])
    grp = np.where(V.r.isin(CENTRE), 'Old kingdom\nof Toledo', np.where(V.r.isin(SOUTH), 'Andalusia\nand Murcia', 'Rest'))
    tab = pd.crosstab(grp, V.v).reindex(['Old kingdom\nof Toledo', 'Andalusia\nand Murcia', 'Rest']).fillna(0)
    share = tab.div(tab.sum(axis=1), axis=0); left = np.zeros(len(share))
    for v in 'aeiu':
        if v in share:
            ax.barh(range(len(share)), share[v], left=left, color=VCOLS[v], label=f'-{v}que'); left += share[v].values
    ax.legend(fontsize=7, frameon=False, ncol=4, loc='upper center', bbox_to_anchor=(0.5, -0.3 if PRINT else -0.13), handlelength=1.2,
              columnspacing=1.0)
    ax.set_yticks(range(len(share))); ax.set_yticklabels([f'{i}\n(n = {int(n)})' for i, n in zip(share.index, tab.sum(axis=1))], fontsize=7.2)
    ax.invert_yaxis(); ax.set_xlabel('share of names'); ax.set_xlim(0, 1)
    clean(ax)
    ax.set_title('(c) Vowel before -que', fontsize=9, loc='left')
    mb.save(fig, 'fig8_scan.png')


# ------------------------------------------------------------------ supplementary figures

def figS1():
    """All accepted second-level layers of S + I (two-level analysis, k = 18)."""
    t, acc, H, _ = layers_of('SI', 'two', 18)
    sel = acc[acc.level.astype(int) == 2]
    cx, cy, _ = grid_of('SI')
    ncol = 7; nrow = int(np.ceil(len(sel) / ncol))
    fig, axs = plt.subplots(nrow, ncol, figsize=(1.95 * ncol, 2.15 * nrow)); axs = np.atleast_1d(axs).ravel()
    for ax, (_, r) in zip(axs, sel.iterrows()):
        layer_panel(ax, cx, cy, H[r.row], pretty_list(r.top, 5) + (f'\n→ {r.recovers}' if r.recovers else ''), fs=6)
    for ax in axs[len(sel):]:
        ax.axis('off')
    plt.tight_layout(); mb.save(fig, 'figS1_level2.png')


def figS2():
    """Accepted first-level layers of S at k = 18 under the two selection rules."""
    cx, cy, _ = grid_of('S')
    tabs = {v: layers_of('S', v, 18) for v in ['sig', 'eff']}
    ncol = 7
    rows = {v: int(np.ceil(len(tabs[v][1]) / ncol)) for v in tabs}
    nrow = sum(rows.values())
    fig, axs = plt.subplots(nrow, ncol, figsize=(1.95 * ncol, 2.2 * nrow))
    r0 = 0
    for v in ['sig', 'eff']:
        _, acc, H, _ = tabs[v]
        for j in range(rows[v] * ncol):
            ax = axs[r0 + j // ncol, j % ncol]
            if j < len(acc):
                r = acc.iloc[j]
                layer_panel(ax, cx, cy, H[r.row], pretty_list(r.top, 4) + (f'\n→ {r.recovers}' if r.recovers else ''), fs=6)
            else:
                ax.axis('off')
        axs[r0, 0].text(-0.08, 0.5, f'selection by\n{VLAB[v].split("by ")[1]}', transform=axs[r0, 0].transAxes,
                        rotation=90, ha='right', va='center', fontsize=8)
        r0 += rows[v]
    plt.tight_layout(); mb.save(fig, 'figS2_selection.png')


def figS3():
    """Frontier φ(λ) in three corpora at the 50 % threshold, with the Duero and the strongholds."""
    A = os.path.join(R, 'arabic')
    fig, ax = plt.subplots(figsize=(7.5, 3.6))
    for key, c in [('SI', '#111111'), ('I', '#a8322a'), ('G', '#5a8f4d')]:
        fr = pd.read_csv(os.path.join(A, f'frontier_{key}_50.csv'))
        ax.plot(fr.lon, fr.phi, '-', c=c, lw=1.6, label=f'frontier, {key}')
    ax.plot(fr.lon, fr.duero, '-', c='#3f6e8c', lw=2.4, alpha=0.6, label='Duero')
    ax.plot(fr.lon, fr.forts, '--', c='#1f6f5c', lw=1.3, label='Middle March strongholds')
    ax.set_xlabel('longitude'); ax.set_ylabel('latitude'); ax.legend(fontsize=7.5, frameon=False, ncol=2)
    clean(ax); ax.set_title('Frontier of Arabic place names (threshold at 50 %)', fontsize=9, loc='left')
    mb.save(fig, 'figS3_frontiers.png')


def figS4():
    """-ena share of -ena + -ana names by province (S + I, Basque provinces omitted) against years of Islamic rule."""
    p = pd.read_csv(os.path.join(R, 'ena', 'provinces_SI.csv'), index_col=0)
    p = p[~p.index.isin({'ES-SS', 'ES-BI', 'ES-VI', 'ES-NA'})]
    fig, ax = plt.subplots(figsize=(6.5, 4.2))
    sc = ax.scatter(p.years, p.share, s=np.sqrt(p.ena + p.ana) * 6, c=p.lat, cmap='viridis', edgecolor='white', lw=0.4)
    for r, row in p.iterrows():
        ax.text(row.years + 6, row.share, r[3:], fontsize=6.5)
    cb = fig.colorbar(sc, ax=ax, pad=0.01); cb.set_label('mean latitude', fontsize=7.5)
    ax.set_xlabel('years of Islamic rule (year of the Christian conquest − 714)'); ax.set_ylabel('-ena / (-ena + -ana)')
    clean(ax); ax.set_title('-ena by province (S + I; marker area ∝ number of names)', fontsize=9, loc='left')
    mb.save(fig, 'figS4_ena_provinces.png')


if __name__ == '__main__':
    import sys
    todo = sys.argv[1:] or ['fig1', 'fig2', 'fig3', 'fig4', 'fig5', 'fig6', 'fig7', 'fig8', 'figS1', 'figS2', 'figS3', 'figS4']
    for name in todo:
        globals()[name](); print(name, 'ok', flush=True)
