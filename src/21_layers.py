"""Layers by non-negative matrix factorisation and their evaluation against the reference dialect areas.

Three ways of choosing the elements that enter the factorisation:
  sig  selection by significance (z > 20, ≥ 15 effective distinct words);
  eff  selection by effect size (lower 95 % bound of the enrichment ≥ 2);
  two  eff, followed by a second level: inside every territory dominated by an accepted first-level layer
       (≥ 4,000 records) the whole analysis is repeated with statistics computed within the territory.
For each corpus and number of first-level layers k: accepted layers, the reference areas they recover and the main null
model (200 random contiguous areas of the same sizes per area).
Usage: python 21_layers.py S 6,8,...,24   (corpus and grid of k; default k = 18)
Output: results/layers/
"""
import os, sys, time
import numpy as np, pandas as pd
from lib import paths
from lib import layers as L
from lib.reference import AREAS, area_of, macro

DEST = os.path.join(paths.RESULTS, 'layers')
os.makedirs(DEST, exist_ok=True)


def layer_table(grid, levels, lab, valid):
    lat, lon = grid.latlon()
    names = list(AREAS)
    real = np.array([np.array([x == a for x in lab]) for a in names])
    rows, Hall, Hacc = [], [], []
    for lv in levels:
        for j in range(len(lv['H'])):
            h = lv['H'][j]
            acc = bool(lv['multi'][j] and lv['stable'][j])
            rec, pur, auc = L.score_masks(h[None, :], real, valid)
            b = int(np.argmax(pur))
            w = h / (h.sum() + 1e-12)
            core = int((h >= 0.5 * h.max()).sum()) if h.max() > 0 else 0
            rows.append(dict(level=lv['level'], parent=lv.get('parent', ''), layer=j, accepted=acc,
                             multi=bool(lv['multi'][j]), stable=bool(lv['stable'][j]),
                             r_2nd_lowest=round(float(np.sort(lv['fit']['r'][:, j])[1]), 3), core_cells=core,
                             lat=round(float((w * lat).sum()), 2), lon=round(float((w * lon).sum()), 2),
                             best_area=names[b], purity=round(float(pur[b]), 3), auc=round(float(auc[b]), 3),
                             recovers=names[b] if rec[b] else '', top=lv['tops'][j]))
            Hall.append(h)
            if acc:
                Hacc.append(h)
    nl = grid.land.sum()
    return pd.DataFrame(rows), (np.array(Hacc) if Hacc else np.zeros((0, nl))), np.array(Hall)


def run(key, kgrid):
    t0 = time.time()
    d = paths.load_corpus(key)
    grid = L.Grid(d)
    ei = L.element_index(d.tok.values)
    st, eloc = L.element_stats(d, ei)
    st.to_csv(os.path.join(DEST, f'elements_{key}.csv'), index=False, encoding='utf-8')
    lab, regions = L.cell_labels(d, grid, area_of)
    valid = L.valid_cells(lab)
    sizes = {a: int((lab == a).sum()) for a in AREAS if (lab == a).sum() >= 10}
    pseudo = L.pseudo_grown(grid, lab, sizes, R=200, seed=0)
    sel = {'sig': L.select_z(st), 'eff': L.select_lcb(st)}
    diag = pd.DataFrame([(v, macro(st.set_index('element').loc[m, 'topreg'])) for v, ms in sel.items() for m in ms],
                        columns=['selection', 'focus_region']).value_counts().unstack(0, fill_value=0)
    diag.to_csv(os.path.join(DEST, f'selection_by_region_{key}.csv'), encoding='utf-8')
    np.savez_compressed(os.path.join(DEST, f'{key}_grid.npz'), cx=grid.cx, cy=grid.cy, x0=grid.x0, y0=grid.y0,
                        W=grid.W, H=grid.H, lab=np.array([str(a) for a in lab]), regions=np.array([str(r) for r in regions]))
    print(f'[{key}] records {len(d)}, cells {grid.land.sum()}, areas {len(sizes)}, elements {len(st)}, '
          f'selected sig {len(sel["sig"])}, eff {len(sel["eff"])}  ({time.time() - t0:.0f} s)', flush=True)
    summ_path = os.path.join(DEST, f'summary_{key}.csv')
    summ = pd.read_csv(summ_path).to_dict('records') if os.path.exists(summ_path) else []
    done = {(r['variant'], r['k']) for r in summ}
    for k in kgrid:
        if all((v, k) in done for v in ('sig', 'eff', 'two')):
            continue
        L1 = {v: L.level1(grid, eloc, sel[v], len(d), k) for v in ('sig', 'eff')}
        L2s, tinfo = [], []
        for j, cells in L.territories(grid, L1['eff']).items():
            r2 = L.level2(d, grid, ei, cells)
            info = dict(parent=j, parent_top=L1['eff']['tops'][j], n_cells=len(cells),
                        nrec=(r2 or {}).get('nrec', int(np.isin(grid.rec_cell, cells).sum())),
                        n_elements=(r2 or {}).get('nelem', 0), k2=(r2 or {}).get('k', 0))
            if r2 and r2.get('fit') is not None:
                r2['parent'] = L1['eff']['tops'][j]
                L2s.append(r2)
                info['accepted2'] = int((r2['multi'] & r2['stable']).sum())
            tinfo.append(info)
        pd.DataFrame(tinfo).to_csv(os.path.join(DEST, f'{key}_two_k{k}_territories.csv'), index=False, encoding='utf-8')
        for v, levels in [('sig', [L1['sig']]), ('eff', [L1['eff']]), ('two', [L1['eff']] + L2s)]:
            tab, Hacc, Hall = layer_table(grid, levels, lab, valid)
            per, obs, counts = L.evaluate(Hacc, lab, pseudo)
            tab.to_csv(os.path.join(DEST, f'{key}_{v}_k{k}_layers.csv'), index=False, encoding='utf-8')
            per.to_csv(os.path.join(DEST, f'{key}_{v}_k{k}_areas.csv'), index=False, encoding='utf-8')
            np.save(os.path.join(DEST, f'{key}_{v}_k{k}_null.npy'), counts)
            np.savez_compressed(os.path.join(DEST, f'{key}_{v}_k{k}_H.npz'), H=Hacc.astype(np.float32),
                                Hall=Hall.astype(np.float32))
            p = (1 + (counts >= obs).sum()) / (1 + len(counts))
            row = dict(corpus=key, variant=v, k=k, n_layers=len(tab), n_accepted=int(tab.accepted.sum()),
                       n_areas=len(sizes), recovered=obs, null_mean=round(counts.mean(), 2), p=round(p, 4),
                       areas=' | '.join(per[per.recovered].area))
            summ = [r for r in summ if not (r['variant'] == v and r['k'] == k)] + [row]
            print(f'[{key}] k={k} {v}: layers {len(tab)}, accepted {row["n_accepted"]}, recovered {obs}/{len(sizes)} '
                  f'(chance {row["null_mean"]}, p={row["p"]})  ({time.time() - t0:.0f} s)', flush=True)
        pd.DataFrame(summ).sort_values(['variant', 'k']).to_csv(summ_path, index=False, encoding='utf-8')


if __name__ == '__main__':
    key = sys.argv[1]
    kgrid = [int(x) for x in sys.argv[2].split(',')] if len(sys.argv) > 2 else [18]
    run(key, kgrid)
