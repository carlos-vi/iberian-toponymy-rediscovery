"""Robustness of the layer evaluation (computed from the saved layer maps; no refitting).

1. Recovery thresholds: purity 0.4/0.5/0.6 × AUC 0.85/0.90/0.95.
2. Null models: random contiguous areas grown from a seed (main), compact areas (the m cells nearest to a seed) and
   random contiguous unions of provinces.
3. Alternative reference partitions: Navarre separate from the Basque area; Salamanca in Leonese; both; a single
   Andalusian area; a single Asturleonese area.
4. Removing accepted layers whose core has fewer than 10 or 20 cells.
5. Replication across corpora at k = 18: areas recovered in ≥ 2 of S, M, I, G, and the number expected by chance
   (per-area null rates, corpora treated as independent; and the ceiling that holds under any dependence).
Output: results/layers/robustness_*.csv
"""
import os, itertools
import numpy as np, pandas as pd
from lib import paths
from lib import layers as L
from lib.reference import AREAS, partition, area_lookup

D = os.path.join(paths.RESULTS, 'layers')
K = list(range(6, 25, 2))
VARIANTS = ['sig', 'eff', 'two']


class G:
    """Minimal grid rebuilt from the saved file."""
    def __init__(self, key):
        g = np.load(os.path.join(D, f'{key}_grid.npz'), allow_pickle=True)
        self.cx, self.cy, self.W, self.H = g['cx'], g['cy'], int(g['W']), int(g['H'])
        self.regions = np.array([None if r == 'None' else r for r in g['regions']], dtype=object)


def labels(grid, variant='main'):
    areas, excl = partition(variant)
    look = area_lookup(areas, excl)
    lab = np.array([look(r) for r in grid.regions], dtype=object)
    sizes = {a: int((lab == a).sum()) for a in areas if (lab == a).sum() >= 10}
    return lab, sizes


def H_of(key, v, k):
    return np.load(os.path.join(D, f'{key}_{v}_k{k}_H.npz'))['H']


def mean_over_k(key, v, lab, pseudo, pur=0.5, auc=0.90, min_core=0, ks=K):
    obs, nul = [], []
    for k in ks:
        H = H_of(key, v, k)
        if min_core and len(H):
            core = (H >= 0.5 * H.max(1, keepdims=True)).sum(1)
            H = H[core >= min_core]
        _, o, counts = L.evaluate(H, lab, pseudo, pur, auc)
        obs.append(o); nul.append(counts)
    nm = np.array(nul).mean(0)
    m = float(np.mean(obs))
    return dict(recovered=round(m, 2), chance=round(float(nm.mean()), 2), excess=round(m - float(nm.mean()), 2),
                p=round((1 + (nm >= m).sum()) / (1 + len(nm)), 4))


def main_metric():
    """Main metric: reference areas recovered, averaged over k = 6, 8, …, 24, against the grown-area null averaged
    the same way (replicate by replicate)."""
    rows = []
    for key in ['S', 'SI']:
        s = pd.read_csv(os.path.join(D, f'summary_{key}.csv'))
        for v in VARIANTS:
            sv = s[(s.variant == v) & s.k.isin(K)].sort_values('k')
            nm = np.mean([np.load(os.path.join(D, f'{key}_{v}_k{k}_null.npy')) for k in sv.k], axis=0)
            m = float(sv.recovered.mean())
            rows.append(dict(corpus=key, variant=v, k_values=len(sv), mean_recovered=round(m, 2),
                             chance=round(float(nm.mean()), 2), excess=round(m - float(nm.mean()), 2),
                             p=round((1 + (nm >= m).sum()) / (1 + len(nm)), 4),
                             mean_accepted_layers=round(float(sv.n_accepted.mean()), 1),
                             recovered_k18=int(sv[sv.k == 18].recovered.iloc[0]), min_recovered=int(sv.recovered.min()),
                             max_recovered=int(sv.recovered.max())))
    t = pd.DataFrame(rows)
    t.to_csv(os.path.join(D, 'main_metric.csv'), index=False)
    print(t.to_string(index=False))


def main():
    main_metric()
    rows = []
    for key in ['S', 'SI']:
        grid = G(key)
        lab, sizes = labels(grid)
        real = L.real_masks(lab, sizes)
        nulls = {'grown': L.pseudo_grown(grid, lab, sizes, R=200, seed=0),
                 'compact': L.pseudo_compact(grid, lab, sizes, R=200, seed=0),
                 'provinces': L.pseudo_provinces(grid, lab, grid.regions, sizes, real, R=200, seed=0)}
        for v in VARIANTS:
            for (pur, auc) in itertools.product([0.4, 0.5, 0.6], [0.85, 0.90, 0.95]):
                rows.append(dict(corpus=key, variant=v, check='thresholds', setting=f'purity {pur}, AUC {auc}',
                                 **mean_over_k(key, v, lab, nulls['grown'], pur, auc)))
            for nn in ['compact', 'provinces']:
                rows.append(dict(corpus=key, variant=v, check='null', setting=nn, **mean_over_k(key, v, lab, nulls[nn])))
            for mc in [10, 20]:
                rows.append(dict(corpus=key, variant=v, check='small layers removed', setting=f'core ≥ {mc} cells',
                                 **mean_over_k(key, v, lab, nulls['grown'], min_core=mc)))
            print(key, v, 'thresholds, nulls and small layers done', flush=True)
        for part in ['navarre', 'salamanca', 'navarre+salamanca', 'one-andalusia', 'asturleonese']:
            lab2, sizes2 = labels(grid, part)
            ps = L.pseudo_grown(grid, lab2, sizes2, R=200, seed=0)
            for v in VARIANTS:
                rows.append(dict(corpus=key, variant=v, check='partition', setting=part, **mean_over_k(key, v, lab2, ps)))
            print(key, part, 'done', flush=True)
    t = pd.DataFrame(rows)
    t.to_csv(os.path.join(D, 'robustness_layers.csv'), index=False, encoding='utf-8')
    pd.set_option('display.width', 220); print(t.to_string(index=False))
    # replication across corpora at k = 18
    names = list(AREAS)
    rep = []
    for v in VARIANTS:
        rec = {}; rate = {}
        for key in ['S', 'M', 'I', 'G']:
            a = pd.read_csv(os.path.join(D, f'{key}_{v}_k18_areas.csv')).set_index('area')
            rec[key] = a.recovered; rate[key] = a.null_rate
        obs = sum(sum(bool(rec[c].get(n, False)) for c in rec) >= 2 for n in names)
        exp = 0.0; ceil = 0.0
        for n in names:
            ps = [float(rate[c][n]) for c in rate if n in rate[c].index]
            # whatever the dependence between corpora, P(recovered in >= 2) <= (sum of the rates) / 2
            ceil += min(1.0, sum(ps) / 2) if len(ps) > 1 else 0.0
            dist = np.zeros(len(ps) + 1); dist[0] = 1
            for p in ps:
                dist[1:] = dist[1:] * (1 - p) + dist[:-1] * p; dist[0] *= (1 - p)
            exp += dist[2:].sum()
        rep.append(dict(variant=v, areas_in_2plus_corpora=obs, expected_by_chance=round(exp, 2),
                        max_by_chance=round(ceil, 2)))
    r = pd.DataFrame(rep); r.to_csv(os.path.join(D, 'robustness_replication.csv'), index=False, encoding='utf-8')
    print(r.to_string(index=False))


if __name__ == '__main__':
    main()
