"""Spatial scan of every name element of a corpus (see lib/scan.py).
Usage: python 30_scan.py SIo      Output: results/scan/scan_{corpus}.csv
"""
import os, sys, time
import numpy as np, pandas as pd
from scipy import sparse
from lib import paths
from lib import scan as SC


def run(key):
    t0 = time.time()
    rng = np.random.default_rng(0)
    d = paths.load_corpus(key, cols=['x', 'y', 'region', 'tok'])
    units, words, cx, cy, creg = SC.build(d)
    B = len(units)
    base = np.bincount(units.c.values, minlength=len(cx)).astype(float)
    A, meta = SC.windows(cx, cy, base, B)
    bw = np.asarray(A.T @ base).ravel()
    print(f'[{key}] records {len(d)}, units {B}, words {len(words)}, cells {len(cx)}, windows {A.shape[1]}', flush=True)
    wm = [SC.word_elements(w) for w in words]
    pm = pd.DataFrame({'wid': np.repeat(np.arange(len(words)), [len(m) for m in wm]), 'm': [m for ms in wm for m in ms]})
    nwords = pm.groupby('m').wid.nunique()
    pm = pm[pm.m.isin(nwords[nwords >= 10].index)]
    um = units[['wid', 'c']].reset_index().merge(pm, on='wid')
    ucount = um.groupby('m').size()
    um = um[um.m.isin(ucount[ucount >= 10].index)]
    mnames, mid = np.unique(um.m.values, return_inverse=True)
    um = um.assign(mid=mid)
    Cm = sparse.csr_matrix((np.ones(len(um)), (um.mid.values, um.c.values)), shape=(len(mnames), len(cx)))
    Ctot = np.asarray(Cm.sum(1)).ravel()
    print(f'[{key}] elements tested {len(mnames)}  ({time.time() - t0:.0f} s)', flush=True)
    um_sorted = um.sort_values('mid')
    starts = np.searchsorted(um_sorted.mid.values, np.arange(len(mnames) + 1))
    Acsc = A.tocsc()
    res = []
    for a in range(0, len(mnames), 400):
        b_ = min(a + 400, len(mnames))
        cw = (Cm[a:b_] @ A).toarray()
        for i in range(b_ - a):
            m = a + i
            C = Ctot[m]
            v = SC.llr(cw[i], bw, C, B)
            v[cw[i] < 5] = 0
            order = np.argsort(-v)
            sub = um_sorted.iloc[starts[m]:starts[m + 1]]
            best = None
            for j in order[:200]:
                if v[j] <= 0:
                    break
                cells = Acsc.indices[Acsc.indptr[j]:Acsc.indptr[j + 1]]
                inside = sub[np.isin(sub.c.values, cells)]
                if inside.wid.nunique() >= 5:
                    best = (j, inside); break
            if best is None:
                res.append((mnames[m], C, sub.wid.nunique(), np.nan, np.nan, np.nan, 0, 0, 0.0, 0.0, 0, '', ''))
                continue
            j, inside = best
            c = cw[i, j]; bb = bw[j]
            top = words[inside.wid.value_counts().index[:8]]
            cc = int(meta[j, 0])
            res.append((mnames[m], C, sub.wid.nunique(), cc, cy[cc] / 110.57, cx[cc] / (111.32 * np.cos(np.radians(40))),
                        int(meta[j, 1]), int(c), round(c * 1.0 / bb / (C / B), 2), float(v[j]), inside.wid.nunique(),
                        creg[cc], ' '.join(top)))
    t = pd.DataFrame(res, columns=['element', 'units', 'words', 'cell', 'lat', 'lon', 'radius', 'c', 'rr_in', 'llr',
                                   'words_in', 'region', 'examples'])
    rep = SC.c_bins(t.units.values)
    need = sorted(set(rep.values()))
    print(f'[{key}] null for {len(need)} values of C  ({time.time() - t0:.0f} s)', flush=True)
    nulls = SC.null_table(need, base / B, A, bw, B, rng)
    t['p'] = [(1 + (SC.NREP - np.searchsorted(nulls[rep[C]], x, side='left'))) / (SC.NREP + 1) if x > 0 else 1.0
              for C, x in zip(t.units.values, t.llr.fillna(0).values)]
    p = t.p.values; o = np.argsort(p); q = np.empty_like(p)
    q[o] = np.minimum.accumulate((p[o] * len(p) / np.arange(1, len(p) + 1))[::-1])[::-1]
    t['q'] = np.minimum(q, 1)
    t = t.sort_values('llr', ascending=False)
    t.to_csv(paths.results('scan', f'scan_{key}.csv'), index=False, encoding='utf-8')
    print(f'[{key}] significant at q < 0.05: {(t.q < 0.05).sum()} of {len(t)}  ({time.time() - t0:.0f} s)', flush=True)


if __name__ == '__main__':
    for key in sys.argv[1:]:
        run(key)
