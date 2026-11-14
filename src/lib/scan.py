"""Kulldorff's Bernoulli spatial scan, run separately for every name element.

Units: distinct pairs (word, 20 km cell). Elements: word endings of 2–5 letters and beginnings of 3–5 letters carried by
≥ 10 distinct words and ≥ 10 units. Windows: circles of 30, 60, 100 and 150 km around every cell with units, except
those holding more than a quarter of all units. An element's cluster is the window with the highest likelihood ratio
(excess only) among those containing ≥ 5 distinct words that carry the element. Null: 999 multinomial placements of
the element's units (exact C up to 60; ±7 % bins above), which require ≥ 5 units per window (slightly conservative).
False discovery rate: Benjamini–Hochberg across all elements of a corpus.
"""
import numpy as np, pandas as pd
from scipy import sparse
from scipy.special import xlogy
from scipy.spatial import cKDTree

CELL = 20.0
RADII = [30, 60, 100, 150]
MAXFRAC = 0.25
NREP = 999


def word_elements(w):
    out = []
    for L in (3, 4, 5):
        if len(w) > L:
            out.append('^' + w[:L])
    for L in (2, 3, 4, 5):
        if len(w) > L + 1:
            out.append(w[-L:] + '$')
    return out


def llr(c, b, C, B):
    """Bernoulli log-likelihood ratio (Kulldorff 1997), for an excess inside the window only."""
    c = np.asarray(c, float); b = np.asarray(b, float)
    oc, ob = C - c, B - b
    l1 = xlogy(c, c / np.maximum(b, 1e-12)) + xlogy(b - c, (b - c) / np.maximum(b, 1e-12))
    l2 = xlogy(oc, oc / np.maximum(ob, 1e-12)) + xlogy(ob - oc, (ob - oc) / np.maximum(ob, 1e-12))
    l0 = xlogy(C, C / B) + xlogy(B - C, (B - C) / B)
    v = l1 + l2 - l0
    return np.where((c / np.maximum(b, 1e-12)) > (oc / np.maximum(ob, 1e-12)), v, 0.0)


def build(d):
    """Units (word, cell) of a corpus, the cell centres (km) and the majority region of each cell."""
    x0, y0 = d.x.min(), d.y.min()
    gx = ((d.x - x0) // CELL).astype(int).values; gy = ((d.y - y0) // CELL).astype(int).values
    ucells, cell = np.unique(gy * 10000 + gx, return_inverse=True)
    reg = pd.DataFrame({'c': cell, 'r': d.region.values}).groupby('c').r.agg(lambda s: s.value_counts().index[0])
    cx = (ucells % 10000 + 0.5) * CELL + x0; cy = (ucells // 10000 + 0.5) * CELL + y0
    units = pd.DataFrame({'w': [t for tk in d.tok for t in tk.split()],
                          'c': np.repeat(cell, [len(tk.split()) for tk in d.tok])}).drop_duplicates().reset_index(drop=True)
    words, wid = np.unique(units.w.values, return_inverse=True)
    units['wid'] = wid
    return units, words, cx, cy, reg.reindex(np.arange(len(ucells))).values


def windows(cx, cy, base, B):
    tree = cKDTree(np.c_[cx, cy])
    rows, cols, meta = [], [], []
    j = 0
    for r in RADII:
        for c, nb in enumerate(tree.query_ball_point(np.c_[cx, cy], r)):
            if base[nb].sum() > MAXFRAC * B:
                continue
            rows += nb; cols += [j] * len(nb); meta.append((c, r)); j += 1
    A = sparse.csr_matrix((np.ones(len(rows)), (rows, cols)), shape=(len(cx), j))
    return A, np.array(meta)


def null_table(Cvals, pcell, A, bw, B, rng, nrep=NREP):
    """Sorted null maxima of the likelihood ratio for each number of units C."""
    out = {}
    AT = A.T.tocsr()
    for C in Cvals:
        cnt = rng.multinomial(C, pcell, size=nrep).astype(float)
        cw = np.asarray(AT @ cnt.T).T
        v = llr(cw, bw[None, :], C, B)
        v[cw < 5] = 0
        out[C] = np.sort(v.max(1))
    return out


def c_bins(Cs):
    """C exact up to 60; above, geometric bins of ±7 % represented by their central C."""
    rep = {}
    edges = [61.0]
    while edges[-1] < max(Cs.max(), 62) * 1.2:
        edges.append(edges[-1] * 1.15)
    for C in np.unique(Cs):
        if C <= 60:
            rep[C] = int(C)
        else:
            i = np.searchsorted(edges, C, side='right') - 1
            rep[C] = int(round(np.sqrt(edges[i] * edges[i + 1])))
    return rep


def element_units(units, words, element):
    """Units carried by the words that contain `element` (same rule as word_elements)."""
    pre = element.startswith('^'); core = element.strip('^$')
    sel = np.array([(w.startswith(core) and len(w) > len(core)) if pre else (w.endswith(core) and len(w) > len(core) + 1)
                    for w in words])
    return units[sel[units.wid.values]]
