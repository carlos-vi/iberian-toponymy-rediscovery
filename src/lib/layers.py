"""Name elements, their spatial clustering, layers (NMF) and their evaluation against reference areas.

- element_index(): every word-initial string of 3–5 letters and word-final string of 2–5 letters that occurs in at
  least 40 records, with the 'lexical key' of the word that carries it (records, keys).
- element_stats(): clustering 'between different words' among the 25 nearest neighbours: enrichment E, z score and the
  lower 95 % confidence bound of E (LCB), on a whole corpus or on a subset of records.
- select_z() (selection by significance) and select_lcb() (selection by effect size).
- Element-by-cell matrix, NMF with stability across random starts (Hungarian matching), accepted layers, territories
  and the second, within-territory level of the hierarchical variant.
- Evaluation: recovery rule for reference areas and three null models of random areas.
"""
import re
import numpy as np, pandas as pd
from scipy.ndimage import gaussian_filter
from scipy.spatial import cKDTree
from scipy.optimize import linear_sum_assignment
from scipy.stats import chi2, rankdata
from sklearn.decomposition import NMF

KNN = 25          # neighbours
MINN = 40         # minimum number of records carrying an element
CELL = 20.0       # km
SIGMA = 1.25      # Gaussian smoothing, in cells (≈ 25 km)
DIM = re.compile(r"(ill|it|ic|ej|uel|in|et|ot|on|uc|ich|iñ|inh|zuel)$")


def lemma(w):
    """Approximate lemma (no plural, no final vowel, no diminutive) used as the lexical key of initial elements."""
    w = re.sub(r"(es|s)$", "", w) if len(w) > 5 else w
    w = re.sub(r"[aeo]$", "", w)
    w2 = DIM.sub("", w)
    return w2 if len(w2) >= 3 else w


# ---------------------------------------------------------------- element index

class ElementIndex:
    """element -> (record indices, integer lexical keys)."""

    def __init__(self, names, starts, recs, keys):
        self.names, self.starts, self.recs, self.keys = names, starts, recs, keys
        self.pos = {m: i for i, m in enumerate(names)}

    def __contains__(self, m):
        return m in self.pos

    def __getitem__(self, m):
        i = self.pos[m]
        return self.recs[self.starts[i]:self.starts[i + 1]], self.keys[self.starts[i]:self.starts[i + 1]]

    def items(self):
        for i, m in enumerate(self.names):
            yield m, (self.recs[self.starts[i]:self.starts[i + 1]], self.keys[self.starts[i]:self.starts[i + 1]])


def element_index(tok_series, minn=MINN):
    toks = [tk.split() for tk in tok_series]
    rec = np.repeat(np.arange(len(toks)), [len(t) for t in toks]).astype(np.int32)
    flat = [t for ts in toks for t in ts]
    posn = np.concatenate([np.arange(len(t)) for t in toks]).astype(np.int16) if flat else np.zeros(0, np.int16)
    uniq, tok_id = np.unique(np.array(flat, dtype=object), return_inverse=True)
    mor, key, tid = [], [], []
    for j, t in enumerate(uniq):
        lem = lemma(t)
        for L in (3, 4, 5):
            if len(t) > L:
                mor.append('^' + t[:L]); key.append('P' + lem); tid.append(j)
        for L in (2, 3, 4, 5):
            if len(t) > L + 1:
                mor.append(t[-L:] + '$'); key.append('S' + t[:4]); tid.append(j)
    mor_names, mor_id = np.unique(np.array(mor, dtype=object), return_inverse=True)
    _, key_id = np.unique(np.array(key, dtype=object), return_inverse=True)
    tid = np.array(tid, np.int32)
    del mor, key
    mor_s, key_s = mor_id.astype(np.int32), key_id.astype(np.int32)
    first = np.searchsorted(tid, np.arange(len(uniq)))
    count = np.searchsorted(tid, np.arange(len(uniq)), side='right') - first
    occ_n = count[tok_id]
    o_rec = np.repeat(rec, occ_n)
    o_pos = np.repeat(posn, occ_n)
    sel = np.repeat(first[tok_id] - (np.cumsum(occ_n) - occ_n), occ_n)
    sel += np.arange(len(sel))
    o_mor, o_key = mor_s[sel], key_s[sel]
    del sel, mor_s, key_s
    # one entry per (element, record): the first word of the record that carries the element
    o = np.lexsort((o_pos, o_rec, o_mor))
    o_mor, o_rec, o_key = o_mor[o], o_rec[o], o_key[o]
    keep = np.ones(len(o_mor), bool)
    keep[1:] = (o_mor[1:] != o_mor[:-1]) | (o_rec[1:] != o_rec[:-1])
    o_mor, o_rec, o_key = o_mor[keep], o_rec[keep], o_key[keep]
    cnt = np.bincount(o_mor, minlength=len(mor_names))
    good = np.nonzero(cnt >= minn)[0]
    gm = np.isin(o_mor, good)
    o_mor, o_rec, o_key = o_mor[gm], o_rec[gm], o_key[gm]
    remap = -np.ones(len(mor_names), np.int64); remap[good] = np.arange(len(good))
    starts = np.concatenate([[0], np.cumsum(np.bincount(remap[o_mor], minlength=len(good)))])
    return ElementIndex(list(mor_names[good]), starts, o_rec.astype(np.int64), o_key.astype(np.int64))


def element_stats(d, ei, sub=None):
    """Clustering statistics of every element on the records `sub` (indices of d; None = all).

    Returns (table, local index). A neighbour counts as a hit only if it carries the element in a word with a different
    lexical key, so that a repeated name ('X de Cameros') does not count as a regional element. The expected count
    uses the corpus-wide rate; zr uses the rate of the neighbour's province (descriptive).
    """
    if sub is None:
        sub = np.arange(len(d))
    sub = np.asarray(sub)
    n = len(sub)
    loc = -np.ones(len(d), np.int64); loc[sub] = np.arange(n)
    ds = d.iloc[sub]
    xyv = ds[['x', 'y']].values
    _, nb = cKDTree(xyv).query(xyv, k=KNN + 1)
    nb = np.ascontiguousarray(nb[:, 1:]).astype(np.int32)
    reg_codes, reg = np.unique(ds.region.values, return_inverse=True)
    reg = reg.astype(np.int32)
    nreg = len(reg_codes)
    reg_n = np.bincount(reg, minlength=nreg).astype(float)
    R = reg[nb]
    keyof = np.full(n, -1, np.int64)
    rows, names, starts, recs_l, keys_l = [], [], [0], [], []
    xs, ys = ds.x.values, ds.y.values
    for m, (I0, k0) in ei.items():
        li = loc[I0]
        ok = li >= 0
        if ok.sum() < MINN:
            continue
        I = li[ok]; keys = k0[ok]
        nm = len(I)
        uk, kid = np.unique(keys, return_inverse=True)
        kc = np.bincount(kid).astype(float)
        keyof[I] = kid
        kn = keyof[nb[I]]
        hits = ((kn >= 0) & (kn != kid[:, None])).sum()
        pdiff = (nm - kc[kid]) / (n - 1)
        exp = KNN * pdiff.sum()
        var = KNN * (pdiff * (1 - pdiff)).sum()
        z = (hits - exp) / np.sqrt(max(var, 1e-9))
        cnt = np.bincount(reg[I], minlength=nreg).astype(float)
        kr = np.zeros((len(uk), nreg)); np.add.at(kr, (kid, reg[I]), 1.0)
        pr = (cnt[R[I]] - kr[kid[:, None], R[I]]) / reg_n[R[I]]
        er = pr.sum()
        zr = (hits - er) / np.sqrt(max((pr * (1 - pr)).sum(), 1e-9))
        keyof[I] = -1
        top = cnt.argmax()
        neff = np.exp(-(kc / nm * np.log(kc / nm)).sum())
        lcb = (0.5 * chi2.ppf(0.025, 2 * hits) / exp) if (hits > 0 and exp > 0) else 0.0
        rows.append((m, nm, len(uk), neff, kc.max() / nm, hits / max(exp, 1e-9), z, hits / max(er, 1e-9), zr,
                     reg_codes[top], cnt[top] / nm, (cnt[top] / nm) / (reg_n[top] / n), xs[I].mean(), ys[I].mean(),
                     int(hits), exp, lcb))
        names.append(m); starts.append(starts[-1] + nm); recs_l.append(I); keys_l.append(keys)
    st = pd.DataFrame(rows, columns=['element', 'n', 'ntypes', 'neff', 'domshare', 'E', 'z', 'Er', 'zr', 'topreg',
                                     'topfrac', 'toplift', 'mx', 'my', 'hits', 'exp', 'lcb'])
    st = st.sort_values('z', ascending=False).reset_index(drop=True)
    ei_loc = ElementIndex(names, np.array(starts), np.concatenate(recs_l) if recs_l else np.zeros(0, np.int64),
                          np.concatenate(keys_l) if keys_l else np.zeros(0, np.int64))
    return st, ei_loc


# ---------------------------------------------------------------- selection

def _unnest(s, score, ratio=1.5, order=None):
    """Among nested elements on the same side (e.g. -eta and -ieta) keep the shorter unless the longer scores ≥ 1.5×."""
    keep = []
    sc = dict(zip(s.element, s[score]))
    for m in (order if order is not None else sorted(s.element, key=len)):
        core = m.strip('^$')
        redundant = False
        for k in keep:
            kc = k.strip('^$')
            if ((m[0] == '^') == (k[0] == '^')) and ((m[0] == '^' and core.startswith(kc)) or (m[-1] == '$' and core.endswith(kc))):
                if sc[m] < ratio * sc[k]:
                    redundant = True; break
        if not redundant:
            keep.append(m)
    return keep


def select_z(st, zmin=20, neffmin=15):
    """Selection by significance: z > 20 and at least 15 effective distinct words."""
    s = st[(st.z > zmin) & (st.neff >= neffmin)]
    return _unnest(s, 'z', order=list(s.sort_values('element', key=lambda c: c.str.len()).element))


def select_lcb(st, lcbmin=2.0, nmin=40, neffmin=15, cap=None):
    """Selection by effect size: lower 95 % bound of the enrichment ≥ 2 (optionally the `cap` best)."""
    s = st[(st.n >= nmin) & (st.neff >= neffmin) & (st.lcb >= lcbmin)]
    keep = _unnest(s, 'lcb')
    if cap is not None and len(keep) > cap:
        lc = dict(zip(s.element, s.lcb))
        keep = sorted(keep, key=lambda m: -lc[m])[:cap]
    return keep


# ---------------------------------------------------------------- element-by-cell matrix

class Grid:
    """Fixed 20 km grid of a corpus (origin at its minimum x, y)."""

    def __init__(self, d):
        self.x0, self.y0 = d.x.min(), d.y.min()
        self.gx = ((d.x.values - self.x0) // CELL).astype(int)
        self.gy = ((d.y.values - self.y0) // CELL).astype(int)
        self.W, self.H = self.gx.max() + 1, self.gy.max() + 1
        tot = np.zeros((self.H, self.W)); np.add.at(tot, (self.gy, self.gx), 1)
        self.land = tot > 0
        self.cell_of = -np.ones((self.H, self.W), np.int64)
        ys, xs = np.nonzero(self.land)
        self.cell_of[ys, xs] = np.arange(len(ys))
        self.cy, self.cx = ys, xs
        self.rec_cell = self.cell_of[self.gy, self.gx]

    def centres_km(self):
        return (self.cx + 0.5) * CELL + self.x0, (self.cy + 0.5) * CELL + self.y0

    def latlon(self):
        x, y = self.centres_km()
        return y / 110.57, x / (111.32 * np.cos(np.radians(40.0)))


def build_X(grid, sub, ei_loc, elements):
    """Smoothed log-ratio of each element's local rate to its overall rate (negative values set to 0; rows sum to 1),
    on the cells that contain records of `sub`."""
    gx, gy = grid.gx[sub], grid.gy[sub]
    H, W = grid.H, grid.W
    tot = np.zeros((H, W)); np.add.at(tot, (gy, gx), 1)
    tots = gaussian_filter(tot, SIGMA)
    landm = tot > 0
    n = len(sub)
    X = np.zeros((len(elements), landm.sum()))
    for i, m in enumerate(elements):
        I = ei_loc[m][0]
        c = np.zeros((H, W)); np.add.at(c, (gy[I], gx[I]), 1)
        cs = gaussian_filter(c, SIGMA)
        p = len(I) / n
        r = (cs + 0.5 * p) / (tots * p + 0.5 * p)
        X[i] = np.maximum(np.log(r[landm]), 0)
    X = X / (X.sum(1, keepdims=True) + 1e-12)
    return X, grid.cell_of[landm]


# ---------------------------------------------------------------- NMF, stability, accepted layers

def nmf(X, k, init='nndsvda', seed=0, max_iter=1500):
    mdl = NMF(k, init=init, random_state=seed, max_iter=max_iter)
    W = mdl.fit_transform(X)
    return W, mdl.components_


def matched_r(Ha, Hb):
    """Correlation of each layer of Ha with its partner in Hb (Hungarian matching on map correlations)."""
    k = len(Ha)
    C = np.nan_to_num(np.corrcoef(Ha, Hb)[:k, k:])
    r, c = linear_sum_assignment(-C)
    out = np.zeros(k); out[r] = C[r, c]
    return out


def fit_k(X, k, seeds=(1, 2, 3, 4)):
    W, H = nmf(X, k)
    rs = np.array([matched_r(H, nmf(X, k, init='random', seed=s)[1]) for s in seeds])
    return dict(k=k, W=W, H=H, r=rs, stab=rs.mean())


def accepted(fit, elements, share=0.35, nmulti=3, rmin=0.8, nstable=3):
    """A layer is accepted if ≥ 3 elements have ≥ 35 % of their weight on it and it is matched with r ≥ 0.8 in ≥ 3 of the
    4 random restarts."""
    W = fit['W']
    Wn = W / (W.sum(1, keepdims=True) + 1e-12)
    multi = (Wn >= share).sum(0) >= nmulti
    stable = (fit['r'] >= rmin).sum(0) >= nstable
    tops = []
    for j in range(fit['k']):
        order = np.argsort(-W[:, j])
        tops.append(' '.join([elements[i] for i in order[:14] if Wn[i, j] > share][:12]))
    return multi, stable, tops


def level1(grid, ei_loc, elements, nrec, k):
    X, cols = build_X(grid, np.arange(nrec), ei_loc, elements)
    f = fit_k(X, k)
    multi, stable, tops = accepted(f, elements)
    Hg = np.zeros((k, grid.land.sum())); Hg[:, cols] = f['H']
    return dict(fit=f, H=Hg, multi=multi, stable=stable, tops=tops, elements=elements, level=1, parent='')


def territories(grid, L1, umin=0.25):
    """Each cell is assigned to its dominant accepted layer (normalised loading ≥ 0.25, else the nearest assigned cell)."""
    acc = np.nonzero(L1['multi'] & L1['stable'])[0]
    if len(acc) == 0:
        return {}
    Hs = L1['H'][acc]
    U = Hs / (Hs.max(1, keepdims=True) + 1e-12)
    dom = U.argmax(0); val = U.max(0)
    ok = val >= umin
    if (~ok).any() and ok.any():
        pts = np.c_[grid.cx, grid.cy]
        _, nn = cKDTree(pts[ok]).query(pts[~ok])
        dom[~ok] = dom[ok][nn]
    return {int(acc[j]): np.nonzero(dom == j)[0] for j in range(len(acc))}


def level2(d, grid, ei, cells, kmax=10, cap=150, nrec_min=4000):
    """Second level: the whole analysis repeated inside a territory with ≥ 4,000 records (statistics computed within
    the territory, the 150 elements with the highest LCB, k chosen in 2–10 as the largest k with stability ≥ 0.9 × max)."""
    sub = np.nonzero(np.isin(grid.rec_cell, cells))[0]
    if len(sub) < nrec_min:
        return None
    st, eloc = element_stats(d, ei, sub)
    elements = select_lcb(st, cap=cap)
    if len(elements) < 10:
        return dict(nrec=len(sub), nelem=len(elements), fit=None)
    X, cols = build_X(grid, sub, eloc, elements)
    fits = {k: fit_k(X, k) for k in range(2, min(kmax, len(elements) - 1) + 1)}
    stab = {k: f['stab'] for k, f in fits.items()}
    best = max(stab.values())
    kk = max(k for k, v in stab.items() if v >= 0.9 * best)
    f = fits[kk]
    multi, stable, tops = accepted(f, elements)
    Hg = np.zeros((kk, grid.land.sum())); Hg[:, cols] = f['H']
    return dict(nrec=len(sub), nelem=len(elements), k=kk, stab=stab, fit=f, H=Hg, multi=multi, stable=stable,
                tops=tops, elements=elements, level=2)


# ---------------------------------------------------------------- evaluation

def cell_labels(d, grid, lookup):
    """Reference area of each cell: the majority region of its records, mapped with `lookup`."""
    t = pd.DataFrame({'c': grid.rec_cell, 'r': d.region.values})
    maj = t.groupby('c').r.agg(lambda s: s.value_counts().index[0])
    regions = maj.reindex(np.arange(grid.land.sum())).values
    return np.array([lookup(r) for r in regions], dtype=object), regions


def valid_cells(lab):
    return np.array([(a is not None) and (a != 'EXCL') for a in lab])


def _neighbours(grid, valid):
    idx = -np.ones((grid.H, grid.W), np.int64)
    idx[grid.cy[valid], grid.cx[valid]] = np.nonzero(valid)[0]
    nbrs = []
    for c in range(len(valid)):
        y, x = grid.cy[c], grid.cx[c]
        nb = []
        for dy, dx in ((1, 0), (-1, 0), (0, 1), (0, -1)):
            yy, xx = y + dy, x + dx
            if 0 <= yy < grid.H and 0 <= xx < grid.W and idx[yy, xx] >= 0:
                nb.append(idx[yy, xx])
        nbrs.append(nb)
    return nbrs


def pseudo_grown(grid, lab, sizes, R=200, seed=0):
    """Main null: random contiguous regions of each size grown from a random seed (4-neighbour region growing)."""
    rng = np.random.default_rng(seed)
    valid = valid_cells(lab)
    nbrs = _neighbours(grid, valid)
    vcells = np.nonzero(valid)[0]
    out = {}
    for name, m in sizes.items():
        masks = []
        for _ in range(R):
            for _try in range(50):
                s = rng.choice(vcells)
                region = {s}; frontier = list(set(nbrs[s])); fset = set(frontier)
                while len(region) < m and frontier:
                    j = rng.integers(len(frontier))
                    c = frontier[j]; frontier[j] = frontier[-1]; frontier.pop(); fset.discard(c)
                    region.add(c)
                    for q in nbrs[c]:
                        if q not in region and q not in fset:
                            frontier.append(q); fset.add(q)
                if len(region) == m:
                    break
            mk = np.zeros(len(lab), bool); mk[list(region)] = True
            masks.append(mk)
        out[name] = np.array(masks)
    return out


def pseudo_compact(grid, lab, sizes, R=200, seed=0):
    """Compact null: the m valid cells nearest to a random seed."""
    rng = np.random.default_rng(seed)
    valid = valid_cells(lab)
    vidx = np.nonzero(valid)[0]
    pts = np.c_[grid.cx, grid.cy].astype(float)
    tree = cKDTree(pts[vidx])
    out = {}
    for name, m in sizes.items():
        masks = np.zeros((R, len(lab)), bool)
        for r in range(R):
            s = rng.choice(vidx)
            _, nn = tree.query(pts[s], k=m)
            masks[r, vidx[np.atleast_1d(nn)]] = True
        out[name] = masks
    return out


def pseudo_provinces(grid, lab, regions, sizes, real, R=200, seed=0):
    """Province null: random contiguous unions of provinces grown until they reach the size of the area; unions that
    reproduce the real area's set of provinces are discarded."""
    rng = np.random.default_rng(seed)
    valid = valid_cells(lab)
    nbrs = _neighbours(grid, valid)
    reg = np.array([r if v else None for r, v in zip(regions, valid)], dtype=object)
    provs = sorted({r for r in reg if r is not None})
    cells_of = {p: np.nonzero(reg == p)[0] for p in provs}
    adj = {p: set() for p in provs}
    for c in np.nonzero(valid)[0]:
        for q in nbrs[c]:
            if reg[q] != reg[c]:
                adj[reg[c]].add(reg[q])
    out = {}
    for name, m in sizes.items():
        real_set = {reg[c] for c in np.nonzero(real[name])[0]}
        masks = []
        tries = 0
        while len(masks) < R and tries < 50 * R:
            tries += 1
            s = provs[rng.integers(len(provs))]
            chosen = {s}; size = len(cells_of[s]); frontier = set(adj[s])
            while size < m and frontier:
                p = sorted(frontier)[rng.integers(len(frontier))]
                frontier.discard(p); chosen.add(p); size += len(cells_of[p])
                frontier |= {q for q in adj[p] if q not in chosen}
            if chosen == real_set:
                continue
            mk = np.zeros(len(lab), bool)
            for p in chosen:
                mk[cells_of[p]] = True
            masks.append(mk)
        out[name] = np.array(masks)
    return out


def score_masks(Hacc, masks, valid, pur_min=0.5, auc_min=0.90):
    """For each mask: is there an accepted layer with ≥ pur_min of its core (H ≥ 0.5 max) inside the mask and an AUC
    ≥ auc_min for separating the mask's cells from all other valid cells? Returns (recovered, purity, AUC)."""
    if len(Hacc) == 0:
        z = np.zeros(len(masks))
        return z.astype(bool), z, z
    Hv = Hacc[:, valid]
    M = masks[:, valid].astype(float)
    core = (Hv >= 0.5 * Hacc.max(1, keepdims=True)).astype(float)
    pur = (M @ core.T) / np.maximum(core.sum(1), 1)[None, :]
    ranks = np.vstack([rankdata(h) for h in Hv])
    n1 = M.sum(1, keepdims=True); n0 = M.shape[1] - n1
    auc = ((M @ ranks.T) - n1 * (n1 + 1) / 2) / np.maximum(n1 * n0, 1)
    ok = (pur >= pur_min) & (auc >= auc_min)
    j = np.where(ok.any(1), np.argmax(np.where(ok, auc, -1), 1), np.argmax(pur, 1))
    rows = np.arange(len(masks))
    return ok.any(1), pur[rows, j], auc[rows, j]


def real_masks(lab, names):
    return {a: np.array([x == a for x in lab]) for a in names}


def evaluate(Hacc, lab, pseudo, pur_min=0.5, auc_min=0.90):
    """Areas recovered and the null distribution of the number recovered (one pseudo-area per real area per replicate)."""
    valid = valid_cells(lab)
    names = list(pseudo)
    real = np.array([np.array([x == a for x in lab]) for a in names])
    rec, pur, auc = score_masks(Hacc, real, valid, pur_min, auc_min)
    per = pd.DataFrame(dict(area=names, n_cells=real.sum(1), recovered=rec, purity=pur.round(3), auc=auc.round(3)))
    R = min(len(v) for v in pseudo.values())
    nullrec = np.zeros((len(names), R), bool)
    for i, a in enumerate(names):
        nullrec[i] = score_masks(Hacc, pseudo[a][:R], valid, pur_min, auc_min)[0]
    per['null_rate'] = nullrec.mean(1).round(3)
    return per, int(rec.sum()), nullrec.sum(0)
