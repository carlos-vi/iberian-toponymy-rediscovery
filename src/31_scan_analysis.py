"""Analysis of the scans.

1. Validation (S+I, all words): share of the elements of known strata (lib.reference.TARGETS) whose cluster is
   significant (q < 0.05) and centred in the expected regions, against 30 random elements of similar size per element.
2. Where the clusters of selected elements fall (S+I opaque names; G, M and I opaque names as replications).
3. The strongest clusters among opaque names (nested elements on the same side with centres < 100 km merged), classified
   by hand: K = stratum or naming system no longer productive; D = morphology or vocabulary of a living language
   (Basque included); U = undocumented. Ranking by concentration and the rarity of -eque among endings of its size.
4. Secondary clusters (the best window sharing no cell with the main one; conservative p) and an exact
   without-replacement null for the -eque cluster.
Output: results/scan/*.csv
"""
import os
import numpy as np, pandas as pd
from lib import paths
from lib import scan as SC
from lib.reference import TARGETS

D = os.path.join(paths.RESULTS, 'scan')
TOLEDO = (39.857, -4.024)
AND = {'ES-' + c for c in 'AL CA CO GR H J MA SE'.split()}
CHECKS = [('ena$', AND), ('iello$', {'ES-O', 'ES-LE', 'ES-ZA', 'ES-HU', 'ES-Z'}), ('^beni', {'ES-V', 'ES-A', 'ES-CS', 'ES-PM'}),
          ('^guad', AND | {'ES-' + c for c in 'BA CC AB CR CU GU TO MU M'.split()} | {f'PT-{x}' for x in '09 03 08 19 14 18 16'.split()}),
          ('eque$', None), ('aque$', None), ('que$', None)]
CLASS = {   # hand classification of the 50 strongest opaque clusters (examples in the scan table)
    'koa$': ('Basque', 'D'), 'ac$': ('Gallo-Romance -ACUM', 'K'), 'rri$': ('Basque', 'D'), 'aga$': ('Basque', 'D'),
    'gi$': ('Basque', 'D'), 'ie$': ('Occitan', 'D'), 'er$': ('Catalan', 'D'), 'eta$': ('Basque', 'D'), 'ea$': ('Basque', 'D'),
    'tza$': ('Basque', 'D'), 'ena$': ('Basque', 'D'), 'ão$': ('Portuguese', 'D'), 'iro$': ('Galician', 'D'),
    'di$': ('Basque', 'D'), '^bini': ('Arabic', 'K'), '^etx': ('Basque', 'D'), 'ira$': ('Galician', 'D'),
    'eira$': ('Portuguese', 'D'), '^ola': ('Basque', 'D'), 'de$': ('Germanic genitive (-sindi > -sende)', 'K'),
    '^err': ('Basque', 'D'), 'll$': ('Catalan', 'D'), '^larr': ('Basque', 'D'), '^vill': ('Castilian repopulation villa-', 'K'),
    'txe$': ('Basque', 'D'), 'eko$': ('Basque', 'D'), 'ux$': ('Occitan', 'D'), 'as$': ('Galician', 'D'), 'ki$': ('Basque', 'D'),
    'ako$': ('Basque', 'D'), 'ndi$': ('Basque', 'D'), 'ru$': ('Asturian', 'D'), 'ido$': ('Galician', 'D'),
    'ain$': ('Basque-Navarrese -ain', 'K'), 'ães$': ('Portuguese -ANES (Latin genitive of personal names)', 'K'),
    'txu$': ('Basque', 'D'), 'illa$': ('Castilian', 'D'), '^goi': ('Basque', 'D'), 'llu$': ('Asturian (-iellu)', 'D'),
    'ha$': ('Portuguese', 'D'), 'mil$': ('Germanic personal names (-mir)', 'K'), 'rrota$': ('Basque', 'D'),
    '^agi': ('Basque', 'D'), 'ka$': ('Basque', 'D'), '^lla': ('Asturian', 'D'), 'ès$': ('Occitan', 'D'),
    'ões$': ('Portuguese', 'D'), '^vila': ('Galician', 'D'), 'ou$': ('Occitan', 'D'), '^nava': ('Pre-Roman nava', 'K'),
}


def km(lat1, lon1, lat2, lon2):
    return np.hypot((lat1 - lat2) * 110.57, (lon1 - lon2) * 111.32 * np.cos(np.radians(40)))


def load(key):
    t = pd.read_csv(os.path.join(D, f'scan_{key}.csv'), keep_default_na=False)
    for c in ['lat', 'lon', 'q', 'p', 'llr', 'units', 'rr_in', 'words_in']:
        t[c] = pd.to_numeric(t[c], errors='coerce')
    return t


def summary():
    rows = []
    for key in ['SI', 'SIo', 'Go', 'Mo', 'Io']:
        t = load(key)
        rows.append(dict(corpus=key, tested=len(t), significant=int((t.q < 0.05).sum()), share=round((t.q < 0.05).mean(), 3)))
    pd.DataFrame(rows).to_csv(os.path.join(D, 'scan_summary.csv'), index=False)
    print(pd.DataFrame(rows).to_string(index=False))


def validation():
    t = load('SI').set_index('element')
    rng = np.random.default_rng(0)
    tg = {m for ms, _, _ in TARGETS for m in ms}
    pool = t[~t.index.isin(tg)]
    rows, null = [], []
    for ms, exp, stratum in TARGETS:
        for m in ms:
            if m not in t.index:
                rows.append(dict(element=m, stratum=stratum, present=False, hit=False)); continue
            r = t.loc[m]
            hit = bool((r.q < 0.05) and (r.region in exp))
            rows.append(dict(element=m, stratum=stratum, present=True, hit=hit, q=r.q, region=r.region, lat=r.lat,
                             lon=r.lon, radius=r.radius, examples=r.examples))
            cand = pool[(pool.units >= 0.67 * r.units) & (pool.units <= 1.5 * r.units)]
            if len(cand):
                s = cand.iloc[rng.integers(0, len(cand), 30)]
                null += list((s.q < 0.05) & s.region.isin(exp))
    v = pd.DataFrame(rows); v.to_csv(os.path.join(D, 'validation.csv'), index=False, encoding='utf-8')
    pres = v[v.present]
    s = dict(targets=len(v), present=len(pres), hits=int(pres.hit.sum()), recall_present=round(pres.hit.mean(), 3),
             recall_all=round(v.hit.sum() / len(v), 3), random_recall=round(float(np.mean(null)), 3), random_n=len(null))
    pd.DataFrame([s]).to_csv(os.path.join(D, 'validation_summary.csv'), index=False)
    by = v.groupby('stratum').agg(n=('element', 'size'), present=('present', 'sum'), hits=('hit', 'sum'))
    by.to_csv(os.path.join(D, 'validation_by_stratum.csv'))
    print(s); print(by.to_string())


def locations():
    rows = []
    for key in ['SIo', 'Go', 'Mo', 'Io']:
        t = load(key); t['rank'] = np.arange(1, len(t) + 1)
        for m, regs in CHECKS:
            r = t[t.element == m]
            if not len(r):
                rows.append(dict(corpus=key, element=m, present=False)); continue
            r = r.iloc[0]
            dist = km(r.lat, r.lon, *TOLEDO)
            ok = bool(r.q < 0.05 and (r.region in regs if regs is not None else dist < 120))
            rows.append(dict(corpus=key, element=m, present=True, rank=int(r['rank']), of=len(t), q=r.q, region=r.region,
                             lat=round(r.lat, 2), lon=round(r.lon, 2), radius=r.radius, km_from_Toledo=round(dist),
                             in_expected_place=ok, examples=r.examples))
    p = pd.DataFrame(rows); p.to_csv(os.path.join(D, 'cluster_locations.csv'), index=False, encoding='utf-8')
    print(p[['corpus', 'element', 'q', 'region', 'lat', 'lon', 'radius', 'km_from_Toledo', 'in_expected_place']].to_string(index=False))


def merge_nested(t, maxkm=100):
    keep = []
    for _, r in t.sort_values('llr', ascending=False).iterrows():
        core = r.element.strip('^$'); pre = r.element.startswith('^')
        if not any(k.element.startswith('^') == pre and (core in k.element.strip('^$') or k.element.strip('^$') in core)
                   and km(r.lat, r.lon, k.lat, k.lon) < maxkm for k in keep):
            keep.append(r)
    return pd.DataFrame(keep)


def leads():
    t = load('SIo')
    m = merge_nested(t[t.q < 0.05]).reset_index(drop=True)
    m['rank_merged'] = np.arange(1, len(m) + 1)
    top = m.head(50).copy()
    top['family'] = [CLASS.get(x, ('?', 'U'))[0] for x in top.element]
    top['class'] = [CLASS.get(x, ('?', 'U'))[1] for x in top.element]
    top.to_csv(os.path.join(D, 'top50_classified.csv'), index=False, encoding='utf-8')
    summ = top.groupby('class').size().reindex(['K', 'D', 'U'], fill_value=0)
    fam = top.groupby('family').size().sort_values(ascending=False)
    pd.concat([summ.rename('n').to_frame().assign(kind='class'), fam.rename('n').to_frame().assign(kind='family')]).to_csv(
        os.path.join(D, 'top50_summary.csv'))
    print('classes:', summ.to_dict())
    # ranking by concentration and rarity of -eque among endings of its size
    s = t[(t.q < 0.05) & (t.words_in >= 8)].sort_values('rr_in', ascending=False).reset_index(drop=True)
    s['rank_rr'] = np.arange(1, len(s) + 1)
    s.head(100).to_csv(os.path.join(D, 'concentration_top100.csv'), index=False, encoding='utf-8')
    rr_eque = float(t[t.element == 'eque$'].rr_in.iloc[0])
    e = t[(~t.element.str.startswith('^')) & (t.element.str.len().between(4, 5)) & t.units.between(20, 60)]
    sharp = e[(e.q < 0.05) & (e.rr_in >= rr_eque)].sort_values('rr_in', ascending=False)
    sharp.to_csv(os.path.join(D, 'endings_as_concentrated_as_eque.csv'), index=False, encoding='utf-8')
    rar = dict(endings_3_4_letters_20_60_units=len(e), significant_share=round((e.q < 0.05).mean(), 3),
               as_concentrated_as_eque=len(sharp), rr_eque=rr_eque,
               eque_rank_by_concentration=int(s[s.element == 'eque$'].rank_rr.iloc[0]), significant_with_8_words=len(s))
    pd.DataFrame([rar]).to_csv(os.path.join(D, 'eque_rarity.csv'), index=False)
    print(rar)


def secondary(key, elements):
    rng = np.random.default_rng(1)
    d = paths.load_corpus(key, cols=['x', 'y', 'region', 'tok'])
    units, words, cx, cy, creg = SC.build(d)
    B = len(units)
    base = np.bincount(units.c.values, minlength=len(cx)).astype(float)
    A, meta = SC.windows(cx, cy, base, B)
    bw = np.asarray(A.T @ base).ravel()
    Acsc = A.tocsc()
    rows = []
    for m in elements:
        sub = SC.element_units(units, words, m)
        C = len(sub)
        cw = np.asarray(A.T @ np.bincount(sub.c.values, minlength=len(cx)).astype(float)).ravel()
        v = SC.llr(cw, bw, C, B); v[cw < 5] = 0
        found = []
        for j in np.argsort(-v)[:2000]:
            if v[j] <= 0:
                break
            cells = Acsc.indices[Acsc.indptr[j]:Acsc.indptr[j + 1]]
            if any(len(np.intersect1d(cells, f[1])) for f in found):
                continue
            inside = sub[np.isin(sub.c.values, cells)]
            if inside.wid.nunique() >= 5:
                found.append((j, cells, inside))
            if len(found) == 2:
                break
        nl = SC.null_table([C], base / B, A, bw, B, rng)[C]
        for rank, (j, cells, inside) in enumerate(found, 1):
            cc = int(meta[j, 0])
            p = (1 + (SC.NREP - np.searchsorted(nl, v[j], side='left'))) / (SC.NREP + 1)
            rows.append(dict(corpus=key, element=m, cluster=rank, units=C, lat=round(cy[cc] / 110.57, 2),
                             lon=round(cx[cc] / (111.32 * np.cos(np.radians(40))), 2), radius=int(meta[j, 1]),
                             c=int(cw[j]), llr=round(float(v[j]), 2), p=p, region=creg[cc],
                             examples=' '.join(words[inside.wid.value_counts().index[:10]])))
        if key == 'SIo' and m == 'eque$':
            # exact null without replacement for the main -eque cluster
            obs = float(v.max()); nrep = 4999; mx = np.zeros(nrep); rng2 = np.random.default_rng(2)
            AT = A.T.tocsr()
            for r in range(nrep):
                pick = rng2.choice(B, C, replace=False)
                cnt = np.bincount(units.c.values[pick], minlength=len(cx)).astype(float)
                cwr = AT @ cnt
                vv = SC.llr(cwr, bw, C, B); vv[cwr < 5] = 0
                mx[r] = vv.max()
            pd.DataFrame([dict(element='eque$', units=C, llr_observed=round(obs, 2), null_max=round(mx.max(), 2),
                               replicates=nrep, p=(1 + (mx >= obs).sum()) / (1 + nrep))]).to_csv(
                os.path.join(D, 'eque_exact_null.csv'), index=False)
    return pd.DataFrame(rows)


def main():
    summary(); validation(); locations(); leads()
    s = pd.concat([secondary('SIo', ['ena$', 'aque$', 'eque$', 'que$', 'ana$']), secondary('Io', ['ena$', 'aque$'])])
    s.to_csv(os.path.join(D, 'secondary_clusters.csv'), index=False, encoding='utf-8')
    print(s.to_string(index=False))


if __name__ == '__main__':
    main()
