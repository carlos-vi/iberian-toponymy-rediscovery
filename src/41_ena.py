"""The suffix -ena (Marchena, Lucena) in the south.

Types = distinct (opaque word, province) pairs. Lift = share of the suffix's types that fall in core Baetica (Seville,
Córdoba, Huelva, Cádiz, Málaga, Jaén) / share of all types there; 95 % bootstrap intervals over the suffix's types.
Computed with all provinces and without the Basque provinces (where -ena is the Basque genitive -(r)ena of farmhouse
names). Share of -ena among -ena + -ana by groups of provinces; Spearman correlation with the years of Islamic rule,
raw and controlling for latitude (partial rank correlation); repeated names.
Output: results/ena/*.csv
"""
import os
import numpy as np, pandas as pd
from scipy.stats import spearmanr, fisher_exact, rankdata
from lib import paths
from lib.reference import CONQ

D = os.path.join(paths.RESULTS, 'ena'); os.makedirs(D, exist_ok=True)
BAET = {'ES-' + c for c in 'SE CO H CA MA J'.split()}
CORE = {'ES-' + c for c in 'SE CO H CA'.split()}
NASR = {'ES-' + c for c in 'GR AL MA'.split()}
VASC = {'ES-' + c for c in 'SS BI VI NA'.split()}
GROUPS = {'Core Baetica (SE CO H CA)': CORE, 'Nasrid kingdom (GR AL MA)': NASR, 'Jaén (J)': {'ES-J'},
          'Ebro valley and Teruel (Z HU TE)': {'ES-Z', 'ES-HU', 'ES-TE'}, 'Murcia (MU)': {'ES-MU'},
          'Never conquered north (C LU OR PO O S)': {'ES-' + c for c in 'C LU OR PO O S'.split()},
          'Duero basin (LE ZA VA P BU SG AV SO SA)': {'ES-' + c for c in 'LE ZA VA P BU SG AV SO SA'.split()}}
CORPORA = ['So', 'Go', 'Mo', 'Io', 'SIo']
rng = np.random.default_rng(0)


def types(key):
    d = paths.load_corpus(key, cols=['x', 'y', 'lat', 'region', 'tok'])
    d = d[~d.region.str.startswith('FR')]
    T = pd.DataFrame([(t, r, x, y, la) for tk, r, x, y, la in zip(d.tok, d.region, d.x, d.y, d.lat) for t in tk.split()],
                     columns=['t', 'r', 'x', 'y', 'lat'])
    return T.groupby(['t', 'r'], as_index=False).agg(x=('x', 'mean'), y=('y', 'mean'), lat=('lat', 'mean'))


def with_suffix(T, suf):
    return T[T.t.str.endswith(suf) & (T.t.str.len() >= len(suf) + 2)]


def lift(T, suf, B=2000):
    s = with_suffix(T, suf)
    base = T.r.isin(BAET).mean()
    inb = s.r.isin(BAET).values
    bt = [inb[rng.integers(0, len(inb), len(inb))].mean() / base for _ in range(B)]
    lo, hi = np.percentile(bt, [2.5, 97.5])
    return dict(n=len(s), in_baetica=int(inb.sum()), lift=inb.mean() / base, lo=lo, hi=hi)


def partial_rank(y, x, z):
    ry, rx, rz = rankdata(y), rankdata(x), rankdata(z)
    res = lambda a, b: a - np.polyval(np.polyfit(b, a, 1), b)
    return float(np.corrcoef(res(ry, rz), res(rx, rz))[0, 1])


def province_table(T):
    T = T[T.r.str.startswith('ES-')]
    lat = T.groupby('r').lat.mean()
    e = with_suffix(T, 'ena').groupby('r').size(); a = with_suffix(T, 'ana').groupby('r').size()
    p = pd.DataFrame({'ena': e, 'ana': a}).fillna(0).astype(int)
    p['share'] = p.ena / (p.ena + p.ana)
    p['years'] = [CONQ.get(r, np.nan) - 714 for r in p.index]
    p['lat'] = lat.reindex(p.index)
    return p[(p.ena + p.ana) >= 5].dropna()


def main():
    store, rows = {}, []
    for key in CORPORA:
        T = types(key); store[key] = T
        e, a = lift(T, 'ena'), lift(T, 'ana')
        rows.append(dict(corpus=key, ena_types=e['n'], ena_in_baetica=e['in_baetica'], ena_lift=round(e['lift'], 2),
                         ena_lo=round(e['lo'], 2), ena_hi=round(e['hi'], 2), ana_types=a['n'], ana_in_baetica=a['in_baetica'],
                         ana_lift=round(a['lift'], 2), ana_lo=round(a['lo'], 2), ana_hi=round(a['hi'], 2)))
    lifts = pd.DataFrame(rows); lifts.to_csv(os.path.join(D, 'lift_all_provinces.csv'), index=False)
    pd.set_option('display.width', 250); print(lifts.to_string(index=False))
    # Basque homograph
    I = with_suffix(store['Io'], 'ena')
    basque = dict(ena_types_I=len(I), in_bizkaia_gipuzkoa=int(I.r.isin({'ES-BI', 'ES-SS'}).sum()),
                  examples=' '.join(I[I.r == 'ES-BI'].t.head(12)))
    pd.DataFrame([basque]).to_csv(os.path.join(D, 'basque_homograph.csv'), index=False); print(basque)
    # by province in S+I (all provinces): years of rule, Fisher core Baetica vs Nasrid
    p = province_table(store['SIo'])
    p.sort_values('share', ascending=False).to_csv(os.path.join(D, 'provinces_SI.csv'))
    core = p.loc[p.index.isin(CORE)].sum(); nas = p.loc[p.index.isin(NASR)].sum()
    odds, pf = fisher_exact([[core.ena, core.ana], [nas.ena, nas.ana]])
    pd.DataFrame([dict(core_ena=int(core.ena), core_ana=int(core.ana), nasrid_ena=int(nas.ena), nasrid_ana=int(nas.ana),
                       fisher_odds=round(odds, 3), fisher_p=round(pf, 4))]).to_csv(os.path.join(D, 'core_vs_nasrid.csv'), index=False)
    # without the Basque provinces
    rows = []
    for key in CORPORA:
        T = store[key]; T = T[~T.r.isin(VASC)]
        e, a = lift(T, 'ena'), lift(T, 'ana')
        r = dict(corpus=key, ena_lift=round(e['lift'], 2), ena_lo=round(e['lo'], 2), ena_hi=round(e['hi'], 2),
                 ana_lift=round(a['lift'], 2), ana_lo=round(a['lo'], 2), ana_hi=round(a['hi'], 2))
        es, an = with_suffix(T, 'ena'), with_suffix(T, 'ana')
        for g, st in GROUPS.items():
            ne, na = int(es.r.isin(st).sum()), int(an.r.isin(st).sum())
            r[f'{g} ena'] = ne; r[f'{g} ana'] = na; r[f'{g} share'] = round(ne / max(ne + na, 1), 2)
        pt = province_table(T)
        rho, pv = spearmanr(pt.years, pt.share)
        r.update(spearman_years=round(rho, 3), spearman_p=round(pv, 4), n_provinces=len(pt),
                 spearman_years_latitude=round(spearmanr(pt.years, pt.lat)[0], 3),
                 partial_years_given_latitude=round(partial_rank(pt.share, pt.years, pt.lat), 3))
        rows.append(r)
    nb = pd.DataFrame(rows); nb.to_csv(os.path.join(D, 'without_basque.csv'), index=False)
    print(nb.T.to_string())
    # repeated names (S+I, without Basque provinces)
    T = store['SIo']; T = T[~T.r.isin(VASC) & T.r.str.startswith('ES-')]
    es = with_suffix(T, 'ena')
    rep = es.groupby('t').r.nunique().sort_values(ascending=False)
    rep.rename('provinces').to_csv(os.path.join(D, 'repeated_names_SI.csv'))
    pd.DataFrame([dict(types=len(es), distinct_words=es.t.nunique())]).to_csv(os.path.join(D, 'repeated_summary_SI.csv'), index=False)
    print('types', len(es), 'distinct words', es.t.nunique()); print(rep.head(8).to_string())
    store['SIo'].to_csv(os.path.join(D, 'types_SIo.csv'), index=False, encoding='utf-8')


if __name__ == '__main__':
    main()
