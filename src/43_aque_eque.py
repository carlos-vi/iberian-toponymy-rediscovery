"""Opaque place names in -aque/-eque in the old kingdom of Toledo.

Types = distinct (opaque word, province) pairs, mainland Spain and Portugal. Centre = Toledo, Ciudad Real, Cuenca,
Albacete, Guadalajara and Madrid; south = Andalusia and Murcia.
1. Lift: share of the -aque/-eque types (words of ≥ 6 letters) that fall in the centre / share of all types there.
2. Vowel before -que: a/e against i/u, centre against south, one-sided Fisher exact test.
Both are also computed on distinct words (a word counts once in an area however many provinces it occurs in).
Output: results/aque_eque/*.csv
"""
import os, re
import numpy as np, pandas as pd
from scipy.stats import fisher_exact
from lib import paths

D = os.path.join(paths.RESULTS, 'aque_eque'); os.makedirs(D, exist_ok=True)
CENTER = {'ES-' + c for c in 'TO CR CU AB GU M'.split()}
SOUTH = {'ES-' + c for c in 'AL CA CO GR H J MA SE MU'.split()}
VQUE = re.compile(r'([aeiou])que$')
CORPORA = ['So', 'Go', 'Mo', 'Io', 'SIo']


def types(key):
    d = paths.load_corpus(key, cols=['lat', 'lon', 'region', 'tok', 'name'])
    d = d[~d.region.str.startswith('FR')]
    rows = [(t, r, la, lo, n) for tk, r, la, lo, n in zip(d.tok, d.region, d.lat, d.lon, d.name) for t in tk.split()]
    T = pd.DataFrame(rows, columns=['t', 'r', 'lat', 'lon', 'name'])
    return T.groupby(['t', 'r'], as_index=False).agg(lat=('lat', 'mean'), lon=('lon', 'mean'), records=('name', 'size'),
                                                     example=('name', 'first'))


def fisher(c_ae, c_iu, s_ae, s_iu):
    odds, p = fisher_exact([[c_ae, c_iu], [s_ae, s_iu]], alternative='greater')
    return (round(odds, 2) if np.isfinite(odds) else odds), p


def analyse(T):
    T = T.assign(v=T.t.str.extract(VQUE)[0])
    q = T[T.v.isin(['a', 'e']) & (T.t.str.len() >= 6)]
    V = T[T.v.isin(list('aeiu')) & (T.t.str.len() >= 6)]
    # types
    base = T.r.isin(CENTER).mean()
    lift = q.r.isin(CENTER).mean() / base if len(q) else np.nan
    n = {f'{a}_{g}': int((V.r.isin(A) & V.v.isin(vs)).sum())
         for a, A in [('centre', CENTER), ('south', SOUTH)] for g, vs in [('ae', ['a', 'e']), ('iu', ['i', 'u'])]}
    odds, p = fisher(n['centre_ae'], n['centre_iu'], n['south_ae'], n['south_iu'])
    # distinct words
    W = T.assign(c=T.r.isin(CENTER)).groupby('t').agg(c=('c', 'any'), v=('v', 'first'))
    qw = W[W.v.isin(['a', 'e']) & (W.index.str.len() >= 6)]
    lift_w = qw.c.mean() / W.c.mean() if len(qw) else np.nan
    nw = {f'{a}_{g}': V[V.r.isin(A) & V.v.isin(vs)].t.nunique()
          for a, A in [('centre', CENTER), ('south', SOUTH)] for g, vs in [('ae', ['a', 'e']), ('iu', ['i', 'u'])]}
    odds_w, p_w = fisher(nw['centre_ae'], nw['centre_iu'], nw['south_ae'], nw['south_iu'])
    res = dict(types_aque_eque=len(q), in_centre=int(q.r.isin(CENTER).sum()), centre_share_all=round(base, 4),
               lift=round(lift, 2), **n, fisher_odds=odds, fisher_p=round(p, 5),
               words_aque_eque=len(qw), words_in_centre=int(qw.c.sum()), lift_words=round(lift_w, 2),
               **{f'words_{k}': v for k, v in nw.items()}, fisher_odds_words=odds_w, fisher_p_words=round(p_w, 5))
    return res, q, V


def main():
    rows, keep = [], {}
    for key in CORPORA:
        r, q, V = analyse(types(key))
        rows.append(dict(corpus=key, **r)); keep[key] = (q, V)
    res = pd.DataFrame(rows)
    res.to_csv(os.path.join(D, 'tests.csv'), index=False)
    pd.set_option('display.width', 250); pd.set_option('display.max_colwidth', 60)
    print(res.T.to_string())
    q, V = keep['SIo']
    V.sort_values(['v', 'r']).to_csv(os.path.join(D, 'vque_SIo.csv'), index=False, encoding='utf-8')
    q.sort_values(['r', 't']).to_csv(os.path.join(D, 'aque_eque_names_SIo.csv'), index=False, encoding='utf-8')
    print(q.sort_values('r')[['t', 'r', 'records', 'example']].to_string(index=False))


if __name__ == '__main__':
    main()
