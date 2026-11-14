"""The Castilian wedge: two isoglosses across northern Spain.

1. -iello / (-iello + -illo): last word of each name (first form), in 0.5° bands of longitude within 41.5–43.6° N;
   segments: west of 5° W and the Castilian strip 4.5–1.5° W; provinces west of 5° W.
2. Latin F- kept (f-) or turned into h-: eleven matched lexical pairs on all words (≥ 3 letters, generic words
   included) of full names, counted as distinct (word, 20 km cell) types; Basque provinces excluded (Basque has an
   unrelated initial h-). Segments in the band 41.5–43.6° N: west of 5.5° W, 4.5–2.0° W and east of 1.0° W; and Spain
   south of 40° N outside the Catalan-speaking provinces. Each pair is also dropped in turn.
Output: results/wedge/*.csv
"""
import os, re
import numpy as np, pandas as pd
from lib import paths
from lib.normalize import norm, raw_last, TOK

D = os.path.join(paths.RESULTS, 'wedge'); os.makedirs(D, exist_ok=True)
OLD = re.compile(r"(?:iello|iella|iellos|iellas|iellu|iecho|ieḷḷu|ieyu)$")
NEW = re.compile(r"(?:illo|illa|illos|illas)$")
VASC = {'ES-SS', 'ES-BI', 'ES-VI', 'ES-NA'}
CATAL = {'ES-CS', 'ES-V', 'ES-A', 'ES-PM'}
FONT_EXCL = {'fonte', 'fontes', 'font', 'fonts'}
PAIRS = [('FONTE', r'^font[aeiou]\w+', r'^hont'), ('FAUCE', r'^fo[zc](es)?$', r'^ho[zc](es)?$'),
         ('FURNU', r'^forn', r'^horn'), ('FERRU', r'^ferr', r'^herr'), ('FICU', r'^figu', r'^higu'),
         ('FENU', r'^fen(ar|al|ares|ales)$', r'^hen(ar|al|ares|ales)$'), ('FOVEA', r'^fo[yij][ao]s?$', r'^ho[yi][ao]s?$'),
         ('FURCA', r'^forc', r'^horc'), ('FILICTU', r'^fe[lt][e]?ch|^felg', r'^hel[e]?ch|^helg'),
         ('FUNDU', r'^fond', r'^hond'), ('FENUCULU', r'^fen[oe]ll|^fenoj|^finoj|^funch', r'^hinoj')]
PAIRS = [(n, re.compile(f), re.compile(h)) for n, f, h in PAIRS]
LONB = np.arange(-9.5, 3.51, 0.5)
SEGS = {'west': lambda t: t.lat.between(41.5, 43.6) & (t.lon < -5.5),
        'centre': lambda t: t.lat.between(41.5, 43.6) & t.lon.between(-4.5, -2.0),
        'east': lambda t: t.lat.between(41.5, 43.6) & (t.lon > -1.0),
        'east_spain': lambda t: t.lat.between(41.5, 43.6) & (t.lon > -1.0) & ~t.region.str.startswith('FR'),
        'south': lambda t: (t.lat < 40) & t.region.str.startswith('ES-') & ~t.region.isin(CATAL)}


def iello(key):
    d = paths.load_corpus(key, cols=['lat', 'lon', 'name', 'region'])
    w = d[d.lat.between(41.5, 43.6) & d.lon.between(-7.5, 1.0) & ~d.region.str.startswith('FR')].copy()
    w['last'] = w.name.map(raw_last)
    w['old'] = w['last'].map(lambda s: bool(OLD.search(s))); w['new'] = w['last'].map(lambda s: bool(NEW.search(s)))
    g = w.groupby(pd.cut(w.lon, np.arange(-7.5, 1.01, 0.5)), observed=False).agg(old=('old', 'sum'), new=('new', 'sum'))
    g['share'] = g.old / (g.old + g.new).replace(0, np.nan)
    g['lon'] = [iv.mid for iv in g.index]; g['corpus'] = key
    seg = lambda s: (int(s.old.sum()), int(s.new.sum()), s.old.sum() / max(s.old.sum() + s.new.sum(), 1))
    west, mid = seg(w[w.lon < -5]), seg(w[w.lon.between(-4.5, -1.5)])
    prov = w[(w.lon < -5) & w.region.str.startswith('ES')].groupby('region')[['old', 'new']].sum()
    prov = prov[(prov.old + prov.new) > 0].assign(corpus=key, share=lambda x: (x.old / (x.old + x.new)).round(3))
    return g.reset_index(drop=True), west, mid, prov


def classify(word, pairs):
    for n, f, h in pairs:
        if n == 'FONTE' and word in FONT_EXCL:
            continue
        if f.search(word):
            return n, 'f'
        if h.search(word):
            return n, 'h'
    return None


def fh_types(key):
    d = paths.load_corpus(key, cols=['lat', 'lon', 'x', 'y', 'name', 'region'])
    d = d[~d.region.isin(VASC)]
    rows = []
    for lat, lon, x, y, name, reg in zip(d.lat, d.lon, d.x, d.y, d.name, d.region):
        for form in re.split(r"\s*/\s*", str(name)):
            for w in TOK.findall(norm(form)):
                if len(w) >= 3:
                    c = classify(w, PAIRS)
                    if c:
                        rows.append((w, int(x // 20), int(y // 20), lat, lon, reg, c[0], c[1]))
    t = pd.DataFrame(rows, columns=['w', 'cx', 'cy', 'lat', 'lon', 'region', 'pair', 'fh'])
    return t.drop_duplicates(['w', 'cx', 'cy']).reset_index(drop=True)


def segments(t, drop=None):
    if drop:
        t = t[t.pair != drop]
    out = {}
    for s, f in SEGS.items():
        x = t[f(t)]
        out[s] = (len(x), (x.fh == 'f').mean() if len(x) else np.nan)
    return out


def main():
    rows, segs, provs = [], [], []
    for key in ['S', 'M', 'G', 'I']:
        g, west, mid, prov = iello(key)
        rows.append(g); provs.append(prov)
        segs.append(dict(corpus=key, west_old=west[0], west_new=west[1], west_share=round(west[2], 3),
                         castile_old=mid[0], castile_new=mid[1], castile_share=round(mid[2], 4)))
    pd.concat(rows).to_csv(os.path.join(D, 'iello_transect.csv'), index=False)
    pd.DataFrame(segs).to_csv(os.path.join(D, 'iello_segments.csv'), index=False)
    pd.concat(provs).to_csv(os.path.join(D, 'iello_west_by_province.csv'))
    print(pd.DataFrame(segs).to_string(index=False)); print(pd.concat(provs).to_string())
    res, loo, words = [], [], []
    for key in ['S', 'M', 'G', 'I']:
        t = fh_types(key)
        t.to_csv(os.path.join(D, f'fh_types_{key}.csv'), index=False, encoding='utf-8')
        sg = segments(t)
        res.append(dict(corpus=key, **{f'{s}_n': n for s, (n, _) in sg.items()},
                        **{f'{s}_share_f': round(v, 3) for s, (_, v) in sg.items()}))
        for n, _, _ in PAIRS:
            sd = segments(t, drop=n)
            loo.append(dict(corpus=key, dropped=n, **{f'{s}_share_f': round(v, 3) for s, (_, v) in sd.items()},
                            **{f'{s}_n': m for s, (m, _) in sd.items()}))
        band = t[t.lat.between(41.5, 43.6)]
        g = band.groupby(pd.cut(band.lon, LONB), observed=False).fh.agg(n='size', f=lambda s: (s == 'f').sum())
        g['share'] = g.f / g.n.replace(0, np.nan); g['lon'] = [iv.mid for iv in g.index]; g['corpus'] = key
        g.reset_index(drop=True).to_csv(os.path.join(D, f'fh_transect_{key}.csv'), index=False)
        for s, f in SEGS.items():
            x = t[f(t)]
            for fh in 'fh':
                for w, n in x[x.fh == fh].w.value_counts().head(12).items():
                    words.append(dict(corpus=key, segment=s, form=fh, word=w, types=n))
        pc = t.assign(segment=np.select([SEGS[s](t) for s in ['west', 'centre', 'east', 'south']],
                                        ['west', 'centre', 'east', 'south'], 'other'))
        pc.groupby(['pair', 'segment', 'fh']).size().unstack(['segment', 'fh'], fill_value=0).to_csv(
            os.path.join(D, f'fh_pairs_{key}.csv'))
    r = pd.DataFrame(res); r.to_csv(os.path.join(D, 'fh_segments.csv'), index=False)
    l = pd.DataFrame(loo); l.to_csv(os.path.join(D, 'fh_leave_one_pair_out.csv'), index=False)
    pd.DataFrame(words).to_csv(os.path.join(D, 'fh_frequent_words.csv'), index=False, encoding='utf-8')
    worst = l.groupby('corpus').agg(west_min=('west_share_f', 'min'), east_min=('east_share_f', 'min'),
                                    centre_max=('centre_share_f', 'max'), south_max=('south_share_f', 'max'))
    worst.to_csv(os.path.join(D, 'fh_leave_one_pair_out_extremes.csv'))
    pd.set_option('display.width', 250)
    print(r.T.to_string()); print(worst.to_string())


if __name__ == '__main__':
    main()
