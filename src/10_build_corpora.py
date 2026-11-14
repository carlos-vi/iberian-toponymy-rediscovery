"""Build the corpora from the raw sources.

S   OpenStreetMap settlements (place = city … farm, except locality) in Spain, Portugal and Andorra, plus GeoNames
    populated places (class P) of south-western France (42.3–45.3° N, 2.0° W – 3.3° E).
M   OpenStreetMap localities (place = locality; field names).
G   GeoNames, classes P and A, Spain, Portugal and Andorra.
I   Official populated places of Spain (IGN/INE: population nuclei and other populated places).
SI  S and I merged.
Every name is reduced to its words (lib.normalize.tokens); records without words are dropped and identical word sets in
the same cell of a 2 km grid are collapsed. Region (province, Portuguese district, French department): from the nearest
GeoNames P/A record (OSM, GeoNames) or from the INE province code (I). Opaque versions (suffix 'o') keep only the words
that are not words of a living language (lib.opacity).
Output: data/derived/corpus_{S,M,G,I,SI}{,o}.csv and results/data/corpus_density.csv
"""
import numpy as np, pandas as pd
from scipy.spatial import cKDTree
from lib import paths
from lib.normalize import tokens, xy
from lib import opacity
from lib.reference import INE2GN, MACRO, AREA_KM2, macro

GN_COLS = ['gid', 'name', 'ascii', 'alt', 'lat', 'lon', 'fclass', 'fcode', 'cc', 'cc2', 'a1', 'a2', 'a3', 'a4', 'pop',
           'elev', 'dem', 'tz', 'mod']


def load_geonames(cc):
    d = pd.read_csv(paths.raw('geonames', f'{cc}.txt'), sep='\t', header=None, names=GN_COLS, dtype=str, quoting=3,
                    keep_default_na=False, encoding='utf-8')
    d['lat'] = d.lat.astype(float); d['lon'] = d.lon.astype(float)
    return d


def dedup(d):
    d = d.assign(cx=(d.x // 2).astype(int), cy=(d.y // 2).astype(int))
    return d.drop_duplicates(['tok', 'cx', 'cy']).drop(columns=['cx', 'cy']).reset_index(drop=True)


def osm_geonames_corpora():
    gn = pd.concat([load_geonames(c) for c in ['ES', 'PT', 'AD', 'FR']], ignore_index=True)

    def region(r):
        if r.cc == 'ES': return 'ES-' + r.a2 if r.a2 else ''
        if r.cc == 'PT': return 'PT-' + r.a1 if r.a1 else ''
        if r.cc == 'AD': return 'AD'
        return 'FR-' + r.a2 if r.a2 else ''
    gp = gn[gn.fclass.isin(['P', 'A'])].copy()
    gp['region'] = gp.apply(region, axis=1)
    gp = gp[(gp.region != '') & ~gp.region.str.endswith('-00')]
    gx, gy = xy(gp.lat.values, gp.lon.values)
    isfr = gp.region.str.startswith('FR-').values
    trees = {f: (cKDTree(np.c_[gx[isfr == f], gy[isfr == f]]), gp.region.values[isfr == f]) for f in (False, True)}
    osm = pd.read_csv(paths.raw('osm', 'osm_places_iberia.tsv'), sep='\t', quoting=3, dtype=str, keep_default_na=False,
                      encoding='utf-8')
    osm.columns = ['type', 'id', 'lat', 'lon', 'name', 'place', 'name_es', 'name_pt']
    osm['lat'] = osm.lat.astype(float); osm['lon'] = osm.lon.astype(float)
    osm = osm[(osm.lat > 35.5) & (osm.lat < 44.2) & (osm.lon > -10) & (osm.lon < 4.6)]   # mainland and Balearics
    osm['src'] = 'osm'
    fr = gn[(gn.cc == 'FR') & (gn.fclass == 'P') & gn.lat.between(42.3, 45.3) & gn.lon.between(-2.0, 3.3)].copy()
    fr['place'] = 'gn_fr'; fr['src'] = 'gnfr'
    g = gn[gn.cc.isin(['ES', 'PT', 'AD']) & gn.fclass.isin(['P', 'A']) & gn.lat.between(35.5, 44.2) & gn.lon.between(-10, 4.6)].copy()
    g['place'] = g.fcode; g['src'] = 'gn'
    cols = ['lat', 'lon', 'name', 'place', 'src']
    corpora = {'S': pd.concat([osm[osm.place != 'locality'][cols], fr[cols]], ignore_index=True),
               'M': osm[osm.place == 'locality'][cols].reset_index(drop=True),
               'G': g[cols].reset_index(drop=True)}
    out = {}
    for k, d in corpora.items():
        x, y = xy(d.lat.values, d.lon.values)
        d['x'] = x; d['y'] = y
        reg = np.empty(len(d), dtype=object)
        for f in (False, True):
            sel = (d.src.values == 'gnfr') == f
            if sel.any():
                tr, ro = trees[f]
                _, idx = tr.query(np.c_[x[sel], y[sel]])
                reg[sel] = ro[idx]
        d['region'] = reg
        d['tok'] = d.name.map(lambda n: ' '.join(tokens(n)))
        out[k] = dedup(d[d.tok != ''])
        print(k, len(d), '->', len(out[k]), 'records', flush=True)
    return out


def ign_corpus():
    raw = pd.concat([pd.read_csv(paths.raw('ign', f'ign_{c}.csv'), keep_default_na=False, encoding='utf-8')
                     .assign(src='ign_' + c) for c in ['nuc', 'opob']], ignore_index=True)
    raw = raw[raw.cpro.isin(INE2GN.keys())]
    raw = raw[(raw.latitud > 35.5) & (raw.latitud < 44.2) & (raw.longitud > -10) & (raw.longitud < 4.6)]
    I = pd.DataFrame(dict(lat=raw.latitud.values, lon=raw.longitud.values, name=raw.nombre.values, place=raw.tipo.values,
                          src=raw.src.values))
    I['x'], I['y'] = xy(I.lat.values, I.lon.values)
    I['region'] = ['ES-' + INE2GN[c] for c in raw.cpro.values]
    I['tok'] = I.name.map(lambda n: ' '.join(tokens(n)))
    n0 = len(I)
    I = dedup(I[I.tok != ''])
    print('I', n0, '->', len(I), 'records', flush=True)
    return I


def opaque(d, cache):
    types = sorted({t for tk in d.tok for t in tk.split()} - set(cache))
    cache.update({t: opacity.transparent(t) for t in types})
    op = d.tok.map(lambda tk: ' '.join(t for t in tk.split() if not cache[t]))
    do = d.assign(tok=op)
    return do[do.tok != ''].reset_index(drop=True)


def main():
    C = osm_geonames_corpora()
    C['I'] = ign_corpus()
    C['SI'] = dedup(pd.concat([C['S'], C['I']], ignore_index=True))
    print('SI', len(C['SI']), 'records', flush=True)
    for k, d in C.items():
        d.to_csv(paths.derived(f'corpus_{k}.csv'), index=False, encoding='utf-8')
    cache = {}
    for k in ['S', 'M', 'G', 'I', 'SI']:
        do = opaque(C[k], cache)
        do.to_csv(paths.derived(f'corpus_{k}o.csv'), index=False, encoding='utf-8')
        print(f'{k}o: {len(do)} of {len(C[k])} records keep at least one opaque word', flush=True)
    rows = {k: C[k].region.map(macro).value_counts() for k in ['S', 'I', 'SI', 'M', 'G']}
    t = pd.DataFrame(rows).fillna(0).astype(int)
    for k in list(rows):
        t[f'{k}_per_1000km2'] = [round(t.loc[g, k] / AREA_KM2[g] * 1000, 1) if g in AREA_KM2 else np.nan for g in t.index]
    t = t.sort_values('S', ascending=False)
    t.loc['TOTAL', ['S', 'I', 'SI', 'M', 'G']] = [len(C[k]) for k in ['S', 'I', 'SI', 'M', 'G']]
    t.to_csv(paths.results('data', 'corpus_density.csv'), encoding='utf-8')
    pd.set_option('display.width', 250)
    print(t.to_string())


if __name__ == '__main__':
    main()
