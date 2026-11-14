"""Arabic place names and the northern frontier of al-Andalus.

Share of distinct names with an Arabic marker (lib.reference.STRICT, minus the loanwords in LEXAR; opaque words only)
in 15 km cells, numerator and denominator smoothed (σ = 2.3 cells), mainland Spain without the Basque provinces.
Levels: median share in the never-conquered north (Galicia, Asturias, Cantabria) and in the provinces conquered in or
after 1164. Frontier φ(λ): going north from 38.6° N at each longitude λ (−6.5° to −2.2°, steps of 0.1°), the first
latitude where the share falls below the threshold and stays below it for the next 50 km. Threshold at 25, 50 and 75 %
of the way from the northern to the southern level. Distances to the Duero and to a line through eight strongholds of
the Umayyad Middle March, by segment. Province-level association with the years of Islamic rule.
Output: results/arabic/*.csv
"""
import os, warnings
import numpy as np, pandas as pd
from scipy.ndimage import gaussian_filter
from scipy.stats import rankdata
from lib import paths
from lib.reference import STRICT, LEXAR, VASC, RIVER_LL, FORTS, CONQ

D = os.path.join(paths.RESULTS, 'arabic'); os.makedirs(D, exist_ok=True)
CELL = 15.0; SIG = 2.3
NORTH = {'ES-' + c for c in 'C LU OR PO O S'.split()}
SOUTH = {r for r, y in CONQ.items() if r.startswith('ES-') and y >= 1164}
LAMBDAS = np.round(np.arange(-6.5, -2.19, 0.1), 2)
SEGMENTS = {'west (6.5–5.1° W)': LAMBDAS <= -5.1, 'centre (5.0–3.5° W)': (LAMBDAS > -5.1) & (LAMBDAS <= -3.5),
            'east (3.0–2.2° W)': LAMBDAS >= -3.0}
C40 = np.cos(np.radians(40.0))


def surface(key):
    full = paths.load_corpus(key, cols=['lat', 'lon', 'x', 'y', 'region', 'tok'])
    opq = {t for tk in paths.load_corpus(key + 'o', cols=['tok']).tok for t in tk.split()}
    full = full[full.region.str.startswith('ES-') & ~full.region.isin(VASC)]
    rec = [(x, y, r, t, (t in opq) and bool(STRICT.match(t)) and not LEXAR.match(t))
           for x, y, r, tk in zip(full.x, full.y, full.region, full.tok) for t in tk.split()]
    d = pd.DataFrame(rec, columns=['x', 'y', 'r', 't', 'ar'])
    x0, y0 = d.x.min(), d.y.min()
    d['gx'] = ((d.x - x0) // CELL).astype(int); d['gy'] = ((d.y - y0) // CELL).astype(int)
    u = d.drop_duplicates(['gx', 'gy', 't'])
    W, H = u.gx.max() + 1, u.gy.max() + 1
    tot = np.zeros((H, W)); ar = np.zeros((H, W))
    np.add.at(tot, (u.gy, u.gx), 1); np.add.at(ar, (u.gy, u.gx), u.ar.astype(float))
    st, sa = gaussian_filter(tot, SIG), gaussian_filter(ar, SIG)
    rate = np.where(st > 3, sa / np.maximum(st, 1e-9) * 1000, np.nan)
    reg = d.groupby(['gy', 'gx']).r.agg(lambda s: s.value_counts().index[0])
    R = np.full((H, W), None, dtype=object)
    for (gy, gx), r in reg.items():
        R[gy, gx] = r
    lat = (np.arange(H) + 0.5) * CELL / 110.57 + y0 / 110.57
    lon = ((np.arange(W) + 0.5) * CELL + x0) / (111.32 * C40)
    north = np.nanmedian([rate[i, j] for i, j in zip(*np.nonzero(np.isin(R, list(NORTH))))])
    south = np.nanmedian([rate[i, j] for i, j in zip(*np.nonzero(np.isin(R, list(SOUTH))))])
    return dict(rate=rate, R=R, lat=lat, lon=lon, x0=x0, y0=y0, tot=tot, north=north, south=south)


def interp(points, lam):
    p = sorted(points, key=lambda t: t[1])
    return np.interp(lam, [q[1] for q in p], [q[0] for q in p])


def frontier(s, L):
    phi = []
    for lam in LAMBDAS:
        j = int(np.argmin(np.abs(s['lon'] - lam)))
        col = s['rate'][:, j]
        rows = [i for i in range(len(s['lat'])) if 38.6 <= s['lat'][i] <= 43.4 and not np.isnan(col[i])]
        f = np.nan
        for a, i in enumerate(rows):
            if col[i] < L:
                ahead = [k for k in rows[a + 1:] if s['lat'][k] <= s['lat'][i] + 50 / 110.57]
                if all(col[k] < L for k in ahead):
                    f = s['lat'][i]; break
        phi.append(f)
    return np.array(phi)


def dose_response(key):
    """Share of distinct Arabic-marked names by province against years of Islamic rule (and latitude)."""
    full = paths.load_corpus(key, cols=['lat', 'region', 'tok'])
    op = paths.load_corpus(key + 'o', cols=['region', 'tok'])
    rows = []
    for p in sorted(set(full.region) & set(CONQ)):
        if p in VASC:
            continue
        fp = full[full.region == p]
        allt = {t for tk in fp.tok for t in tk.split()}
        if len(allt) < 300:
            continue
        art = {t for tk in op[op.region == p].tok for t in tk.split() if STRICT.match(t) and not LEXAR.match(t)}
        rows.append((p, len(allt), len(art), CONQ[p] - 714, fp.lat.mean()))
    t = pd.DataFrame(rows, columns=['province', 'names', 'arabic', 'years', 'lat']).set_index('province')
    t['lograte'] = np.log((t.arabic + 1) / (t.names + 2))
    r = np.corrcoef(t.years, t.lograte)[0, 1]
    res = lambda a, b: a - np.polyval(np.polyfit(b, a, 1), b)
    pr = np.corrcoef(res(t.years.values, t.lat.values), res(t.lograte.values, t.lat.values))[0, 1]
    return t, dict(corpus=key, provinces=len(t), r_years=round(r, 3), r_years_latitude=round(np.corrcoef(t.years, t.lat)[0, 1], 3),
                   partial_r_given_latitude=round(pr, 3))


def main():
    duero = [(la, lo) for _, la, lo in RIVER_LL]
    forts = [(la, lo) for _, la, lo in FORTS]
    dl, fl = interp(duero, LAMBDAS), interp(forts, LAMBDAS)
    rows, seg_rows, prof = [], [], []
    for key in ['SI', 'I', 'G']:
        s = surface(key)
        np.savez_compressed(os.path.join(D, f'surface_{key}.npz'), rate=s['rate'], tot=s['tot'], lat=s['lat'],
                            lon=s['lon'], x0=s['x0'], y0=s['y0'], cell=CELL, north=s['north'], south=s['south'])
        for frac in [0.25, 0.5, 0.75]:
            L = s['north'] + frac * (s['south'] - s['north'])
            phi = frontier(s, L)
            pd.DataFrame(dict(lon=LAMBDAS, phi=phi, duero=dl, forts=fl)).to_csv(
                os.path.join(D, f'frontier_{key}_{int(frac * 100)}.csv'), index=False)
            west = (LAMBDAS >= -6.5) & (LAMBDAS <= -3.5); east = (LAMBDAS >= -3.0) & (LAMBDAS <= -2.2)
            if frac == 0.5:
                rows.append(dict(corpus=key, north_per_mille=round(s['north'], 2), south_per_mille=round(s['south'], 2),
                                 threshold=round(L, 2),
                                 west_centre_km_south_of_Duero=round(np.nanmean((dl - phi)[west]) * 110.57, 1),
                                 east_km_north_of_Duero=round(np.nanmean((phi - dl)[east]) * 110.57, 1),
                                 mean_km_to_strongholds=round(np.nanmean(np.abs(phi - fl)) * 110.57, 1),
                                 mean_km_to_Duero=round(np.nanmean(np.abs(phi - dl)) * 110.57, 1)))
            for sn, m in SEGMENTS.items():
                seg_rows.append(dict(corpus=key, threshold_fraction=frac, segment=sn,
                                     km_to_strongholds=round(np.nanmean(np.abs(phi - fl)[m]) * 110.57, 1),
                                     km_to_Duero=round(np.nanmean(np.abs(phi - dl)[m]) * 110.57, 1)))
        lat = s['lat']
        for name, (lo0, lo1) in {'west (6.5–5.1° W)': (-6.6, -5.1), 'west and centre (6.5–3.5° W)': (-6.6, -3.5),
                                 'east (3.0–2.2° W)': (-3.0, -2.2)}.items():
            js = np.nonzero((s['lon'] >= lo0) & (s['lon'] <= lo1))[0]
            with np.errstate(all='ignore'), warnings.catch_warnings():
                warnings.simplefilter('ignore', RuntimeWarning)
                pr = np.nanmean(np.where(s['tot'][:, js] > 0, s['rate'][:, js], np.nan), axis=1)
            for la, v in zip(lat, pr):
                if 38.6 <= la <= 43.4:
                    prof.append(dict(corpus=key, band=name, lat=round(la, 3), per_mille=round(v, 2) if np.isfinite(v) else np.nan))
    pd.DataFrame(rows).to_csv(os.path.join(D, 'frontier_summary.csv'), index=False)
    pd.DataFrame(seg_rows).to_csv(os.path.join(D, 'frontier_segments.csv'), index=False)
    pd.DataFrame(prof).to_csv(os.path.join(D, 'profiles.csv'), index=False)
    pd.set_option('display.width', 250)
    print(pd.DataFrame(rows).T.to_string()); print(pd.DataFrame(seg_rows).to_string(index=False))
    dr, tabs = [], []
    for key in ['S', 'G', 'M']:
        t, s = dose_response(key); dr.append(s); tabs.append(t.assign(corpus=key))
    pd.DataFrame(dr).to_csv(os.path.join(D, 'dose_response.csv'), index=False)
    pd.concat(tabs).to_csv(os.path.join(D, 'dose_response_provinces.csv'))
    print(pd.DataFrame(dr).to_string(index=False))


if __name__ == '__main__':
    main()
