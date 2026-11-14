"""Download the official populated places of Spain from the OGC Features API of the Instituto Geográfico Nacional
(CC BY 4.0): collections 'nuc' (population nuclei of the INE) and 'opob' (other populated places).

Province by province (INE code cpro), following the 'next' links, without geometries (the properties include latitude and
longitude). Each province is checked against the number of records announced by the service (numberMatched).
Output: data/raw/ign/ign_{nuc,opob}.csv and data/raw/ign/download.log. The snapshot used in the paper was downloaded
on 24 October 2026.
Usage: python 03_download_ign.py [--refresh]
"""
import json, os, sys, time, urllib.request
import pandas as pd
from lib import paths
from lib.reference import INE2GN

BASE = 'https://api-features.ign.es/collections'
COLS = ['id', 'nombre', 'tipo', 'cpro', 'codine', 'habitantes', 'ine', 'nucrural', 'aislada', 'capital', 'latitud',
        'longitud', 'altitud']


def fetch(url):
    req = urllib.request.Request(url, headers={'User-Agent': 'iberian-toponymy-research/1.0'})
    for attempt in range(5):
        try:
            with urllib.request.urlopen(req, timeout=300) as r:
                return json.load(r)
        except Exception as e:
            print('retry', attempt, url[-80:], e, flush=True)
            time.sleep(10 * (attempt + 1))
    raise SystemExit('download failed: ' + url)


def province(coll, cpro):
    url = f'{BASE}/{coll}/items?f=json&limit=10000&skipGeometry=true&cpro={cpro}'
    rows, matched = [], None
    while url:
        d = fetch(url)
        matched = d.get('numberMatched', matched)
        for f in d['features']:
            p = f['properties']
            rows.append([f.get('id')] + [p.get(c) for c in COLS[1:]])
        url = next((l['href'] for l in d.get('links', []) if l.get('rel') == 'next'), None)
        if d.get('numberReturned', 0) == 0:
            url = None
        time.sleep(0.5)
    return rows, matched


def main(refresh=False):
    dest = paths.raw('ign')
    os.makedirs(dest, exist_ok=True)
    if all(os.path.exists(os.path.join(dest, f'ign_{c}.csv')) for c in ['nuc', 'opob']) and not refresh:
        print('exists, skipped:', dest); return
    log = open(os.path.join(dest, 'download.log'), 'w', encoding='utf-8')
    for coll in ['nuc', 'opob']:
        allrows = []
        for cpro in sorted(INE2GN):
            rows, matched = province(coll, cpro)
            msg = f'{coll} cpro={cpro:02d} {INE2GN[cpro]:>2}: downloaded {len(rows)} (distinct ids {len({r[0] for r in rows})}) of {matched}'
            print(msg, flush=True); log.write(msg + '\n'); log.flush()
            allrows += rows
        d = pd.DataFrame(allrows, columns=COLS).drop_duplicates('id')
        d.to_csv(os.path.join(dest, f'ign_{coll}.csv'), index=False, encoding='utf-8')
        log.write(f'{coll}: total {len(d)}\n')
    log.close()


if __name__ == '__main__':
    main('--refresh' in sys.argv)
