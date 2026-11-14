"""Download the Natural Earth 1:10m vector layers used for the base maps (public domain): land, rivers, lakes,
first-order administrative boundaries and international boundaries.
Output: data/raw/naturalearth/*.geojson
Usage: python 05_download_naturalearth.py [--refresh]
"""
import os, sys, urllib.request
from lib import paths

URL = 'https://raw.githubusercontent.com/nvkelso/natural-earth-vector/master/geojson/{name}.geojson'
NAMES = ['ne_10m_land', 'ne_10m_rivers_lake_centerlines', 'ne_10m_rivers_europe', 'ne_10m_lakes',
         'ne_10m_admin_1_states_provinces_lines', 'ne_10m_admin_0_boundary_lines_land']


def main(refresh=False):
    os.makedirs(paths.raw('naturalearth'), exist_ok=True)
    for name in NAMES:
        out = paths.raw('naturalearth', f'{name}.geojson')
        if os.path.exists(out) and not refresh:
            print('exists, skipped:', out); continue
        req = urllib.request.Request(URL.format(name=name), headers={'User-Agent': 'iberian-toponymy-research/1.0'})
        with urllib.request.urlopen(req, timeout=600) as r:
            open(out, 'wb').write(r.read())
        print('saved', out)


if __name__ == '__main__':
    main('--refresh' in sys.argv)
