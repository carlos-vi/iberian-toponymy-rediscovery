"""Download the GeoNames country extracts for Spain, Portugal, Andorra and France (CC BY 4.0).

Output: data/raw/geonames/{ES,PT,AD,FR}.txt (tab-separated, the 'geoname' table of GeoNames). The snapshot used in the
paper was downloaded on 19 October 2026.
Usage: python 02_download_geonames.py [--refresh]
"""
import io, os, sys, urllib.request, zipfile
from lib import paths

URL = 'https://download.geonames.org/export/dump/{cc}.zip'


def main(refresh=False):
    os.makedirs(paths.raw('geonames'), exist_ok=True)
    for cc in ['ES', 'PT', 'AD', 'FR']:
        out = paths.raw('geonames', f'{cc}.txt')
        if os.path.exists(out) and not refresh:
            print('exists, skipped:', out); continue
        req = urllib.request.Request(URL.format(cc=cc), headers={'User-Agent': 'iberian-toponymy-research/1.0'})
        with urllib.request.urlopen(req, timeout=600) as r:
            z = zipfile.ZipFile(io.BytesIO(r.read()))
        open(out, 'wb').write(z.read(f'{cc}.txt'))
        print('saved', out)


if __name__ == '__main__':
    main('--refresh' in sys.argv)
