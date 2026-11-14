"""Download OpenStreetMap populated places and localities of Spain, Portugal and Andorra (Overpass API, ODbL).

Query: data/raw/osm/places_iberia.overpassql (place = city, town, village, hamlet, isolated_dwelling, locality, farm;
named features; output as TSV with type, id, latitude, longitude, name, place, name:es, name:pt).
Output: data/raw/osm/osm_places_iberia.tsv. The public Overpass endpoint is often at capacity: the script retries with
exponential back-off. OpenStreetMap changes continuously; the snapshot used in the paper was downloaded on 20 October 2026.
Usage: python 01_download_osm.py [--refresh]
"""
import os, sys, time, urllib.request, urllib.parse
from lib import paths

ENDPOINT = 'https://overpass-api.de/api/interpreter'
QUERY = """[out:csv(::type,::id,::lat,::lon,name,place,"name:es","name:pt";true;"\\t")][timeout:900][maxsize:1073741824];
(area["ISO3166-1"="ES"][admin_level=2];area["ISO3166-1"="PT"][admin_level=2];area["ISO3166-1"="AD"][admin_level=2];)->.a;
nwr(area.a)["place"~"^(city|town|village|hamlet|isolated_dwelling|locality|farm)$"]["name"];
out center;
"""


def main(refresh=False):
    os.makedirs(paths.raw('osm'), exist_ok=True)
    qpath = paths.raw('osm', 'places_iberia.overpassql')
    open(qpath, 'w', encoding='utf-8').write(QUERY)
    out = paths.raw('osm', 'osm_places_iberia.tsv')
    if os.path.exists(out) and not refresh:
        print('exists, skipped:', out); return
    data = urllib.parse.urlencode({'data': QUERY}).encode()
    wait = 30
    for attempt in range(8):
        try:
            req = urllib.request.Request(ENDPOINT, data=data, headers={'User-Agent': 'iberian-toponymy-research/1.0'})
            with urllib.request.urlopen(req, timeout=1200) as r:
                body = r.read()
            if body.startswith(b'@type'):
                open(out, 'wb').write(body)
                print('saved', out, len(body) // 1_000_000, 'MB'); return
            print('unexpected answer, retrying:', body[:200])
        except Exception as e:
            print('attempt', attempt, e)
        time.sleep(wait); wait *= 2
    raise SystemExit('Overpass download failed')


if __name__ == '__main__':
    main('--refresh' in sys.argv)
