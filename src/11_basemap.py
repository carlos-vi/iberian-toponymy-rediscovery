"""Base-map layers for the figures from Natural Earth 1:10m (land, rivers, lakes, provincial and national boundaries),
clipped to the study area and simplified (Douglas–Peucker).
Output: data/derived/basemap.json
"""
import json
import numpy as np
from lib import paths

BBOX = (-10.0, 35.6, 4.6, 44.2)   # lon0, lat0, lon1, lat1 (with a margin of 0.5° when clipping)


def simplify(pts, tol):
    if len(pts) < 3:
        return pts
    a, b = np.array(pts[0]), np.array(pts[-1]); P = np.array(pts)
    ab = b - a; L = np.hypot(*ab)
    d = np.abs(ab[0] * (P[:, 1] - a[1]) - ab[1] * (P[:, 0] - a[0])) / L if L > 0 else np.hypot(*(P - a).T)
    i = int(np.argmax(d))
    if d[i] > tol:
        return simplify(pts[:i + 1], tol)[:-1] + simplify(pts[i:], tol)
    return [pts[0], pts[-1]]


def inside(c):
    lon0, lat0, lon1, lat1 = BBOX
    return lon0 - .5 <= c[0] <= lon1 + .5 and lat0 - .5 <= c[1] <= lat1 + .5


def lines(name, tol=.01, keep=lambda p: True):
    gj = json.load(open(paths.raw('naturalearth', f'{name}.geojson'), encoding='utf-8'))
    out = []
    for f in gj['features']:
        p = f['properties']
        if not keep(p):
            continue
        g = f['geometry']
        if g['type'] == 'LineString': parts = [g['coordinates']]
        elif g['type'] in ('MultiLineString', 'Polygon'): parts = g['coordinates']
        elif g['type'] == 'MultiPolygon': parts = [r for poly in g['coordinates'] for r in poly]
        else: parts = []
        for part in parts:
            seg = []
            for c in part:
                if inside(c):
                    seg.append(c)
                elif seg:
                    if len(seg) > 1: out.append((p, seg))
                    seg = []
            if len(seg) > 1:
                out.append((p, seg))
    return [(p, [[round(x, 3), round(y, 3)] for x, y in simplify([tuple(c[:2]) for c in seg], tol)]) for p, seg in out]


def main():
    land = [s for _, s in lines('ne_10m_land', .01)]
    rivers = []
    for fn in ('ne_10m_rivers_lake_centerlines', 'ne_10m_rivers_europe'):
        for p, s in lines(fn, .008, lambda p: p.get('featurecla', '') != 'Lake Centerline'):
            rivers.append(dict(n=p.get('name_es') or p.get('name') or '', r=int(p.get('scalerank') or 9), c=s))
    lakes = [s for _, s in lines('ne_10m_lakes', .005)]
    prov = [s for _, s in lines('ne_10m_admin_1_states_provinces_lines', .01,
                                lambda p: (p.get('ADM0_NAME') or p.get('adm0_name')) in ('Spain', 'Portugal'))]
    border = [s for _, s in lines('ne_10m_admin_0_boundary_lines_land', .01)]
    json.dump(dict(bbox=BBOX, land=land, rivers=rivers, lakes=lakes, prov=prov, border=border),
              open(paths.derived('basemap.json'), 'w', encoding='utf-8'), ensure_ascii=False, separators=(',', ':'))
    print('land', len(land), 'rivers', len(rivers), 'provincial lines', len(prov), 'borders', len(border))


if __name__ == '__main__':
    main()
