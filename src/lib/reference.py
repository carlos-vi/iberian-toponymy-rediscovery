"""Reference information taken from the literature: dialect areas, element–region associations of known strata,
dates of the Christian conquest, Arabic markers, the Duero and the strongholds of the Umayyad Middle March."""
import re

# INE province code (cpro) -> province code used by GeoNames (the 'ES-xx' region labels of the corpora)
INE2GN = {1: 'VI', 2: 'AB', 3: 'A', 4: 'AL', 5: 'AV', 6: 'BA', 7: 'PM', 8: 'B', 9: 'BU', 10: 'CC', 11: 'CA', 12: 'CS',
          13: 'CR', 14: 'CO', 15: 'C', 16: 'CU', 17: 'GI', 18: 'GR', 19: 'GU', 20: 'SS', 21: 'H', 22: 'HU', 23: 'J',
          24: 'LE', 25: 'L', 26: 'LO', 27: 'LU', 28: 'M', 29: 'MA', 30: 'MU', 31: 'NA', 32: 'OR', 33: 'O', 34: 'P',
          36: 'PO', 37: 'SA', 39: 'S', 40: 'SG', 41: 'SE', 42: 'SO', 43: 'T', 44: 'TE', 45: 'TO', 46: 'V', 47: 'VA',
          48: 'BI', 49: 'ZA', 50: 'Z'}

_es = lambda s: ['ES-' + x for x in s.split()]
_pt = lambda s: ['PT-' + x for x in s.split()]
_fr = lambda s: ['FR-' + x for x in s.split()]

# ------------------------------------------------------------------ reference dialect areas (Table S1)
# Zamora Vicente (1967) and Alvar (1996) for Spanish; Veny (1982) Catalan; Cintra (1971) Galician-Portuguese;
# Bec (1963) Occitan; Michelena (1973) for the Basque toponymic area. Mixed provinces are excluded.
AREAS = {
    'Galician': _es('C LU OR PO'),
    'Asturian': _es('O'),
    'Leonese': _es('LE ZA'),
    'Cantabrian': _es('S'),
    'Basque': _es('SS BI VI NA'),
    'Upper Aragonese': _es('HU'),
    'Ebro Aragon and Teruel': _es('Z TE'),
    'Northern Castilian': _es('BU P VA SG AV SO LO'),
    'Southern Castilian': _es('TO CR CU AB GU M'),
    'Extremaduran': _es('CC BA'),
    'Western Andalusian': _es('SE H CA CO MA'),
    'Eastern Andalusian': _es('GR AL J'),
    'Murcian': _es('MU'),
    'Catalan': _es('B GI L T') + ['AD', 'FR-66'],
    'Valencian': _es('CS V A'),
    'Balearic': _es('PM'),
    'Northern Portuguese': _pt('20 04 17 21 05 22 11'),
    'Central-Southern Portuguese': _pt('07 06 13 14 18 16 19 08 03 09'),
    'Gascon': _fr('40 32 65'),
    'Languedocien': _fr('81 82 11 34 12 46'),
    'Northern Occitan': _fr('19 15'),
}
EXCLUDED = set(_es('SA') + _pt('02') + _fr('64 33 47 31 09 24 16 17 30 43 48 63'))


def partition(variant='main'):
    """Reference partition and excluded regions; the alternatives are used only for robustness checks."""
    areas = {a: list(r) for a, r in AREAS.items()}
    excl = set(EXCLUDED)
    if variant in ('navarre', 'navarre+salamanca'):
        areas['Basque'] = _es('SS BI VI'); areas['Navarrese'] = _es('NA')
    if variant in ('salamanca', 'navarre+salamanca'):
        areas['Leonese'] = _es('LE ZA SA'); excl.discard('ES-SA')
    if variant == 'one-andalusia':
        areas['Andalusian'] = areas.pop('Western Andalusian') + areas.pop('Eastern Andalusian')
    if variant == 'asturleonese':
        areas['Asturleonese'] = areas.pop('Asturian') + areas.pop('Leonese')
    return areas, excl


def area_lookup(areas, excluded):
    reg2 = {r: a for a, rs in areas.items() for r in rs}
    return lambda region: reg2.get(region, 'EXCL' if region in excluded else None)


area_of = area_lookup(AREAS, EXCLUDED)

# ------------------------------------------------------------------ element–region associations of known strata
_AND = 'AL CA CO GR H J MA SE'.split(); _EXT = 'BA CC'.split(); _CLM = 'AB CR CU GU TO'.split()
_VAL = 'A V CS'.split(); _CAT = 'B GI L T'.split(); _GAL = 'C LU OR PO'.split(); _VASC = 'SS BI VI NA'.split()
_CYL = 'AV BU LE P SA SG SO VA ZA'.split()
_ES = lambda xs: {'ES-' + x for x in xs}
_PT_N = {'PT-' + x for x in '20 04 17 21 05 22 02 11 07'.split()}
_PT_S = {'PT-' + x for x in '09 03 08 19 14 18 16'.split()}
_GASC = {'FR-' + x for x in '64 65 40 32 31 33 47'.split()}
_ALLFR = {f'FR-{i:02d}' for i in range(1, 96)}
TARGETS = [  # (elements, expected regions, stratum)
    (['^guad'], _ES(_AND + _EXT + _CLM + ['MU', 'M']) | _PT_S, 'Arabic'),
    (['^beni', '^bini', '^ben'], _ES(_VAL + ['PM', 'MU', 'AL', 'GR', 'MA']), 'Arabic'),
    (['^alc', '^alm', '^alb'], _ES(_AND + ['MU'] + _VAL + _CLM + ['Z', 'TE']) | _PT_S, 'Arabic'),
    (['eta$', 'aga$', 'egi$', 'egui$', 'ola$', 'tza$', 'za$', 'uri$', 'iri$', 'berri$', '^etxe', '^eche'],
     _ES(_VASC) | {'FR-64'}, 'Basque'),
    (['ues$', 'oz$'], _ES(['NA', 'HU']) | _GASC, 'Aquitanian-Pyrenean'),
    (['os$'], _GASC, 'Aquitanian (Gascony)'),
    (['ac$'], _ALLFR, 'Gallo-Romance -ACUM'),
    (['an$'], _GASC, 'Latin -ANUM (Gascony)'),
    (['ià$', 'à$'], _ES(_CAT + ['PM']), 'Latin fundus (Catalan)'),
    (['ell$', 'ella$'], _ES(_CAT + ['PM'] + _VAL), 'Catalan Romance'),
    (['iz$'], _ES(_GAL) | _PT_N, 'Genitive -ICI'),
    (['ães$', 'ãe$'], _PT_N, 'Germanic/Latin genitive'),
    (['mil$', 'mir$', 'ufe$', 'ofe$', 'ife$', 'ide$'], _ES(_GAL + ['LE', 'ZA', 'O']) | _PT_N, 'Germanic personal names'),
    (['obre$', 'bre$'], _ES(_GAL), 'Celtic -briga'),
    (['uelo$', 'uela$'], _ES(_CYL + _CLM + _AND + _EXT + ['M', 'LO']), 'Castilian Romance'),
    (['eiro$', 'eira$', 'iño$'], _ES(_GAL) | _PT_N | _PT_S, 'Galician-Portuguese'),
]

# ------------------------------------------------------------------ Islamic rule
# year of the definitive Christian conquest (approximate, by province or district); years of rule = year − 714
CONQ = {}
def _setd(codes, year, pre='ES-'):
    for c in codes.split():
        CONQ[pre + c] = year
_setd('C LU OR PO O S BI SS VI', 750); _setd('NA', 800); _setd('GI', 785); _setd('B', 801)
_setd('LE', 850); _setd('BU', 860); _setd('P', 870); _setd('ZA', 893); _setd('VA', 900); _setd('LO', 923)
_setd('20 04 17', 868, 'PT-'); _setd('21 05', 900, 'PT-'); _setd('22', 1058, 'PT-'); _setd('02 07', 1064, 'PT-'); _setd('11', 1199, 'PT-')
_setd('SA', 1100); _setd('AV SG', 1088); _setd('SO', 1100); _setd('M', 1083); _setd('TO GU', 1085)
_setd('HU', 1096); _setd('Z', 1118); _setd('L', 1149); _setd('T', 1148)
_setd('13 18 14', 1147, 'PT-'); _setd('06', 1165, 'PT-'); _setd('16', 1230, 'PT-')
_setd('CU', 1177); _setd('TE', 1171); _setd('CR', 1212); _setd('CC', 1229); _setd('BA', 1230); _setd('AB', 1242)
_setd('19', 1217, 'PT-'); _setd('08', 1165, 'PT-'); _setd('03', 1232, 'PT-'); _setd('PM', 1229)
_setd('CS', 1233); _setd('V', 1238); _setd('MU', 1243); _setd('A', 1266)
_setd('CO', 1236); _setd('J', 1246); _setd('SE', 1248); _setd('H CA', 1262); _setd('09', 1249, 'PT-')
_setd('MA', 1487); _setd('AL', 1489); _setd('GR', 1492)

# Arabic markers (applied to opaque words only) and frequent Arabic loanwords of Castilian that are excluded
STRICT = re.compile(r"^(al[bcdfghjklmnpqrstvxzç]|beni|bini|guad|uad|calat|medin|rafal|rahal|alqu|almun|alcaz|rabit|rapit|aljam)")
LEXAR = re.compile(r"^(alde|alcor|alcub|alberg|almenar|alfoz|alfar|alcant|alber[ct]|alcanta|albarr|alcorn|alt[aeiou]|alto|alb[aeiou]?$)")
VASC = {'ES-BI', 'ES-SS', 'ES-VI', 'ES-NA'}

# the Duero (polyline through riverside towns, east to west) and the strongholds of the Middle March
RIVER_LL = [('Soria', 41.764, -2.465), ('Almazán', 41.486, -2.531), ('San Esteban de Gormaz', 41.574, -3.204),
            ('Aranda de Duero', 41.670, -3.689), ('Roa', 41.697, -3.926), ('Peñafiel', 41.602, -4.114),
            ('Tudela de Duero', 41.585, -4.581), ('Tordesillas', 41.502, -5.001), ('Toro', 41.523, -5.393),
            ('Zamora', 41.503, -5.745), ('Fermoselle', 41.317, -6.396)]
FORTS = [('Coria', 40.00, -6.54), ('Talavera', 39.96, -4.83), ('Madrid', 40.42, -3.70), ('Talamanca', 40.74, -3.51),
         ('Guadalajara', 40.63, -3.17), ('Gormaz', 41.52, -3.01), ('Atienza', 41.20, -2.87), ('Medinaceli', 41.17, -2.43)]

# ------------------------------------------------------------------ regions for the density table
MACRO = {'Galicia': 'C LU OR PO', 'Asturias': 'O', 'Cantabria': 'S', 'Basque Country and Navarre': 'BI SS VI NA',
         'Castile and León': 'AV BU LE P SA SG SO VA ZA', 'La Rioja': 'LO', 'Aragon': 'HU Z TE', 'Catalonia': 'B GI L T',
         'Valencia': 'A V CS', 'Balearic Islands': 'PM', 'Murcia': 'MU', 'Castile-La Mancha': 'AB CR CU GU TO',
         'Madrid': 'M', 'Extremadura': 'BA CC', 'Andalusia': 'AL CA CO GR H J MA SE'}
AREA_KM2 = {'Galicia': 29575, 'Asturias': 10604, 'Cantabria': 5321, 'Basque Country and Navarre': 17634,
            'Castile and León': 94226, 'La Rioja': 5045, 'Aragon': 47720, 'Catalonia': 32108, 'Valencia': 23255,
            'Balearic Islands': 4992, 'Murcia': 11313, 'Castile-La Mancha': 79461, 'Madrid': 8028, 'Extremadura': 41634,
            'Andalusia': 87599, 'Portugal': 89015}


def macro(region):
    if region.startswith('PT'):
        return 'Portugal'
    if region.startswith('FR') or region == 'AD':
        return 'France and Andorra'
    c = region[3:]
    for g, cs in MACRO.items():
        if c in cs.split():
            return g
    return 'other'
