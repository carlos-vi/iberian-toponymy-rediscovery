"""Normalisation of place names into 'words' (articles, saints' names and generic words removed)."""
import re
import unicodedata
import numpy as np

STOP = set("""
de del da do das dos la las el els les lo los o a os as al ao aos y e i en em na no nas nos d l s
sa ses es sant santa san santo sao são sta sto santiago santi sainte saint st ste the of du des et le aux en sur sous
nuestra nossa senhora señora virgen mare deu nostra
""".split())
GENERIC = set("""
casa casas caseta caserio caserío casal casais casar cortijo cortijos cortijada mas masia masía can cal cases
monte montes cerro cerros loma lomas llano llanos cuesta cuestas fuente fuentes fonte fontes font fonts pozo pozos poço
arroyo arroyos barranco barrancos riu rio río rios ribeira ribera rivera valle valles vale vales vall valls val
sierra serra serras alto alta altos altas bajo baja baixo baixa baixos puente ponte pont molino molinos moinho moinhos
molí moli camino caminho cami camí peña peñas pena penas penya penyes collado coll puerto porto port prado prados prat prats
campo campos camp camps huerta huertas horta hortas hort quinta quintas granja granjas venta ventas lugar barrio bairro
barri aldea aldeia aldeas aldeias torre torres castillo castelo castell iglesia igreja ermita cruz creu pico picos puig
tossal cabezo cabeço umbria umbría solana dehesa hoya hoyo vega vegas nava plana pla partida paraje finca cova cueva cuevas
raso lomo ladera majada corral corrales barco playa praia platja cala cabo punta isla ilha illa estacion estación
urbanizacion urbanización poligono polígono colonia poblado nuevo nueva nuevos novo nova novas viejo vieja velho velha
grande grandes chico chica pequeño pequeña arriba abajo cima fondo fundo mayor menor medio meio baix dalt
cabana cabanas cabaña cabañas borda bordas pajar pajares era eras mina minas tejar tejera alberca balsa casilla casillas
granja quinta herdade monte barranc font hostal hotel restaurante bar camping estacio carretera ctra cooperativa
lieu dit hameau ferme maison bois moulin chateau château pont mont croix grange bastide borde
""".split())
SAINT = {'san', 'santa', 'santo', 'sant', 'sao', 'são', 'sta', 'sto', 'santi', 'saint', 'sainte', 'st', 'ste', 'sancti'}
KEEP = {'ñ', 'ç', 'ã', 'õ', 'à', 'è', 'ò'}
TOK = re.compile(r"[a-zñçãõàèò]+")


def norm(s):
    """Lower case; strip diacritics except ñ, ç, ã, õ and the grave accents (regional spelling features)."""
    s = unicodedata.normalize('NFC', str(s).lower())
    out = []
    for ch in s:
        if ch in KEEP:
            out.append(ch)
        else:
            d = unicodedata.normalize('NFD', ch)
            out.append(''.join(c for c in d if unicodedata.category(c) != 'Mn'))
    return ''.join(out)


def tokens(name):
    """Distinct words (≥ 4 letters) of a name; bilingual forms 'A / B' contribute both; saints' names are dropped."""
    toks = []
    for part in re.split(r"\s+/\s+|/", name):
        prev = ''
        for t in TOK.findall(norm(part)):
            if prev in SAINT:
                prev = t
                continue
            prev = t
            if len(t) >= 4 and t not in STOP and t not in GENERIC:
                toks.append(t)
    return sorted(set(toks))


def raw_last(name):
    """Last word (≥ 4 letters, accents kept) of the first form of a name."""
    s = unicodedata.normalize('NFC', str(name).lower())
    s = re.split(r"\s+/\s+|/", s)[0]
    ws = [w for w in re.findall(r"[^\W\d_]+", s) if len(w) >= 4]
    return ws[-1] if ws else ''


def xy(lat, lon):
    """Equirectangular projection centred on 40° N, in km (adequate at peninsular scale)."""
    return np.asarray(lon) * 111.32 * np.cos(np.radians(40.0)), np.asarray(lat) * 110.57
