"""Lexical opacity: a word is transparent if it is (or decomposes into) a word of a living language of the area.

Lexicon: 50,000-word frequency lists (OpenSubtitles, hermitdave/FrequencyWords) for Spanish, Portuguese, Galician,
Catalan, Basque and French, normalised like the place names. A word is transparent if
  - it is in the lexicon, or is after removing a plural, diminutive or augmentative ending, or
  - it splits into two frequent words (≥ 4 letters each), optionally joined by a connector (de, del, la…), or
  - it is a frequent toponymic compound head (valde-, villa-, fuente-…) followed by a lexical word.
"""
from functools import lru_cache
from . import paths
from .normalize import norm

LANGS = ['es', 'pt', 'gl', 'ca', 'eu', 'fr']
CONNECT = ['de', 'del', 'la', 'el', 'los', 'las', 'les', 'els', 'da', 'do', 'das', 'dos', 'y', 'i', 'e', 'a', 'o', 'al', 'ao']
SUFFIXES = ['s', 'es', 'illo', 'illa', 'illos', 'illas', 'ito', 'ita', 'ico', 'ica', 'ejo', 'eja', 'uelo', 'uela', 'ino', 'ina',
            'iño', 'iña', 'inho', 'inha', 'ete', 'eta', 'et', 'ell', 'ella', 'on', 'ona', 'ones', 'azo', 'al', 'ar', 'era', 'ero',
            'eira', 'eiro', 'edo', 'eda', 'ela', 'elo', 'ais', 'iello', 'iella', 'ín', 'in', 'ot', 'ota', 'ote', 'uca', 'uco']
HEADS = ('valde', 'val', 'villa', 'vila', 'fuente', 'fuen', 'fonte', 'font', 'castro', 'castell', 'torre', 'quintana',
         'quintanilla', 'casa', 'mas', 'puig', 'pena', 'peña', 'monte', 'mont', 'santa', 'san', 'sant', 'cabeza', 'cabezo', 'ribera')
LEX = None
SEG = None


def load_lex(minlen=3, top=None):
    lex = set()
    for lang in LANGS:
        with open(paths.raw('lexicons', f'{lang}_50k.txt'), encoding='utf-8') as f:
            for r, line in enumerate(f):
                if top and r >= top:
                    break
                w = norm(line.split(' ')[0])
                if len(w) >= minlen and w.isalpha():
                    lex.add(w)
    return lex


def init():
    global LEX, SEG
    if LEX is None:
        LEX = load_lex()
        SEG = load_lex(minlen=4, top=15000)


def in_lex(w):
    if w in LEX:
        return True
    for s in SUFFIXES:
        if w.endswith(s) and len(w) - len(s) >= 3:
            b = w[:-len(s)]
            if b in LEX or b + 'a' in LEX or b + 'o' in LEX or b + 'e' in LEX:
                return True
    return False


def seg_word(w):
    if w in SEG:
        return True
    return any(w.endswith(s) and len(w) - len(s) >= 4 and w[:-len(s)] in SEG for s in ('s', 'es'))


@lru_cache(maxsize=None)
def segmentable(w):
    for i in range(4, len(w) - 3):
        head, rest = w[:i], w[i:]
        if not seg_word(head):
            continue
        if seg_word(rest):
            return True
        for c in CONNECT:
            if rest.startswith(c) and len(rest) - len(c) >= 4 and seg_word(rest[len(c):]):
                return True
    for p in HEADS:
        if w.startswith(p) and len(w) - len(p) >= 3:
            rest = w[len(p):]
            for c in [''] + CONNECT:
                r2 = rest[len(c):] if rest.startswith(c) else None
                if r2 and len(r2) >= 3 and in_lex(r2):
                    return True
    return False


def transparent(w):
    init()
    return in_lex(w) or segmentable(w)
