"""Folder layout of the project and loaders for the derived corpora."""
import os
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
RAW = os.path.join(ROOT, 'data', 'raw')
DERIVED = os.path.join(ROOT, 'data', 'derived')
RESULTS = os.path.join(ROOT, 'results')
FIGURES = os.path.join(ROOT, 'figures')
for _p in (RAW, DERIVED, RESULTS, FIGURES):
    os.makedirs(_p, exist_ok=True)


def raw(*parts):
    return os.path.join(RAW, *parts)


def derived(*parts):
    return os.path.join(DERIVED, *parts)


def results(*parts):
    path = os.path.join(RESULTS, *parts)
    os.makedirs(os.path.dirname(path), exist_ok=True)
    return path


def figure(name):
    return os.path.join(FIGURES, name)


def load_corpus(key, cols=None):
    """Corpus `key` (S, M, G, I, SI, or the opaque versions So, Mo, Go, Io, SIo)."""
    return pd.read_csv(derived(f'corpus_{key}.csv'), keep_default_na=False, encoding='utf-8', usecols=cols)
