"""Download the 50,000-word frequency lists (OpenSubtitles 2018, hermitdave/FrequencyWords, CC BY-SA 4.0) used to
decide whether a word of a place name is a word of a living language: Spanish, Portuguese, Galician, Catalan, Basque
and French.
Output: data/raw/lexicons/{lang}_50k.txt
Usage: python 04_download_lexicons.py [--refresh]
"""
import os, sys, urllib.request
from lib import paths

URL = 'https://raw.githubusercontent.com/hermitdave/FrequencyWords/master/content/2018/{lang}/{lang}_50k.txt'


def main(refresh=False):
    os.makedirs(paths.raw('lexicons'), exist_ok=True)
    for lang in ['es', 'pt', 'gl', 'ca', 'eu', 'fr']:
        out = paths.raw('lexicons', f'{lang}_50k.txt')
        if os.path.exists(out) and not refresh:
            print('exists, skipped:', out); continue
        req = urllib.request.Request(URL.format(lang=lang), headers={'User-Agent': 'iberian-toponymy-research/1.0'})
        with urllib.request.urlopen(req, timeout=300) as r:
            open(out, 'wb').write(r.read())
        print('saved', out)


if __name__ == '__main__':
    main('--refresh' in sys.argv)
