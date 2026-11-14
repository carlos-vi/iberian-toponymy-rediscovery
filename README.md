# Rediscovery as validation

Code and data for the article

> Vidal-Perez, C. *Rediscovery as validation: unsupervised analysis of Iberian place names recovers known dialect areas, isoglosses, and historical frontiers.* Manuscript submitted for publication.

The analysis takes about 840,000 georeferenced place-name records of the Iberian Peninsula and south-western France from three public sources (OpenStreetMap, GeoNames and the official gazetteer of Spain) and asks whether an unsupervised method recovers what philologists and historians already know. It

1. extracts recurrent word beginnings and endings without any list of suffixes;
2. factorises their spatial co-occurrence into layers and scores the layers against 21 reference dialect areas, with explicit chance levels;
3. scans every element for local excesses (Kulldorff's spatial scan statistic);
4. measures four known results: the Castilian wedge, the suffix *-ena*, the northern edge of Arabic toponymy, and the names in *-aque/-eque* around Toledo.

This repository holds everything needed to reproduce the article: the code and the snapshot of the raw data it was run on. Every table, figure and number is generated from them by one command.

## Contents

```
src/                  the analysis; run_all.py runs every step in order
  lib/                shared modules: paths, normalisation, opacity filter, reference areas, layers, scan, maps
  01–05_download_*.py downloads of the raw data
  10_build_corpora.py corpora S, M, G, I and S + I, and their opaque versions
  11_basemap.py       base map for the figures
  21–22_*.py          layers, their evaluation against the reference areas, robustness checks
  30–31_*.py          spatial scan and its analysis
  40–43_*.py          targeted measurements: wedge, -ena, Arabic frontier, -aque/-eque
  50_figures.py       figures of the article and of the supplement
  60_numbers.py       every number quoted in the article, in one file
  70_supplement.py    supplementary tables
data/raw/             snapshot of the raw data used in the article (October 2026) and SOURCES.md,
                      which gives the source, date and licence of each file
supplement/           template of the supplementary material and the script that renders it
requirements.txt      Python packages, with the versions used
```

The Natural Earth layers used only to draw the base maps are not stored here; `run_all.py` downloads them (about 47 MB) on the first run.

## Reproducing the analysis

Requirements: Python 3.11 and the packages in `requirements.txt`. About 2 GB of memory and 35–60 minutes on a laptop.

```bash
pip install -r requirements.txt
python src/run_all.py
```

The run writes:

| Folder | Contents |
|---|---|
| `data/derived/` | the corpora built from the raw data, and the base map |
| `results/` | every table behind the text and the figures (CSV), and `numbers.json` with the numbers quoted in the article |
| `figures/` | the eight figures of the article and the four of the supplement |
| `supplement/` | the supplementary material (`supplementary.md`, `.html` and, if Chrome, Chromium or Edge is installed, `.pdf`) |
| `logs/` | one log per step |

All random procedures use fixed seeds, so a run from the snapshot reproduces every number of the article exactly. Other options:

```bash
python src/run_all.py --from 30             # resume from a given step
python src/run_all.py --download --refresh  # download current data first; results will differ slightly
```

## Where each result comes from

| Article | Script | Output |
|---|---|---|
| Table 1, Figure 1 | `10_build_corpora.py` | `results/data/` |
| Table 2, Figures 2–4 | `21_layers.py`, `22_layers_robustness.py` | `results/layers/` |
| Section 5.2, Figure 5 | `40_wedge.py` | `results/wedge/` |
| Section 5.3, Figure 6 | `41_ena.py`, `31_scan_analysis.py` | `results/ena/`, `results/scan/` |
| Section 5.4, Figure 7 | `42_arabic_frontier.py` | `results/arabic/` |
| Section 5.5, Figure 8 | `30_scan.py`, `31_scan_analysis.py`, `43_aque_eque.py` | `results/scan/`, `results/aque_eque/` |
| Supplementary tables | `70_supplement.py` | `supplement/supplementary.md` |

## Licences

The code is released under the MIT licence (see `LICENSE`). The data keep the licences of their sources:

- OpenStreetMap: © OpenStreetMap contributors, ODbL 1.0;
- GeoNames: CC BY 4.0;
- Instituto Geográfico Nacional (Spain): CC BY 4.0;
- FrequencyWords word lists (H. Dave): CC BY-SA 4.0;
- Natural Earth: public domain.

Corpora derived from OpenStreetMap data (S, M and S + I) are subject to the ODbL. See `data/raw/SOURCES.md`.

## Contact

Carlos Vidal-Perez · [ORCID 0009-0006-2500-3672](https://orcid.org/0009-0006-2500-3672) · carlosvidalpe97@gmail.com
