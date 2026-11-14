# Raw data: sources, dates and licences

All raw files in this folder were downloaded by the scripts in `src/01–05_download_*.py`. The files used in the article are kept as a snapshot because three of the sources (OpenStreetMap, GeoNames, IGN) change over time. Running `python src/run_all.py --download --refresh` replaces the snapshot with current data. The results will then differ slightly from the published ones.

| Folder | Source | What it contains | Downloaded | Licence |
|---|---|---|---|---|
| `osm/` | OpenStreetMap, through the public Overpass API (`https://overpass-api.de/api/interpreter`). Query in `osm/places_iberia.overpassql` | Named features of Spain, Portugal and Andorra tagged `place` = city, town, village, hamlet, isolated_dwelling, locality or farm. One row per feature with type, id, latitude, longitude, `name`, `place`, `name:es` and `name:pt` (`osm_places_iberia.tsv`, 43 MB) | 20 Oct 2026 | ODbL 1.0, © OpenStreetMap contributors |
| `geonames/` | GeoNames country dumps (`https://download.geonames.org/export/dump/{ES,PT,AD,FR}.zip`) | Full gazetteers of Spain, Portugal, Andorra and France (`*.txt`, tab-separated, standard GeoNames columns) | 19 Oct 2026 | CC BY 4.0, GeoNames |
| `ign/` | Instituto Geográfico Nacional (Spain), OGC API – Features (`https://api-features.ign.es/collections/{nuc,opob}/items`) | Population nuclei of the INE (`ign_nuc.csv`, 35,811 records) and other populated places (`ign_opob.csv`, 143,467 records) of 48 provinces. Mainland Spain and the Balearic Islands only: the Canary Islands, Ceuta and Melilla are not requested. `download.log` compares the records received with the number announced by the service for each province | 24 Oct 2026 | CC BY 4.0, © Instituto Geográfico Nacional |
| `lexicons/` | FrequencyWords, 2018 edition (H. Dave; word lists derived from the OpenSubtitles corpus), `https://github.com/hermitdave/FrequencyWords` | The 50,000 most frequent words of Spanish, Portuguese, Galician, Catalan, Basque and French (`{es,pt,gl,ca,eu,fr}_50k.txt`) | 20 Oct 2026 | Lists: CC BY-SA 4.0; code: MIT (see the repository) |
| `naturalearth/` | Natural Earth 1:10m vectors, GeoJSON edition (`https://github.com/nvkelso/natural-earth-vector`) | Land, rivers (global and European supplement), lakes, first-order administrative lines and land boundaries. Used only to draw the maps | 20 Oct 2026 | Public domain |

Further notes:

- **French data.** France enters the study only through GeoNames populated places (class P) of the south-western box 42.3–45.3° N, 2.0° W – 3.3° E. These records are added to corpus S so that the Gascon and Occitan reference areas can be scored.
- **Processing.** Normalisation, deduplication and the opacity filter are applied by `src/10_build_corpora.py`. It writes its output to `data/derived/`. Nothing in `data/raw/` is modified by the analysis.
- **Size.** About 150 MB uncompressed.
