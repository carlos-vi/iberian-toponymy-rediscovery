# Supplementary material — Rediscovery as validation: unsupervised analysis of Iberian place names recovers known dialect areas, isoglosses, and historical frontiers

This document gives the details of the data and methods and the full results behind every number in the article.

The code and the snapshot of the raw data are in the public repository https://github.com/carlos-vi/iberian-toponymy-rediscovery, and the file and folder names used below refer to it. The repository stores `src/`, `data/raw/` and `supplement/`. The folders `data/derived/`, `results/`, `figures/` and `logs/` are not stored: they are written when the code is run (§S9). The tables of this document are generated automatically from the files in `results/` by `src/70_supplement.py`; the script that produces each file is given in §S9.

---

## S1. Data

**Sources.** Each corpus is downloaded by its own script (see `data/raw/SOURCES.md` for URLs, dates and licences):

- OpenStreetMap: features of Spain, Portugal and Andorra with a `place` tag, through the Overpass API;
- GeoNames: dumps for Spain, Portugal, Andorra and France;
- Instituto Geográfico Nacional: collections *nuc* and *opob* of its OGC API Features, province by province;
- FrequencyWords: lists of the 50,000 most frequent words of six languages;
- Natural Earth: vectors used for the maps.

**Corpora.**

- **S:** OSM `place` ∈ {city, town, village, hamlet, isolated_dwelling, farm}, plus GeoNames class P in the window 42.3–45.3° N, 2.0° W – 3.3° E of France.
- **M:** OSM `place = locality`.
- **G:** GeoNames classes P and A of Spain, Portugal and Andorra.
- **I:** IGN *nuc* and *opob*.
- **S + I:** S and I merged.

All corpora are clipped to 35.5–44.2° N and 10° W – 4.6° E.

**Normalisation** (`src/lib/normalize.py`). The steps are:

- names are lower-cased and Unicode-decomposed, and diacritics are removed except ñ, ç, ã, õ and the grave accents (à, è, ò);
- alternative forms separated by "/" (bilingual names) are split, and each form is tokenised into words of letters;
- articles, prepositions and conjunctions of the six languages are removed;
- the word following a saint's title (*san*, *santa*, *santo*, *sant*, *são*, *saint*…) is removed with the title;
- a closed list of generic words for settlements, buildings and landforms is removed (*casa*, *casas*, *cortijo*, *caserío*, *mas*, *masía*, *quinta*, *monte*, *cerro*, *fuente*, *río*, *puente*, *barrio*, *urbanización*, *polígono*, *lugar*…; the full list is in the code);
- words of fewer than four letters are dropped.

A record keeps the set of its remaining words (`tok`). Records without words are dropped, and records with the same word set in the same 2 km cell are collapsed.

**Regions.** Each record of S, M and G takes the province (Spain), district (Portugal) or department (France) of the nearest GeoNames record of class P or A. Records of I take the province of their INE code.

**Opacity** (`src/lib/opacity.py`). A word is transparent if any of the following holds:

- it is among the 50,000 most frequent words of Spanish, Portuguese, Galician, Catalan, Basque or French;
- removing a plural, diminutive, augmentative or other frequent ending (*-s*, *-es*, *-illo*, *-ito*, *-ico*, *-ejo*, *-uelo*, *-iño*, *-inho*, *-ete*, *-et*, *-ell*, *-on*, *-al*, *-ar*, *-era*, *-ero*, *-eira*, *-edo*…; the full list is in the code) gives such a word, possibly after restoring a final vowel;
- it splits into two words of at least four letters from the 15,000 most frequent of each list, possibly joined by a connector (*de*, *del*, *la*, *dos*…);
- it is a frequent toponymic head (*valde-*, *villa-*, *fuente-*, *castro-*, *torre-*, *quintana-*…) followed by a transparent word.

All other words are opaque. The opaque version of a corpus keeps, in each record, only its opaque words, and drops records left without words.

**Table S1.** Size of the corpora.

{{table corpus_sizes}}

**Table S2.** Records and density (records per 1,000 km², in brackets) by region. Region areas are official land areas.

{{table density}}

## S2. Reference areas

The 21 reference areas follow the dialect classifications of the handbooks cited in the article (Zamora Vicente 1967; Alvar 1996; Veny 1982; Cintra 1971; Bec 1963) and Michelena (1973) for the Basque toponymic area. Each is a union of provinces, Portuguese districts and French departments.

- **Transition units.** Provinces where two areas meet, and French departments at the edge of the window, are not scored. Their cells are excluded from both the areas and the background of the AUC.
- **Cell labels.** Each 20 km cell takes the region of the majority of its records.
- **Size limit.** Areas with fewer than 10 cells in a corpus are not scored in that corpus. This is why I and M score fewer areas than S, which also covers France.

**Table S3.** Reference areas.

{{table areas}}

## S3. Layers

**Element statistics** (`src/lib/layers.py`). For an element carried by *n* records, each record's 25 nearest neighbours are examined. A neighbour is a **hit** if it carries the element in a word with a different lexical key: the first four letters of the word for endings, and for beginnings an approximate lemma (the word without plural, final vowel and diminutive ending). The expected number of hits is the sum, over the element's records, of 25 × (*n* − *n*~key~)/(*N* − 1). Here *N* is the number of records and *n*~key~ the number of the element's records that share the record's key. The statistics are:

- *E* = hits / expected;
- *z* = (hits − expected) / √(binomial variance);
- LCB = ½ χ²~0.025~(2 · hits) / expected, the lower limit of the exact 95 % Poisson interval for *E*;
- *n*~eff~ = exp(entropy of the distribution of the element's records over its lexical keys).

Nested elements on the same side of the word (e.g. *-eta* and *-ieta*) are pruned: the longer one is kept only if its score is at least 1.5 times that of the shorter one.

**Matrix.** For each selected element, counts per 20 km cell and the counts of all records are smoothed with a Gaussian kernel (σ = 1.25 cells = 25 km). The entry for each cell is log[(element + ½*p*) / (all · *p* + ½*p*)], where *p* is the element's overall rate; negative values are set to zero and each row is normalised to sum 1.

**Factorisation.** scikit-learn NMF (initialisation `nndsvda`, up to 1,500 iterations), refitted from four random initialisations. Layers are matched to the main fit by the Hungarian algorithm on map correlations.

**Acceptance.** A layer is accepted if two conditions hold:

- at least three elements have ≥ 35 % of their (row-normalised) weight on it;
- its match correlation is ≥ 0.8 in at least three of the four refits.

**Two-level analysis.**

- **Territories.** Each cell is assigned to the accepted first-level layer with the highest normalised loading, provided it is ≥ 0.25 of that layer's maximum; other cells take the layer of the nearest assigned cell.
- **Second level.** In each territory with ≥ 4,000 records, the element statistics are recomputed on the territory's records. The 150 elements with the highest LCB that pass LCB ≥ 2 and *n*~eff~ ≥ 15 are factorised for *k* = 2–10. The chosen *k* is the largest whose mean match correlation reaches 90 % of the best.

**Evaluation.**

- **Recovery rule.** An area is recovered if an accepted layer has at least half of its core (cells with loading ≥ ½ of its maximum) inside the area and an AUC ≥ 0.90. The AUC is computed with the layer's loading as score, separating the area's cells from all other scored cells.
- **Null models.** Each null replaces every area by a random region of the same number of scored cells, with 200 replicates; the null count is the number of areas "recovered" in a replicate. The three nulls are:
  - **grown** (main): 4-neighbour region growing from a random seed cell;
  - **compact**: the *m* scored cells nearest to a random seed;
  - **provinces**: random contiguous unions of provinces, grown until they reach the size of the area, excluding unions that reproduce the real area.
- **Main metric.** The mean over *k* = 6, 8, …, 24 of the number of areas recovered, compared with the replicate-wise mean of the null counts over the same *k*.

**Table S4.** Elements selected in S by each rule, by region of maximum enrichment.

{{table selection_regions}}

**Table S5.** Main metric for S and S + I.

{{table main_metric}}

**Table S6.** Areas recovered in S for each *k*: recovered (chance; accepted layers).

{{table summary_S}}

**Table S7.** Areas recovered in S + I for each *k*: recovered (chance; accepted layers).

{{table summary_SI}}

**Table S8.** M, I and G at *k* = 18.

{{table other_corpora}}

**Table S9.** Robustness of the main metric: areas recovered / chance (p), mean over *k*.

{{table robustness}}

**Table S10.** Replication at *k* = 18 across S, M, I and G. Expected values assume independent corpora, with per-area chance rates from the grown null. Because the corpora overlap, the true chance level is higher; the ceiling holds under any dependence between corpora, since an area recovered by chance in two corpora needs two chance recoveries (for each area, half the sum of its chance rates).

{{table replication}}

**Table S11.** Territories of the two-level analysis in S + I (*k* = 18).

{{table territories}}

**Table S12.** Accepted layers of the two-level analysis in S + I (*k* = 18).

{{table layers_SI}}

![Figure S1](../figures/figS1_level2.png)

**Figure S1.** All accepted second-level layers of S + I (two-level analysis, *k* = 18).

![Figure S2](../figures/figS2_selection.png)

**Figure S2.** Accepted first-level layers of S at *k* = 18 under selection by significance (top) and by effect size (bottom).

## S4. Spatial scan

**Units and elements** (`src/lib/scan.py`, `src/30_scan.py`).

- **Units.** A unit is a distinct pair (word, 20 km cell).
- **Elements tested.** Endings of 2–5 letters (the word must be at least two letters longer than the ending) and beginnings of 3–5 letters (the word must be longer than the beginning), carried by ≥ 10 distinct words and ≥ 10 units.
- **Windows.** Circles of radius 30, 60, 100 and 150 km around each cell that holds units. Windows holding more than 25 % of all units are dropped.

**Statistic.** Let *B* be the total number of units, *C* the element's units, and *b*, *c* their numbers inside a window. The statistic is the Bernoulli log-likelihood ratio for an excess inside the window:

*c* log(*c*/*b*) + (*b* − *c*) log((*b* − *c*)/*b*) + (*C* − *c*) log((*C* − *c*)/(*B* − *b*)) + (*B* − *b* − *C* + *c*) log((*B* − *b* − *C* + *c*)/(*B* − *b*)) − *C* log(*C*/*B*) − (*B* − *C*) log((*B* − *C*)/*B*)

It is set to 0 when the rate inside the window is not higher than outside, or when *c* < 5. The element's cluster is the window with the highest ratio that contains ≥ 5 distinct words carrying the element.

**Null and false discovery rate.**

- **Null distribution.** For each element, *C* units are placed at random among cells with probability proportional to their units (multinomial), 999 times. The maximum ratio over windows is recorded, with windows of fewer than five units set to 0.
- **Binning.** Null tables are computed for each exact *C* up to 60 and for geometric bins of ±7 % above.
- **p and q.** p = (1 + number of null maxima ≥ observed) / 1,000. Benjamini–Hochberg q-values are computed across all elements of a corpus.

**Second clusters and exact null.**

- **Second clusters.** The second cluster is the best window that shares no cell with the first and holds ≥ 5 distinct words. Its p-value uses the null of the maximum, which is conservative.
- **Exact null for *-eque*.** For the main *-eque* cluster of S + I (opaque names), the *C* units are drawn without replacement from the observed units, 4,999 times.

**Table S13.** Elements tested and significant clusters.

{{table scan_summary}}

**Table S14.** Elements of known strata used to validate the scan. Clusters must be significant and centred in one of the expected provinces; the provinces are listed in `src/lib/reference.py`.

{{table targets}}

**Table S15.** Validation (S + I, all words), by stratum.

{{table validation}}

**Table S16.** Validation, element by element.

{{table validation_elements}}

**Table S17.** Clusters of selected elements in S + I and in G, M and I (opaque names). The expected place is Andalusia for *-ena*, Asturias–León or Aragon for *-iello*, Valencia or the Balearic Islands for *beni-* and the south for *guad-*; for *-eque*, *-aque* and *-que* it is a centre within 120 km of Toledo.

{{table cluster_locations}}

**Table S18.** The 50 strongest clusters among opaque names (S + I), after merging nested elements, with their hand classification. Classes:

- K: stratum or naming system no longer productive;
- D: morphology or vocabulary of a living language;
- U: undocumented.

{{table top50}}

**Table S19.** Second clusters.

{{table secondary}}

**Table S20.** How rare is the concentration of *-eque*? Endings with 20–60 units, and the exact null.

{{table eque_scan}}

## S5. The Castilian wedge

***-iello*** (`src/40_wedge.py`).

- **Word.** The last word of at least four letters of each name (first form, accents kept).
- **Matching.** Endings *-iello*, *-iella*, *-iellos*, *-iellas*, *-iellu*, *-iecho*, *-ieḷḷu* and *-ieyu* are counted against *-illo*, *-illa*, *-illos* and *-illas*.
- **Band.** 41.5–43.6° N in Spain, in bins of 0.5° of longitude.

**F- and h-.**

- **Unit.** Distinct pairs (word, 20 km cell), counted on all words of all forms of each name, generic words included.
- **Pairs.** The eleven matched pairs are given in Table S21. *fonte*, *fontes*, *font* and *fonts*, the simple words for 'spring', are not counted, so that FONTE is represented only by derivatives.
- **Exclusions.** The Basque provinces are excluded.
- **Segments.** West: 41.5–43.6° N, west of 5.5° W. Centre: 4.5–2.0° W. East: east of 1.0° W (also Spain only). South: Spain south of 40° N outside Castellón, Valencia, Alicante and the Balearic Islands.

**Table S21.** Matched pairs for Latin F-.

{{table fh_pairs}}

**Table S22.** *-iello* / (*-iello* + *-illo*): west of 5° W and in the Castilian strip (4.5–1.5° W).

{{table iello}}

**Table S23.** *-iello* west of 5° W by province.

{{table iello_prov}}

**Table S24.** Share of F- forms (number of types) by segment.

{{table fh}}

**Table S25.** Dropping each pair in turn: most extreme values of the share of F- forms.

{{table fh_loo}}

## S6. *-ena*

**Definitions** (`src/41_ena.py`).

- **Types.** Distinct pairs (opaque word, province), Spain and Portugal.
- **Suffix.** A word carries the suffix if it ends in *-ena* (or *-ana*) and has at least two letters before it.
- **Lift.** Core Baetica: Seville, Córdoba, Huelva, Cádiz, Málaga and Jaén. Bootstrap intervals resample the suffix's types (2,000 replicates).
- **Correlations.** Rank correlations use provinces with at least five *-ena* + *-ana* types. The partial correlation regresses the ranks of both variables on the rank of the mean latitude.

**Table S26.** Lift of *-ena* and *-ana* in core Baetica.

{{table ena_lift}}

**Table S27.** Share of *-ena* among *-ena* + *-ana* types by groups of provinces (Basque provinces omitted), and its association with the years of Islamic rule.

{{table ena_groups}}

**Table S28.** By province (S + I, all provinces).

{{table ena_provinces}}

**Table S29.** The Basque homograph and other quantities.

{{table ena_misc}}

![Figure S4](../figures/figS4_ena_provinces.png)

**Figure S4.** Share of *-ena* by province against the years of Islamic rule (S + I; Basque provinces omitted).

## S7. Arabic toponymy and the frontier

**Markers** (`src/lib/reference.py`). The markers are applied to opaque words only, and the counts are of distinct names (cell × word).

**Surface and levels.**

- **Surface.** 15 km cells; numerator and denominator smoothed with σ = 2.3 cells; cells whose smoothed denominator is below 3 are left empty.
- **Northern level.** The median over the cells of A Coruña, Lugo, Ourense, Pontevedra, Asturias and Cantabria.
- **Southern level.** The median over the provinces conquered in or after 1164 (Table S31).

**Frontier.**

- **Rule.** At each longitude from 6.5° to 2.2° W, in steps of 0.1°, the frontier is the first latitude from 38.6° N northwards where the surface is below the threshold and stays below it in all cells of the next 50 km.
- **Comparison lines.** The Duero is a polyline through riverside towns from Soria to Fermoselle. The strongholds are given in the article.
- **Distances.** Distances are in latitude (km) at equal longitude.

**Table S30.** Arabic markers and exclusions.

{{table markers}}

**Table S31.** Approximate year of the definitive Christian conquest used for the years of Islamic rule (year − 714).

{{table conquest}}

**Table S32.** The frontier at the 50 % threshold in three corpora.

{{table frontier}}

**Table S33.** Mean distance (km) of the frontier to the strongholds / to the Duero, by segment and threshold.

{{table frontier_segments}}

**Table S34.** Arabic names by province: correlation of the log share of distinct names with an Arabic marker with the years of Islamic rule, and partial correlation given the mean latitude. Provinces with ≥ 300 distinct words; Basque provinces omitted.

{{table dose}}

![Figure S3](../figures/figS3_frontiers.png)

**Figure S3.** The frontier at the 50 % threshold in S + I, I and G.

## S8. *-aque/-eque*

**Definitions** (`src/43_aque_eque.py`).

- **Types.** Distinct pairs (opaque word, province) of Spain and Portugal; words of at least six letters.
- **Areas.** Centre: Toledo, Ciudad Real, Cuenca, Albacete, Guadalajara and Madrid. South: Andalusia and Murcia.
- **Tests.** One-sided Fisher exact test for a higher proportion of *a/e* before *-que* in the centre.

**Table S35.** Lift and vowel contrast in five corpora, on types and on distinct words.

{{table eque}}

**Table S36.** Opaque names in *-aque* and *-eque* (S + I).

{{table eque_names}}

## S9. Reproduction

**Running the analysis.**

- **Repository.** https://github.com/carlos-vi/iberian-toponymy-rediscovery. It holds the code (`src/`), the snapshot of the raw data (`data/raw/`) and the template of this document (`supplement/`).
- **Requirements.** Python 3.11 with the packages in `requirements.txt`. The analysis needs about 2 GB of memory and runs in 35–60 minutes on a laptop.
- **Run.** `python src/run_all.py` runs every step from the raw snapshot in `data/raw/`; with `--download` it first downloads the data again. The Natural Earth layers of the base maps are not part of the snapshot and are downloaded on the first run. The run creates `data/derived/`, `results/`, `figures/` and `logs/` (one log per step) and renders this document in `supplement/`.

**Table S37.** Scripts and their outputs.

| Step | Script | Output |
|---|---|---|
| Download | `01–05_download_*.py` | `data/raw/` |
| Corpora | `10_build_corpora.py` | `data/derived/corpus_*.csv`, `results/data/corpus_density.csv` |
| Base map | `11_basemap.py` | `data/derived/basemap.json` |
| Layers | `21_layers.py` | `results/layers/` (Tables 2, S4–S8, S11–S12; Figures 2–4, S1–S2) |
| Robustness | `22_layers_robustness.py` | `results/layers/main_metric.csv`, `robustness_*.csv` (Tables S5, S9–S10) |
| Scan | `30_scan.py` | `results/scan/scan_*.csv` |
| Scan analysis | `31_scan_analysis.py` | `results/scan/` (Tables S13–S20; Figure 8) |
| Wedge | `40_wedge.py` | `results/wedge/` (Tables S21–S25; Figure 5) |
| *-ena* | `41_ena.py` | `results/ena/` (Tables S26–S29; Figures 6, S4) |
| Arabic names | `42_arabic_frontier.py` | `results/arabic/` (Tables S30–S34; Figures 7, S3) |
| *-aque/-eque* | `43_aque_eque.py` | `results/aque_eque/` (Tables S35–S36; Figure 8) |
| Figures | `50_figures.py` | `figures/` |
| Numbers | `60_numbers.py` | `results/numbers.json` |
| Supplement | `70_supplement.py` | `supplement/supplementary.md` |
| Supplement (HTML, PDF) | `supplement/build.py` | `supplement/supplementary.html`, `.pdf` |
