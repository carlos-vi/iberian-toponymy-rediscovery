"""Collects every number quoted in the manuscript into one file, results/numbers.json, so that each figure in the
text can be traced to the table that produced it (key = path of the table inside results/).
"""
import os, json
import numpy as np, pandas as pd
from lib import paths

R = paths.RESULTS
TABLES = ['data/corpus_density.csv', 'data/corpus_sizes.csv',
          'layers/main_metric.csv', 'layers/robustness_layers.csv', 'layers/robustness_replication.csv',
          'layers/selection_by_region_S.csv', 'layers/selection_by_region_SI.csv',
          'scan/scan_summary.csv', 'scan/validation_summary.csv', 'scan/validation_by_stratum.csv',
          'scan/cluster_locations.csv', 'scan/top50_summary.csv', 'scan/eque_rarity.csv', 'scan/eque_exact_null.csv',
          'scan/secondary_clusters.csv',
          'wedge/iello_segments.csv', 'wedge/iello_west_by_province.csv', 'wedge/fh_segments.csv',
          'wedge/fh_leave_one_pair_out_extremes.csv',
          'ena/lift_all_provinces.csv', 'ena/without_basque.csv', 'ena/basque_homograph.csv', 'ena/core_vs_nasrid.csv',
          'ena/repeated_summary_SI.csv',
          'arabic/frontier_summary.csv', 'arabic/frontier_segments.csv', 'arabic/dose_response.csv',
          'aque_eque/tests.csv']


def corpus_sizes():
    rows = []
    for key in ['S', 'M', 'G', 'I', 'SI']:
        d = paths.load_corpus(key, cols=['region', 'tok'])
        o = paths.load_corpus(key + 'o', cols=['tok'])
        words = {t for tk in d.tok for t in tk.split()}
        rows.append(dict(corpus=key, records=len(d), distinct_words=len(words), opaque_records=len(o),
                         spain=int(d.region.str.startswith('ES').sum()), portugal=int(d.region.str.startswith('PT').sum()),
                         france=int(d.region.str.startswith('FR').sum()), andorra=int((d.region == 'AD').sum())))
    t = pd.DataFrame(rows)
    t.to_csv(paths.results('data', 'corpus_sizes.csv'), index=False)
    print(t.to_string(index=False))


def main():
    corpus_sizes()
    out = {}
    for rel in TABLES:
        p = os.path.join(R, rel)
        if not os.path.exists(p):
            print('missing:', rel); continue
        t = pd.read_csv(p, keep_default_na=False)
        out[rel] = json.loads(t.to_json(orient='records', force_ascii=False))
    for key in ['S', 'SI', 'M', 'I', 'G']:
        p = os.path.join(R, 'layers', f'summary_{key}.csv')
        if os.path.exists(p):
            out[f'layers/summary_{key}.csv'] = json.loads(pd.read_csv(p).to_json(orient='records', force_ascii=False))
    json.dump(out, open(os.path.join(R, 'numbers.json'), 'w', encoding='utf-8'), ensure_ascii=False, indent=1)
    print('tables collected:', len(out))


if __name__ == '__main__':
    main()
