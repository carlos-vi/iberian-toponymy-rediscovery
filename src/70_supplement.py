"""Builds supplement/supplementary.md from supplement/supplementary_template.md: every line of the form
{{table NAME}} is replaced by a Markdown table drawn from the files in results/ (builders below).
"""
import os, re
import numpy as np, pandas as pd
from lib import paths
from lib.reference import AREAS, EXCLUDED, TARGETS
from lib.mapbase import pretty

R = paths.RESULTS
MAN = os.path.join(paths.ROOT, 'supplement')
VNAME = {'sig': 'significance', 'eff': 'effect size', 'two': 'two-level'}


def read(rel, **kw):
    return pd.read_csv(os.path.join(R, rel), keep_default_na=False, **kw)


def md(df, floatfmt='{:.2f}'):
    def cell(x):
        if isinstance(x, (float, np.floating)):
            if np.isnan(x):
                return '–'
            return f'{int(x):,}' if float(x).is_integer() and abs(x) >= 1000 else floatfmt.format(x)
        if isinstance(x, (int, np.integer)):
            return f'{x:,}' if abs(x) >= 10000 else str(x)
        return str(x).replace('|', '/').replace('\n', ' ')
    cols = list(df.columns)
    lines = ['| ' + ' | '.join(map(str, cols)) + ' |', '|' + '---|' * len(cols)]
    lines += ['| ' + ' | '.join(cell(v) for v in row) + ' |' for row in df.itertuples(index=False)]
    return '\n'.join(lines)


def as_int(x):
    """Integer display ('–' for missing), with thousands separators."""
    try:
        v = float(x)
    except (TypeError, ValueError):
        return '–'
    return '–' if np.isnan(v) else f'{int(round(v)):,}'


def pv(p):
    p = float(p)
    return '< 0.001' if p < 0.001 else f'{p:.3f}'


# ------------------------------------------------------------------ data
def t_corpus_sizes():
    t = read('data/corpus_sizes.csv')
    return md(t.rename(columns={'corpus': 'Corpus', 'records': 'Records', 'distinct_words': 'Distinct words',
                                'opaque_records': 'Records with ≥ 1 opaque word', 'spain': 'Spain', 'portugal': 'Portugal',
                                'france': 'France', 'andorra': 'Andorra'}))


def t_density():
    t = read('data/corpus_density.csv')
    t = t[~t.region.isin(['other', 'TOTAL', 'France and Andorra'])]
    out = pd.DataFrame({'Region': t.region})
    for k in ['S', 'I', 'SI', 'M', 'G']:
        out[k] = [f'{int(float(a)):,} ({float(b):.0f})' if b != '' else '–' for a, b in zip(t[k], t[f'{k}_per_1000km2'])]
    return md(out)


def t_areas():
    rows = [(i + 1, a, ' '.join(r[3:] if r.startswith('ES-') else r for r in rs))
            for i, (a, rs) in enumerate(AREAS.items())]
    t = pd.DataFrame(rows, columns=['#', 'Reference area',
                                    'Units (Spanish provinces by ISO code; PT- Portuguese districts; FR- French departments; AD Andorra)'])
    return md(t) + '\n\nNot scored (transition areas): ' + ', '.join(sorted(EXCLUDED)) + '.'


def t_targets():
    rows = [(s, ' '.join(pretty(m) for m in ms), len(exp)) for ms, exp, s in TARGETS]
    return md(pd.DataFrame(rows, columns=['Stratum', 'Elements', 'Number of expected provinces']))


# ------------------------------------------------------------------ layers
def t_main_metric():
    t = read('layers/main_metric.csv')
    t['variant'] = t.variant.map(VNAME); t['p'] = t.p.map(pv); t['corpus'] = t.corpus.replace({'SI': 'S + I'})
    return md(t[['corpus', 'variant', 'mean_recovered', 'chance', 'excess', 'p', 'mean_accepted_layers', 'recovered_k18',
                 'min_recovered', 'max_recovered']].rename(columns={
        'corpus': 'Corpus', 'variant': 'Selection', 'mean_recovered': 'Areas recovered (mean over k)', 'chance': 'Chance',
        'excess': 'Excess', 'p': 'p', 'mean_accepted_layers': 'Accepted layers (mean)', 'recovered_k18': 'Recovered at k = 18',
        'min_recovered': 'Min', 'max_recovered': 'Max'}))


def t_summary(key):
    s = read(f'layers/summary_{key}.csv')
    s['cell'] = [f'{r} ({c:.2f}; {n})' for r, c, n in zip(s.recovered, s.null_mean, s.n_accepted)]
    t = s.pivot(index='k', columns='variant', values='cell')[['sig', 'eff', 'two']].reset_index()
    return md(t.rename(columns={'k': 'k', **{v: VNAME[v] for v in VNAME}}))


def t_other_corpora():
    rows = []
    for key in ['M', 'I', 'G']:
        s = read(f'layers/summary_{key}.csv')
        for v in ['sig', 'eff', 'two']:
            r = s[(s.variant == v) & (s.k == 18)].iloc[0]
            rows.append(dict(Corpus=key, Selection=VNAME[v], **{'Areas scored': r.n_areas, 'Accepted layers': r.n_accepted,
                             'Recovered': r.recovered, 'Chance': r.null_mean, 'p': pv(r.p), 'Areas': r.areas}))
    return md(pd.DataFrame(rows))


def t_robustness():
    t = read('layers/robustness_layers.csv')
    t['cell'] = [f'{a:.1f} / {c:.1f} ({pv(p)})' for a, c, p in zip(t.recovered, t.chance, t.p)]
    t['col'] = t.corpus.replace({'SI': 'S + I'}) + ' ' + t.variant.map(VNAME)
    out = t.pivot_table(index=['check', 'setting'], columns='col', values='cell', aggfunc='first', sort=False).reset_index()
    order = [f'{k} {VNAME[v]}' for k in ['S', 'S + I'] for v in ['sig', 'eff', 'two']]
    return md(out[['check', 'setting'] + [c for c in order if c in out.columns]].rename(columns={'check': 'Check', 'setting': 'Setting'}))


def t_replication():
    t = read('layers/robustness_replication.csv')
    t['variant'] = t.variant.map(VNAME)
    return md(t.rename(columns={'variant': 'Selection', 'areas_in_2plus_corpora': 'Areas recovered in ≥ 2 of S, M, I, G',
                                'expected_by_chance': 'Expected by chance (independent corpora)',
                                'max_by_chance': 'Ceiling of chance (any dependence)'}))


def t_territories():
    t = read('layers/SI_two_k18_territories.csv')
    t['parent_top'] = t.parent_top.map(lambda s: ' '.join(pretty(x) for x in s.split()[:5]))
    cols = ['parent_top', 'n_cells', 'nrec', 'n_elements', 'k2'] + (['accepted2'] if 'accepted2' in t else [])
    for c in cols[1:]:
        t[c] = t[c].map(as_int)
    t = t[cols].rename(columns={'parent_top': 'First-level layer (leading elements)', 'n_cells': 'Cells', 'nrec': 'Records',
                                'n_elements': 'Elements', 'k2': 'k (second level)', 'accepted2': 'Accepted'})
    return md(t)


def t_layers_SI():
    t = read('layers/SI_two_k18_layers.csv')
    t = t[t.accepted.astype(str) == 'True']
    t['top'] = t.top.map(lambda s: ' '.join(pretty(x) for x in s.split()[:8]))
    t['recovers'] = t.recovers.replace('', '–')
    return md(t[['level', 'top', 'core_cells', 'lat', 'lon', 'best_area', 'purity', 'auc', 'recovers']].rename(columns={
        'level': 'Level', 'top': 'Leading elements', 'core_cells': 'Core cells', 'lat': 'Lat', 'lon': 'Lon',
        'best_area': 'Best-matching area', 'purity': 'Purity', 'auc': 'AUC', 'recovers': 'Recovers'}))


def t_selection_regions():
    t = read('layers/selection_by_region_S.csv').rename(columns={'focus_region': 'Region of maximum enrichment',
                                                                   'sig': 'significance', 'eff': 'effect size'})
    return md(t)


# ------------------------------------------------------------------ scan
def t_scan_summary():
    t = read('scan/scan_summary.csv')
    return md(t.rename(columns={'corpus': 'Corpus', 'tested': 'Elements tested', 'significant': 'Significant (q < 0.05)',
                                'share': 'Share'}))


def t_validation():
    t = read('scan/validation_by_stratum.csv')
    s = read('scan/validation_summary.csv').iloc[0]
    return md(t.rename(columns={'stratum': 'Stratum', 'n': 'Elements', 'present': 'Tested', 'hits': 'Recovered'})) + (
        f'\n\nTotal: {s.hits} of {s.present} tested elements ({float(s.recall_present):.0%}); {s.hits} of {s.targets} listed '
        f'({float(s.recall_all):.0%}). Random elements of similar size: {float(s.random_recall):.1%} (n = {s.random_n}).')


def t_validation_elements():
    t = read('scan/validation.csv')
    t['element'] = t.element.map(pretty)
    t['q'] = pd.to_numeric(t.q, errors='coerce').map(lambda x: '–' if np.isnan(x) else pv(x))
    t['radius'] = t.radius.map(as_int)
    return md(t[['stratum', 'element', 'present', 'hit', 'q', 'region', 'radius', 'examples']].rename(columns={
        'stratum': 'Stratum', 'element': 'Element', 'present': 'Tested', 'hit': 'Recovered', 'q': 'q', 'region': 'Centre',
        'radius': 'Radius (km)', 'examples': 'Examples inside the cluster'}))


def t_top50():
    t = read('scan/top50_classified.csv')
    t['element'] = t.element.map(pretty)
    t['examples'] = t.examples.map(lambda s: ' '.join(s.split()[:5]))
    t['radius'] = t.radius.map(as_int); t['words_in'] = t.words_in.map(as_int)
    return md(t[['rank_merged', 'element', 'region', 'radius', 'words_in', 'rr_in', 'family', 'class', 'examples']].rename(columns={
        'rank_merged': 'Rank', 'element': 'Element', 'region': 'Centre', 'radius': 'Radius (km)', 'words_in': 'Words inside',
        'rr_in': 'Relative risk', 'family': 'Attribution', 'class': 'Class', 'examples': 'Examples'}))


def t_cluster_locations():
    t = read('scan/cluster_locations.csv')
    t = t[t.present.astype(str) == 'True'].copy()
    t['element'] = t.element.map(pretty)
    t['q'] = pd.to_numeric(t.q).map(pv)
    t['radius'] = t.radius.map(as_int); t['km_from_Toledo'] = t.km_from_Toledo.map(as_int)
    t['corpus'] = t.corpus.str[:-1].replace({'SI': 'S + I'})
    return md(t[['corpus', 'element', 'q', 'region', 'lat', 'lon', 'radius', 'km_from_Toledo', 'in_expected_place']].rename(columns={
        'corpus': 'Corpus', 'element': 'Element', 'q': 'q', 'region': 'Centre', 'lat': 'Lat', 'lon': 'Lon',
        'radius': 'Radius (km)', 'km_from_Toledo': 'km from Toledo', 'in_expected_place': 'Where expected'}))


def t_secondary():
    t = read('scan/secondary_clusters.csv')
    t['element'] = t.element.map(pretty); t['p'] = t.p.map(pv)
    t['examples'] = t.examples.map(lambda s: ' '.join(s.split()[:6]))
    t['radius'] = t.radius.map(as_int); t['c'] = t.c.map(as_int)
    t['corpus'] = t.corpus.str[:-1].replace({'SI': 'S + I'})
    return md(t[['corpus', 'element', 'cluster', 'region', 'lat', 'lon', 'radius', 'c', 'llr', 'p', 'examples']].rename(columns={
        'corpus': 'Corpus', 'element': 'Element', 'cluster': 'Cluster', 'region': 'Centre', 'lat': 'Lat', 'lon': 'Lon',
        'radius': 'Radius (km)', 'c': 'Units inside', 'llr': 'LLR', 'p': 'p', 'examples': 'Examples'}))


def t_eque_scan():
    r = read('scan/eque_rarity.csv').iloc[0]; e = read('scan/eque_exact_null.csv').iloc[0]
    rows = [('Endings of 3–4 letters with 20–60 units', f'{int(r.endings_3_4_letters_20_60_units)}'),
            ('… with a significant cluster (share)', f'{float(r.significant_share):.2f}'),
            ('… at least as concentrated as -eque (-eque included)', f'{int(r.as_concentrated_as_eque)}'),
            ('Relative risk inside the -eque cluster', f'{float(r.rr_eque):.2f}'),
            ('Rank of -eque by concentration among significant clusters with ≥ 8 words',
             f'{int(r.eque_rank_by_concentration)} of {int(r.significant_with_8_words):,}'),
            ('Exact null (sampling without replacement): replicates', f'{int(e.replicates):,}'),
            ('Observed LLR / largest null LLR', f'{float(e.llr_observed):.2f} / {float(e.null_max):.2f}'),
            ('Exact p', f'{float(e.p):.4f}')]
    return md(pd.DataFrame(rows, columns=['Quantity', 'Value']))


# ------------------------------------------------------------------ wedge
def t_iello():
    t = read('wedge/iello_segments.csv')
    return md(t.rename(columns={'corpus': 'Corpus', 'west_old': '-iello west of 5° W', 'west_new': '-illo west of 5° W',
                                'west_share': 'Share', 'castile_old': '-iello 4.5–1.5° W', 'castile_new': '-illo 4.5–1.5° W',
                                'castile_share': 'Share '}), floatfmt='{:.3f}')


def t_iello_prov():
    t = read('wedge/iello_west_by_province.csv')
    return md(t.rename(columns={'region': 'Province', 'old': '-iello', 'new': '-illo', 'corpus': 'Corpus', 'share': 'Share'}),
              floatfmt='{:.3f}')


def t_fh():
    t = read('wedge/fh_segments.csv')
    segs = ['west', 'centre', 'east', 'east_spain', 'south']
    out = pd.DataFrame({'Corpus': t.corpus})
    for s in segs:
        out[s.replace('_', ' ')] = [f'{float(v):.2f} ({int(n)})' for v, n in zip(t[f'{s}_share_f'], t[f'{s}_n'])]
    return md(out)


def t_fh_loo():
    t = read('wedge/fh_leave_one_pair_out_extremes.csv')
    return md(t.rename(columns={'corpus': 'Corpus', 'west_min': 'West, lowest', 'east_min': 'East, lowest',
                                'centre_max': 'Centre, highest', 'south_max': 'South, highest'}))


# ------------------------------------------------------------------ -ena
def t_ena_lift():
    a = read('ena/lift_all_provinces.csv'); b = read('ena/without_basque.csv')
    out = pd.DataFrame({'Corpus': a.corpus.str[:-1],
                        '-ena types (in core Baetica)': [f'{n} ({k})' for n, k in zip(a.ena_types, a.ena_in_baetica)],
                        '-ena lift, all provinces [95 % CI]': [f'{x:.2f} [{l:.2f}–{h:.2f}]' for x, l, h in zip(a.ena_lift, a.ena_lo, a.ena_hi)],
                        '-ena lift, without Basque provinces': [f'{x:.2f} [{l:.2f}–{h:.2f}]' for x, l, h in zip(b.ena_lift, b.ena_lo, b.ena_hi)],
                        '-ana types (in core Baetica)': [f'{n} ({k})' for n, k in zip(a.ana_types, a.ana_in_baetica)],
                        '-ana lift, all provinces': [f'{x:.2f} [{l:.2f}–{h:.2f}]' for x, l, h in zip(a.ana_lift, a.ana_lo, a.ana_hi)]})
    return md(out)


def t_ena_groups():
    b = read('ena/without_basque.csv')
    groups = [c[:-len(' share')] for c in b.columns if c.endswith(' share')]
    rows = []
    for g in groups:
        rows.append([g] + [f'{float(r[f"{g} share"]):.2f} ({int(r[f"{g} ena"])}/{int(r[f"{g} ena"]) + int(r[f"{g} ana"])})'
                           for _, r in b.iterrows()])
    t = pd.DataFrame(rows, columns=['Provinces'] + [c[:-1] for c in b.corpus])
    t2 = b[['corpus', 'n_provinces', 'spearman_years', 'spearman_p', 'spearman_years_latitude', 'partial_years_given_latitude']]
    t2 = t2.assign(corpus=t2.corpus.str[:-1]).rename(columns={
        'corpus': 'Corpus', 'n_provinces': 'Provinces', 'spearman_years': 'Spearman ρ (years of rule)', 'spearman_p': 'p',
        'spearman_years_latitude': 'ρ (years, latitude)', 'partial_years_given_latitude': 'Partial ρ given latitude'})
    return md(t) + '\n\n' + md(t2, floatfmt='{:.3f}')


def t_ena_provinces():
    t = read('ena/provinces_SI.csv')
    t = t.rename(columns={t.columns[0]: 'Province'})
    return md(t[['Province', 'ena', 'ana', 'share', 'years', 'lat']].rename(columns={
        'ena': '-ena', 'ana': '-ana', 'share': 'Share', 'years': 'Years of rule', 'lat': 'Mean latitude'}))


def t_ena_misc():
    b = read('ena/basque_homograph.csv').iloc[0]; c = read('ena/core_vs_nasrid.csv').iloc[0]
    r = read('ena/repeated_summary_SI.csv').iloc[0]
    rows = [('-ena types in I', b.ena_types_I), ('… in Bizkaia and Gipuzkoa', b.in_bizkaia_gipuzkoa),
            ('Examples (Bizkaia)', b.examples),
            ('Core Baetica: -ena / -ana (S + I)', f'{c.core_ena} / {c.core_ana}'),
            ('Nasrid provinces: -ena / -ana (S + I)', f'{c.nasrid_ena} / {c.nasrid_ana}'),
            ('Fisher exact test, core Baetica vs Nasrid: odds ratio, p', f'{c.fisher_odds}, {pv(c.fisher_p)}'),
            ('-ena types without Basque provinces (S + I) / distinct words', f'{r.types} / {r.distinct_words}')]
    return md(pd.DataFrame(rows, columns=['Quantity', 'Value']))


# ------------------------------------------------------------------ Arabic names
def t_frontier():
    t = read('arabic/frontier_summary.csv')
    return md(t.rename(columns={'corpus': 'Corpus', 'north_per_mille': 'North (‰)', 'south_per_mille': 'South (‰)',
                                'threshold': 'Threshold at 50 % (‰)', 'west_centre_km_south_of_Duero': 'km south of the Duero, 6.5–3.5° W',
                                'east_km_north_of_Duero': 'km north of the Duero, 3.0–2.2° W',
                                'mean_km_to_strongholds': 'Mean km to strongholds', 'mean_km_to_Duero': 'Mean km to Duero'}),
              floatfmt='{:.1f}')


def t_frontier_segments():
    t = read('arabic/frontier_segments.csv')
    t['cell'] = [f'{a:.0f} / {b:.0f}' for a, b in zip(t.km_to_strongholds, t.km_to_Duero)]
    t['thr'] = (t.threshold_fraction * 100).astype(int).astype(str) + ' %'
    out = t.pivot_table(index=['corpus', 'thr'], columns='segment', values='cell', aggfunc='first', sort=False).reset_index()
    return md(out.rename(columns={'corpus': 'Corpus', 'thr': 'Threshold'}))


def t_dose():
    t = read('arabic/dose_response.csv')
    return md(t.rename(columns={'corpus': 'Corpus', 'provinces': 'Provinces', 'r_years': 'r (years of rule)',
                                'r_years_latitude': 'r (years, latitude)', 'partial_r_given_latitude': 'Partial r given latitude'}),
              floatfmt='{:.3f}')


# ------------------------------------------------------------------ -aque/-eque
def t_eque():
    t = read('aque_eque/tests.csv')
    out = pd.DataFrame({'Corpus': t.corpus.str[:-1],
                        '-aque/-eque types (centre)': [f'{a} ({b})' for a, b in zip(t.types_aque_eque, t.in_centre)],
                        'Lift (types)': t.lift,
                        'Centre a·e / i·u': [f'{a} / {b}' for a, b in zip(t.centre_ae, t.centre_iu)],
                        'South a·e / i·u': [f'{a} / {b}' for a, b in zip(t.south_ae, t.south_iu)],
                        'Fisher p (types)': t.fisher_p.map(pv),
                        'Words (centre)': [f'{a} ({b})' for a, b in zip(t.words_aque_eque, t.words_in_centre)],
                        'Lift (words)': t.lift_words,
                        'Fisher p (words)': t.fisher_p_words.map(pv)})
    return md(out)


def t_eque_names():
    t = read('aque_eque/aque_eque_names_SIo.csv')
    t = t.sort_values(['r', 't'])
    return md(t[['t', 'r', 'records', 'example']].rename(columns={'t': 'Word', 'r': 'Province', 'records': 'Records',
                                                                   'example': 'Example'}))


def t_conquest():
    from lib.reference import CONQ
    g = {}
    for r, y in CONQ.items():
        g.setdefault(y, []).append(r)
    rows = [(y, y - 714, ' '.join(sorted(rs))) for y, rs in sorted(g.items())]
    return md(pd.DataFrame(rows, columns=['Year of definitive Christian conquest (approx.)', 'Years of Islamic rule',
                                          'Provinces / districts']))


def t_fh_pairs():
    import importlib
    W = importlib.import_module('40_wedge')
    show = lambda pat: ' or '.join(f'`{x}`' for x in pat.split('|'))
    rows = [(n, show(f.pattern), show(h.pattern)) for n, f, h in W.PAIRS]
    return md(pd.DataFrame(rows, columns=['Latin root', 'F- form (regular expression)', 'h- form (regular expression)']))


def t_markers():
    from lib.reference import STRICT, LEXAR
    alts = lambda pat: ', '.join(f'`{x}`' for x in pat[2:-1].split('|'))   # '^(a|b|c)' -> `a`, `b`, `c`
    return md(pd.DataFrame([('Arabic markers (word beginnings, opaque words only)', alts(STRICT.pattern)),
                            ('Excluded loanwords and lookalikes', alts(LEXAR.pattern))],
                           columns=['List', 'Alternatives (regular expressions anchored at the start of the word)']))


BUILDERS = {k[2:]: v for k, v in globals().items() if k.startswith('t_') and callable(v)}
BUILDERS['summary_S'] = lambda: t_summary('S'); BUILDERS['summary_SI'] = lambda: t_summary('SI')


def main():
    tpl = open(os.path.join(MAN, 'supplementary_template.md'), encoding='utf-8').read()
    def sub(m):
        name = m.group(1)
        if name not in BUILDERS:
            raise SystemExit(f'unknown table: {name}')
        return BUILDERS[name]()
    out = re.sub(r'^\{\{table (\w+)\}\}[ \t]*$', sub, tpl, flags=re.M)
    open(os.path.join(MAN, 'supplementary.md'), 'w', encoding='utf-8').write(out)
    print('supplementary.md written,', len(re.findall(r'^\| ', out, re.M)), 'table rows')


if __name__ == '__main__':
    main()
