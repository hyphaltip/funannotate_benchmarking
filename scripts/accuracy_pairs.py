#!/usr/bin/env python3
"""Paired accuracy comparison vs RefSeq between cells, from compare_predictions.py TSVs.

Excludes non-RefSeq references (accession GCA_*). For each pair, uses genomes
compared in both cells, reports median A, median B, median paired difference
(B-A), n B>A / n B<A, and the Wilcoxon signed-rank p-value.

Usage: accuracy_pairs.py samples.csv tsv1,tsv2,... cellA:cellB [...]
"""
import csv, sys, statistics as S, collections
from scipy.stats import wilcoxon

samples = {r['accession']: r for r in csv.DictReader(open(sys.argv[1]))}
rows = []
for f in sys.argv[2].split(','):
    rows += list(csv.DictReader(open(f), delimiter='\t'))
pairs = [p.split(':') for p in sys.argv[3:]]

METRICS = [
    ('gene_ratio', lambda r: float(r['query_gene_count']) / float(r['ref_gene_count'])),
    ('exon_sn', 'gffcompare_exon_level_sensitivity'), ('exon_pr', 'gffcompare_exon_level_precision'),
    ('intron_sn', 'gffcompare_intron_level_sensitivity'), ('intron_pr', 'gffcompare_intron_level_precision'),
    ('ichain_sn', 'gffcompare_intron_chain_level_sensitivity'), ('ichain_pr', 'gffcompare_intron_chain_level_precision'),
    ('tx_sn', 'gffcompare_transcript_level_sensitivity'), ('tx_pr', 'gffcompare_transcript_level_precision'),
    ('locus_sn', 'gffcompare_locus_level_sensitivity'), ('locus_pr', 'gffcompare_locus_level_precision'),
    ('missed_loci_pct', 'gffcompare_missed_loci_pct'),
    ('ref_not_recovered_pct', lambda r: 100 * float(r['overlap_ref_genes_not_recovered']) / float(r['overlap_ref_genes'])),
]

def val(r, m):
    try:
        return m(r) if callable(m) else float(r[m])
    except (ValueError, KeyError, ZeroDivisionError):
        return None

by = collections.defaultdict(dict)
for r in rows:
    if r.get('compare_status') != 'compared' or r['accession'].startswith('GCA_'):
        continue
    by[r['cell']][r['genome']] = r

subset = None
if len(sys.argv) > 3 and sys.argv[-1].startswith('--rnaseq='):
    pass

print('\t'.join(['cellA', 'cellB', 'metric', 'n', 'medA', 'medB', 'med_diff_B-A', 'B>A', 'B<A', 'p']))
for a, b in pairs:
    g = sorted(set(by[a]) & set(by[b]))
    for name, m in METRICS:
        xs, ys = [], []
        for x in g:
            va, vb = val(by[a][x], m), val(by[b][x], m)
            if va is not None and vb is not None:
                xs.append(va); ys.append(vb)
        if not xs:
            continue
        d = [y - x for x, y in zip(xs, ys)]
        try:
            p = wilcoxon(xs, ys).pvalue if len(xs) >= 6 and any(d) else float('nan')
        except ValueError:
            p = float('nan')
        print('\t'.join(map(str, [a, b, name, len(xs), f'{S.median(xs):.3f}', f'{S.median(ys):.3f}',
              f'{S.median(d):+.3f}', sum(v > 0 for v in d), sum(v < 0 for v in d), f'{p:.2g}'])))
