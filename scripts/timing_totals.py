#!/usr/bin/env python3
"""Per-species end-to-end timing (RNASEQ_PREPARE + TRAIN + PREDICT), paired across cells.

RNASEQ_PREPARE is keyed by species tag (e.g. Aspergillus_fumigatus) and runs
once per species; TRAIN and PREDICT are keyed by genome (Aspergillus_fumigatus_Af293).
Where Trinity runs (RNASEQ_PREPARE vs inside TRAIN) differs between cells, so only
the sum is comparable. A species is kept for a pair only when every genome of that
species has a final PREDICT row in both cells, and the species' TRAIN/RNASEQ rows
exist in both or in neither (same work profile).

Usage: timing_totals.py tasks.tsv samples.csv cellA:cellB [...]
"""
import csv, sys, statistics, collections
from scipy.stats import wilcoxon

tasks, samples = sys.argv[1], sys.argv[2]
pairs = [p.split(':') for p in sys.argv[3:]]


# genome tags in traces are "<Genus>_<species>_<strain...>"; species tag is the first two tokens
def sp_of(tag):
    return '_'.join(tag.split('_')[:2])

agg = collections.defaultdict(lambda: collections.defaultdict(lambda: [0.0, 0.0]))  # cell -> sp -> [wall,cpu]
have = collections.defaultdict(lambda: collections.defaultdict(set))  # cell -> sp -> stages seen
pred_genomes = collections.defaultdict(lambda: collections.defaultdict(set))
for r in csv.DictReader(open(tasks), delimiter='\t'):
    if r['role'] != 'final' or 'TRAIN_PREDICT' not in r['process']:
        continue
    st = r['process'].split(':')[-1]
    sp = sp_of(r['genome'])
    a = agg[r['cell']][sp]
    a[0] += float(r['realtime_s']) / 3600
    a[1] += float(r['cpu_hours'] or 0)
    have[r['cell']][sp].add(st if st != 'RNASEQ_PREPARE' else 'RNASEQ_PREPARE')
    if st == 'FUNANNOTATE_PREDICT':
        pred_genomes[r['cell']][sp].add(r['genome'])

print('\t'.join(['cellA', 'cellB', 'n_species', 'n_with_rnaseq', 'medA_wall_h', 'medB_wall_h', 'med_ratio_wall',
                 'p_wall', 'medA_cpu_h', 'medB_cpu_h', 'med_ratio_cpu', 'p_cpu', 'sumA_cpu_h', 'sumB_cpu_h']))
for a, b in pairs:
    keep = []
    for sp in set(pred_genomes[a]) & set(pred_genomes[b]):
        if pred_genomes[a][sp] != pred_genomes[b][sp]:
            continue
        train_a = 'FUNANNOTATE_TRAIN' in have[a][sp] or 'RNASEQ_PREPARE' in have[a][sp]
        train_b = 'FUNANNOTATE_TRAIN' in have[b][sp] or 'RNASEQ_PREPARE' in have[b][sp]
        if train_a != train_b:
            continue
        keep.append((sp, train_a))
    aw = [agg[a][s][0] for s, _ in keep]; bw = [agg[b][s][0] for s, _ in keep]
    ac = [agg[a][s][1] for s, _ in keep]; bc = [agg[b][s][1] for s, _ in keep]
    if not keep:
        print(a, b, 0); continue
    rw = [y / x for x, y in zip(aw, bw) if x > 0]; rc = [y / x for x, y in zip(ac, bc) if x > 0]
    pw = wilcoxon(aw, bw).pvalue if len(keep) >= 6 else float('nan')
    pc = wilcoxon(ac, bc).pvalue if len(keep) >= 6 else float('nan')
    print('\t'.join(map(str, [a, b, len(keep), sum(t for _, t in keep),
          f'{statistics.median(aw):.2f}', f'{statistics.median(bw):.2f}', f'{statistics.median(rw):.2f}', f'{pw:.2g}',
          f'{statistics.median(ac):.2f}', f'{statistics.median(bc):.2f}', f'{statistics.median(rc):.2f}', f'{pc:.2g}',
          f'{sum(ac):.0f}', f'{sum(bc):.0f}'])))
