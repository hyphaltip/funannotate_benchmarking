#!/usr/bin/env python3
"""Paired per-genome timing comparison between cells, per stage.

Reads a collect_metrics.py --tasks-out TSV. Uses only role=='final' rows
(last COMPLETED attempt per task). For each (cellA, cellB, stage) pair, keeps
genomes that have a final row in both cells, and reports n, median of each,
median of per-genome ratio B/A, and the Wilcoxon signed-rank p-value.
"""
import csv, sys, statistics, collections, itertools
from scipy.stats import wilcoxon

tasks = sys.argv[1]
pairs = [p.split(':') for p in sys.argv[2:]]
STAGES = ['RNASEQ_PREPARE', 'FUNANNOTATE_TRAIN', 'FUNANNOTATE_PREDICT']

d = collections.defaultdict(dict)   # (cell, stage) -> genome -> (wall_h, cpu_h, cpus)
for r in csv.DictReader(open(tasks), delimiter='\t'):
    if r['role'] != 'final':
        continue
    st = r['process'].split(':')[-1]
    if st not in STAGES:
        continue
    d[(r['cell'], st)][r['genome']] = (float(r['realtime_s']) / 3600, float(r['cpu_hours'] or 0), r['cpus'])

print('\t'.join(['cellA', 'cellB', 'stage', 'n', 'medA_wall_h', 'medB_wall_h', 'med_ratio_wall_B/A',
                 'p_wall', 'medA_cpu_h', 'medB_cpu_h', 'med_ratio_cpu_B/A', 'p_cpu', 'sumA_wall_h', 'sumB_wall_h']))
for a, b in pairs:
    for st in STAGES:
        A, B = d.get((a, st), {}), d.get((b, st), {})
        g = sorted(set(A) & set(B))
        if not g:
            print(f'{a}\t{b}\t{st}\t0')
            continue
        aw = [A[x][0] for x in g]; bw = [B[x][0] for x in g]
        ac = [A[x][1] for x in g]; bc = [B[x][1] for x in g]
        rw = [y / x for x, y in zip(aw, bw) if x > 0]
        rc = [y / x for x, y in zip(ac, bc) if x > 0]
        def p(x, y):
            try:
                return f'{wilcoxon(x, y).pvalue:.2g}' if len(x) >= 6 else 'NA'
            except ValueError:
                return 'NA'
        print('\t'.join(map(str, [a, b, st, len(g),
              f'{statistics.median(aw):.2f}', f'{statistics.median(bw):.2f}',
              f'{statistics.median(rw):.2f}' if rw else 'NA', p(aw, bw),
              f'{statistics.median(ac):.2f}', f'{statistics.median(bc):.2f}',
              f'{statistics.median(rc):.2f}' if rc else 'NA', p(ac, bc),
              f'{sum(aw):.1f}', f'{sum(bw):.1f}'])))
