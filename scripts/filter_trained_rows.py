# Drop rows for RNA-seq genomes whose cell has no PASA training models, so a
# pair compares only genomes trained in both cells. RNA-seq genome = has a
# training dir in v1.9.0-rc1_container. Writes filtered TSV to stdout.
import csv, os, sys
def trained(cell, g):
    f = f'runs/{cell}/genome_annotation_training/{g}/training/funannotate_train.pasa.gff3'
    if not os.path.exists(f) or os.path.getsize(f) == 0:
        return False
    with open(f) as fh:
        return any('\tmRNA\t' in l for l in fh)
rnaseq = set(os.listdir('runs/v1.9.0-rc1_container/genome_annotation_training'))
r = csv.DictReader(open(sys.argv[1]), delimiter='\t')
w = csv.DictWriter(sys.stdout, r.fieldnames, delimiter='\t', lineterminator='\n')
w.writeheader()
drop = 0
for row in r:
    if row['genome'] in rnaseq and not trained(row['cell'], row['genome']):
        drop += 1
        print('drop', row['cell'], row['genome'], file=sys.stderr)
        continue
    w.writerow(row)
print('dropped', drop, file=sys.stderr)
