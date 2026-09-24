#!/usr/bin/bash -l
#SBATCH -N 1 -n 1 -c 2 --mem 16gb --time 12:00:00
#SBATCH --partition=epyc
#SBATCH --job-name=fabench_compare
#SBATCH --output=/bigdata/stajichlab/shared/projects/BFD/Funannotate_benchmarking/runs/%x.%j.log

# gffcompare + bedtools scoring vs RefSeq (scripts/compare_predictions.py) for
# the named cells under one cells dir. Re-run after BUSCO finishes to pick up
# the BUSCO columns.
#
# Submit:
#   sbatch launch/compare_cells.sbatch <cells_dir> <out.tsv> <cell> [<cell> ...]

set -euo pipefail

BENCH_ROOT="/bigdata/stajichlab/shared/projects/BFD/Funannotate_benchmarking"
cd "$BENCH_ROOT"

CELLS_DIR="${1:?usage: compare_cells.sbatch <cells_dir> <out.tsv> <cell>...}"
OUT="${2:?}"
shift 2

module load gffcompare/0.12.10 bedtools/2.30.0

python3 scripts/compare_predictions.py \
    --cells-dir "$CELLS_DIR" --cells "$@" --samples samples.csv --out "$OUT"
