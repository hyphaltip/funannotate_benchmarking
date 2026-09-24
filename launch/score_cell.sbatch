#!/usr/bin/bash -l
#SBATCH -N 1 -n 1 -c 16 --mem 32gb --time 12:00:00
#SBATCH --partition=epyc
#SBATCH --job-name=fabench_score
#SBATCH --output=/bigdata/stajichlab/shared/projects/BFD/Funannotate_benchmarking/runs/%x.%j.log

# BUSCO (protein mode) on every delivered predict_results/*.proteins.fa of ONE
# cell. Idempotent: run_busco_completeness.py skips genomes that already have a
# short_summary. gffcompare scoring (compare_predictions.py) is fast and is run
# separately across all cells once BUSCO is done.
#
# Submit:
#   sbatch launch/score_cell.sbatch <cell> [cells_dir]
#   sbatch launch/score_cell.sbatch v1.9.0-rc1_container.sif20260921 archive

set -euo pipefail

BENCH_ROOT="/bigdata/stajichlab/shared/projects/BFD/Funannotate_benchmarking"
cd "$BENCH_ROOT"

CELL="${1:?usage: sbatch launch/score_cell.sbatch <cell> [cells_dir]}"
CELLS_DIR="${2:-runs}"

module load busco/6.0.0

python3 scripts/run_busco_completeness.py \
    --cells-dir "$CELLS_DIR" --cells "$CELL" --samples samples.csv \
    --cpus "${SLURM_CPUS_PER_TASK:-16}"
