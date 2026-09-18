#!/usr/bin/bash -l
#SBATCH -N 1 -n 1 -c 1 --mem 1gb --time 00:10:00
#SBATCH --partition=short
#SBATCH --job-name=fabench_launch_all
#SBATCH --output=/bigdata/stajichlab/shared/projects/BFD/Funannotate_benchmarking/runs/%x.%j.log

# Submit launch/run_cell.sbatch for every cell in conf/cells.tsv. This job
# itself does no pipeline work -- it just fans out one sbatch submission per
# cell, so each cell's nextflow head process runs as its own independent
# SLURM job (its own -resume state, its own failure blast radius).
#
# Submit:
#   sbatch launch/run_all_cells.sbatch
#
# To (re)launch a subset, use launch/run_cell.sbatch directly instead:
#   sbatch launch/run_cell.sbatch v1.8.17_container

set -euo pipefail

BENCH_ROOT="/bigdata/stajichlab/shared/projects/BFD/Funannotate_benchmarking"
cd "$BENCH_ROOT"

CELLS_TSV="conf/cells.tsv"
# Skip comment lines (leading #, e.g. the retired-beta.10 note) and the
# header row itself, wherever it falls -- not just "line 1".
awk -F'\t' '$1!="" && $1!~/^#/ && $1!="cell" {print $1}' "$CELLS_TSV" | while read -r CELL; do
    [ -n "$CELL" ] || continue
    echo "[INFO] submitting cell: $CELL"
    sbatch --job-name="fabench_${CELL}" launch/run_cell.sbatch "$CELL"
done
