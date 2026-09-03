#!/usr/bin/bash -l
#SBATCH -N 1 -n 1 -c 2 --mem 8gb --time 01:00:00
#SBATCH --partition=short
#SBATCH --job-name=fabench_validate
#SBATCH --output=/bigdata/stajichlab/shared/projects/BFD/Funannotate_benchmarking/runs/%x.%j.log

# Validation gate (DESIGN.md "Validation gate"): diff a cell's ALREADY-COMPLETED
# genome_annotation/ output against Fungi_BFD's known-good results for the same
# genome tag. This does not launch the pipeline -- it only inspects output that
# launch/run_cell.sbatch has already produced, so run it once enough genomes
# have finished predict in a cell (see scripts/validate_harness.py --n-genomes).
#
# Submit:
#   sbatch launch/validate_harness.sbatch <cell>
#   sbatch launch/validate_harness.sbatch v1.8.17_conda --genomes Aspergillus_fumigatus_Af293

set -euo pipefail

BENCH_ROOT="/bigdata/stajichlab/shared/projects/BFD/Funannotate_benchmarking"
cd "$BENCH_ROOT"

CELL="${1:-${CELL:-}}"
[ -n "$CELL" ] || { echo "[ERROR] usage: sbatch launch/validate_harness.sbatch <cell> [extra validate_harness.py args]" >&2; exit 1; }
shift || true

module load miniconda3 gffcompare

python3 scripts/validate_harness.py --cell "$CELL" "$@"
