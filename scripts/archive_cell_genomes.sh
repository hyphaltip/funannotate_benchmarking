#!/usr/bin/bash
# Move selected genomes' outputs out of a benchmark cell so the next
# `launch/run_cell.sbatch <cell>` (-resume) redoes only those genomes.
#
# Moves, per genome, into archive/<cell>.<tag>/ (outside runs/, so
# collect_metrics.py / compare_predictions.py do not see it as a cell):
#   runs/<cell>/genome_annotation/<genome>
#   runs/<cell>/genome_annotation_training/<genome>
#   runs/<cell>/rnaseq_data/<Genus_species>.{trinity-GG.fasta,stringtie.gtf,junctions.bed}
# rnaseq_data/ is RNASEQ_PREPARE's storeDir: while those files exist it never
# re-runs Trinity, regardless of -resume or changed reads.
#
# Usage: scripts/archive_cell_genomes.sh <cell> <tag> <genome_dir_name>...
#   scripts/archive_cell_genomes.sh v1.9.0-rc1_conda_rust badreads20260924 \
#       Rhizopus_microsporus_ATCC_52813 Batrachochytrium_dendrobatidis_JAM81
set -euo pipefail

BENCH_ROOT="/bigdata/stajichlab/shared/projects/BFD/Funannotate_benchmarking"
CELL="${1:?usage: archive_cell_genomes.sh <cell> <tag> <genome>...}"
TAG="${2:?tag required}"
shift 2
[ "$#" -gt 0 ] || { echo "[ERROR] no genomes given" >&2; exit 1; }

SRC="$BENCH_ROOT/runs/$CELL"
DST="$BENCH_ROOT/archive/$CELL.$TAG"
[ -d "$SRC" ] || { echo "[ERROR] no such cell: $SRC" >&2; exit 1; }
mkdir -p "$DST/genome_annotation" "$DST/genome_annotation_training" "$DST/rnaseq_data"

for G in "$@"; do
    # species tag = first two underscore-separated tokens (Genus_species)
    SP=$(echo "$G" | cut -d_ -f1,2)
    for d in genome_annotation genome_annotation_training; do
        if [ -e "$SRC/$d/$G" ]; then
            mv "$SRC/$d/$G" "$DST/$d/"
            echo "[INFO] $CELL: moved $d/$G"
        fi
    done
    for ext in trinity-GG.fasta stringtie.gtf junctions.bed; do
        if [ -e "$SRC/rnaseq_data/$SP.$ext" ]; then
            mv "$SRC/rnaseq_data/$SP.$ext" "$DST/rnaseq_data/"
            echo "[INFO] $CELL: moved rnaseq_data/$SP.$ext"
        fi
    done
done
echo "[INFO] archived to $DST"
