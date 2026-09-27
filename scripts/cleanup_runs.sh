#!/usr/bin/bash
# Delete benchmark run data whose results are already extracted into results/*.tsv.
#
# Dry run by default: prints what would be removed. Pass --apply to delete.
#
# Tier A -- safe now (2026-09-25). Retired/diagnostic cells, pilot, and archive/.
#   Before deleting a retired cell, a small bundle is written to
#   snapshots/retired_cells_20260925/<cell>.tar.zst with each genome's
#   predict_results/*.{gff3,proteins.fa}, busco_completeness short_summary,
#   logs/nextflow traces and samples.csv. Scored tables:
#     results/gene_content_comparison.retired_cells_20260925.tsv
#     results/metrics*.retired_cells_20260925.tsv
#     results/gene_content_comparison.archive_sif20260921.tsv   (archive/*.sif20260921)
#     results/gene_content_comparison.runs_20260924.tsv         (archive/v1.9.0-rc1_conda.trinityfail20260924)
#   archive/*.badreads20260924 hold outputs trained on the wrong RNA-seq reads;
#   they are not scored and are not wanted.
#   runs/pilot is NOT removed: runs/v1.8.17_conda/lib/augustus/3.5 is a symlink
#   into runs/pilot/lib. Only runs/pilot/work goes.
#   Also: runs/v1.8.17_conda/*.stale_* (set aside 2026-09-19, not read by the
#   current run) and the top-level work/ (GeneMark sidecar Nextflow work dir;
#   runs/genemark_sidecar/output/ has .gtf/.mod for all 65 genomes and is kept).
#
# Tier B -- only after the 5 active cells are finished AND scored
#   (BUSCO + compare_predictions.py + collect_metrics.py). Removes each active
#   cell's Nextflow work/ dir. After this, -resume on that cell reruns
#   everything. genome_annotation/, genome_annotation_training/, rnaseq_data/
#   and logs/ are kept.
#
# Tier C -- rebuildable intermediates inside genome_annotation_training/<g>/training/
#   of the 5 active cells. Removed per genome:
#     trinity_gg/  getBestModel/  genome.fasta.gmap/  normalize/
#     pasa/*.transdecoder_dir/   (TransDecoder scratch; pasa/ outputs are kept)
#     hisat2.genome.*.ht2  hisat2.coordSorted.bam
#     trinity.fasta.clean*  trinity.fasta.cln  trinity.fasta.cidx
#   Kept: trinity.fasta, trinity/transcript.alignments.{bam,gff3} (targets of the
#   funannotate_train.* symlinks), funannotate_train.*, pasa.step1.gff3, the rest
#   of pasa/, kallisto.tsv, Trinity-gg.log and logfiles/.
#   A genome is touched only if (a) it has funannotate_train.pasa.gff3 (failed
#   trainings are kept whole for debugging), (b) its predict_results/*.proteins.fa
#   is newer than that file (PREDICT ran after the current TRAIN), and (c) no running or
#   pending SLURM task of that cell names it (a TRAIN retry can resume from
#   trinity_gg/ checkpoints).
#
# Usage:
#   scripts/cleanup_runs.sh A            # dry run, tier A
#   scripts/cleanup_runs.sh A --apply
#   scripts/cleanup_runs.sh B --apply
#   scripts/cleanup_runs.sh C --apply
set -euo pipefail

BENCH_ROOT="/bigdata/stajichlab/shared/projects/BFD/Funannotate_benchmarking"
cd "$BENCH_ROOT"

TIER="${1:?usage: cleanup_runs.sh A|B|C [--apply]}"
APPLY=0
[ "${2:-}" = "--apply" ] && APPLY=1

RETIRED_CELLS=(
    v1.9.0-beta10_conda v1.9.0-beta10_conda_rust
    v1.9.0-beta10_container v1.9.0-beta10_container_rust
    v1.9.0-beta11_conda v1.9.0-beta11_conda_rust v1.9.0-beta11_container_rust
    v1.9.0-beta12_container v1.9.0-beta12_container_rust
    v1.9.0-beta13_blat_container_rust v1.9.0-beta13_gmap_conda_rust
    v1.9.0-beta13_gmap_container_rust v1.9.0-beta13_gmap_pasaT1_container_rust
    v1.8.17_container v1.8.17_pasaiso_conda
    trinity_repro
)
ACTIVE_CELLS=(
    v1.8.17_conda v1.9.0-rc1_conda v1.9.0-rc1_conda_rust
    v1.9.0-rc1_container v1.9.0-rc1_container_rust
    v1.9.0-rc3_container_rust
)
SNAP="snapshots/retired_cells_20260925"

remove() {
    local p="$1"
    [ -e "$p" ] || { echo "[skip] $p (absent)"; return; }
    if [ "$APPLY" -eq 1 ]; then
        echo "[rm]   $p"
        rm -rf -- "$p"
    else
        echo "[dry]  would remove $p"
    fi
}

bundle_cell() {
    local cell="$1" out="$SNAP/$1.tar.zst"
    [ -s "$out" ] && { echo "[keep] $out exists"; return; }
    if [ "$APPLY" -eq 0 ]; then
        echo "[dry]  would write $out"
        return
    fi
    mkdir -p "$SNAP"
    ( cd runs && find "$cell" \
        \( -path "*/predict_results/*.gff3" -o -path "*/predict_results/*.proteins.fa" \
           -o -name "short_summary*.txt" -o -path "$cell/logs/nextflow/*trace*.txt" \
           -o -path "$cell/samples.csv" \) -print0 \
      | tar --null -T - -cf - ) | zstd -T0 -q -o "$out"
    # refuse to delete the cell if the bundle did not come out
    zstd -tq "$out" || { echo "[ERROR] bad bundle $out" >&2; exit 1; }
    echo "[tar]  $out ($(du -h "$out" | cut -f1))"
}

case "$TIER" in
A)
    for c in "${RETIRED_CELLS[@]}"; do
        [ -d "runs/$c" ] || continue
        bundle_cell "$c"
        remove "runs/$c"
    done
    for a in archive/*/; do
        remove "${a%/}"
    done
    remove runs/pilot/work
    for d in runs/v1.8.17_conda/*.stale_*; do
        remove "$d"
    done
    n_gm=$(ls runs/genemark_sidecar/output/*.genemark.gtf 2>/dev/null | wc -l)
    if [ "$n_gm" -ge 65 ]; then
        remove work
    else
        echo "[keep] work (only $n_gm GeneMark gtf files in runs/genemark_sidecar/output)"
    fi
    remove tmp_asmtest
    ;;
B)
    for c in "${ACTIVE_CELLS[@]}"; do
        if squeue -h -u "$USER" -o %j | grep -qx "fabench_$c"; then
            echo "[ERROR] fabench_$c is still running; not touching runs/$c/work" >&2
            continue
        fi
        remove "runs/$c/work"
    done
    ;;
C)
    # "<cell> <genome>" for every running/pending Nextflow task of this user
    busy=$(for j in $(squeue -h -u "$USER" -o %i); do
               wd=$(scontrol show job "$j" | grep -o 'WorkDir=[^ ]*' | cut -d= -f2)
               [ -f "$wd/.command.run" ] || continue
               cell=$(sed -n 's|.*/runs/\([^/]*\)/work/.*|\1|p' <<<"$wd")
               g=$(grep -o 'FUNANNOTATE_[A-Z]* ([^)]*)' "$wd/.command.run" | head -1 | sed 's/.*(\(.*\))/\1/' || true)
               [ -n "$cell" ] && [ -n "$g" ] && echo "$cell $g"
           done | sort -u)
    for c in "${ACTIVE_CELLS[@]}"; do
        for t in "runs/$c"/genome_annotation_training/*/training; do
            g=$(basename "$(dirname "$t")")
            if grep -qxF "$c $g" <<<"$busy"; then
                echo "[keep] $c/$g (task running or pending)"; continue
            fi
            # keep failed trainings whole (debug material)
            pasa="$t/funannotate_train.pasa.gff3"
            if [ ! -s "$pasa" ]; then
                echo "[keep] $c/$g (no PASA training models)"; continue
            fi
            # PREDICT must have finished after the current training output
            prot=$(ls "runs/$c/genome_annotation/$g/predict_results/"*.proteins.fa 2>/dev/null | head -1 || true)
            if [ -z "$prot" ] || [ ! "$prot" -nt "$pasa" ]; then
                echo "[keep] $c/$g (no PREDICT output newer than funannotate_train.pasa.gff3)"; continue
            fi
            for p in "$t/trinity_gg" "$t/getBestModel" "$t/genome.fasta.gmap" "$t/normalize" \
                     "$t"/pasa/*.transdecoder_dir "$t"/hisat2.genome.*.ht2 "$t/hisat2.coordSorted.bam" \
                     "$t"/trinity.fasta.clean* "$t/trinity.fasta.cln" "$t/trinity.fasta.cidx"; do
                [ -e "$p" ] && remove "$p"
            done
        done
    done
    ;;
*)
    echo "unknown tier: $TIER" >&2
    exit 1
    ;;
esac
[ "$APPLY" -eq 1 ] || echo "(dry run -- add --apply to delete)"
