#!/usr/bin/bash -l
#SBATCH -N 1 -n 1 -c 16 --mem 64gb --time 8:00:00
#SBATCH --partition=epyc
#SBATCH --job-name=fabench_td_test
#SBATCH --output=/bigdata/stajichlab/shared/projects/BFD/Funannotate_benchmarking/runs/%x.%j.log

# Controlled TransDecoder 5.7.1 vs 6.0.0 test on ONE identical transcript set,
# plus two leniency arms for 6.0.0.
#
# WHY: beta.12 training models lose ~21 points of exon-structure retention
# versus 1.8.17 (73.0% -> 51.7% mean across six genomes), and stage-splitting
# put the loss at the TransDecoder step, not PASA. The funannotate call
# (train.py) and PASA's parameter construction (pasa_asmbls_to_training_set.dbi)
# are byte-for-byte identical between the two versions, and both TransDecoder
# CLIs expose the same options with the same defaults -- so the change is
# internal to TransDecoder 6.0.0. See docs/beta12-exon-collapse-and-input-integrity.md.
#
# Input is Botrytis cinerea Trinity transcripts, subsampled. Not PASA assemblies
# (funannotate deletes $TRAINDIR/pasa after training) but the same class of
# input: assembled transcripts in, ORF choice out.
#
# Arms:
#   A  TD 5.7.1  defaults                                   <- the 1.8.17 behaviour
#   B  TD 6.0.0  defaults                                   <- the beta.12 behaviour
#   C  TD 6.0.0  --retain_long_orfs_mode strict --retain_long_orfs_length 300
#   D  TD 6.0.0  -m 50   (min protein length, default 100)
#
# C and D are the leniency question: C bypasses the dynamic 1%FDR threshold,
# D lowers the length floor. NOTE both are currently UNREACHABLE through
# funannotate -- PASA builds its TransDecoder command from only -t/-G/-v/-T/-m/-S,
# and funannotate passes none of them. Reaching C needs a PASA patch; reaching
# D needs a one-line funannotate change.
set -euo pipefail

BENCH=/bigdata/stajichlab/shared/projects/BFD/Funannotate_benchmarking
SIF=/bigdata/stajichlab/shared/singularity_cache/funannotate-1.9.0-beta.12-norust.sif
TD5=/bigdata/stajichlab/shared/condaenv/funannotate-1.8.17/opt/transdecoder
TD6=/pixi/.pixi/envs/base/opt/transdecoder
SRC="$BENCH/runs/v1.9.0-beta12_container/genome_annotation_training/Botrytis_cinerea_B05.10/training/trinity.fasta.clean"
OUT="$BENCH/results/transdecoder_version_test"
WORK="${SCRATCH:?SCRATCH must be set by SLURM}/td_test"
NSEQ=8000

source /etc/profile.d/modules.sh
module load singularity

[ -s "$SRC" ] || { echo "[ERROR] $SRC missing" >&2; exit 1; }
mkdir -p "$OUT" "$WORK"
cd "$WORK"

echo "[INFO] host=$(hostname) work=$WORK"
echo "[INFO] === subsampling $NSEQ transcripts from $SRC ==="
awk -v n=$NSEQ '/^>/{c++} c>n{exit} {print}' "$SRC" > "$WORK/input.fasta"
echo "[INFO] input transcripts: $(grep -c '^>' "$WORK/input.fasta")"

run_arm () {
    local arm="$1"; shift
    local bin_longorfs="$1"; shift
    local bin_predict="$1"; shift
    local use_sif="$1"; shift
    local d="$WORK/$arm"
    mkdir -p "$d"; cp "$WORK/input.fasta" "$d/"
    echo "[INFO] === arm $arm : extra opts: $* ==="
    if [ "$use_sif" = "sif" ]; then
        singularity exec "$SIF" bash -c "cd $d && perl $bin_longorfs -t input.fasta -O $d/td >/dev/null 2>&1 && perl $bin_predict -t input.fasta -O $d/td $* >/dev/null 2>&1" || {
            echo "[WARN] arm $arm failed"; return 0; }
    else
        ( cd "$d" && perl "$bin_longorfs" -t input.fasta -O "$d/td" >/dev/null 2>&1 && \
                     perl "$bin_predict"  -t input.fasta -O "$d/td" "$@" >/dev/null 2>&1 ) || {
            echo "[WARN] arm $arm failed"; return 0; }
    fi
    local gff="$d/input.fasta.transdecoder.gff3"
    [ -s "$gff" ] || gff=$(find "$d" -name "*.transdecoder.gff3" | head -1)
    if [ -s "$gff" ]; then
        cp "$gff" "$OUT/${arm}.transdecoder.gff3"
        local cds=$(awk -F'\t' '$3=="CDS"' "$gff" | wc -l)
        local mrna=$(awk -F'\t' '$3=="mRNA"' "$gff" | wc -l)
        local len=$(awk -F'\t' '$3=="CDS"{s+=$5-$4+1} END{print s}' "$gff")
        echo "[RESULT] $arm  mRNA=$mrna  CDS=$cds  total_CDS_bp=$len"
        [ "$mrna" -gt 0 ] && echo "[RESULT] $arm  mean_CDS_per_mRNA=$(awk -v a=$cds -v b=$mrna 'BEGIN{printf "%.3f",a/b}')  mean_CDS_bp=$(awk -v a=$len -v b=$mrna 'BEGIN{printf "%.1f",a/b}')"
    else
        echo "[WARN] arm $arm produced no gff3"
    fi
}

run_arm A_td571_default   "$TD5/TransDecoder.LongOrfs"     "$TD5/TransDecoder.Predict"     local
run_arm B_td600_default   "$TD6/util/TransDecoder.LongOrfs" "$TD6/util/TransDecoder.Predict" sif
run_arm C_td600_strict300 "$TD6/util/TransDecoder.LongOrfs" "$TD6/util/TransDecoder.Predict" sif --retain_long_orfs_mode strict --retain_long_orfs_length 300
run_arm D_td600_m50       "$TD6/util/TransDecoder.LongOrfs" "$TD6/util/TransDecoder.Predict" sif -m 50

echo "[INFO] === summary ==="
grep "^\[RESULT\]" "$BENCH/runs/fabench_td_test.${SLURM_JOB_ID}.log" 2>/dev/null || true
echo "[INFO] outputs in $OUT"
rm -rf "$WORK"
