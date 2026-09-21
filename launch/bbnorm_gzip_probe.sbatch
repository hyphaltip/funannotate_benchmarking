#!/usr/bin/bash -l
#SBATCH -N 1 -n 1 -c 4 --mem 32gb --time 1:00:00
#SBATCH --partition=short
#SBATCH --job-name=fabench_bbnorm_probe
#SBATCH --output=/bigdata/stajichlab/shared/projects/BFD/Funannotate_benchmarking/runs/%x.%j.log

# Small controlled probe: does bbnorm's assertion failure come from MULTI-MEMBER
# gzip input, or from the read data itself?
#
# Two bbnorm jobs on BYTE-IDENTICAL read content, differing only in gzip framing:
#   MULTI  : cat a.gz b.gz > in.gz          (what rebuild_rnaseq.sbatch produced)
#   SINGLE : zcat a.gz b.gz | gzip > in.gz  (one gzip member)
#
# WHY: jobs 28956332 and 28958145 both died inside BBTools paired-read streaming
# with off-by-one assertions (java -ea):
#   AssertionError: List size mismatch: 716 vs 717 / 5015!=5016 ...
#   RuntimeException: BBNorm terminated in an error state; the output may be corrupt.
# The reads are NOT unpaired: input R1 = R2 = 49,119,772 (delta 0) and repair.sh
# retained 100% with zero singletons. So the earlier "R1/R2 out of step" reading
# of the PairStreamer trace was wrong. The remaining structural difference is that
# the input was built with `cat SRR..._1.fastq.gz SRR..._2.fastq.gz`, i.e. a
# multi-member gzip, which zcat and repair.sh handle but BBTools' parallel reader
# may not.
#
# 2M pairs, minutes of runtime. Purpose is to decide between:
#   (1) re-run with -da (suppress assertions)   -- only if data is at fault
#   (2) normalize per-accession, never concatenating -- if framing is at fault
#   (3) drop bbnorm for this genome              -- last resort, costs comparability
set -euo pipefail

BENCH=/bigdata/stajichlab/shared/projects/BFD/Funannotate_benchmarking
AUX=/bigdata/stajichlab/shared/condaenv/nf_funannotate1-aux
STAGE="$BENCH/runs/rnaseq_rebuild/Rhodotorula_toruloides/raw_staging"
WORK="${SCRATCH:?SCRATCH must be set by SLURM}/bbnorm_probe"
NPAIR=1000000   # per half; 2M pairs total per arm

export PATH="$AUX/bin:$PATH"
mkdir -p "$WORK"; cd "$WORK"

for f in "$STAGE/all_R1.fastq.gz" "$STAGE/all_R2.fastq.gz"; do
    [ -s "$f" ] || { echo "[ERROR] missing $f" >&2; exit 1; }
done

echo "[INFO] host=$(hostname) work=$WORK"
echo "[INFO] === building halves ($NPAIR pairs each, single-member gzip) ==="
for m in 1 2; do
    zcat "$STAGE/all_R${m}.fastq.gz" | head -n $((NPAIR*4))          | gzip > "$WORK/a_R${m}.gz"
    zcat "$STAGE/all_R${m}.fastq.gz" | head -n $((NPAIR*8)) | tail -n $((NPAIR*4)) | gzip > "$WORK/b_R${m}.gz"
done
ls -l "$WORK"/a_R*.gz "$WORK"/b_R*.gz | awk '{print "  "$NF"  "$5}'

echo "[INFO] === arm MULTI: cat a.gz b.gz (multi-member) ==="
for m in 1 2; do cat "$WORK/a_R${m}.gz" "$WORK/b_R${m}.gz" > "$WORK/multi_R${m}.gz"; done
echo "[INFO] === arm SINGLE: zcat a.gz b.gz | gzip (one member) ==="
for m in 1 2; do zcat "$WORK/a_R${m}.gz" "$WORK/b_R${m}.gz" | gzip > "$WORK/single_R${m}.gz"; done

# Sanity: identical read content in both arms.
for arm in multi single; do
    n=$(zcat "$WORK/${arm}_R1.gz" | awk 'NR%4==1' | wc -l)
    s=$(zcat "$WORK/${arm}_R1.gz" | md5sum | cut -d' ' -f1)
    echo "[INFO] $arm: R1 reads=$n  content_md5=$s"
done

run_arm () {
    local arm="$1"
    echo "[INFO] ================ bbnorm on $arm ================"
    if bbnorm.sh -Xmx24g in="$WORK/${arm}_R1.gz" in2="$WORK/${arm}_R2.gz" \
            out="$WORK/${arm}_norm_R1.gz" out2="$WORK/${arm}_norm_R2.gz" \
            target=30 threads=4 > "$WORK/${arm}.bbnorm.log" 2>&1; then
        local n=$(zcat "$WORK/${arm}_norm_R1.gz" 2>/dev/null | awk 'NR%4==1' | wc -l)
        echo "[RESULT] $arm  SUCCESS  normalized_pairs=$n"
    else
        echo "[RESULT] $arm  FAILED (exit $?)"
        grep -mE "AssertionError|RuntimeException" -m 3 "$WORK/${arm}.bbnorm.log" 2>/dev/null | head -3 || true
        grep -E "AssertionError|RuntimeException" "$WORK/${arm}.bbnorm.log" | head -3 || true
    fi
}

run_arm multi
run_arm single

echo "[INFO] ================ verdict ================"
echo "[INFO] MULTI fails + SINGLE succeeds -> gzip framing is the cause (option 2)"
echo "[INFO] both fail                      -> the read data is at fault (option 1 or 3)"
echo "[INFO] both succeed                   -> failure is scale-dependent, not framing"
rm -rf "$WORK"
