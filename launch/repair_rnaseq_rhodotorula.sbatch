#!/usr/bin/bash -l
#SBATCH -N 1 -n 1 -c 16 --mem 250gb --time 12:00:00
#SBATCH --partition=epyc
#SBATCH --job-name=fabench_rnaseq_repair
#SBATCH --output=/bigdata/stajichlab/shared/projects/BFD/Funannotate_benchmarking/runs/%x.%j.log

# Repair R1/R2 pairing in the already-downloaded Rhodotorula toruloides reads,
# then normalize what is still paired.
#
# WHY: job 28956332 downloaded SRR29721378 + SRR29721379 cleanly (24,996,630 and
# 24,123,142 pairs; mate1 == mate2 within each accession) but bbnorm died with
#   java.lang.AssertionError: List size mismatch: 6 vs 5
#     at stream.PairStreamer.nextList(PairStreamer.java:88)
#     at jgi.BBMerge.makeInsertHistogram(BBMerge.java:1339)
# PairStreamer is BBTools' pairing layer, so R1 and R2 are out of step somewhere
# past the head of the file. parallel-fastq-dump splits each accession into 15
# blocks and reassembles them; that reassembly is the known desync source (the
# same one modules/local/sra_fetch.nf's rename_headers override works around).
#
# Rather than re-downloading, this repairs the concatenated pair that was
# rescued from job 28956332's node-local scratch and keeps whatever survives
# as a matched pair. Singletons are written out and counted, not discarded
# silently, so the loss is on the record.
#
# ecc=t is dropped: error correction is what pulled BBMerge (and therefore
# PairStreamer) into the run, and it is not needed for training-read
# normalization.
set -euo pipefail

BENCH=/bigdata/stajichlab/shared/projects/BFD/Funannotate_benchmarking
AUX=/bigdata/stajichlab/shared/condaenv/nf_funannotate1-aux
TAG=Rhodotorula_toruloides
OUT="$BENCH/runs/rnaseq_rebuild/$TAG"
STAGE="$OUT/raw_staging"
WORK="${SCRATCH:?SCRATCH must be set by SLURM}/rnaseq_repair_$TAG"

export PATH="$AUX/bin:$PATH"
mkdir -p "$OUT" "$WORK"

for f in "$STAGE/all_R1.fastq.gz" "$STAGE/all_R2.fastq.gz"; do
    [ -s "$f" ] || { echo "[ERROR] $f missing or empty" >&2; exit 1; }
done
echo "[INFO] host=$(hostname) work=$WORK stage=$STAGE"
ls -l "$STAGE"

# Copy to node-local NVMe first; repair.sh is I/O heavy and /bigdata is shared.
echo "[INFO] === staging to node-local scratch ==="
cp "$STAGE/all_R1.fastq.gz" "$STAGE/all_R2.fastq.gz" "$WORK/"

echo "[INFO] === counting input reads (both mates) ==="
N1=$(zcat "$WORK/all_R1.fastq.gz" | awk 'NR%4==1' | wc -l)
N2=$(zcat "$WORK/all_R2.fastq.gz" | awk 'NR%4==1' | wc -l)
echo "[INFO] input R1=$N1  R2=$N2  delta=$((N1-N2))"

echo "[INFO] === repair.sh (restore pairing order) ==="
repair.sh -Xmx200g in="$WORK/all_R1.fastq.gz" in2="$WORK/all_R2.fastq.gz" \
    out="$WORK/fixed_R1.fastq.gz" out2="$WORK/fixed_R2.fastq.gz" \
    outs="$WORK/singletons.fastq.gz" repair overwrite=t

P1=$(zcat "$WORK/fixed_R1.fastq.gz" | awk 'NR%4==1' | wc -l)
P2=$(zcat "$WORK/fixed_R2.fastq.gz" | awk 'NR%4==1' | wc -l)
SG=$(zcat "$WORK/singletons.fastq.gz" | awk 'NR%4==1' | wc -l)
echo "[INFO] repaired pairs R1=$P1 R2=$P2  singletons=$SG"
[ "$P1" -eq "$P2" ] || { echo "[ERROR] repair left R1 != R2 ($P1 vs $P2)" >&2; exit 1; }
echo "[INFO] retained $P1 of $N1 input pairs ($(awk -v a=$P1 -v b=$N1 'BEGIN{printf "%.2f", 100*a/b}')%)"

echo "[INFO] === bbnorm (target=30, no ecc) ==="
bbnorm.sh -Xmx200g in="$WORK/fixed_R1.fastq.gz" in2="$WORK/fixed_R2.fastq.gz" \
    out="$WORK/norm_R1.fastq.gz" out2="$WORK/norm_R2.fastq.gz" target=30 threads=16

echo "[INFO] === fastp ==="
fastp --in1 "$WORK/norm_R1.fastq.gz" --in2 "$WORK/norm_R2.fastq.gz" \
      --out1 "$OUT/${TAG}_norm_R1.fastq.gz" --out2 "$OUT/${TAG}_norm_R2.fastq.gz" \
      --cut_right --cut_right_window_size 4 --cut_right_mean_quality 5 \
      --length_required 25 --thread 16 \
      --json "$OUT/${TAG}.fastp.json" --html "$OUT/${TAG}.fastp.html"

# SE placeholder: this rebuild is paired-only by design.
: > "$OUT/${TAG}_norm_SE.fastq.gz"

NP=$(zcat "$OUT/${TAG}_norm_R1.fastq.gz" | awk 'NR%4==1' | wc -l)
NQ=$(zcat "$OUT/${TAG}_norm_R2.fastq.gz" | awk 'NR%4==1' | wc -l)
echo "[INFO] FINAL normalized pairs: R1=$NP R2=$NQ"
ls -l "$OUT"/${TAG}_norm_*.fastq.gz

[ "$NP" -eq "$NQ" ] || { echo "[ERROR] final R1 != R2 ($NP vs $NQ)" >&2; exit 1; }
if [ "$NP" -lt 100000 ]; then
    echo "[ERROR] only $NP normalized pairs -- refusing to call this a success" >&2
    exit 1
fi
echo "[INFO] repair+normalize complete for $TAG"
rm -rf "$WORK"
