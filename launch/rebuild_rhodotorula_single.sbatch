#!/usr/bin/bash -l
#SBATCH -N 1 -n 1 -c 16 --mem 128gb --time 8:00:00
#SBATCH --partition=epyc
#SBATCH --job-name=fabench_rhodo_single
#SBATCH --output=/bigdata/stajichlab/shared/projects/BFD/Funannotate_benchmarking/runs/%x.%j.log

# Rhodotorula toruloides RNA-seq from ONE accession, no concatenation.
#
# ACCESSION: SRR29721378 only (SRP518254 / PRJNA1132283, strain NP11, 24,996,630
# pairs). That is ample -- Malassezia trained cleanly on 2.2M normalized pairs.
#
# WHY SINGLE-MEMBER GZIP MATTERS HERE. Two previous jobs died inside BBTools
# paired-read streaming with off-by-one assertions (java -ea):
#   AssertionError: List size mismatch: 716 vs 717 / 5015!=5016 ...
#   RuntimeException: BBNorm terminated in an error state; the output may be corrupt.
# The reads were never unpaired -- R1 = R2 = 49,119,772 (delta 0), and repair.sh
# retained 100% with zero singletons. The remaining structural difference was
# multi-member gzip.
#
# Dropping the `cat` of two accessions is NOT sufficient on its own:
# parallel-fastq-dump splits each accession into N blocks and merges them with
# shutil.copyfileobj in binary mode (parallel-fastq-dump:78-83), so EVERY output
# it produces is already multi-member gzip, single accession or not. So we
# explicitly re-compress to a single gzip member before BBTools touches the data.
#
# Reads are taken from the rescued staging copy (first 24,996,630 pairs of
# all_R*.fastq.gz are SRR29721378, which was concatenated first) -- verified
# inline below -- so no re-download is needed.
set -euo pipefail

BENCH=/bigdata/stajichlab/shared/projects/BFD/Funannotate_benchmarking
AUX=/bigdata/stajichlab/shared/condaenv/nf_funannotate1-aux
TAG=Rhodotorula_toruloides
ACC=SRR29721378
NPAIR=24996630
OUT="$BENCH/runs/rnaseq_rebuild/$TAG"
STAGE="$OUT/raw_staging"
WORK="${SCRATCH:?SCRATCH must be set by SLURM}/rhodo_single"

export PATH="$AUX/bin:$PATH"

# zcat | head closes the pipe early, killing zcat with SIGPIPE; under
# `set -o pipefail` that 141 aborts the script (this killed jobs 28960921 and
# 28960888 within seconds). Disable pipefail for the truncating pipeline only,
# in a subshell so the strictness stays on everywhere else. Downstream read-count
# assertions catch a genuinely truncated or corrupt file.
head_gz() { ( set +o pipefail; zcat "$1" | head -n "$2" ); }
mkdir -p "$OUT" "$WORK"; cd "$WORK"
for f in "$STAGE/all_R1.fastq.gz" "$STAGE/all_R2.fastq.gz"; do
    [ -s "$f" ] || { echo "[ERROR] missing $f" >&2; exit 1; }
done
echo "[INFO] host=$(hostname) work=$WORK acc=$ACC"

echo "[INFO] === verifying the accession boundary ==="
LAST=$(head_gz "$STAGE/all_R1.fastq.gz" $((NPAIR*4-3)) | tail -n 1 | cut -d. -f1 | tr -d '@')
NEXT=$(head_gz "$STAGE/all_R1.fastq.gz" $((NPAIR*4+1)) | tail -n 1 | cut -d. -f1 | tr -d '@')
echo "[INFO] read $NPAIR is from: $LAST   next read is from: $NEXT"
[ "$LAST" = "$ACC" ] || { echo "[ERROR] read $NPAIR is $LAST, expected $ACC" >&2; exit 1; }
[ "$NEXT" != "$ACC" ] || { echo "[ERROR] read $((NPAIR+1)) is still $ACC -- boundary wrong" >&2; exit 1; }
echo "[INFO] boundary confirmed"

echo "[INFO] === extracting $ACC and re-compressing as a SINGLE gzip member ==="
for m in 1 2; do
    head_gz "$STAGE/all_R${m}.fastq.gz" $((NPAIR*4)) | pigz -p 8 > "$WORK/${ACC}_R${m}.fastq.gz"
done
# pigz writes ONE gzip member per input stream, which is the whole point of
# re-compressing here; verify read counts and member count explicitly.
for m in 1 2; do
    c=$(zcat "$WORK/${ACC}_R${m}.fastq.gz" | awk 'NR%4==1' | wc -l)
    mem=$(python3 - "$WORK/${ACC}_R${m}.fastq.gz" <<'PYEOF'
import sys, zlib
n = 0
with open(sys.argv[1], "rb") as fh:
    d = zlib.decompressobj(31)
    while True:
        chunk = fh.read(1 << 22)
        if not chunk:
            break
        while chunk:
            d.decompress(chunk)
            if not d.eof:
                break
            n += 1
            chunk = d.unused_data
            d = zlib.decompressobj(31)
print(n)
PYEOF
)
    echo "[INFO] mate$m: $c reads, $mem gzip member(s), $(stat -c %s "$WORK/${ACC}_R${m}.fastq.gz") bytes"
    [ "$c" -eq "$NPAIR" ] || { echo "[ERROR] mate$m has $c reads, expected $NPAIR" >&2; exit 1; }
    [ "$mem" -eq 1 ] || { echo "[ERROR] mate$m has $mem gzip members, expected 1" >&2; exit 1; }
done

echo "[INFO] === repair.sh: verify pairing before normalizing ==="
repair.sh -Xmx100g in="$WORK/${ACC}_R1.fastq.gz" in2="$WORK/${ACC}_R2.fastq.gz" \
    out="$WORK/fixed_R1.fastq.gz" out2="$WORK/fixed_R2.fastq.gz" \
    outs="$WORK/singletons.fastq.gz" repair overwrite=t
P1=$(zcat "$WORK/fixed_R1.fastq.gz" | awk 'NR%4==1' | wc -l)
P2=$(zcat "$WORK/fixed_R2.fastq.gz" | awk 'NR%4==1' | wc -l)
SG=$(zcat "$WORK/singletons.fastq.gz" | awk 'NR%4==1' | wc -l)
echo "[INFO] repaired: R1=$P1 R2=$P2 singletons=$SG (of $NPAIR input pairs)"
[ "$P1" -eq "$P2" ] || { echo "[ERROR] repair left R1 != R2" >&2; exit 1; }

echo "[INFO] === bbnorm (target=30, no ecc, assertions ON) ==="
bbnorm.sh -Xmx100g in="$WORK/fixed_R1.fastq.gz" in2="$WORK/fixed_R2.fastq.gz" \
    out="$WORK/norm_R1.fastq.gz" out2="$WORK/norm_R2.fastq.gz" target=30 threads=16

echo "[INFO] === fastp ==="
fastp --in1 "$WORK/norm_R1.fastq.gz" --in2 "$WORK/norm_R2.fastq.gz" \
      --out1 "$OUT/${TAG}_norm_R1.fastq.gz" --out2 "$OUT/${TAG}_norm_R2.fastq.gz" \
      --cut_right --cut_right_window_size 4 --cut_right_mean_quality 5 \
      --length_required 25 --thread 16 \
      --json "$OUT/${TAG}.fastp.json" --html "$OUT/${TAG}.fastp.html"
: > "$OUT/${TAG}_norm_SE.fastq.gz"

NP=$(zcat "$OUT/${TAG}_norm_R1.fastq.gz" | awk 'NR%4==1' | wc -l)
NQ=$(zcat "$OUT/${TAG}_norm_R2.fastq.gz" | awk 'NR%4==1' | wc -l)
echo "[INFO] FINAL normalized pairs: R1=$NP R2=$NQ"
ls -l "$OUT"/${TAG}_norm_*.fastq.gz
[ "$NP" -eq "$NQ" ] || { echo "[ERROR] final R1 != R2 ($NP vs $NQ)" >&2; exit 1; }
[ "$NP" -ge 100000 ] || { echo "[ERROR] only $NP pairs -- refusing to call this success" >&2; exit 1; }
echo "[INFO] complete for $TAG from $ACC"
rm -rf "$WORK"
