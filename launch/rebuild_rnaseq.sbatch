#!/usr/bin/bash -l
#SBATCH -N 1 -n 1 -c 16 --mem 96gb --time 12:00:00
#SBATCH --partition=epyc
#SBATCH --job-name=fabench_rnaseq_rebuild
#SBATCH --output=/bigdata/stajichlab/shared/projects/BFD/Funannotate_benchmarking/runs/%x.%j.log

# Parameterized RNA-seq rebuild. Generalized from
# launch/rebuild_rnaseq_malassezia.sbatch, which was written for one genome.
#
# WHY this exists: modules/local/sra_fetch.nf truncates all three outputs to
# empty up front ("Output files must always exist (storeDir requirement)") and
# exits 0 on download failure. storeDir then caches the empty result as a
# success and never retries. Rebuilding out-of-band is the only way to replace
# a bad library. See docs/malassezia-evm-undercall.md.
#
# Mirrors SRA_FETCH's chain: parallel-fastq-dump -> bbnorm (target=30 ecc=t)
# -> fastp. Writes to a STAGING dir, NOT the shared production store.
#
# Usage:
#   sbatch --export=ALL,TAG=Genus_species,ACCS="SRRxxx SRRyyy" launch/rebuild_rnaseq.sbatch
set -euo pipefail

BENCH=/bigdata/stajichlab/shared/projects/BFD/Funannotate_benchmarking
AUX=/bigdata/stajichlab/shared/condaenv/nf_funannotate1-aux
: "${TAG:?TAG must be set, e.g. TAG=Rhodotorula_toruloides}"
: "${ACCS:?ACCS must be set, e.g. ACCS=\"SRR29721378 SRR29721379\"}"
MIN_PAIRS="${MIN_PAIRS:-100000}"
OUT="$BENCH/runs/rnaseq_rebuild/$TAG"
WORK="${SCRATCH:?SCRATCH must be set by SLURM}/rnaseq_rebuild_$TAG"

export PATH="$AUX/bin:$PATH"
mkdir -p "$OUT" "$WORK"
cd "$WORK"

echo "[INFO] host=$(hostname) scratch=$WORK out=$OUT"
echo "[INFO] tag=$TAG accessions: $ACCS"

for ACC in $ACCS; do
    echo "[INFO] === downloading $ACC ==="
    parallel-fastq-dump --sra-id "$ACC" --threads 16 --outdir "$WORK" \
        --split-files --gzip
    for m in 1 2; do
        f="$WORK/${ACC}_${m}.fastq.gz"
        [ -s "$f" ] || { echo "[ERROR] $f missing or empty after download" >&2; exit 1; }
        n=$(zcat "$f" | awk 'NR%4==1' | wc -l)
        echo "[INFO] $ACC mate$m: $n reads"
        [ "$n" -gt 0 ] || { echo "[ERROR] $ACC mate$m has zero reads" >&2; exit 1; }
    done
done

echo "[INFO] === concatenating ==="
cat $(for A in $ACCS; do echo "$WORK/${A}_1.fastq.gz"; done) > "$WORK/all_R1.fastq.gz"
cat $(for A in $ACCS; do echo "$WORK/${A}_2.fastq.gz"; done) > "$WORK/all_R2.fastq.gz"
echo "[INFO] combined R1 pairs: $(zcat "$WORK/all_R1.fastq.gz" | awk 'NR%4==1' | wc -l)"

echo "[INFO] === bbnorm (target=30 ecc=t) ==="
bbnorm.sh -Xmx80g in="$WORK/all_R1.fastq.gz" in2="$WORK/all_R2.fastq.gz" \
    out="$WORK/norm_R1.fastq.gz" out2="$WORK/norm_R2.fastq.gz" target=30 ecc=t threads=16

echo "[INFO] === fastp ==="
fastp --in1 "$WORK/norm_R1.fastq.gz" --in2 "$WORK/norm_R2.fastq.gz" \
      --out1 "$OUT/${TAG}_norm_R1.fastq.gz" --out2 "$OUT/${TAG}_norm_R2.fastq.gz" \
      --cut_right --cut_right_window_size 4 --cut_right_mean_quality 5 \
      --length_required 25 --thread 16 \
      --json "$OUT/${TAG}.fastp.json" --html "$OUT/${TAG}.fastp.html"

# SE placeholder: this rebuild is paired-only by design.
: > "$OUT/${TAG}_norm_SE.fastq.gz"

NP=$(zcat "$OUT/${TAG}_norm_R1.fastq.gz" | awk 'NR%4==1' | wc -l)
echo "[INFO] FINAL normalized pairs: $NP"
echo "[INFO] sizes:"; ls -l "$OUT"/${TAG}_norm_*.fastq.gz

# Fail loudly rather than leaving an empty result that looks like success.
if [ "$NP" -lt "$MIN_PAIRS" ]; then
    echo "[ERROR] only $NP normalized pairs -- refusing to call this a success" >&2
    exit 1
fi
echo "[INFO] rebuild complete for $TAG"
rm -rf "$WORK"
