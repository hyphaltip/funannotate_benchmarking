#!/usr/bin/bash -l
#SBATCH -N 1 -n 1 -c 16 --mem 96gb --time 12:00:00
#SBATCH --partition=epyc
#SBATCH --job-name=fabench_rnaseq_rebuild
#SBATCH --output=/bigdata/stajichlab/shared/projects/BFD/Funannotate_benchmarking/runs/%x.%j.log

# Rebuild Malassezia globosa RNA-seq from SRP487370 (PRJNA1071990).
#
# WHY: the previous Malassezia_globosa_norm_* files in the shared store were
# 0-byte R1/R2 + a 13.8 MB SE file, giving only 185k mapped reads and a 4x gene
# undercall. Root cause: modules/local/sra_fetch.nf truncates all three outputs
# to empty up front ("Output files must always exist (storeDir requirement)")
# and exits 0 on download failure -- storeDir then caches the empty result as a
# success, so it never retries. See docs/malassezia-evm-undercall.md.
#
# Mirrors SRA_FETCH's chain: parallel-fastq-dump -> bbnorm (target=30 ecc=t)
# -> fastp. Writes to a STAGING dir, NOT the shared production store.
set -euo pipefail

BENCH=/bigdata/stajichlab/shared/projects/BFD/Funannotate_benchmarking
AUX=/bigdata/stajichlab/shared/condaenv/nf_funannotate1-aux
TAG=Malassezia_globosa
ACCS="SRR27830013 SRR27830002"
OUT="$BENCH/runs/rnaseq_rebuild/$TAG"
WORK="${SCRATCH:?SCRATCH must be set by SLURM}/rnaseq_rebuild_$TAG"

export PATH="$AUX/bin:$PATH"
mkdir -p "$OUT" "$WORK"
cd "$WORK"

echo "[INFO] host=$(hostname) scratch=$WORK out=$OUT"
echo "[INFO] accessions: $ACCS"

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
if [ "$NP" -lt 100000 ]; then
    echo "[ERROR] only $NP normalized pairs -- refusing to call this a success" >&2
    exit 1
fi
echo "[INFO] rebuild complete for $TAG"
rm -rf "$WORK"
