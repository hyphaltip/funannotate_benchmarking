#!/usr/bin/bash -l
#SBATCH -N 1 -n 1 -c 8 --mem 16gb --time 1:00:00
#SBATCH --partition=stajichlab
#SBATCH --job-name=fabench_verify_rnaseq
#SBATCH --output=/bigdata/stajichlab/shared/projects/BFD/Funannotate_benchmarking/runs/%x.%j.log

# Map a sample of read pairs to the genome with hisat2 and report the overall
# alignment rate. Guard against installing reads from the wrong organism:
# Rhizopus microsporus' old reads mapped 26 of 49.2M.
#
# Usage:
#   sbatch --export=ALL,R1=...,R2=...,GENOME=...fa.gz,NPAIRS=200000 launch/verify_rnaseq_mapping.sbatch
set -euo pipefail
echo "[INFO] start $(date)"
: "${R1:?}" "${R2:?}" "${GENOME:?}"
NPAIRS="${NPAIRS:-200000}"
ENV=/bigdata/stajichlab/shared/condaenv/funannotate-1.9.0-rc.1
export PATH="$ENV/bin:$PATH"
W="${SCRATCH:?}/verify_rnaseq"
mkdir -p "$W"
zcat "$GENOME" > "$W/genome.fa"
hisat2-build -p 8 -q "$W/genome.fa" "$W/idx"
# head closes the pipe early; zcat's SIGPIPE (141) must not trip pipefail.
{ zcat "$R1" || true; } | head -n $((NPAIRS * 4)) > "$W/r1.fq"
{ zcat "$R2" || true; } | head -n $((NPAIRS * 4)) > "$W/r2.fq"
[ -s "$W/r1.fq" ] && [ -s "$W/r2.fq" ] || { echo "[ERROR] no reads sampled" >&2; exit 1; }
hisat2 -p 8 --max-intronlen 3000 -x "$W/idx" -1 "$W/r1.fq" -2 "$W/r2.fq" 2> "$W/summary.txt" > /dev/null
echo "[INFO] R1=$R1"
echo "[INFO] GENOME=$GENOME  sampled pairs=$NPAIRS"
cat "$W/summary.txt"
rm -rf "$W"
