#!/usr/bin/bash -l
#SBATCH -N 1 -n 1 -c 8 --mem 32gb --time 2:00:00
#SBATCH --partition=short
#SBATCH --job-name=fabench_strand_probe
#SBATCH --output=/bigdata/stajichlab/shared/projects/BFD/Funannotate_benchmarking/runs/%x.%j.log

# Can library strandedness be inferred from spliced reads alone, with no
# reference annotation? Measure it rather than assume it.
#
# METHOD (annotation-free, unlike RSeQC infer_experiment.py and
# check_strandedness, which both need a reference gene model / transcriptome):
# for a spliced alignment, hisat2 sets XS:A from the SPLICE MOTIF (GT-AG etc.),
# which gives the transcript's genomic strand. Compare that with the read's own
# alignment strand (FLAG 0x10) for READ 1 only (FLAG 0x40):
#     agree ~= 1.0  -> FR  (second-strand; R1 sense)
#     agree ~= 0.0  -> RF  (first-strand, dUTP; R1 antisense) -- most common
#     agree ~= 0.5  -> unstranded
#
# This probe reports the raw distribution so the separation can be judged, not
# just a verdict. It deliberately does NOT pass --rna-strandness to hisat2:
# doing so would constrain XS and make the test circular.
#
# Known limits this probe is designed to expose:
#   * needs enough SPLICED read-1 alignments; intron-poor genomes (S. cerevisiae
#     ~1.06 CDS/mRNA, Z. rouxii ~1.03) may not supply them
#   * XS:A is aligner-specific (hisat2 yes; STAR only with --outSAMstrandField
#     intronMotif; minimap2 writes ts:A read-relative instead)
set -euo pipefail

BENCH=/bigdata/stajichlab/shared/projects/BFD/Funannotate_benchmarking
SIF=/bigdata/stajichlab/shared/singularity_cache/funannotate-1.9.0-beta.12-norust.sif
CELL="$BENCH/runs/v1.9.0-beta12_container"
GENOME="$CELL/genome_annotation_training/Botrytis_cinerea_B05.10/training/genome.fasta"
R1="$CELL/rnaseq_reads/Botrytis_cinerea_norm_R1.fastq.gz"
R2="$CELL/rnaseq_reads/Botrytis_cinerea_norm_R2.fastq.gz"
NREAD=400000        # read pairs to sample
OUT="$BENCH/results/strandedness_probe"
WORK="${SCRATCH:?SCRATCH must be set by SLURM}/strand_probe"

source /etc/profile.d/modules.sh; module load singularity

# zcat | head closes the pipe early, killing zcat with SIGPIPE; under
# `set -o pipefail` that 141 aborts the script (this killed jobs 28960921 and
# 28960888 within seconds). Disable pipefail for the truncating pipeline only,
# in a subshell so the strictness stays on everywhere else. Downstream read-count
# assertions catch a genuinely truncated or corrupt file.
head_gz() { ( set +o pipefail; zcat "$1" | head -n "$2" ); }
mkdir -p "$OUT" "$WORK"; cd "$WORK"
for f in "$GENOME" "$R1" "$R2"; do [ -s "$f" ] || { echo "[ERROR] missing $f" >&2; exit 1; }; done
echo "[INFO] host=$(hostname) sampling $NREAD pairs"

head_gz "$R1" $((NREAD*4)) > "$WORK/s_R1.fq"
head_gz "$R2" $((NREAD*4)) > "$WORK/s_R2.fq"
echo "[INFO] sampled $(( $(wc -l < "$WORK/s_R1.fq") / 4 )) pairs"

singularity exec -B /bigdata:/bigdata "$SIF" bash -c "
    set -euo pipefail
    cd $WORK
    hisat2-build -p 8 $GENOME idx > hisat2-build.log 2>&1
    hisat2 -p 8 -x idx -1 s_R1.fq -2 s_R2.fq --no-unal -S aln.sam 2> hisat2.log
    tail -6 hisat2.log
    samtools view -b aln.sam | samtools sort -@ 4 -o aln.bam
    samtools index aln.bam
"

echo "[INFO] === strandedness statistic ==="
singularity exec -B /bigdata:/bigdata "$SIF" samtools view -F 0x904 "$WORK/aln.bam" \
 | awk '
   {
     cig=$6; if (cig !~ /N/) next;          # spliced alignments only
     flag=$2
     if (and(flag,0x1) && !and(flag,0x40)) next;   # read 1 only (or SE)
     xs=""
     for (i=12;i<=NF;i++) if ($i ~ /^XS:A:/) { xs=substr($i,6,1); break }
     if (xs != "+" && xs != "-") { noxs++; next }
     rstrand = and(flag,0x10) ? "-" : "+"
     tot++
     if (rstrand == xs) agree++; else disagree++
   }
   END {
     printf "[RESULT] spliced read-1 alignments with XS:A : %d\n", tot
     printf "[RESULT] no XS:A (non-canonical, excluded)   : %d\n", noxs+0
     if (tot < 1000) {
       printf "[RESULT] VERDICT: INSUFFICIENT DATA (<1000) -- must fall back to --stranded no\n"
       exit 0
     }
     f = agree/tot
     printf "[RESULT] read-strand == XS strand           : %d (%.4f)\n", agree, f
     printf "[RESULT] read-strand != XS strand           : %d (%.4f)\n", disagree, 1-f
     if (f >= 0.90)      printf "[RESULT] VERDICT: FR (second-strand)  confidence=%.4f\n", f
     else if (f <= 0.10) printf "[RESULT] VERDICT: RF (first-strand)   confidence=%.4f\n", 1-f
     else if (f >= 0.40 && f <= 0.60) printf "[RESULT] VERDICT: UNSTRANDED (no)  f=%.4f\n", f
     else                printf "[RESULT] VERDICT: AMBIGUOUS f=%.4f -- must fall back to --stranded no\n", f
   }' | tee "$OUT/botrytis_strandedness.txt"

echo "[INFO] hisat2 summary:"; tail -4 "$WORK/hisat2.log"
rm -rf "$WORK"
