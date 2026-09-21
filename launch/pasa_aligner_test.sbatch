#!/usr/bin/bash -l
#SBATCH -N 1 -n 1 -c 8 --mem 64gb --time 24:00:00
#SBATCH --partition=epyc
#SBATCH --job-name=fabench_aligner_test
#SBATCH --output=/bigdata/stajichlab/shared/projects/BFD/Funannotate_benchmarking/runs/%x.%j.log

# Does PASA's aligner choice explain beta.12's shorter ORFs?
#
# ESTABLISHED SO FAR (Botrytis cinerea B05.10, funannotate_train.pasa.gff3):
#   mean exon/gene  1.8.17 = 3.12   beta12 = 3.12   <- IDENTICAL transcript structure
#   mean CDS/gene   1.8.17 = 2.69   beta12 = 1.96   <- ORFs cover less of it
#   CDS/exon ratio  1.8.17 = 0.863  beta12 = 0.630
# beta12 also has 15,075 three_prime_UTR features for 7,505 genes (1.8.17: 10,256),
# i.e. ORFs terminate early and the remaining exons become long multi-exon 3'UTRs.
#
# TransDecoder is NOT the cause: 5.7.1 and 6.0.0 give byte-identical output on the
# same transcripts, and util/cdna_alignment_orf_to_genome_orf.pl is byte-identical
# (md5 03d0f69760953885cb662416ca454909). TransDecoder is deterministic, so a
# shorter ORF at the same exon count means the assembled transcript SEQUENCE
# differs -- consistent with shifted splice boundaries causing a frameshift and a
# premature stop.
#
# The two arms used different PASA aligners:
#   1.8.17 -> --ALIGNERS gmap           (funannotate strips minimap2; harness passed "minimap2 gmap")
#   beta12 -> --ALIGNERS minimap2,blat  (funannotate forces minimap2 in; harness passed nothing)
# See funannotate train.py:425-430 (1.8.17) vs train.py:530-541 (beta.12, commit c4175c0).
#
# This job re-runs ONLY the PASA step on beta.12's own existing training inputs
# (same Trinity, same BAMs, same genome), varying only --ALIGNERS. Everything else
# is held constant, so any difference is the aligner.
#
#   arm gmap      : --ALIGNERS gmap            (reproduces 1.8.17's aligner set)
#   arm mm2blat   : --ALIGNERS minimap2,blat   (reproduces beta.12's, as a control)
#
# beta.12 hard-codes `filtaligners = ['minimap2']`, so gmap-only is unreachable via
# the CLI. We bind-mount a patched copy of the funannotate package to reach it.
set -euo pipefail

BENCH=/bigdata/stajichlab/shared/projects/BFD/Funannotate_benchmarking
SIF=/bigdata/stajichlab/shared/singularity_cache/funannotate-1.9.0-beta.12-norust.sif
LIVE=/rhome/jstajich/projects/funannotate/funannotate-live/funannotate
SRCCELL="$BENCH/runs/v1.9.0-beta12_container"
GENOME="$SRCCELL/input_clean_genomes/GCF_000143535.2.masked.fasta.gz"
SRCTRAIN="$SRCCELL/genome_annotation_training/Botrytis_cinerea_B05.10"
OUT="$BENCH/results/aligner_test"
WORK="${SCRATCH:?SCRATCH must be set by SLURM}/aligner_test"

source /etc/profile.d/modules.sh
module load singularity

mkdir -p "$OUT" "$WORK"

# ── patched funannotate package: filtaligners = [] so --aligners is honoured ──
PKG="$WORK/funannotate_patched"
rm -rf "$PKG"; cp -r "$LIVE" "$PKG"
sed -i "s/^\( *\)filtaligners = \['minimap2'\]/\1filtaligners = []/" "$PKG/train.py"
grep -n "filtaligners = " "$PKG/train.py" | head -3
grep -q "filtaligners = \[\]" "$PKG/train.py" || { echo "[ERROR] patch did not apply" >&2; exit 1; }
echo "[INFO] patch applied"

R1="$SRCCELL/rnaseq_reads/Botrytis_cinerea_norm_R1.fastq.gz"
R2="$SRCCELL/rnaseq_reads/Botrytis_cinerea_norm_R2.fastq.gz"
for f in "$R1" "$R2"; do [ -s "$f" ] || { echo "[ERROR] reads missing: $f" >&2; exit 1; }; done
echo "[INFO] reads: $(stat -Lc %s "$R1") / $(stat -Lc %s "$R2") bytes"

GEN="$WORK/genome.fa"
zcat "$GENOME" > "$GEN"
echo "[INFO] genome: $(grep -c '^>' "$GEN") contigs"

run_arm () {
    local arm="$1"; local aligners="$2"
    local o="$WORK/$arm"
    echo "[INFO] ================ arm $arm : --aligners $aligners ================"
    rm -rf "$o"; mkdir -p "$o/training"
    # Reuse beta12's own Trinity + alignments; drop PASA/kallisto products so PASA re-runs.
    for f in trinity.fasta trinity.fasta.clean trinity.alignments.bam trinity.alignments.gff3 \
             transcript.alignments.bam transcript.alignments.gff3; do
        [ -e "$SRCTRAIN/training/$f" ] && cp -L "$SRCTRAIN/training/$f" "$o/training/"
    done
    ls -l "$o/training/" | head
    singularity exec \
        --bind "$PKG":/pixi/.pixi/envs/base/lib/python3.8/site-packages/funannotate \
        --bind "$WORK":"$WORK" --bind "$BENCH":"$BENCH" \
        "$SIF" bash -c "
            export AUGUSTUS_CONFIG_PATH=$SRCCELL/lib/augustus/3.5/config
            export FUNANNOTATE_DB=/bigdata/stajichlab/shared/lib/funannotate_db
            export TMPDIR=$WORK/tmp_$arm; mkdir -p \$TMPDIR
            cd $o && funannotate train -i $GEN -o $o \
                --trinity $o/training/trinity.fasta \
                --left_norm $R1 --right_norm $R2 \
                --species 'Botrytis cinerea' --strain B05.10 \
                --cpus 8 --memory 48G --header_length 24 \
                --jaccard_clip --no-progress --max_intronlen 3000 \
                --aligners $aligners --pasa_db sqlite
        " > "$OUT/${arm}.train.log" 2>&1 || echo "[WARN] arm $arm exited non-zero (see log)"

    grep -oE '\-\-ALIGNERS [a-z0-9,]+' "$OUT/${arm}.train.log" | sort -u | sed "s/^/[INFO] $arm PASA got /" || true
    local G="$o/training/pasa.step1.gff3"
    if [ -s "$G" ]; then
        cp "$G" "$OUT/${arm}.pasa.step1.gff3"
        /usr/bin/python3.12 - "$G" "$arm" <<'PY'
import sys,re,collections
p,arm=sys.argv[1],sys.argv[2]
ex=collections.Counter(); cds=collections.Counter()
for ln in open(p):
    if ln.startswith("#"): continue
    f=ln.rstrip().split("\t")
    if len(f)<9: continue
    m=re.search(r'Parent=([^;,]+)',f[8])
    if not m: continue
    if f[2]=="exon": ex[m.group(1)]+=1
    elif f[2]=="CDS": cds[m.group(1)]+=1
n=len(ex) or 1
hist=collections.Counter(cds.values())
multi=sum(v for k,v in hist.items() if k>=2)
print(f"[RESULT] {arm}  models={len(ex)}  exon/gene={sum(ex.values())/n:.2f}  "
      f"CDS/gene={sum(cds.values())/max(len(cds),1):.2f}  "
      f"CDS_over_exon={sum(cds.values())/max(sum(ex.values()),1):.3f}  "
      f"multiCDS={multi} ({100*multi/max(len(cds),1):.1f}%)")
PY
    else
        echo "[RESULT] $arm  NO pasa.step1.gff3 produced"
    fi
}

run_arm gmap    gmap
run_arm mm2blat "minimap2 blat"

echo "[INFO] ================ summary ================"
echo "[INFO] reference values from the real runs (Botrytis):"
echo "[INFO]   1.8.17 (--ALIGNERS gmap)          exon/gene=3.12 CDS/gene=2.69 CDS_over_exon=0.863"
echo "[INFO]   beta12 (--ALIGNERS minimap2,blat) exon/gene=3.12 CDS/gene=1.96 CDS_over_exon=0.630"
echo "[INFO] note: these reference values are post-getBestModel; the arms above are pre-getBestModel"
echo "[INFO]       (pasa.step1.gff3), so compare arms against EACH OTHER, not against the reference."
echo "[INFO] outputs in $OUT"
