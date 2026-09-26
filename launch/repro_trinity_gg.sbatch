#!/usr/bin/bash -l
#SBATCH -N 1 -n 1 -c 16 --mem 64gb --time 2:00:00
#SBATCH --partition=stajichlab
# c01-c03 (abu_dhabi Opteron) lack BMI2; Trinity bamsifter SIGILLs (ret 132) there. Job 29069427.
#SBATCH --constraint=milan
#SBATCH --job-name=fabench_repro_trinity
#SBATCH --output=/bigdata/stajichlab/shared/projects/BFD/Funannotate_benchmarking/runs/%x.%j.log

# Re-run RNASEQ_PREPARE's `funannotate train --stop_after_trinity` for one
# species in a given conda env, and KEEP Trinity-gg.log and the funannotate
# logs (the pipeline runs this on node-local $SCRATCH, so a failure's Trinity
# log is lost when the job ends).
#
# Usage:
#   sbatch --export=ALL,ENVNAME=funannotate-1.9.0-rc.1,TAG=Batrachochytrium_dendrobatidis,\
#       OUT=Batrachochytrium_dendrobatidis_JAM81,SPECIES="Batrachochytrium dendrobatidis",STRAIN=JAM81,\
#       GENOME=runs/v1.9.0-rc1_conda/input_clean_genomes/GCF_000203795.2.masked.fasta.gz \
#       launch/repro_trinity_gg.sbatch
set -euo pipefail
: "${ENVNAME:?}" "${TAG:?}" "${OUT:?}" "${SPECIES:?}" "${STRAIN:?}" "${GENOME:?}"
BENCH=/bigdata/stajichlab/shared/projects/BFD/Funannotate_benchmarking
READS=/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/rnaseq_reads
KEEP="$BENCH/runs/trinity_repro/$ENVNAME.$OUT.${SLURM_JOB_ID}"
mkdir -p "$KEEP"
cd "$BENCH"
GENOME=$(readlink -f "$GENOME")
W="${SCRATCH:?}/repro"
mkdir -p "$W/trinity_work"
cd "$W"
zcat "$GENOME" > genome.fa
ln -s "$READS/${TAG}_norm_R1.fastq.gz" .
ln -s "$READS/${TAG}_norm_R2.fastq.gz" .
module load miniconda3 2>/dev/null || true
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "/bigdata/stajichlab/shared/condaenv/$ENVNAME"
export FUNANNOTATE_DB=/bigdata/stajichlab/shared/lib/funannotate_db
export TMPDIR="$W/trinity_work"
set +e
funannotate train -i genome.fa -o "$W/$OUT" \
    --left_norm "${TAG}_norm_R1.fastq.gz" --right_norm "${TAG}_norm_R2.fastq.gz" --aligners minimap2 \
    --species "$SPECIES" --strain "$STRAIN" --cpus 16 --memory 52G --header_length 24 \
    --jaccard_clip --no-progress --min_coverage 4 --max_intronlen 3000 \
    --stop_after_trinity --no_trimmomatic
RC=$?
set -e
echo "[INFO] funannotate train exit=$RC"
cp -p "$W/$OUT"/training/Trinity-gg.log "$W/$OUT"/logfiles/*.log "$KEEP"/ 2>/dev/null || true
ls -la "$W/$OUT/training" > "$KEEP/training.ls" 2>/dev/null || true
echo "[INFO] logs kept in $KEEP"; ls -la "$KEEP"
grep -nE "rror|died|Died|exception|Exception" "$KEEP"/Trinity-gg.log | head -20 || true
tail -30 "$KEEP"/Trinity-gg.log || true
