#!/usr/bin/bash -l
#SBATCH -N 1 -n 1 -c 2 --mem 8gb --time 7-00:00:00
#SBATCH --partition=stajichlab
#SBATCH --job-name=fabench_genemark_sidecar
#SBATCH --output=/bigdata/stajichlab/shared/projects/BFD/Funannotate_benchmarking/runs/%x.%j.log

# Run GeneMark exactly ONCE per genome (nf_funannotate1's genemark_sidecar.nf,
# always container-mode / braker3, always fresh ES self-training -- no RNA-seq
# hints, so the result is identical regardless of which cell's own training
# produced what), shared by all 6 comparison cells via --genemark_sidecar_dir.
# See DESIGN.md "GeneMark sidecar" and README.md "The 6 cells" for why: GeneMark
# itself is not part of what's being compared, and container cells were
# previously getting a DIFFERENT GeneMark (braker3-bundled) than conda cells
# (host-licensed genemarkESET) -- a real confound for the conda-vs-container
# question, on top of 6x redundant compute.
#
# Submit ONCE, before (or concurrently with -- it only needs to finish before
# a given genome reaches FUNANNOTATE_PREDICT in each cell) launching the cells:
#   sbatch launch/run_genemark_sidecar.sbatch
#   sbatch launch/run_genemark_sidecar.sbatch --n_test 5     # extra args pass through
#
# Any cell's samples.csv works as the --samples input (all cells share the
# same genome set / masked-genome paths); v1.8.17_conda's is used below simply
# because it's always populated first.

set -euo pipefail

BENCH_ROOT="/bigdata/stajichlab/shared/projects/BFD/Funannotate_benchmarking"
cd "$BENCH_ROOT"

NF_FUNANNOTATE1="${NF_FUNANNOTATE1:-/rhome/jstajich/projects/nf/nf_funannotate1}"
SAMPLES="${SAMPLES:-runs/v1.8.17_conda/samples.csv}"
TARGET="${TARGET:-runs/genemark_sidecar/output}"

[ -f "$SAMPLES" ] || { echo "[ERROR] $SAMPLES not found -- run scripts/build_run_samples.py first" >&2; exit 1; }

module load nextflow

mkdir -p logs/nextflow "$TARGET"

NXF_OPTS="-Xms512m -Xmx4g" \
nextflow run "${NF_FUNANNOTATE1}/genemark_sidecar.nf" \
    -c "${NF_FUNANNOTATE1}/nextflow.config" \
    -profile genemark_sidecar,slurm,singularity \
    --samples "$SAMPLES" \
    --target "$TARGET" \
    -work-dir work/genemark_sidecar \
    -resume \
    "$@"
