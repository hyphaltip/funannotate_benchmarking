#!/usr/bin/bash -l
#SBATCH -N 1 -n 1 -c 2 --mem 64gb --time 02:00:00
#SBATCH --partition=short
#SBATCH --job-name=pull_funannotate_1.8.17
#SBATCH --output=/bigdata/stajichlab/shared/projects/BFD/Funannotate_benchmarking/runs/pilot/logs/pull_1.8.17.%j.log
#SBATCH --error=/bigdata/stajichlab/shared/projects/BFD/Funannotate_benchmarking/runs/pilot/logs/pull_1.8.17.%j.err

# Pull the funannotate 1.8.17 Docker Hub image (ghcr has no 1.8.17 tag) into
# the shared singularity cache as funannotate-1.8.17.sif for the
# v1.8.17_container benchmark cell.
# Submit: sbatch scripts/pull_funannotate_image.sh

set -euo pipefail

SIF_DIR=/bigdata/stajichlab/shared/lib/singularity_cache
export TMPDIR="${SCRATCH:?}"

source /etc/profile.d/modules.sh
module load apptainer

cd "$SIF_DIR"
apptainer pull docker://nextgenusfs/funannotate:v1.8.17
mv -f funannotate_v1.8.17.sif "$SIF_DIR/funannotate-1.8.17.sif"
touch "$SIF_DIR/funannotate-1.8.17.done"
echo "PULL_DONE $(date -Is)"
