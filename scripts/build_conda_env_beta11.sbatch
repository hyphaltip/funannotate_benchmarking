#!/usr/bin/bash -l
#SBATCH -N 1 -n 1 -c 4 --mem 16gb --time 02:00:00
#SBATCH --partition=short
#SBATCH --job-name=build_funannotate_beta11
#SBATCH --output=/bigdata/stajichlab/shared/projects/BFD/Funannotate_benchmarking/runs/pilot/logs/build_beta11.%j.log
#SBATCH --error=/bigdata/stajichlab/shared/projects/BFD/Funannotate_benchmarking/runs/pilot/logs/build_beta11.%j.err

# Build funannotate-1.9.0-beta.11 (+ -rust) conda envs to sit alongside the
# beta.10 pair, mirroring what's already been done on the container side
# (funannotate-1.9.0-beta.11.sif was pulled 2026-09-06 from
# ghcr.io/nextgenusfs/funannotate:1.9.0-beta.11 -- same recipe: pip install
# git+https://github.com/nextgenusfs/funannotate.git@v1.9.0-beta.11, commit
# 17c58f5da60da7b3778518fc3ff5902d50e5168e, confirmed via the .sif's
# funannotate-1.9.0.dist-info/direct_url.json).
#
# Approach: clone the existing beta.10 env(s) verbatim (same conda/bioconda
# deps, same augustus/exonerate/etc build) via `conda create --clone`, then
# force-reinstall ONLY the funannotate Python package at the beta.11 git tag
# inside the clone. This is the same "swap the pip package, keep everything
# else" recipe already used to go from 1.8.17 -> 1.9.0-beta.10 (see
# conda-meta/history in each env -- the bioconda/conda-forge package set is
# otherwise identical across funannotate point releases).
#
# Submit: sbatch scripts/build_conda_env_beta11.sbatch

set -euo pipefail

CONDA_ENVS_ROOT="${CONDA_ENVS_ROOT:-/bigdata/stajichlab/shared/condaenv}"
TAG="v1.9.0-beta.11"

source /etc/profile.d/modules.sh
module load miniconda3

build_one() {
    local src="$1" dest="$2"
    if [ -d "$dest" ]; then
        echo "[SKIP] $dest already exists"
        return
    fi
    echo "[INFO] cloning $src -> $dest"
    conda create -y -p "$dest" --clone "$src"
    echo "[INFO] reinstalling funannotate @ ${TAG} into $dest"
    "$dest/bin/pip" install --no-deps --force-reinstall \
        "git+https://github.com/nextgenusfs/funannotate.git@${TAG}"
    "$dest/bin/funannotate" check || true
    echo "[DONE] $dest"
}

build_one "${CONDA_ENVS_ROOT}/funannotate-1.9.0-beta.10"       "${CONDA_ENVS_ROOT}/funannotate-1.9.0-beta.11"
build_one "${CONDA_ENVS_ROOT}/funannotate-1.9.0-beta.10-rust"  "${CONDA_ENVS_ROOT}/funannotate-1.9.0-beta.11-rust"

echo "BUILD_DONE $(date -Is)"
