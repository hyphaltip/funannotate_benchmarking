#!/usr/bin/bash -l
#SBATCH -N 1 -n 1 -c 2 --mem 8gb --time 7-00:00:00
#SBATCH --partition=stajichlab
#SBATCH --job-name=fabench_cell
#SBATCH --output=/bigdata/stajichlab/shared/projects/BFD/Funannotate_benchmarking/runs/%x.%j.log

# Launch ONE benchmark cell (see conf/cells.tsv / DESIGN.md "6 cells per genome")
# of the vendored nf_funannotate1 pipeline against this project's per-cell
# samplesheet under runs/<cell>/.
#
# Submit:
#   sbatch launch/run_cell.sbatch <cell>
#   CELL=v1.8.17_conda sbatch launch/run_cell.sbatch
#
# <cell> must match a row in conf/cells.tsv, e.g.:
#   v1.8.17_conda  v1.8.17_container  v1.9.0-beta10_conda  v1.9.0-beta10_conda_rust
#   v1.9.0-beta10_container  v1.9.0-beta10_container_rust
#
# Extra args after <cell> are passed straight through to nextflow, e.g.:
#   sbatch launch/run_cell.sbatch v1.8.17_conda --n_test 2

set -euo pipefail

BENCH_ROOT="/bigdata/stajichlab/shared/projects/BFD/Funannotate_benchmarking"
cd "$BENCH_ROOT"

CELL="${1:-${CELL:-}}"
[ -n "$CELL" ] || { echo "[ERROR] usage: sbatch launch/run_cell.sbatch <cell>  (see conf/cells.tsv)" >&2; exit 1; }
shift || true

CELLS_TSV="conf/cells.tsv"
ROW=$(awk -F'\t' -v c="$CELL" 'NR>1 && $1==c {print; exit}' "$CELLS_TSV")
[ -n "$ROW" ] || { echo "[ERROR] cell '$CELL' not found in $CELLS_TSV" >&2; exit 1; }

VERSION=$(echo "$ROW" | cut -f2)
PROVISIONING=$(echo "$ROW" | cut -f3)
EVM_BACKEND=$(echo "$ROW" | cut -f4)
CONDA_ENV=$(echo "$ROW" | cut -f5)
CONTAINER_FUNANNOTATE=$(echo "$ROW" | cut -f6)

CELL_DIR="runs/${CELL}"
[ -d "$CELL_DIR" ] || { echo "[ERROR] $CELL_DIR does not exist -- run scripts/build_run_samples.py first" >&2; exit 1; }

# conf/benchmark.yaml pipeline: block.
NF_FUNANNOTATE1="${NF_FUNANNOTATE1:-/rhome/jstajich/projects/nf/nf_funannotate1}"
export CONDA_ENVS_ROOT="${CONDA_ENVS_ROOT:-/bigdata/stajichlab/shared/condaenv}"

# GeneMark is held constant across all 6 cells via the shared sidecar
# (launch/run_genemark_sidecar.sbatch) rather than retrained independently per
# cell -- see DESIGN.md "GeneMark sidecar". Must be an ABSOLUTE path: this
# script cd's into $CELL_DIR (runs/<cell>/) below before invoking nextflow.
GENEMARK_SIDECAR_DIR="${GENEMARK_SIDECAR_DIR:-${BENCH_ROOT}/runs/genemark_sidecar/output}"
if [ ! -d "$GENEMARK_SIDECAR_DIR" ] || [ -z "$(ls -A "$GENEMARK_SIDECAR_DIR" 2>/dev/null)" ]; then
    echo "[WARN] $GENEMARK_SIDECAR_DIR is empty/missing -- run launch/run_genemark_sidecar.sbatch first, or every genome in this cell will fall back to --auto-skip-genemark" >&2
fi

module load nextflow

EXECUTOR="${EXECUTOR:-slurm}"
# ucr_hpcc always stays in the profile list ahead of the provisioning axis --
# it carries the SLURM SCRATCH/TMPDIR safety net (see nf_funannotate1
# run_annotate.sh), independent of conda/container choice.
#
# conf/cells.tsv's provisioning column uses "container", but nf_funannotate1's
# actual Nextflow profile for that axis is named "singularity" (see its
# nextflow.config profiles{} block) -- there is no "container" profile, so
# passing the column value straight through silently drops apptainer.enabled
# and conf/provision_singularity.config (Nextflow only warns on an unknown
# profile name, it doesn't fail).
case "$PROVISIONING" in
    conda)     PROFILE_AXIS="conda" ;;
    container) PROFILE_AXIS="singularity" ;;
    *)         echo "[ERROR] unknown provisioning '$PROVISIONING' for cell '$CELL'" >&2; exit 1 ;;
esac
PROFILES="annotate,${EXECUTOR},ucr_hpcc,${PROFILE_AXIS}"

EXTRA_ARGS=(--genemark_sidecar_dir "$GENEMARK_SIDECAR_DIR")
case "$PROVISIONING" in
    conda)
        [ -n "$CONDA_ENV" ] || { echo "[ERROR] cell '$CELL' has no conda_env set in $CELLS_TSV" >&2; exit 1; }
        EXTRA_ARGS+=(--conda_env "$CONDA_ENV")
        ;;
    container)
        [ -n "$CONTAINER_FUNANNOTATE" ] || { echo "[ERROR] cell '$CELL' has no container_funannotate set in $CELLS_TSV" >&2; exit 1; }
        EXTRA_ARGS+=(--container_funannotate "$CONTAINER_FUNANNOTATE")
        ;;
esac

echo "[INFO] cell=${CELL} version=${VERSION} provisioning=${PROVISIONING} evm_backend=${EVM_BACKEND}"
echo "[INFO] -profile ${PROFILES} ${EXTRA_ARGS[*]}"

cd "$CELL_DIR"
mkdir -p logs/nextflow

# runs/<cell>/ is itself the nf_funannotate1 launch dir (samples.csv, lib/,
# tables/ were populated there by scripts/build_run_samples.py), matching
# nf_funannotate1's own run_annotate.sh convention.
NXF_OPTS="-Xms512m -Xmx4g" \
nextflow run "${NF_FUNANNOTATE1}" \
    -profile "${PROFILES}" \
    "${EXTRA_ARGS[@]}" \
    -resume \
    "$@"
