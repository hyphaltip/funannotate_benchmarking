# Funannotate Benchmarking

Reproducible benchmark comparing **funannotate 1.8.17 vs 1.9.0-beta.10** (with and
without the Rust-optimized EVM/PASA/Trinity build) across **conda vs container**
execution, on a structured random sample of 65 fungal genomes.

Authoritative design: **[DESIGN.md](DESIGN.md)** (supersedes [PLAN.md](PLAN.md)).

## The three questions

1. Does the 1.9.0 Rust reimplementation of Trinity/EVM/PASA improve performance
   over the perl stack (measured inside 1.9.0, rust on vs off)?
2. What do the 1.8.17 → 1.9.0 code changes do to runtime, resource use, and
   gene-content predictions (1.8.17 vs 1.9.0 rust-off)?
3. What does conda vs container execution cost or change (same version)?

## Repository layout

```
Funannotate_benchmarking/
  DESIGN.md / PLAN.md        design (authoritative) + original request (historical)
  samples.csv                append-only dataset-selection record (generated)
  conf/
    benchmark.yaml           dataset-selection tunables (N, quotas, weights)
    cells.tsv                the 6 run cells → nf_funannotate1 invocations
    subphylum_map.tsv        fungal class → subphylum (quota bookkeeping)
    busco_lineage_map.tsv    taxon → BUSCO lineage used for predict/compare
  launch/
    run_cell.sbatch           launch/relaunch ONE cell (reads conf/cells.tsv)
    run_all_cells.sbatch      submit run_cell.sbatch for every cell in conf/cells.tsv
    validate_harness.sbatch   sbatch wrapper for scripts/validate_harness.py
  scripts/
    select_genomes.py        append-only sampler (taxonomic quotas, ≤2/genus, ≤1/species)
    fetch_genomes.py         masked-genome + RefSeq-GFF provisioning (symlink or download)
    build_run_samples.py     benchmark samples.csv → per-cell nf_funannotate1 samples.csv
    validate_harness.py      gate: nf_funannotate1 vs known-good Fungi_BFD outputs
    collect_metrics.py       nextflow trace/report → results/metrics.tsv
    compare_predictions.py   BUSCO + gffcompare + bedtools vs RefSeq GFF
    summarize_report.py      results → results/report.md + figures/tables (not yet written)
  input_clean_genomes/       (generated) symlinks to masked genomes, keyed by accession
  runs/
    v1.8.17_conda/  v1.8.17_container/
    v1.9.0-beta10_conda/  v1.9.0-beta10_conda_rust/
    v1.9.0-beta10_container/  v1.9.0-beta10_container_rust/
    genemark_sidecar/        (unused — sidecar sharing not implemented; see "The 6 cells")
  results/                   (generated) metrics.tsv, gene_content_comparison.tsv, report.*
```

`nf_funannotate1` itself is **not vendored into this repo** — every cell runs
directly against the dev checkout path in `conf/benchmark.yaml`
(`pipeline.nf_funannotate1`); see "Running the pipeline" below.

## Workflow

```bash
# 1. Dataset selection (append-only) — runs against the local Fungi_BFD pool + NCBI
python3 scripts/select_genomes.py --out samples.csv          # creates samples.csv
python3 scripts/select_genomes.py                            # re-run → only tops up missing quotas

# 2. Provision inputs: symlink already-masked genomes from Fungi_BFD_runs,
#    or download+mask anything missing; also fetch RefSeq GFFs for comparison
python3 scripts/fetch_genomes.py --samples samples.csv

# 3. Build per-cell nf_funannotate1 samplesheets (GENOME → masked path)
python3 scripts/build_run_samples.py --samples samples.csv

# 4. Launch the 6 cells (per genome, fixed SLURM resource class; image/env
#    per cell comes from conf/cells.tsv — see "The 6 cells" below)
sbatch launch/run_all_cells.sbatch            # all 6, one sbatch job per cell
sbatch launch/run_cell.sbatch v1.8.17_conda   # relaunch/retry a single cell (-resume)

# 5. Validation gate — once a cell has predict output for a few genomes,
#    diff it against Fungi_BFD's known-good results for the same genome tag
sbatch launch/validate_harness.sbatch v1.8.17_conda

# 6. Collect + compare + report (can run incrementally as cells finish)
python3 scripts/collect_metrics.py     --cells-dir runs
python3 scripts/compare_predictions.py --cells-dir runs
python3 scripts/summarize_report.py    --cells-dir runs   # not yet written
```

Steps 1-3 have already run once (dataset selected; `runs/<cell>/` samplesheets
built). All 6 cells in step 4 have already been launched at least once
(`runs/<cell>/launch_*.log`). **Current status (2026-09-03):** `v1.8.17_conda`
is progressing cleanly; the other 5 have each hit one of two recurring
failures — see "Known issues" below. Re-launching a cell with
`launch/run_cell.sbatch <cell>` is safe (`-resume`) once the underlying fix
for its failure lands.

### Known issues (as of 2026-09-03)

- **`SRA_FETCH` chunked-download race** (hit `v1.8.17_container`,
  `v1.9.0-beta10_conda`, `v1.9.0-beta10_conda_rust`, a different genome each
  time): the fastq-dump/header-fixup pipeline reports
  `cannot open '<species>_R1.fastq.gz': No such file or directory` after an
  otherwise-normal chunked download. Root cause not yet confirmed — needs a
  dedicated repro/debug pass. **Not fixed.**
- **Container disk-full** (`v1.9.0-beta10_container`,
  `v1.9.0-beta10_container_rust`, both `[FAILED] failed=41`): apptainer fell
  back to full sandbox extraction of the `.sif` into node-local `/tmp`
  (`squashfuse not found`) and concurrent tasks exhausted disk. Likely fixed
  by nf_funannotate1's singularity-axis `process.shell = ['/bin/bash']` change
  (a login shell was re-sourcing `/etc/profile` inside `.command.run` and
  wiping the squashfuse-bearing PATH before the container command ran) — **not
  yet re-verified with a relaunch.**
- **GeneMark sidecar not implemented** (see "The 6 cells" below): each cell
  currently runs its own `GENEMARK_RUN` independently rather than sharing one
  cached-per-genome result.

## Design constraints honored by `select_genomes.py`

- Append-only: existing `samples.csv` rows are never removed; only unsatisfied
  quotas get new candidates on re-runs.
- N = 65; floors: ≥10 Ascomycota (≥5 Pezizomycotina), ≥5 Basidiomycota
  (≥2 Agaricomycotina, ≥1 Pucciniomycotina, ≥1 Ustilaginomycotina), ≥3
  Mucoromycota, ≥1 Chytridiomycota, ≥1 Blastocladiomycota.
- Constraints: genome size ≤75 Mb, ≤2 genomes/genus, ≤1 genome/species,
  RefSeq (`GCF_`) preferred; remainder filled proportionally with a diversity
  tilt (class-level sqrt weighting).
- Metadata recorded at selection/fetch time: accession, species, genus, phylum,
  subphylum, genome_size_mb, rnaseq_available (+SRA ids), refseq_annotation_release,
  assembly_level, busco_lineage_used.

## The 6 cells

See `conf/cells.tsv` (source of truth; `launch/run_cell.sbatch` reads it
directly). Each cell is a `-profile annotate,slurm,ucr_hpcc,<provisioning>`
invocation of the vendored `nf_funannotate1`, holding SLURM resource class
fixed across cells for a genome (`conf/benchmark.yaml` pipeline: block).

| cell | funannotate | EVM | provisioning |
|---|---|---|---|
| v1.8.17_conda | 1.8.17 | perl | conda env `funannotate-1.8.17` |
| v1.8.17_container | 1.8.17 | perl | local `.sif` (Docker Hub `nextgenusfs/funannotate:v1.8.17` — ghcr has no 1.8.17 tag) |
| v1.9.0-beta10_conda | 1.9.0-beta.10 | perl | conda env `funannotate-1.9.0-beta.10` |
| v1.9.0-beta10_conda_rust | 1.9.0-beta.10 | rust | conda env `funannotate-1.9.0-beta.10-rust` |
| v1.9.0-beta10_container | 1.9.0-beta.10 | perl | local `.sif`, custom no-rust rebuild (`-norust`) |
| v1.9.0-beta10_container_rust | 1.9.0-beta.10 | rust | local `.sif` pulled from `ghcr.io/nextgenusfs/funannotate:1.9.0-beta.10` (rust-enabled by default) |

**GeneMark sidecar sharing (DESIGN.md's `runs/genemark_sidecar/`) is not
implemented** — each cell runs its own `GENEMARK_RUN` independently. This
matters differently for the two axes:
- **conda cells** (`v1.8.17_conda`, `v1.9.0-beta10_conda`,
  `v1.9.0-beta10_conda_rust`): all three already resolve `genemark_path` to
  the *same* host-licensed `/opt/linux/.../genemarkESET/4.72_lic` module —
  same binary every time — so results are already consistent across conda
  cells without a sidecar. The only cost is **redundant compute** (GeneMark
  reruns per cell instead of once per genome), not a correctness gap.
- **container cells**: `genemark_container_mode=true` routes them through the
  bundled GeneMark 4.72 in the public `teambraker/braker3:v3.1.1` image
  instead of the host module — a **different GeneMark source than the conda
  cells**. For question 3 (conda vs. container cost/accuracy) this is a real
  confound, not just a compute-efficiency gap: any conda-vs-container
  gene-content delta could partly reflect this GeneMark source difference
  rather than the funannotate wrapper/environment itself. Resolving it means
  either building the sidecar (one shared result fed to all 6 cells, as
  DESIGN.md specifies) or explicitly documenting/controlling for the
  GeneMark-source difference in `compare_predictions.py`'s analysis.

## Running the pipeline

`nf_funannotate1` is the standalone Nextflow pipeline (github.com/stajichlab/
nf_funannotate1). `scripts/pin_nextflow.py` (vendoring into `nextflow/`) was
never used — `conf/benchmark.yaml`'s `pipeline.nf_funannotate1` points
directly at the dev checkout (`/rhome/jstajich/projects/nf/nf_funannotate1`,
same repo as `/bigdata/stajichlab/jstajich/projects/nf/nf_funannotate1`), and
`launch/run_cell.sbatch` reads that same path — **this is the current
strategy**, not a temporary stand-in. It means fixes landed in that checkout
(even uncommitted ones) take effect on the next `-resume` relaunch with no
re-pin step. All three conda envs (`funannotate-1.8.17`,
`funannotate-1.9.0-beta.10`, `funannotate-1.9.0-beta.10-rust` +
`nf_funannotate1-aux`) and all three container `.sif` images referenced in
`conf/cells.tsv` are already built under `/bigdata/stajichlab/shared/condaenv`
and `/bigdata/stajichlab/shared/lib/singularity_cache` respectively.

Outputs land in `runs/<cell>/genome_annotation/<species>_<strain>/predict_results/`
(`*.gbk`, `*.gff3`) and trace/report files under `runs/<cell>/logs/nextflow/`.

## Status

- [x] Design (DESIGN.md)
- [x] Selection + fetch + harness scripts
- [x] conda env build (all 3 funannotate envs + aux env) — done
- [x] container images (all 3 `.sif`) — done
- [x] nf_funannotate1 standalone testing — confirmed runnable via both conda and singularity axes
- [x] All 6 cells launched at least once — 1 progressing cleanly, 5 blocked on the two "Known issues" above
- [ ] validate `launch/run_cell.sbatch` relaunches clear the container disk-full failures
- [ ] root-cause + fix the `SRA_FETCH` race
- [ ] resolve GeneMark sidecar sharing (or explicitly control for its absence)
- [ ] validation gate pass
- [ ] N=65 × 6 cells run to completion
- [ ] metrics + accuracy comparison + report (`scripts/summarize_report.py` not yet written)
