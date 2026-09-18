# Funannotate Benchmarking

Reproducible benchmark comparing **funannotate 1.8.17 vs 1.9.0-beta.11** (with and
without the Rust-optimized EVM/PASA/Trinity build) across **conda vs container**
execution, on a structured random sample of 65 fungal genomes.

(1.9.0-beta.10 was retired from the active plan 2026-09-07 in favor of
beta.11 — see DESIGN.md's 2026-09-07 update note. Its completed run data
stays under `runs/v1.9.0-beta10_*/` as a historical record.)

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
    run_genemark_sidecar.sbatch  run GeneMark once per genome (shared by all 6 cells)
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
    genemark_sidecar/output/  one <out>.genemark.gtf/.mod per genome, shared by all 6 cells
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

# 4. GeneMark sidecar — run once per genome, shared by all 6 cells (see
#    "The 6 cells" below). Do this before/alongside step 5; run_cell.sbatch
#    warns (doesn't fail) if a cell launches before this has produced results.
sbatch launch/run_genemark_sidecar.sbatch

# 5. Launch the cells (per genome, fixed SLURM resource class; image/env
#    per cell comes from conf/cells.tsv — see "The 6 cells" below; currently
#    5 active cells; a 6th (v1.9.0-beta11_container, norust) is deferred)
sbatch launch/run_all_cells.sbatch            # all 6, one sbatch job per cell
sbatch launch/run_cell.sbatch v1.8.17_conda   # relaunch/retry a single cell (-resume)

# 6. Validation gate — once a cell has predict output for a few genomes,
#    diff it against Fungi_BFD's known-good results for the same genome tag
sbatch launch/validate_harness.sbatch v1.8.17_conda

# 7. Collect + compare + report (can run incrementally as cells finish)
python3 scripts/collect_metrics.py     --cells-dir runs
python3 scripts/compare_predictions.py --cells-dir runs
python3 scripts/summarize_report.py    --cells-dir runs   # not yet written
```

Steps 1-3 have already run once (dataset selected; `runs/<cell>/` samplesheets
built). All 6 original cells (including the now-retired beta.10 ones) were
launched at least once (`runs/<cell>/launch_*.log`), **before** the GeneMark
sidecar (step 4) existed — those launches used each cell's own independent
`GENEMARK_RUN`, not the sidecar. **Current status (2026-09-07):** all cells
are currently running against a pared-down 10-genome smoke-test subset
(`runs/<cell>/samples.csv`, not yet committed) rather than the full N=65, to
shake out per-cell provisioning bugs cheaply first — see "Known issues"
below for what's still open (`v1.8.17_container` Augustus scripts,
possible container disk-full/`swissprot_fungi.faa` recurrence on
`v1.9.0-beta11_container_rust`). `v1.8.17_conda` is the only cell confirmed
clean end-to-end (9/10 genomes on the smoke-test subset; the 10th,
`Malassezia_globosa`, should now be fixed by the RNA-seq pre-seed above).
Re-launching a cell with `launch/run_cell.sbatch <cell>` is safe (`-resume`)
once the underlying fix for its failure lands; a relaunch after step 4 has
run will pick up the shared sidecar GTF instead of
retraining GeneMark in-cell (a resumed genome whose predict already completed
under the old per-cell GeneMark keeps its old result until something
invalidates that cache — full apples-to-apples parity means running the
sidecar before any cell's first launch, not applicable retroactively without
a fresh run).

Every cell must always start from an **already-masked genome** (never
re-run `GENOME_CLEAN`) and **already-downloaded RNA-seq** (never hit
`SRA_FETCH`/`SRA_QUERY_BATCH`). `launch/run_cell.sbatch` enforces this on
every launch by running, before invoking nextflow:

```bash
python3 scripts/preseed_clean_genomes.py --runs-dir runs   # symlinks masked genomes (pre-existing)
python3 scripts/preseed_rnaseq_reads.py  --runs-dir runs   # symlinks already-normalized reads + sra_query cache
```

`preseed_rnaseq_reads.py` symlinks `<species_tag>_norm_R1/_R2/_SE.fastq.gz`
from the shared, already-downloaded pool
(`conf/benchmark.yaml` `fetch.rnaseq_reads_dir`, currently
`Fungi_BFD_runs/rnaseq_reads`) into each cell, and writes a per-species
`rnaseq_reads/sra_query/<tag>.sra_query.csv` from the master SRA manifest
(`pool.rnaseq`). `run_cell.sbatch` passes `--skip_sra_query true`, so
`FETCH_RNASEQ` (subworkflows/local/fetch_rnaseq.nf) reuses both caches and
never makes an NCBI network call. Both scripts are idempotent/rerunnable —
safe to run before every launch, on any subset of genomes.

### Known issues (as of 2026-09-03, RNA-seq fetch race fixed 2026-09-07)

- ~~**`SRA_FETCH` chunked-download race**~~ **fixed 2026-09-07**: the race
  (hit `v1.8.17_container`, `v1.9.0-beta10_conda`, `v1.9.0-beta10_conda_rust`,
  and recurred against `Malassezia_globosa`/SRR12496022 in every conda cell
  through 2026-09-06) was never actually a bug in fastq-dump chunking —
  `Fungi_BFD_runs/rnaseq_reads` already had the correctly-resolved reads for
  every benchmark species sitting unused; each cell was independently
  re-fetching from SRA via its own `SRA_FETCH`/`SRA_QUERY_BATCH` instead of
  reusing them. `scripts/preseed_rnaseq_reads.py` (see above) closes this by
  pre-seeding every cell from that pool before launch, so `SRA_FETCH` never
  runs at all for a covered species.
- ~~**Stale/truncated per-cell RNA-seq reads silently trusted**~~ **fixed
  2026-09-07** (found live, same day as the fix above): the first version of
  `preseed_rnaseq_reads.py` only checked whether a dest file *existed*, so a
  cell that already had a stale/truncated read pair from before this script
  existed (e.g. `v1.8.17_conda`'s and `v1.8.17_container`'s
  `Allomyces_macrogynus_norm_R1/_R2.fastq.gz` were 23-byte corrupt leftovers
  from an aborted fetch, sitting untouched since 2026-09-02) got skipped as
  "already present" instead of overridden — causing a genuine
  `funannotate train` failure (`Trinity de novo assembly failed`, BAM had 0
  reads) on relaunch. Fixed by making the master pool always authoritative:
  the script now replaces any dest that isn't already the correct symlink to
  the pool copy (`relinked` in its output), not just gaps (`linked`).
- **Container `module`/`singularity` not found *inside* the container —
  hits every genome** (seen on `v1.8.17_container`, 2026-09-07 relaunch,
  unrelated to the RNA-seq issue above): `FUNANNOTATE_TRAIN`'s PASA/mysql
  setup first fails to find `/rhome/jstajich/.pasa/pasa_conf/conf.txt`, then
  `.command.sh` shells out to `module load` / `singularity` from *inside*
  the apptainer container (host commands, not available in-container).
  Confirmed failing on every genome in the cell (Aspergillus_fumigatus,
  Malassezia_globosa, ...), not an isolated genome issue — **blocking**;
  don't bother relaunching this cell until it's root-caused. Likely the same
  nested-container assumption as the Augustus-scripts issue below.
  **Deferred 2026-09-08**: `v1.8.17_container` is not being relaunched —
  question 3 (conda vs. container) is covered in the meantime by the
  1.9.0-beta.11 pair (`v1.9.0-beta11_conda_rust` vs.
  `v1.9.0-beta11_container_rust`); the 1.8.17-generation conda-vs-container
  comparison specifically is just not available until this is root-caused.
- **Container disk-full** (seen on the now-retired `v1.9.0-beta10_container`,
  `v1.9.0-beta10_container_rust`, both `[FAILED] failed=41`): apptainer fell
  back to full sandbox extraction of the `.sif` into node-local `/tmp`
  (`squashfuse not found`) and concurrent tasks exhausted disk. Likely fixed
  by nf_funannotate1's singularity-axis `process.shell = ['/bin/bash']` change
  (a login shell was re-sourcing `/etc/profile` inside `.command.run` and
  wiping the squashfuse-bearing PATH before the container command ran) — not
  container-image-specific, so **watch for it on `v1.9.0-beta11_container_rust`
  too** on its first relaunch.
- ~~**Container `swissprot_fungi.faa` missing**~~ **fixed 2026-09-07**
  (seen on the now-retired `v1.9.0-beta10_container`, `_container_rust`,
  2026-09-06; recurred on `v1.9.0-beta11_container_rust` 2026-09-07):
  `FUNANNOTATE_PREDICT` aborted with `<cell>/lib/swissprot_fungi.faa is not a
  valid file, exiting`. Root cause: `params.proteins`/`params.sbt_template`
  default to `<launchDir>/lib/swissprot_fungi.faa` / `lib/template.sbt`
  (`conf/profile_annotate.config`), but — unlike `lib/augustus/`, which
  `SETUP_AUGUSTUS_CONFIG` auto-seeds — nothing populates either file for a
  new cell. Conda profile reads the real path and mostly got away with it;
  container profile bind-mounts `<launchDir>/lib` and hard-fails when the
  file inside is missing. The `_stale_pre_swissprot_fix_2026-09-04` symlink
  was the same fix landed at the project *root* `lib/` instead of any cell's
  own `runs/<cell>/lib/` — never a Nextflow launchDir, so it never took
  effect. `scripts/preseed_funannotate_lib.py` (wired into
  `launch/run_cell.sbatch`) now symlinks both files into every cell from
  `conf/benchmark.yaml` `fetch.funannotate_lib_dir`
  (`Fungi_BFD/lib/`) on every launch.
- **`v1.8.17_container` missing Augustus helper scripts** (seen 2026-09-06,
  still an active cell): `FUNANNOTATE_PREDICT` fails with
  `unable to locate ... join_mult_hints.pl`, `gff2gbSmallDNA.pl` under
  `<cell>/lib/augustus/3.5/scripts/`. **Not yet root-caused.**
- **PASA MySQL checkpoint mismatch** (found 2026-09-14 after a relaunch broke
  `v1.8.17_conda` end-to-end, which had previously been the one cell confirmed
  clean): `FUNANNOTATE_TRAIN` spins up a fresh, non-persistent MariaDB sidecar
  on every task attempt, but PASA's own on-disk checkpoint markers
  (`training/pasa*/__pasa_<genome>_pasa_mysql_chkpts/*.ok`) persist across
  relaunches/`-resume`. Once a genome's PASA run got past `create_db.ok` once,
  any later re-attempt sees that stale checkpoint, skips re-creating the
  database against the new (empty) MariaDB instance, and dies at
  `update_fli_status.dbi` with `Unknown database '<genome>_pasa'` — hit every
  re-attempted genome in `v1.8.17_conda` (Malassezia_globosa,
  Saccharomyces_cerevisiae, Umbelopsis_ramanniana, Aspergillus_fumigatus,
  ...), including ones whose `predict_results/*.gff3` already existed from an
  earlier success (an upstream input change still invalidated nextflow's task
  cache and forced a re-run). **Fixed 2026-09-17**:
  `scripts/preseed_pasa_checkpoints.py` (wired into `launch/run_cell.sbatch`)
  unconditionally clears every genome's `pasa*` checkpoint dir(s) before each
  launch — harmless for a genome nextflow skips via `-resume`, since a cached
  task never looks at that directory again.

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
| v1.9.0-beta11_conda | 1.9.0-beta.11 | perl | conda env `funannotate-1.9.0-beta.11` (built by `scripts/build_conda_env_beta11.sbatch`) |
| v1.9.0-beta11_conda_rust | 1.9.0-beta.11 | rust | conda env `funannotate-1.9.0-beta.11-rust` (built by `scripts/build_conda_env_beta11.sbatch`) |
| v1.9.0-beta11_container_rust | 1.9.0-beta.11 | rust | local `.sif` pulled from `ghcr.io/nextgenusfs/funannotate:1.9.0-beta.11` (rust-enabled by default; perl-EVM `v1.9.0-beta11_container` (norust rebuild) deferred — see DESIGN.md "Execution design") |

**GeneMark sidecar sharing is implemented** (`launch/run_genemark_sidecar.sbatch`
→ `nf_funannotate1`'s `genemark_sidecar.nf`, `-profile genemark_sidecar,slurm,singularity`):
GeneMark runs exactly once per genome — always container-mode, always the
public `teambraker/braker3:v3.1.1` image, always fresh ES self-training with
no RNA-seq hints (so the result cannot depend on any cell's own training) —
and publishes `<out>.genemark.gtf`/`.mod` to `runs/genemark_sidecar/output/`.
Every cell's `launch/run_cell.sbatch` invocation passes
`--genemark_sidecar_dir` pointing at that directory, so **no cell runs its own
`GENEMARK_RUN` any more** — GeneMark is deliberately held out of the
comparison (not a variable under test) and computed 1x per genome instead of
6x. A genome missing from the sidecar output degrades that genome to
`--auto-skip-genemark` in every cell (a nextflow warning, not a hard
failure), rather than silently reverting to the old per-cell behavior.

This also fixes a real confound the previous (unimplemented) state had:
conda cells were resolving GeneMark to the host-licensed
`genemarkESET/4.72_lic` module while container cells used braker3's bundled
GeneMark instead — a different GeneMark *source* between the two axes being
compared for question 3 (conda vs. container). All 6 cells now use the exact
same GeneMark result per genome, so any conda-vs-container delta observed
going forward reflects the funannotate wrapper/environment only.

## Running the pipeline

`nf_funannotate1` is the standalone Nextflow pipeline (github.com/stajichlab/
nf_funannotate1). `scripts/pin_nextflow.py` (vendoring into `nextflow/`) was
never used — `conf/benchmark.yaml`'s `pipeline.nf_funannotate1` points
directly at the dev checkout (`/rhome/jstajich/projects/nf/nf_funannotate1`,
same repo as `/bigdata/stajichlab/jstajich/projects/nf/nf_funannotate1`), and
`launch/run_cell.sbatch` reads that same path — **this is the current
strategy**, not a temporary stand-in. It means fixes landed in that checkout
(even uncommitted ones) take effect on the next `-resume` relaunch with no
re-pin step. All conda envs (`funannotate-1.8.17`,
`funannotate-1.9.0-beta.11`, `funannotate-1.9.0-beta.11-rust` +
`nf_funannotate1-aux`) and container `.sif` images referenced in
`conf/cells.tsv` are already built under `/bigdata/stajichlab/shared/condaenv`
and `/bigdata/stajichlab/shared/lib/singularity_cache` respectively (the
retired beta.10 envs/images are left in place, just no longer referenced by
`conf/cells.tsv`).

Outputs land in `runs/<cell>/genome_annotation/<species>_<strain>/predict_results/`
(`*.gbk`, `*.gff3`) and trace/report files under `runs/<cell>/logs/nextflow/`.

## Status

- [x] Design (DESIGN.md)
- [x] Selection + fetch + harness scripts
- [x] conda env build (all 3 funannotate envs + aux env) — done
- [x] container images (all 3 `.sif`) — done
- [x] nf_funannotate1 standalone testing — confirmed runnable via both conda and singularity axes
- [x] All 6 cells launched at least once — 1 progressing cleanly, 5 blocked on the two "Known issues" above
- [x] GeneMark sidecar sharing (`launch/run_genemark_sidecar.sbatch` + `--genemark_sidecar_dir`) — implemented, not yet run against the full N=65 set
- [ ] validate `launch/run_cell.sbatch` relaunches clear the container disk-full failures
- [ ] root-cause + fix the `SRA_FETCH` race
- [ ] validation gate pass
- [ ] N=65 × 6 cells run to completion
- [ ] metrics + accuracy comparison + report (`scripts/summarize_report.py` not yet written)
