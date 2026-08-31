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
  scripts/
    select_genomes.py        append-only sampler (taxonomic quotas, ≤2/genus, ≤1/species)
    fetch_genomes.py         masked-genome + RefSeq-GFF provisioning (symlink or download)
    build_run_samples.py     benchmark samples.csv → per-cell nf_funannotate1 samples.csv
    validate_harness.py      gate: nf_funannotate1 vs known-good Fungi_BFD outputs
    collect_metrics.py       nextflow trace/report → results/metrics.tsv
    compare_predictions.py   BUSCO + gffcompare + bedtools vs RefSeq GFF
    summarize_report.py      results → results/report.md + figures/tables
  input_clean_genomes/       (generated) symlinks to masked genomes, keyed by accession
  runs/
    v1.8.17_conda/  v1.8.17_container/
    v1.9.0-beta10_conda/  v1.9.0-beta10_conda_rust/
    v1.9.0-beta10_container/  v1.9.0-beta10_container_rust/
    genemark_sidecar/        shared licensed GeneMark output, keyed by genome
  results/                   (generated) metrics.tsv, gene_content_comparison.tsv, report.*
```

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

# 4. Validation gate (2-3 genomes Fungi_BFD already annotated well)
sbatch launch/validate_harness.sbatch

# 5. Launch the 6 cells (per genome, fixed SLURM resource class)
sbatch launch/run_all_cells.sbatch            # or launch/run_cell.sbatch <cell>

# 6. Collect + compare + report (can run incrementally as cells finish)
python3 scripts/collect_metrics.py     --cells-dir runs
python3 scripts/compare_predictions.py --cells-dir runs
python3 scripts/summarize_report.py    --cells-dir runs
```

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

See `conf/cells.tsv`. Each cell is a `-profile annotate,slurm,<provisioning>`
invocation of the vendored `nf_funannotate1`, holding SLURM resource class fixed
across cells for a genome. GeneMark runs once per genome in
`runs/genemark_sidecar/` (licensed step) and is fed identically into all 6 cells.

| cell | funannotate | EVM | provisioning |
|---|---|---|---|
| v1.8.17_conda | 1.8.17 | perl | conda env `funannotate-1.8.17` |
| v1.8.17_container | 1.8.17 | perl | `ghcr.io/nextgenusfs/funannotate:1.8.17` |
| v1.9.0-beta10_conda | 1.9.0-beta.10 | perl | conda env `funannotate-1.9.0-beta.10` |
| v1.9.0-beta10_conda_rust | 1.9.0-beta.10 | rust | conda env `funannotate-1.9.0-beta.10-rust` |
| v1.9.0-beta10_container | 1.9.0-beta.10 | perl | `ghcr.io/nextgenusfs/funannotate:1.9.0-beta.10` |
| v1.9.0-beta10_container_rust | 1.9.0-beta.10 | rust | rust container tag (open item) |

## Running the pipeline

`nf_funannotate1` is the standalone Nextflow pipeline (github.com/stajichlab/
nf_funannotate1). The harness normally runs a **pinned vendored copy** under
`nextflow/` (see `scripts/pin_nextflow.py`); until that copy is pinned, set
`NF_FUNANNOTATE1` (see `conf/benchmark.yaml`) to the dev checkout. conda envs
are built once into `/bigdata/stajichlab/shared/condaenv` (see
nf_funannotate1 `environments/conda/`) — currently in progress.

Outputs land in `runs/<cell>/genome_annotation/<species>_<strain>/predict_results/`
(`*.gbk`, `*.gff3`) and trace/report files under `runs/<cell>/logs/nextflow/`.

## Status

- [x] Design (DESIGN.md)
- [x] Selection + fetch + harness scripts
- [ ] conda env build (nf_funannotate1 environments/conda) — in progress
- [ ] nf_funannotate1 standalone testing — in progress
- [ ] validation gate pass
- [ ] N=65 × 6 cells run
- [ ] metrics + accuracy comparison + report
