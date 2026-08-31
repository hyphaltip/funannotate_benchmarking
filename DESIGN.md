# Funannotate Benchmark Design

Status: approved, pre-implementation
Supersedes: PLAN.md (kept for historical context; this document is authoritative)
Date: 2026-08-30

## Purpose

Build a standalone, reproducible benchmark comparing two funannotate versions
(1.8.17 and 1.9.0-beta.10) across execution environments, to answer three
questions:

1. Does the 1.9.0 rust reimplementation of Trinity/EVM/PASA improve
   performance over 1.8.17? 1.9.0-beta.10 exposes a +/- rust build toggle,
   so this is measured directly within 1.9.0 (rust on vs. off), isolated
   from the 1.8.17-vs-1.9.0 code-change comparison below.
2. What do the 1.8.17 -> 1.9.0 code changes do to runtime, resource use, and
   gene-content predictions, independent of environment and independent of
   the rust toggle (compare 1.8.17 against 1.9.0 with rust off)?
3. What does conda vs. container execution cost or change, independent of
   funannotate version?

This is deliberately separated from Fungi_BFD's production pipeline so the
benchmark isn't entangled with lab-specific infra and can be rerun cleanly as
funannotate evolves.

## Non-goals (deferred)

- Reuse-of-trained-ab-initio-parameters study (speedup/accuracy of skipping
  per-strain Augustus/GeneMark training) — stays in PLAN.md's "Future"
  section, out of scope for this pass.
- Module-load (Fungi_BFD-style) execution as a distinct axis — superseded by
  the conda-vs-container axis below, which achieves the same
  "non-container" comparison without reviving Fungi_BFD's module-system
  nextflow code.

## Prior findings this design accounts for

From `../Fungi_BFD/.living/findings/`:

- `funannotate-genemark-contribution.md` — GeneMark supplied up to 68% of
  gene models on small fragmented genomes in a prior n=3 study; it is
  licensing-restricted and reportedly absent from the rust/container build.
  Drives the GeneMark sidecar requirement below.
- `genome-size-architecture.md` — genome size correlates with taxonomy;
  affects dataset design (size cap interacts with taxonomic balance).
- `rnaseq-annotation-evidence.md` — only ~20% of species have transcript
  assemblies, availability is bimodal, not uniform across taxa. Drives the
  per-genome RNA-seq metadata column.
- `intraspecific-ani-diversity.md` — some nominal "species" are cryptic
  complexes with large within-label divergence; relevant to the ≤2/genus,
  ≤1/species selection rule.
- `annotation-failure-modes.md` — small fragmented genomes fail at ab-initio
  training thresholds, large genomes fail via OOM/timeout. Drives
  first-class completion-status tracking.
- `funannotate-abinitio-reuse-validation.md` — BUSCO-neutral results can
  hide real gene-count/locus differences. Drives the coordinate-aware
  diffing requirement.

## Repository layout

```
Funannotate_benchmarking/
  PLAN.md                    (original request, historical)
  DESIGN.md                  (this document)
  samples.csv                (append-only; dataset selection output)
  scripts/
    select_genomes.py        # append-only sampler, taxonomic quotas
    fetch_genomes.py         # NCBI datasets download, RefSeq-preferred
    validate_harness.py      # nf_funannotate1 sanity check vs Fungi_BFD outputs
    collect_metrics.py       # parse nextflow trace + resource logs -> tidy table
    compare_predictions.py   # BUSCO + gffcompare + bedtools vs RefSeq GFF
    summarize_report.py      # final tables/plots
  input_clean_genomes/       # symlinks from ../Fungi_BFD_runs/input_clean_genomes
  nextflow/                  # vendored/pinned copy of nf_funannotate1
  runs/
    v1.8.17_conda/
    v1.8.17_container/
    v1.9.0-beta10_conda/
    v1.9.0-beta10_conda_rust/
    v1.9.0-beta10_container/
    v1.9.0-beta10_container_rust/
    genemark_sidecar/        # shared GeneMark evidence, keyed by genome, reused by all 4 cells
  results/
    metrics.tsv               # per genome x cell: runtime, peak mem, CPU-hours, status
    gene_content_comparison.tsv
    report.md / report.html
```

## Dataset selection

`scripts/select_genomes.py`:

- **Append-only**: reads existing `samples.csv`, never removes prior rows;
  adds new candidates only for taxa/quotas not yet satisfied.
- **Target N = 65 genomes total.**
- **Taxonomic floor quotas** (unchanged from PLAN.md, must all be met):
  - ≥10 Ascomycota, including ≥5 Pezizomycotina
  - ≥5 Basidiomycota, including ≥2 Agaricomycotina, ≥1 Pucciniomycotina,
    ≥1 Ustilaginomycotina
  - ≥3 Mucoromycota
  - ≥1 Chytridiomycota
  - ≥1 Blastocladiomycota
  - Remaining ~45 genomes fill proportionally across these groups (weighted
    toward Ascomycota/Basidiomycota diversity, the largest and most
    genomically diverse groups) to reach N=65.
- **Constraints**: genome size ≤75Mb; ≤2 genomes per genus; ≤1 genome per
  species; RefSeq (`GCF_`) accessions preferred over GenBank when both
  exist for a species.
- **Metadata columns** recorded at selection/fetch time (not computed later):
  `accession`, `species`, `genus`, `phylum`, `subphylum`, `genome_size_mb`,
  `rnaseq_available` (bool + SRA/source id if known), `refseq_annotation_release`
  (PGAP vs. legacy, from the NCBI assembly report), `assembly_level`,
  `busco_lineage_used`.
- Genomes are pulled from already-masked assemblies via symlink from
  `../Fungi_BFD_runs/input_clean_genomes` where available; otherwise fetched
  fresh via `fetch_genomes.py` and masked using the same method Fungi_BFD
  uses, to keep masking methodology consistent across the dataset.

## Execution design: 6 cells per genome

1.8.17 gets the conda/container axis only (no rust toggle exists for it).
1.9.0-beta.10 gets conda/container **and** the +/- rust build toggle:

| | conda | container |
|---|---|---|
| **1.8.17** | v1.8.17_conda | v1.8.17_container |
| **1.9.0-beta.10, rust off** | v1.9.0-beta10_conda | v1.9.0-beta10_container |
| **1.9.0-beta.10, rust on** | v1.9.0-beta10_conda_rust | v1.9.0-beta10_container_rust |

- All 6 cells driven by the same vendored nf_funannotate1 pipeline, using its
  native conda and container execution profiles — no Fungi_BFD module-load
  code is touched. nf_funannotate1 is currently being tested standalone
  (outside Fungi_BFD) to confirm it can run this way at all; that testing
  feeds directly into the validation gate below.
- Same SLURM resource class (CPUs, memory, partition) is held fixed across
  all 6 cells for a given genome, so only version, rust toggle, and
  environment vary.
- **GeneMark sidecar**: GeneMark is run once per genome, outside all
  container/conda cells (licensed step), cached under
  `runs/genemark_sidecar/`, and fed as evidence identically into all 6
  cells — so no cell is structurally missing GeneMark's contribution.
- **Single run per cell** (no replicates across the full N=65). If runtime
  numbers look implausibly noisy after the fact, a small representative
  subset can be rerun to establish variance — not planned upfront.
- **Completion status is a first-class output**: each cell records
  `success` / `oom` / `timeout` / `training_failed` / other, alongside
  performance metrics. Failed cells are reported, not silently dropped from
  the summary.

## Metrics collection

`scripts/collect_metrics.py` parses Nextflow's `-with-trace` /
`-with-report` output per cell into `results/metrics.tsv`:

- Wall-clock time, CPU-hours, peak memory, per process and pipeline-total.
- Normalized per genome by genome size (Mb) and BUSCO gene count, so
  cross-genome comparisons account for the obvious "bigger genome = slower"
  effect.
- Container pull/build time tracked separately as a one-time cost, not
  amortized into per-genome numbers.

## Gene-content comparison

`scripts/compare_predictions.py`, for genomes with a RefSeq GFF available:

- BUSCO completeness (complete/fragmented/missing) per cell.
- `gffcompare` class-code summary against the RefSeq GFF.
- `bedtools intersect`-based coordinate-aware diffing to catch same-locus
  evidence-source swaps that BUSCO/count deltas alone would hide (the
  failure mode documented in `funannotate-abinitio-reuse-validation.md`).
- Every comparison is annotated with the genome's `refseq_annotation_release`
  metadata — RefSeq GFFs from different NCBI pipeline eras are not treated
  as a single uniform gold standard.

## Validation gate (must pass before trusting the N=65 run)

`scripts/validate_harness.py` runs the vendored nf_funannotate1 on 2-3
genomes Fungi_BFD has already annotated well, and diffs outputs against
Fungi_BFD's known-good results (BUSCO score, gene count, spot-check GFF
coordinates). This gate exists because nf_funannotate1 is newer and less
battle-tested than Fungi_BFD's production pipeline — the benchmark must not
end up measuring wrapper immaturity instead of funannotate itself.

## Open items to resolve during implementation

- Confirm nf_funannotate1's conda profile actually covers all steps needed
  (predict, EVM, PASA, GeneMark integration) at parity with its container
  profile — this was asserted but not yet verified against the current
  nf_funannotate1 repo state.
- Confirm nf_funannotate1 can in fact run standalone (currently being
  tested); if it cannot, the harness base decision (nf_funannotate1 vs.
  Fungi_BFD nextflow) needs revisiting.
- Confirm how the +/- rust build toggle is actually selected in
  nf_funannotate1 for 1.9.0-beta.10 (build flag, container tag, env var).
- 6 cells/genome (vs. the 4 originally budgeted for) increases compute
  further; confirm N=65 still holds or needs trimming again.
- Confirm what GeneMark license/install is available on the cluster for the
  sidecar step.
- Confirm SLURM resource class to hold fixed across cells (partition,
  cpus, mem) before the first real run.
