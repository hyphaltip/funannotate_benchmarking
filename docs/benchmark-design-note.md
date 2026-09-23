---
topic: benchmark-design-note
description: Maps the two manuscript questions (is 1.9.0 better than 1.8.17; do the rust reimplementations help) onto the cells that can actually answer them, with the confounds in each comparison made explicit. Primary axis is workflow completeness and accuracy vs RefSeq; runtime is secondary.
created: 2026-09-19
last_updated: 2026-09-19
status: active
---

# Benchmark design note: what each cell can and cannot answer

Written 2026-09-19 against the stated manuscript premise:

> 1. Is funannotate 1.9.0 faster and better, in measurable ways, than 1.8.17?
> 2. Do the rust reimplementations of the dependent tools give faster and
>    accurate genome annotation?

**Priority for the first pass: completeness and accuracy.** Get the workflow
running to completion on every genome, produce gene sets, and score them against
the RefSeq reference annotations. Runtime is a secondary axis (see
"Why runtime is not yet measurable" below).

## The rust axis is bigger than EVM -- this drives the design

funannotate 1.9.0 does not swap one component. The rust work spans **EVM,
PASA and Trinity**. That matters because the cells differ in how much of that
stack they vary:

| comparison | EVM | PASA | Trinity | salmon | verdict |
|---|---|---|---|---|---|
| beta12_container vs beta12_container_rust | **perl vs rust** | **perl vs rust** | **perl vs rust** | same 2.7.0 | **the full rust-stack comparison** |
| beta11_conda vs beta11_conda_rust | perl vs rust | conda vs source-built | 2.15.2 vs 2.16.1_rust | 1.10.3 vs 2.7.0 | full-stack, but 3 variables at once |
| v1.8.17_conda vs beta11_conda | 1.1.1 vs 1.9.0's perl | conda 2.5.3 both | 2.15.2 both | 1.10.3 both | version, plus a gmap->minimap2 aligner change |

**CORRECTION 2026-09-19.** An earlier version of this note claimed the beta.12
pair varies EVM only. That was WRONG. It came from comparing the perl DRIVER
scripts (`Launch_PASA_pipeline.pl`, `pasa_asmbls_to_training_set.dbi`, `Trinity`,
`salmon_runner.pl`) -- which are indeed byte-identical between the images -- and
missing that those drivers dispatch to rust HELPER BINARIES when present.

Both images contain the rust implementations; the norust image **deactivates**
them by renaming. Full inventory of `*.disabled` in `-norust.sif` (the rust image
has none):

| component | disabled binaries in norust |
|---|---|
| EVM | `evidence_modeler` |
| PASA | `pasa_rust`, `cdbyank_rust`, `slclust_rust`, `faidx_rust` |
| Trinity | `sam_to_read_coords`, `define_coverage_partitions`, `extract_reads_per_partition`, `fragment_coverage_writer` |

`FUNANNOTATE_EVM_ENGINE` (perl vs rust) is the env-level switch, but the
deactivation is broader than EVM.

**Consequence:** the beta.12 pair IS the full rust-stack comparison -- EVM, PASA
and Trinity together -- and is the right basis for manuscript question 2. The
beta.11 conda pair is no longer needed for that purpose (and carries
salmon/Trinity confounds anyway).

**Corollary for methodology:** do NOT seed shared training evidence between the
beta.12 arms for the main run. Since rust PASA and rust Trinity act during
TRAINING, seeding would bypass the components being measured and would also
erase the runtime difference. Seeding is correct ONLY for a deliberate
EVM-isolation experiment (as done for the Chaetomium single-genome test, which
therefore measures EVM alone, not the full stack).

## Primary goal: completeness and accuracy vs RefSeq

### Scoreability

15 of the 16 benchmark genomes have RefSeq (GCF) references already downloaded
under `results/reference_annotations/`. **Allomyces macrogynus is
`GCA_000151295.1` -- GenBank only, no RefSeq annotation**, so it cannot be scored
against a reference. Keep it for gene-count and runtime comparison; exclude it
from accuracy tables, or the denominator is 15, not 16.

### Metrics, in priority order

1. **Completion rate** -- of 16 genomes x N cells, how many produce
   `predict_results/*.proteins.fa` without degrading to ab-initio. This is the
   headline for "does the workflow run", and today it is the honest weak point:
   most cells are partial (see status below).
2. **Gene-count agreement with RefSeq** -- cheap, already computed ad hoc.
   Useful as a sanity signal, NOT as an accuracy claim: two annotations can
   match on count and differ structurally.
3. **Structural concordance vs RefSeq** -- `scripts/compare_predictions.py`
   (exon/CDS boundary agreement). This is the real accuracy measurement and the
   one the manuscript should lead with.
4. **BUSCO completeness** on the predicted proteome, with the genome-mode BUSCO
   as the reference point. The genome-vs-proteome BUSCO gap is a direct
   "you undercalled" signal (see docs/malassezia-evm-undercall.md).

### Known accuracy caveat already in hand

Malassezia globosa CBS 7966 was annotated at 11-31% of its RefSeq gene count by
*every* configuration including the independent production pipeline, traced to a
227k-read library. It has been re-sourced from SRP487370 (2.2M paired reads,
86.4% mapping) and is re-running. Any accuracy table must state the RNA-seq depth
per genome, because depth -- not version, not backend -- was the dominant effect
in the one genome where the pipelines disagreed wildly.

## Secondary goal: runtime

### Why runtime is not yet measurable

`scripts/collect_metrics.py` globs **every** `annotate_trace.*.txt` snapshot per
cell and sums across them. The merged traces are dominated by the 2026-09-07..18
debugging era:

| cell | trace files | COMPLETED rows | FAILED rows |
|---|---|---|---|
| v1.8.17_conda | 20 | 42 | 232 |
| v1.9.0-beta11_conda | 10 | 40 | 162 |
| v1.9.0-beta11_conda_rust | 15 | 27 | 224 |
| v1.9.0-beta11_container_rust | 4 | 42 | 39 |

So `cpu_hours_sum` currently includes every failed MySQL/PASA/Trinity attempt,
unevenly across cells (39 failed rows in one cell, 224 in another). **Any speed
comparison built on this is an artifact of how much debugging each cell needed.**
Two fixes, either acceptable: a timestamp cutoff that only merges traces written
after the fixes landed, or last-attempt-wins keyed on task name. This must be
resolved before a runtime number appears anywhere.

**Fixed 2026-09-23 (last-completion-wins).** `collect_metrics.py` now keeps
only the last COMPLETED row per (cell, task name) in the timing columns. It
drops CACHED rows, which repeat the original run's realtime and were counted
twice. FAILED/ABORTED rows and superseded completions go into separate
`overhead_task_time_hours_sum` / `overhead_cpu_hours_sum` columns. Genome
status is `success` when every task eventually completed. Limit: failed rows
with no exit code (SLURM-side kill) also have no `%cpu`, so
`overhead_cpu_hours_sum` undercounts them; `overhead_task_time_hours_sum` does
not.

The trace also gains `attempt,cpus,memory,queue,hostname,submit,start,
complete,peak_rss` (nf_funannotate1 `profile_annotate.config`,
`profile_genemark_sidecar.config`), for launches from 2026-09-23 on only.
Before this, a TRAIN task that completed on attempt 1 (8 cpus, epyc) and one
that completed on attempt 3 (24 cpus, highmem) could not be told apart, and
`peak_rss_gb_max` was read from `rss` (resident memory, not the peak). Older
rows leave `cpus`/`queues`/`attempts`/`wall_clock_hours` blank and report
`rss_field=rss`.

### What to measure once it is fixed

Per genome, per cell: wall-clock and CPU-hours for `train` and `predict`
separately. Reporting them separately matters because the rust work targets
different stages -- PASA/Trinity in train, EVM in predict -- and an aggregate
would blur exactly the distinction the manuscript is making.

## Cell status (2026-09-19, 16-genome set)

| cell | train | predict | note |
|---|---|---|---|
| v1.8.17_conda | 9 | 11 | running; 6 new genomes + rescued Malassezia in flight |
| v1.8.17_container | 0 | 0 | DEFERRED -- PASA/mysql shells out to host commands from inside the container |
| v1.9.0-beta11_conda | 8 | 10 | complete for the old 10-genome set; 6 new genomes not started |
| v1.9.0-beta11_conda_rust | 6 | 7 | running (from-scratch Trinity after the salmon fix) |
| v1.9.0-beta11_container_rust | 7 | 10 | complete for the old set |
| v1.9.0-beta12_container | 0 | 0 | new; Chaetomium single-genome EVM test running |
| v1.9.0-beta12_container_rust | 0 | 0 | new; Chaetomium single-genome EVM test running |

Gaps that block the primary goal:
- No cell has all 16 genomes complete.
- No beta.12 **conda** cells exist, so conda-vs-container at beta.12 is not
  testable; Q3 (provisioning) currently only has the beta.11 pair.
- `v1.8.17_container` is unfixed, so 1.8.17 has no container arm.

## Recommended sequence

1. Finish `v1.8.17_conda` at 16/16 -- the stable comparator, unaffected by the
   beta.12 rebuild.
2. Run `v1.9.0-beta12_container` (perl) at 16/16 -- the 1.9.0 arm for Q1.
3. Run `v1.9.0-beta12_container_rust` at 16/16, **training independently** (do
   not seed), so the pair also reports a runtime difference. Seeding remains the
   right choice only for a dedicated EVM-isolation experiment like the Chaetomium
   test.
4. Score 1+2 with `compare_predictions.py` + BUSCO across the 15 scoreable
   genomes -- this is the primary manuscript result.
5. Fix `collect_metrics.py`, then regenerate runtime tables from post-fix traces
   only.

## Provenance warning

`funannotate --version` reports **`v1.9.0` for beta.11 AND beta.12** -- the
binary cannot distinguish them, and every beta.11 artifact already on disk
reports the same string. Record the container path/digest or conda env name as
the provenance field; do not rely on the version string in any results table.
