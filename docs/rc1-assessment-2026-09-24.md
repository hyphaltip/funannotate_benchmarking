---
topic: rc1-assessment-2026-09-24
description: Accuracy (gffcompare vs RefSeq, BUSCO) and runtime of funannotate 1.9.0-rc.1 conda and container cells vs 1.8.17, written overnight 2026-09-24 while the container cells re-run on the rebuilt rc.1 SIFs.
created: 2026-09-24
status: interim -- container numbers are from the PREVIOUS rc.1 image; new-image re-run in progress
---

# funannotate 1.9.0-rc.1 assessment (interim, 2026-09-24)

## 1. What ran overnight

| step | detail |
|---|---|
| SIF rebuild | job 29045856, COMPLETED 00:54. Both images: PASA `PASA_transcripts_and_assemblies_to_GFF3.dbi` md5 `c3d2f07b…` (fix present). norust has 9 `*.disabled` rust binaries; `FUNANNOTATE_EVM_ENGINE` perl/rust as expected. Both still report `v1.9.0-beta.13.dev0+g6f3eaad` (same string as the old image). |
| new SIF paths | `funannotate-1.9.0-rc.1.20260924.sif` (sha256 `f348bc40…`), `funannotate-1.9.0-rc.1-norust.20260924.sif` (sha256 `7bbd60fb…`). Hard links. `conf/cells.tsv` column 6 updated; backup `conf/cells.tsv.bak_20260924`. |
| archive | old-image outputs, traces, `rnaseq_data/` and `pasa_train_failed.tsv` moved to `archive/v1.9.0-rc1_container{,_rust}.sif20260921/` (43 GB + 46 GB). |
| re-run | jobs 29048317 (perl) and 29048318 (rust). All RNASEQ_PREPARE tasks submitted fresh (not cached). |
| scoring | BUSCO 6.0.0 (protein mode) on 7 cells (`launch/score_cell.sbatch`), then `compare_predictions.py` (`launch/compare_cells.sbatch`). No BUSCO failures. |

Why `rnaseq_data/` also had to move: `RNASEQ_PREPARE` uses `storeDir rnaseq_data/`. It skips Trinity when the output FASTA exists, independent of `-resume`. Without the move, the new image would have trained on Trinity output from the old image.

**Throughput limit.** At 02:35, 9 Trinity tasks ran and ~70 waited on `AssocGrpMemLimit` (epyc) / `AssocGrpCpuLimit` (preempt). Other jobs from the same account use that allocation. The re-run will not finish by morning.

## 2. Scope and completion

| cell | genomes with predict output | RNA-seq genomes with PASA training models |
|---|---|---|
| v1.8.17_conda | 15 (old 16-genome set; Rhodotorula only has a stale dir) | -- |
| v1.9.0-rc1_conda | 63 / 65 | 41 / 42 |
| v1.9.0-rc1_conda_rust | 64 / 65 | 37 / 42 |
| v1.9.0-rc1_container (old image) | 64 / 65 | 35 / 42 |
| v1.9.0-rc1_container_rust (old image) | 64 / 65 | 39 / 42 |

- Rhizopus microsporus ATCC 52813 fails in all four rc.1 cells: Trinity-GG produced 0 transcripts ("every partition failed"). Not investigated further.
- rc1_conda also lacks Morchella importuna.
- 65 genomes = 43 with RNA-seq, 22 without. Allomyces is GCA (no RefSeq) and is excluded from gffcompare tables.
- **The 1.8.17 comparisons have only 11-15 genomes.** 1.8.17 was never run on the 65-genome set.

## 3. Accuracy vs RefSeq

Paired per genome, Wilcoxon signed-rank. Restricted to genomes with PASA training models in BOTH cells (`results/accuracy_pairs_trained_both.pre_sif20260924.tsv`). "d" = median paired difference, B minus A, in percentage points. Counts are (B better / B worse).

### 3.1 1.8.17 vs 1.9.0-rc.1

| B vs A = 1.8.17_conda | n | exon Pr | intron Sn | intron-chain Sn | transcript Pr | locus Sn | BUSCO C |
|---|---|---|---|---|---|---|---|
| rc1_conda (perl) | 13 | +0.5 (13/0, p=2e-4) | +0.6 (12/1, p=1e-3) | +0.7 (13/0, p=2e-4) | +1.4 (13/0, p=2e-4) | +0.4 (p=0.01) | +0.1 (p=0.23) |
| rc1_conda_rust | 12 | +0.6 (p=0.02) | 0.0 (p=0.86) | 0.0 (p=0.91) | +1.2 (p=0.03) | +0.3 (p=0.27) | -0.4 (p=0.20) |
| rc1_container (old) | 11 | -1.6 (p=0.70) | +1.5 (p=0.15) | **-2.9 (1/10, p=0.04)** | -1.1 (p=0.76) | -2.9 (p=0.97) | -0.4 (p=0.46) |
| rc1_container_rust (old) | 13 | -1.4 (p=0.95) | +1.3 (p=0.41) | **-3.1 (1/12, p=0.008)** | -1.1 (p=1.0) | -3.5 (p=0.68) | -0.6 (p=0.15) |

- rc.1 conda perl is better than 1.8.17 on every structural metric. The gains are small (0.4-1.4 points) and consistent (12-13 of 13 genomes).
- rc.1 conda rust is equal or slightly better, and fewer metrics reach significance.
- The **old-image** rc.1 containers are ~3 points worse than 1.8.17 on intron-chain sensitivity.
- BUSCO completeness does not differ significantly between 1.8.17 and any rc.1 cell (n=15).

### 3.2 Old-image container vs conda, same rc.1 version

| A -> B | n | intron Pr | intron-chain Sn | transcript Sn | locus Sn | missed loci | BUSCO C |
|---|---|---|---|---|---|---|---|
| conda -> container (perl) | 33 | -3.5 (p=0.005) | -3.5 (3/30, p=2e-5) | -3.1 (p=0.004) | -3.1 (p=0.004) | -0.7 (p=0.16) | -0.7 (p=0.013) |
| conda_rust -> container_rust | 35 | -2.7 (p=0.016) | -3.1 (5/30, p=4e-5) | -3.0 (p=0.04) | -2.9 (p=0.04) | -0.8 (p=0.03) | -0.6 (p=0.04) |

- The old-image containers were worse than conda at the same version.
- The largest per-genome drops are in the Saccharomycetes. For K. lactis, S. cerevisiae, Z. rouxii and L. thermotolerans, intron-chain Sn fell from 55-61 to 11-17.
- In S. cerevisiae, 1.1% of old container_rust mRNAs are multi-CDS. The other cells range from 5.3% to 6.3% (1.8.17: 6.1%).
- In the 22 genomes without RNA-seq, conda and container agree to within 0.1 point. The difference is therefore in the RNA-seq/training path.
- **Cause not determined.** Tool versions differ between conda_rust and the SIF only in diamond (2.2.8 vs 2.2.2), hisat2 (2.2.3 vs 2.2.2) and salmon (1.10.3 vs 2.7.0). Trinity and PASA are source-built in both, and their exact commits were not compared. The new-image re-run is the direct test.

### 3.3 perl vs rust (EVM + PASA + Trinity rust helpers)

| pair | n | intron Sn | intron-chain Sn | transcript Sn | locus Sn | missed loci | BUSCO C |
|---|---|---|---|---|---|---|---|
| conda -> conda_rust | 35 | -0.4 (p=0.003) | -0.2 (p=0.03) | +0.1 (p=0.99) | +0.1 (p=0.88) | +0.2 (p=0.03) | -0.2 (p=0.025) |
| container -> container_rust (old) | 33 | -0.1 (p=0.015) | -0.1 (p=0.011) | -0.2 (p=2e-4) | -0.2 (p=2e-4) | +0.1 (p=1e-4) | -0.1 (p=0.004) |

- Rust is consistently a little worse, by 0.1-0.4 points. The difference is statistically detectable but small.

### 3.4 beta.12 -> rc.1 (PASA duplicate-output-loop fix)

On 15 genomes, beta12_container -> rc1_container (old image):
- BUSCO C: 85.1 -> 93.8 (p=0.023)
- missed loci: 16.2% -> 8.9%
- intron Sn: 39.9 -> 79.9

This confirms the PASA fix restored the collapsed gene structures.

## 4. Runtime

Source: `results/metrics_tasks.pre_sif20260924.tsv` (last-completed-attempt rows).

**Stage totals are not comparable across cells.** In rc1_conda (perl), genome-guided Trinity failed in RNASEQ_PREPARE for all 42 RNA-seq species. It failed within about 1 minute in the Trinity clustering step, in the `$SCRATCH` workdir. The Trinity log was on node scratch and is lost, so the root cause is unknown. TRAIN then ran Trinity again on `/bigdata` at 8 CPUs, and succeeded: A. fumigatus took 8.9 h. So only per-species end-to-end totals (RNASEQ_PREPARE + TRAIN + PREDICT) are compared (`results/timing_totals.pre_sif20260924.tsv`). Most old-image rows have no `cpus` value, so allocation parity is not verified.

| A -> B | species | median wall ratio B/A | p | median CPU-h ratio B/A | p | CPU-h sum A -> B |
|---|---|---|---|---|---|---|
| rc1_conda -> rc1_container (old) | 63 | 0.39 | 1e-7 | 0.70 | 0.001 | 1129 -> 740 |
| rc1_conda_rust -> rc1_container_rust (old) | 64 | 0.64 | 4e-7 | 0.93 | 0.02 | 795 -> 683 |
| rc1_conda -> rc1_conda_rust | 63 | 0.64 | 4e-4 | 0.67 | 1e-4 | 1129 -> 772 |
| rc1_container -> rc1_container_rust (old) | 64 | 0.88 | 0.002 | 0.91 | 0.002 | 761 -> 683 |
| 1.8.17_conda -> rc1_conda | 15 | 1.27 | 0.11 | 1.29 | 2e-4 | 179 -> 338 |
| 1.8.17_conda -> rc1_conda_rust | 15 | 0.83 | 0.85 | 0.74 | 0.49 | 179 -> 216 |
| 1.8.17_conda -> rc1_container (old) | 15 | 0.59 | 0.08 | 0.99 | 0.72 | 179 -> 196 |
| 1.8.17_conda -> rc1_container_rust (old) | 15 | 0.57 | 0.03 | 0.76 | 0.64 | 179 -> 158 |

- **Container vs conda.** The container is faster in wall time for both engines. For perl this is mostly the failed conda Trinity step (an env/pipeline defect, not provisioning). For rust, where both cells ran Trinity in RNASEQ_PREPARE, the container used 7% less CPU and 36% less wall time.
- **Rust vs perl, container.** Rust uses ~10% less end-to-end CPU and wall time. The difference is concentrated in PREDICT (EVM): wall 0.80, CPU 0.76.
- **1.8.17 vs rc.1.** Only container_rust is significantly faster in wall time (0.57, p=0.03, n=15). No rc.1 cell differs significantly in CPU-hours, except rc1_conda perl, which is slower because of the Trinity re-run above.

## 5. Holistic assessment (current data)

- **Accuracy.** 1.9.0-rc.1 in conda is slightly better than 1.8.17: +0.4 to +1.4 points on structure, BUSCO unchanged. This rests on 12-13 genomes. The previous rc.1 container image was ~3 points worse than 1.8.17 on intron-chain sensitivity. It must not be the basis of a "1.9.0 is better" claim until the new-image results are in.
- **Speed.** 1.9.0 is not slower in CPU-hours when Trinity runs correctly. The rust container is the fastest configuration in wall time. Rust saves ~10% end-to-end, mostly in EVM, at a cost of ~0.1-0.4 points of accuracy.
- **Robustness.** rc.1 completes 63-64 of 65 genomes. Open failures:
  - Rhizopus microsporus Trinity: all cells
  - conda-perl Trinity in RNASEQ_PREPARE: all species
  - 3-7 RNA-seq genomes per cell without PASA training models
- **Not yet supported by data.**
  - Any 1.8.17 vs 1.9.0 claim on the 65-genome set: 1.8.17 has 15 genomes.
  - Any claim about the new rc.1 image: re-run in progress.

## 6. New-image re-run: partial results (07:30-08:30)

State at 07:30:
- Trinity-GG done for 28 (perl) / 27 (rust) species.
- 0 PASA training models yet; 57 TRAIN tasks queued.
- 23 genomes predicted per cell, all no-RNA-seq genomes.
- No task errors.

Scored those 23 (`results/gene_content_comparison.rc1_container_new_partial.tsv`, `results/metrics_tasks.rc1_container_new_partial.tsv`):

- **Accuracy, new vs old image.** Median paired difference 0.00 on every gffcompare metric and on BUSCO, perl and rust. Expected: the old-image defect (section 3.2) is in the RNA-seq/training path, and these genomes do not use it. **The test that matters (RNA-seq genomes) is still pending.**
- **Rust vs perl, new image, PREDICT.** Wall ratio 0.85 (p=6e-5), CPU ratio 0.86 (p=7e-5). Same direction and size as the old image (0.80 / 0.76).
- **PREDICT time, new vs old image, same genomes.**
  - perl: 0.87 -> 1.25 h (p=1e-4)
  - rust: 0.85 -> 0.99 h (p=0.002)
  - All new rows ran at 8 cpus. The runs overlapped ~17 concurrent Trinity jobs plus other account workloads. Old-image allocation is unknown. **Do not attribute this slowdown to the image** until the complete run is compared under similar load.

## 7. Open items

1. Finish the new-image container re-run. Then run BUSCO + `compare_predictions.py` + `collect_metrics.py` on it and repeat sections 3.2, 3.3 and 4.
2. Decide whether to run v1.8.17_conda on the full 65-genome set, for a like-for-like version comparison.
3. Debug the conda-perl RNASEQ_PREPARE Trinity failure. Re-run one species with the scratch workdir kept.
4. Debug the Rhizopus microsporus Trinity-GG zero-transcript failure.
