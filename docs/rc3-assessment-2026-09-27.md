---
topic: rc3-assessment-2026-09-27
description: Accuracy (gffcompare vs RefSeq, BUSCO) and runtime of funannotate 1.9.0-rc.3 (container, rust EVM) against 1.9.0-rc.1 and 1.8.17 on the 65-genome benchmark set. Supersedes the interim rc1-assessment-2026-09-24.
created: 2026-09-27
status: final -- all 6 cells complete (65/65 genomes), gffcompare + BUSCO 6.0.0 + runtime
---

# funannotate 1.9.0-rc.3 assessment (2026-09-27)

## 1. Summary

- **Accuracy.** On the 39 comparable RNA-seq genomes, rc.3 is more accurate than 1.8.17 on 8 of 9 gffcompare metrics, by +0.1 to +1.5 points (p<0.05). Intron-chain sensitivity does not differ (p=0.21). On the 23 genomes without RNA-seq, rc.3 is better by +0.2 to +0.6 points (p<0.02).
- **rc.1 to rc.3.** The rc.3 PASA fixes change only the RNA-seq genomes. rc.3 predicts about 1.3% more genes than rc.1 relative to RefSeq and misses 0.4 points fewer loci (37 of 39 genomes, p=1e-7). Intron-chain sensitivity drops by 0.4 points (p=0.01). rc.1 predicted fewer genes than 1.8.17 and missed more loci; rc.3 removes both deficits.
- **Gold-standard genomes.** rc.3 is the best of the three versions on nearly every metric for S. cerevisiae, A. fumigatus, B. cinerea and Z. tritici. S. cerevisiae gains the most: intron-chain Sn 53.9 (1.8.17) to 63.3 (rc.3).
- **BUSCO.** Completeness does not separate the versions on RNA-seq genomes (rc.3 vs 1.8.17: -0.5 points, p=0.15). On genomes without RNA-seq, rc.3 is +0.2 points over 1.8.17 (15/4, p=0.003).
- **Runtime.** 1.9 with the rust EVM in a container uses 0.51x the wall time and 0.75x the CPU-hours of 1.8.17 end to end (62 species). The rc.3 fixes do not change TRAIN or PREDICT time relative to rc.1.
- **RNA-seq read check.** rc.3 rejects RNA-seq reads that do not map to the genome. It rejected the P. teres 0-1 reads (2.73% mapped). Ab initio rc.3 then missed 13.6% of RefSeq loci, against 36-40% for 1.8.17 and rc.1, which trained on those reads.

## 2. Cells and inputs

| cell | version | provisioning | EVM | genomes predicted | RNA-seq genomes with PASA models |
|---|---|---|---|---|---|
| v1.8.17_conda | 1.8.17 | conda | perl | 65 | 42 / 42 |
| v1.9.0-rc1_conda | 1.9.0-rc.1 | conda | perl | 65 | 41 / 42 |
| v1.9.0-rc1_conda_rust | 1.9.0-rc.1 | conda | rust | 65 | 38 / 42 |
| v1.9.0-rc1_container | 1.9.0-rc.1 | container | perl | 65 | 41 / 42 |
| v1.9.0-rc1_container_rust | 1.9.0-rc.1 | container | rust | 65 | 41 / 42 |
| v1.9.0-rc3_container_rust | 1.9.0-rc.3 | container | rust | 65 | 40 / 42 |

- 65 genomes: 42 with RNA-seq, 23 without. Allomyces macrogynus is GCA (no RefSeq annotation) and is left out of gffcompare tables.
- All cells use the same masked genomes, the same pinned RNA-seq reads (B. dendrobatidis: JEL423-only set) and the same GeneMark sidecar models.
- rc.3 image: `funannotate-1.9.0-rc.3.20260926.sif` (sha256 `851732c1...`, dated hard link). It reports `funannotate v1.9.0-rc.3`, `FUNANNOTATE_EVM_ENGINE=rust`.
- **rc.3 reused Trinity.** Its `rnaseq_data/` holds copies of the v1.9.0-rc1_container_rust Trinity-GG, StringTie and junction files (149 files, md5-verified). All 42 RNASEQ_PREPARE tasks were skipped. rc.3 therefore tests only the TRAIN (PASA) and PREDICT changes on identical transcript evidence. rc.3 has no Trinity runtime of its own.
- No rc.3 perl or conda cell was run. rc.3 changes PASA, which affects perl and rust the same way. The rc.1 perl-vs-rust and conda-vs-container results stand (section 5).

## 3. rc.3 run issues and workarounds

1. **BUSCO_COMPLETENESS.** This nf_funannotate1 step started running on 2026-09-26 (commit 71a24ce). It passed `-l <busco_lineages>/<lineage>`, which skips the `lineages/` level. Every call failed. After 3 attempts the global `errorStrategy 'finish'` stopped the run from submitting new tasks. The relaunch (job 29117094) passed `--busco_lineages /srv/projects/db/BUSCO/v10/lineages`. Fixed in nf_funannotate1 branch `fix/rc3-rnaseq-gate-busco-lineages` (commit 9547821, not merged).
2. **RNA-seq concordance gate.** rc.3 train samples 200,000 reads and stops with exit 3 if fewer than 10% map. P. teres 0-1 had 2.73%. nf_funannotate1 treated exit 3 as a hard failure. The `training/.pasa_train_failed` marker was written by hand for P. teres, so predict ran ab initio. Same branch fixes the module (exit 3 plus the gate message now degrades to ab initio).
3. **BUSCO version.** The pipeline step scored 41 rc.3 genomes with BUSCO 6.1.0. All other cells use BUSCO 6.0.0 (`launch/score_cell.sbatch`). The 6.1.0 results were moved to `runs/v1.9.0-rc3_container_rust/busco_completeness_pipeline_6.1.0/`, and those genomes were rescored with 6.0.0.

## 4. Accuracy vs RefSeq

gffcompare against the NCBI RefSeq annotation. Paired per genome, Wilcoxon signed-rank. "d" = median paired difference, B minus A, in percentage points. Counts are (B better / B worse). Genomes are kept only if both cells have PASA models (drops Drepanopeziza brunnea, which has no rc.3 or rc.1 PASA models). P. teres 0-1 is left out of the RNA-seq tables (section 4.4).

### 4.1 RNA-seq genomes (n=39)

| metric | rc.1 -> rc.3 d (B+/B-, p) | 1.8.17 -> rc.3 d (B+/B-, p) | 1.8.17 -> rc.1 d (B+/B-, p) |
|---|---|---|---|
| gene count / RefSeq | +0.013 (38/1, 2e-11) | -0.003 (18/20, 0.70) | -0.016 (4/35, 2e-7) |
| exon Pr | +0.1 (24/12, 0.03) | +0.7 (29/9, 1e-4) | +0.4 (30/5, 9e-4) |
| intron Sn | +0.1 (23/13, 0.14) | +0.2 (23/15, 0.05) | 0.0 (19/16, 0.51) |
| intron Pr | +0.4 (28/9, 6e-4) | +0.5 (30/9, 5e-4) | +0.1 (21/13, 0.06) |
| intron-chain Sn | -0.4 (10/26, 0.01) | -0.4 (14/24, 0.21) | 0.0 (19/17, 0.46) |
| intron-chain Pr | +0.6 (30/9, 3e-4) | +1.0 (30/9, 5e-4) | +0.1 (25/9, 0.03) |
| transcript Sn | +0.9 (29/9, 7e-4) | +1.0 (28/11, 7e-4) | +0.1 (22/16, 0.57) |
| transcript Pr | +0.3 (23/13, 0.19) | +1.5 (30/8, 4e-4) | +1.0 (33/6, 3e-4) |
| locus Sn | +0.9 (29/9, 6e-4) | +1.0 (28/10, 9e-4) | +0.1 (22/16, 0.59) |
| missed loci % (lower is better) | -0.4 (37 fewer / 1 more, 1e-7) | -0.1 (22 fewer / 9 more, 0.02) | +0.3 (5 fewer / 32 more, 4e-5) |

### 4.2 Genomes without RNA-seq (n=23)

- rc.1 -> rc.3: no metric differs (all p>0.08). Expected: PASA is not used.
- 1.8.17 -> rc.3: +0.2 to +0.6 points on exon Pr, intron Sn/Pr, intron-chain Sn/Pr, transcript Sn/Pr and locus Sn (14-20 of 23 better, p<0.02). Missed loci do not differ (p=0.24).

### 4.3 Gold-standard genomes

| genome | cell | genes (RefSeq) | intron Sn / Pr | intron-chain Sn / Pr | transcript Sn / Pr | missed loci % |
|---|---|---|---|---|---|---|
| S. cerevisiae S288C | 1.8.17 | 5687 (6459) | 53.5 / 58.4 | 53.9 / 56.6 | 83.4 / 95.0 | 10.6 |
| | rc.1 | 5680 | 57.4 / 70.2 | 57.1 / 69.3 | 84.0 / 95.8 | 10.9 |
| | rc.3 | 5736 | 63.3 / 73.6 | 63.3 / 73.1 | 85.2 / 96.2 | 10.1 |
| A. fumigatus Af293 | 1.8.17 | 10118 (10085) | 76.3 / 83.0 | 53.7 / 55.6 | 52.4 / 52.2 | 7.0 |
| | rc.1 | 9896 | 76.4 / 82.7 | 53.9 / 55.7 | 52.5 / 53.5 | 7.4 |
| | rc.3 | 10023 | 77.3 / 83.4 | 54.0 / 57.5 | 54.3 / 54.6 | 7.1 |
| B. cinerea B05.10 | 1.8.17 | 12424 (11698) | 83.5 / 93.2 | 65.0 / 76.5 | 57.5 / 63.4 | 3.2 |
| | rc.1 | 12139 | 83.3 / 93.2 | 64.5 / 76.5 | 57.3 / 64.6 | 3.3 |
| | rc.3 | 12284 | 84.1 / 93.9 | 64.5 / 77.9 | 57.6 / 64.1 | 3.1 |
| Z. tritici IPO323 | 1.8.17 | 13101 (10941) | 62.1 / 67.0 | 45.6 / 41.1 | 46.2 / 38.6 | 7.6 |
| | rc.1 | 12484 | 62.4 / 67.0 | 46.4 / 41.8 | 46.8 / 41.0 | 7.9 |
| | rc.3 | 12885 | 62.2 / 68.7 | 45.8 / 44.0 | 48.9 / 41.5 | 7.3 |

rc.1 and rc.3 rows are container_rust. All four genomes have RNA-seq and PASA models in all three cells. Exon-level Sn/Pr is not shown: gffcompare requires exact exon ends, so UTR differences dominate it.

### 4.4 Pyrenophora teres 0-1 (rejected RNA-seq)

| cell | genes (RefSeq 11799) | intron-chain Sn / Pr | transcript Sn / Pr | missed loci % |
|---|---|---|---|---|
| 1.8.17 (trained on the reads) | 7898 | 41.7 / 51.1 | 36.0 / 53.8 | 36.1 |
| rc.1 container_rust (trained on the reads) | 7485 | 38.7 / 49.2 | 32.8 / 51.6 | 39.5 |
| rc.3 container_rust (reads rejected, ab initio) | 10437 | 70.6 / 73.6 | 67.4 / 76.2 | 13.6 |

- Training on off-target reads lost about one third of the RefSeq genes.
- This is one genome, not a statistical result. P. teres is excluded from the paired RNA-seq tests because the reads, not the software, drive its result.
- The other 41 RNA-seq genomes passed the rc.3 gate (>=10% mapped).

## 5. Runtime

Source: `results/metrics_tasks.runs_20260927.tsv` (last completed attempt per task).

### 5.1 End to end, 1.8.17 vs 1.9 (Trinity included)

Per-species sum of RNASEQ_PREPARE + TRAIN + PREDICT (`scripts/timing_totals.py`). rc.1 cells are used because rc.3 has no Trinity time of its own.

| A -> B | species | median wall ratio B/A (p) | median CPU-h ratio B/A (p) | CPU-h sum A -> B |
|---|---|---|---|---|
| 1.8.17 -> rc1_container_rust | 62 | 0.51 (4e-9) | 0.75 (2e-4) | 1154 -> 866 |
| 1.8.17 -> rc1_container | 62 | 0.60 (2e-7) | 0.86 (0.01) | 1154 -> 952 |
| 1.8.17 -> rc1_conda_rust | 62 | 0.65 (2e-4) | 0.74 (4e-6) | 1154 -> 807 |
| 1.8.17 -> rc1_conda (perl) | 42 | 1.10 (0.30) | 1.04 (0.39) | 1003 -> 1057 |
| rc1_container -> rc1_container_rust | 65 | 0.86 (3e-10) | 0.88 (1e-8) | 971 -> 881 |
| rc1_conda_rust -> rc1_container_rust | 65 | 0.76 (5e-4) | 1.06 (0.05) | 820 -> 881 |

- rc1_conda (perl) is not faster than 1.8.17. For 13 of its 42 RNA-seq species, Trinity failed in RNASEQ_PREPARE and ran again inside TRAIN at 8 CPUs. This is an env/pipeline defect, not a property of 1.9.
- timing_totals keeps a species only when both cells have the same stage profile; this drops 3 species from the 1.8.17 pairs and all no-RNA-seq species from the rc1_conda pair.

### 5.2 rc.1 -> rc.3 (container_rust), per stage

| stage | n | median wall ratio (p) | median CPU-h ratio (p) |
|---|---|---|---|
| FUNANNOTATE_TRAIN | 42 | 0.95 (0.46) | 0.97 (0.14) |
| FUNANNOTATE_PREDICT | 65 | 0.97 (0.11) | 0.94 (0.004) |

The rc.3 fixes do not change runtime, apart from 6% less PREDICT CPU time.

### 5.3 Caveats

- Traces do not record hostnames. From SLURM accounting, the TRAIN and PREDICT jobs that could be matched ran almost all on `r` (epyc) nodes; 1 of 30 matched 1.8.17 TRAIN jobs ran on a `c` node. Coverage is partial (8 of 42 rc1_container_rust TRAIN jobs matched), so node type is not fully excluded as a confound.
- nf_funannotate1 has kept RNASEQ_PREPARE and TRAIN off nodes c01-c30 since 2026-09-25 (commits bdee611, 390e3c8). rc.1 ran before that.
- Runs overlapped other group workloads; queue wait is not included in any number.

## 6. What the data supports

- 1.9.0-rc.3 is at least as accurate as 1.8.17 on every gffcompare metric except intron-chain sensitivity, where it does not differ significantly. The gains are small (0.2-1.5 points) and consistent across genomes.
- 1.9.0 with the rust EVM in a container halves end-to-end wall time and uses about 25% less CPU than 1.8.17.
- The rc.3 PASA fixes recover the gene-count and missed-loci deficit that rc.1 had against 1.8.17.
- rc.3's RNA-seq gate prevents a large accuracy loss when the reads are off-target (one case in this set).

**Not supported by this data:**
- Any claim about rc.3 conda or rc.3 perl EVM. Not run.
- Any rc.3 end-to-end runtime that includes Trinity. Trinity was reused from rc.1.
- A general claim about the RNA-seq gate from a single genome.

## 7. BUSCO completeness

BUSCO 6.0.0, protein mode, odb10 lineage per genome from `samples.csv` (`launch/score_cell.sbatch`). All 6 x 65 results are 6.0.0. Paired Wilcoxon on complete % (C); same genome filter as section 4, but Allomyces (no RefSeq) is included, so the RNA-seq n is 40. Table: `results/busco_pairs.runs_20260927.tsv`.

| A -> B | RNA-seq: n, median C A -> B, d (B+/B-, p) | no RNA-seq: n, median C A -> B, d (B+/B-, p) |
|---|---|---|
| rc1_container_rust -> rc3_container_rust | 40, 94.8 -> 94.5, -0.10 (13/24, 0.25) | 23, 93.1 -> 93.0, 0.00 (4/6, 0.96) |
| 1.8.17 -> rc3_container_rust | 40, 95.0 -> 94.5, -0.50 (13/27, 0.15) | 23, 92.6 -> 93.0, +0.20 (15/4, 0.003) |
| 1.8.17 -> rc1_container_rust | 40, 95.0 -> 94.8, -0.25 (10/25, 0.046) | 23, 92.6 -> 93.1, +0.20 (17/3, 0.002) |
| 1.8.17 -> rc1_conda | 40, 95.0 -> 95.0, 0.00 (16/19, 0.80) | 23, 92.6 -> 93.1, +0.20 (19/1, 1e-4) |
| rc1_conda -> rc1_container | 40, 95.0 -> 95.0, 0.00 (15/16, 0.52) | 23, 93.1 -> 93.1, 0.00 (6/9, 0.51) |
| rc1_conda_rust -> rc1_container_rust | 37, 94.9 -> 94.9, 0.00 (12/13, 0.11) | 23, 93.0 -> 93.1, 0.00 (5/6, 0.93) |
| rc1_container -> rc1_container_rust | 40, 95.0 -> 94.8, -0.10 (6/22, 0.003) | 23, 93.1 -> 93.1, 0.00 (6/10, 0.12) |
| rc1_conda -> rc1_conda_rust | 37, 95.1 -> 94.9, 0.00 (10/17, 0.02) | 23, 93.1 -> 93.0, 0.00 (5/10, 0.06) |

- BUSCO medians are 93-95% in every cell; differences are at most 0.5 points.
- Rust vs perl costs about 0.1 points of BUSCO C on RNA-seq genomes (22 of 40 worse, p=0.003, container).
- New-image rc.1 container and conda agree (p=0.52). The old-image container deficit in the 09-24 report is gone.

Gold-standard genomes and P. teres, BUSCO C / M (%):

| genome | 1.8.17 | rc.1 container_rust | rc.3 container_rust |
|---|---|---|---|
| S. cerevisiae S288C | 97.5 / 1.9 | 97.6 / 1.7 | 98.0 / 1.3 |
| A. fumigatus Af293 | 95.9 / 1.9 | 95.7 / 2.2 | 96.2 / 2.1 |
| B. cinerea B05.10 | 98.7 / 0.5 | 98.2 / 0.6 | 98.2 / 0.6 |
| Z. tritici IPO323 | 96.6 / 1.9 | 97.2 / 1.4 | 97.0 / 1.4 |
| P. teres 0-1 | 72.1 / 18.0 | 67.2 / 22.1 | 96.7 / 1.6 |

## 8. Files

| file | contents |
|---|---|
| `results/gene_content_comparison.runs_20260927.tsv` | gffcompare + BUSCO, all 6 cells x 65 genomes |
| `results/accuracy_pairs_rnaseq.runs_20260927.tsv` | section 4.1 (RNA-seq genomes, P. teres excluded) plus rc.1 conda/container and perl/rust pairs |
| `results/accuracy_pairs_nornaseq.runs_20260927.tsv` | section 4.2, same pairs |
| `results/busco_pairs.runs_20260927.tsv` | section 7 |
| `results/gold_standard_accuracy.runs_20260927.tsv` | sections 4.3, 4.4 and 7, all 6 cells |
| `results/metrics*.runs_20260927.tsv` | runtime, all 6 cells + GeneMark sidecar |
| `results/*.retired_cells_20260925.tsv` | gffcompare + runtime for the deleted beta10/11/12/13, 1.8.17_container and pasaiso cells |
| `scripts/filter_trained_rows.py` | drops RNA-seq genomes without PASA models before pairing |
| `scripts/cleanup_runs.sh` | tiered deletion of run data (A: retired cells/archive; B: work/; C: training intermediates) |
