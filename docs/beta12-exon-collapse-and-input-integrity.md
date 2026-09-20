---
topic: beta12-exon-collapse-and-input-integrity
description: OPEN. v1.9.0-beta.12 loses roughly one exon per gene on six of fifteen benchmark genomes versus v1.8.17, dropping mean gffcompare locus sensitivity from 55.9 to 48.1. Stage-splitting the training pipeline isolates the cause to the TransDecoder step (5.7.1 -> 6.0.0), not PASA and not the rust reimplementations, which are accuracy-neutral. Separately, the benchmark's SRA_QUERY cache holds transient network failures cached as successful empty results, and Rhodotorula toruloides' RNA-seq was replaced after mapping at 4.5%.
created: 2026-09-20
last_updated: 2026-09-20
status: OPEN -- TransDecoder identified as the stage; controlled same-input test not yet run
---

# v1.9.0-beta.12 exon-structure collapse, and RNA-seq input integrity

> Opened 2026-09-20 from BFD/Funannotate_benchmarking. Two independent threads,
> recorded together because they were found in the same pass and both bear on
> whether the benchmark's accuracy numbers can be published.
>
> Thread 1 (**blocking**): beta.12 predicts genes in the right places but with
> the wrong internal structure on a subset of genomes. Root cause isolated to
> the TransDecoder step.
>
> Thread 2 (**input integrity**): the benchmark's RNA-seq provenance is weaker
> than assumed, and one genome's reads were wrong outright.

---

## Thread 1 — beta.12 loses exon structure

### What the accuracy run showed

`scripts/compare_predictions.py` over 45 (cell, genome) pairs, scored against
RefSeq with gffcompare and bedtools:

| cell | mean locus sensitivity | mean locus precision |
|---|---|---|
| v1.8.17_conda | 55.9 | 59.1 |
| v1.9.0-beta12_container | 48.1 | 54.0 |
| v1.9.0-beta12_container_rust | 47.9 | 54.0 |

The result is bimodal, not a uniform shift. Five genomes drop by roughly half
(Botrytis 67.3 -> 29.2, Fusarium 65.4 -> 38.4, Aspergillus 52.5 -> 26.7,
Chaetomium 35.8 -> 21.0, Batrachochytrium 44.6 -> 24.9) while four improve
(Malassezia, Schizophyllum, Yarrowia, Zygosaccharomyces).

### It is not a scoring or format artifact

This was the first hypothesis and it is **wrong**. Both GFF3 files are
well-formed with identical feature types (`gene`/`mRNA`/`exon`/`CDS`/`tRNA`).
The bedtools coordinate-overlap tally even slightly favours beta.12 on Botrytis:
9,955/10,353 predicted genes (96.2%) overlap a RefSeq gene, versus
11,660/12,424 (93.8%) for 1.8.17.

The genes are in the right places. Their internal structure is wrong.

### The measurement that shows it: CDS per mRNA

| genome | RefSeq | 1.8.17 | beta12 | beta12-rust |
|---|---|---|---|---|
| Aspergillus fumigatus | 3.07 | 3.00 | **1.84** | 1.83 |
| Batrachochytrium dendrobatidis | 4.37 | 4.75 | **2.98** | 2.96 |
| Botrytis cinerea | 2.98 | 2.87 | **1.92** | 1.92 |
| Chaetomium globosum | 3.12 | 2.73 | **1.77** | 1.77 |
| Fusarium verticillioides | 2.75 | 2.74 | **1.97** | 1.97 |
| Malassezia globosa | 1.49 | 1.97 | **1.24** | 1.25 |
| Allomyces macrogynus | 3.36 | 3.57 | 3.91 | 3.93 |
| Cryptococcus neoformans | 6.27 | 6.07 | 6.49 | 6.48 |
| Pneumocystis murina | 6.00 | 7.17 | 7.18 | 7.18 |
| Schizophyllum commune | 5.16 | 5.64 | 6.03 | 6.04 |
| Tilletiopsis washingtonensis | 4.46 | 4.96 | 5.11 | 5.14 |
| Umbelopsis ramanniana | 4.64 | 5.26 | 5.41 | 5.41 |
| Saccharomyces cerevisiae | 1.06 | 1.06 | 1.01 | 1.01 |
| Yarrowia lipolytica | 1.17 | 1.28 | 1.12 | 1.13 |
| Zygosaccharomyces rouxii | 1.03 | 1.15 | 1.02 | 1.02 |

The six collapsed genomes are exactly the six with the large gffcompare drops.
Intron-level sensitivity tracks it (Botrytis 83.5 -> 36.0).

**The rust and norust arms are identical to two decimal places on every row.**
This is a 1.9.0 change present in both, not a rust effect.

### Isolating the stage

`pasa.step1.gff3` (PASA's raw alignment assemblies) and
`funannotate_train.pasa.gff3` (after TransDecoder selects training models) split
the training pipeline at the version boundary.

**Stage 1 -- PASA output, exons per assembly. Essentially unchanged:**

| genome | 1.8.17 | beta12 |
|---|---|---|
| Botrytis cinerea | 3.73 | 3.62 |
| Chaetomium globosum | 3.25 | 3.50 (higher) |
| Fusarium verticillioides | 3.46 | 3.40 |
| Aspergillus fumigatus | 4.61 | 3.90 |
| Cryptococcus neoformans | 7.76 | 7.31 |

**Stage 2 -- after TransDecoder, CDS per training model. Collapses:**

| genome | 1.8.17 | beta12 | change |
|---|---|---|---|
| Botrytis cinerea | 2.69 | 1.96 | -27% |
| Chaetomium globosum | 2.21 | 1.70 | -23% |
| Fusarium verticillioides | 2.62 | 1.94 | -26% |
| Aspergillus fumigatus | 2.70 | 1.84 | -32% |
| Cryptococcus neoformans | 5.72 | 3.27 | -43% |

**Fraction of exon structure retained across the TransDecoder step** (stage 2
divided by stage 1), which normalizes for the differing inputs:

| genome | 1.8.17 (TD 5.7.1) | beta12 (TD 6.0.0) |
|---|---|---|
| Botrytis cinerea | 72.1% | 54.2% |
| Chaetomium globosum | 68.0% | 48.6% |
| Fusarium verticillioides | 75.7% | 57.2% |
| Aspergillus fumigatus | 58.4% | 47.3% |
| Cryptococcus neoformans | 73.6% | 44.8% |
| Schizophyllum commune | 90.0% | 58.1% |
| **mean** | **73.0%** | **51.7%** |

Six of six genomes, no exceptions, mean drop of 21 points.

### Conclusion and what is still unproven

**PASA 2.5.3 -> 2.6.0 is exonerated**: its per-assembly exon structure is
preserved, and for Chaetomium is slightly better. **TransDecoder 5.7.1 -> 6.0.0
is selecting systematically shorter, fewer-exon ORFs from the same input.**

Versions:

| | 1.8.17 conda | beta12 container |
|---|---|---|
| TransDecoder | 5.7.1 (Jul 2023) | 6.0.0 (Mar 2026) |
| PASA | 2.5.3 | 2.6.0_rust |

TransDecoder 6.0.0's changelog records a restructure -- *"phase-specific
executables are now provided under `util/`"* -- which is the same change behind
the PATH probe fix in `modules/local/funannotate_train.nf`, confirming the
container runs the new layout.

**Not yet done:** the controlled same-input experiment -- feed one identical
PASA assembly set through TransDecoder 5.7.1 and 6.0.0 and compare CDS/model.
The retention ratio above already controls for input differences and the effect
is uniform across six genomes spanning three phyla, but the direct test is what
should back a published claim.

### A second, separate defect in the same data

beta.12 produces 30-40% fewer PASA assemblies from equal or slightly larger
Trinity inputs:

| genome | 1.8.17 assemblies | beta12 assemblies |
|---|---|---|
| Botrytis cinerea | 23,402 | 15,266 |
| Aspergillus fumigatus | 26,978 | 16,306 |
| Cryptococcus neoformans | 24,709 | 17,717 |
| Fusarium verticillioides | 28,386 | 20,213 |
| Chaetomium globosum | 27,190 | 24,049 |

This is a PASA-level loss, independent of the TransDecoder structure collapse.
**Not investigated.**

### Why only some genomes fail at prediction

Every genome checked is degraded at the *training* stage, including the two
whose final predictions are fine -- Cryptococcus loses 43% of its training-model
exon structure yet still predicts at 6.49 CDS/mRNA against a RefSeq 6.27.
Augustus/EVM recovers from degraded training on some genomes and not others.
The six visible failures therefore understate how widespread the underlying
problem is.

### Bearing on the manuscript

On present evidence the "is 1.9.0 better than 1.8.17" question answers **no** on
gene structure, but the cause looks like a dependency version bump rather than
funannotate's own code. If TransDecoder 6.0.0's ORF selection is confirmed as
the cause and either reverted or configured, the comparison must be re-run
before any accuracy claim. The rust question is unaffected: **the rust
reimplementations are accuracy-neutral** (mean locus sensitivity 47.9 vs 48.1,
and identical CDS/mRNA to two decimals on all 15 genomes).

---

## Thread 2 — RNA-seq input integrity

### SRA_QUERY caches transient network failures as success

`esearch` fails intermittently and **exits 0**:

```
curl: (56) OpenSSL SSL_read: ... unexpected eof while reading, errno 0
 ERROR:  curl command failed with: 56
rc=0
```

The empty `runinfo` flows into the awk filter in `modules/local/sra_query.nf`,
which writes no rows, leaving a header-only CSV. The module reports
`[INFO] Found 0 paired-end SRA accessions` and `storeDir` caches that
permanently. The module's `set -euo pipefail` does not catch it because the
failure is inside `timeout 300 bash -c "..."` and the pipeline exits 0.

This is the same cache-a-failure-as-success class as the known `SRA_FETCH` bug,
in a second module.

Demonstrated on *S. cerevisiae* (txid559292): BFD's cached CSV has 0
accessions, but three immediate retries all returned 251 runinfo lines and the
filter yields 5 accessions. Its 147 MB of real reads trained cleanly, so the
cached zero was never true.

### Empty cell-local query CSVs are intentional, but load-bearing

`scripts/preseed_rnaseq_reads.py` deliberately writes header-only CSVs as a
sentinel so `skip_sra_query` finds a cache hit. It populates rows from
`pool.rnaseq` = `Fungi_BFD_runs/samples.rnaseq_sra.csv`, **which contains 3 rows,
all for *Albifimbria verrucaria*, a species not in the benchmark.** Hence 15/16
empty in every beta11/beta12 cell -- expected output of that configuration, not
15 failed queries.

The hazard: in `subworkflows/local/fetch_rnaseq.nf:107-112` an empty CSV routes
to `no_data` -> `WRITE_EMPTY_READS`, which writes 0-byte FASTQs. That process has
`storeDir "${launchDir}/rnaseq_reads"`, so preseeded symlinks satisfy its
outputs and it is skipped. It works only by that coincidence. **Remove a symlink
and the species silently becomes ab-initio-only with 0-byte reads and no error.**

### The BFD query cache is not provenance for the BFD reads

BFD has 7,498 populated query CSVs, but they cannot be used to document which
accessions produced the reads:

- All were written `2026-07-16`; the reads beside them were downloaded May-June
  2026. The cache was regenerated weeks later, independently of the download.
- *S. cerevisiae* has 0 cached accessions beside 147 MB of real reads.
- `modules/local/sra_query.nf:56` sorts `-k2,2r` -- **most-recent ReleaseDate
  first**, then takes the top 5. The set is a moving target by construction.

Drift is observed, not merely structural. A live *S. cerevisiae* query returns
four of its top five released **2026-08-30** -- after the cache was written and
months after the reads were downloaded. (For *Chaetomium globosum*, by contrast,
a live query today returns the same five accessions as the cache: drift is
real but not universal.)

### Read status, 16-genome set

14 of 16 have real reads in the BFD store. *Pneumocystis murina* is 0 bytes with
0 accessions -- genuinely no usable RNA-seq, ab-initio only, legitimate.
*Rhodotorula toruloides* was 0 bytes despite 5 cached accessions.

### Recommendation

**Keep using the cached BFD reads; do not let the pipeline re-query SRA.** A
release-date-sorted top-5 means a re-run months from now can train on different
reads, which would make the 1.8.17-vs-1.9.0 comparison uncontrolled.

**But pin them explicitly rather than relying on empty stubs.** Two options,
undecided:

1. Publish the reads without per-accession provenance -- describe them as a
   frozen normalized read set and cite the BFD run. Honest, no re-download, but
   not reproducible by a reviewer.
2. Re-derive provenance by re-downloading from a pinned accession list. Gives a
   citable input set, but means re-fetching and re-training up to 14 genomes and
   changes the training data underneath results already collected.

A cheap partial check before deciding: for a few genomes, map a sample of each
cached read file against the candidate accessions' reads to test whether BFD's
listed accessions match the stored FASTQ. Not run.

### Rhodotorula toruloides: reads replaced

The previous accession set (PRJNA169538: SRR31947201, SRR31947212, SRR8733521,
SRR31947211, SRR31947205) mapped at **4.5%** against the NP11 reference
`GCF_000320785.1`. Candidates tested:

| dataset | mapping rate |
|---|---|
| previous set (PRJNA169538) | 4.5% |
| SRP655405 | 0.46% -- rejected |
| **SRP518254 (PRJNA1132283, strain NP11)** | **96.18%** (93.57% concordant-exactly-once) -- accepted |

The reference assembly was cleared of suspicion: strain-matched reads map at 96%.

Committed in `d8f5e1e`. `samples.csv` updated; stale Rhodotorula artifacts set
aside across 7 cells as `.stale_rnaseq_20260920`; the query CSV rewritten in all
11 cells with the two verified accessions so a clean-sweep run does not hit the
empty-query path.

### bbnorm pairing failure during the rebuild

Both accessions downloaded clean and mate-balanced (SRR29721378: 24,996,630
pairs; SRR29721379: 24,123,142 pairs; mate1 == mate2 within each). bbnorm then
died:

```
java.lang.AssertionError: List size mismatch: 6 vs 5
  at stream.PairStreamer.nextList(PairStreamer.java:88)
  at jgi.BBMerge.makeInsertHistogram(BBMerge.java:1339)
```

`PairStreamer` is BBTools' pairing layer, so R1 and R2 are out of step somewhere
past the head of the file. `parallel-fastq-dump` splits each accession into 15
blocks and reassembles them; that reassembly is the known desync source -- the
same one the pipeline's `rename_headers` blacklist override works around.

`launch/rebuild_rnaseq.sbatch` only counted combined R1, never R2. That gap is
fixed in `launch/repair_rnaseq_rhodotorula.sbatch`, which counts both mates,
runs `repair.sh` keeping singletons on the record, drops `ecc=t` (error
correction is what pulled BBMerge, and therefore PairStreamer, into the run),
and asserts R1 == R2 at the end.

---

## Open items

- [ ] Controlled TransDecoder 5.7.1 vs 6.0.0 test on one identical assembly set
- [ ] Investigate the 30-40% PASA assembly-count loss in beta.12
- [ ] Decide Thread 2 provenance option 1 vs 2
- [ ] Fix `sra_query.nf` to fail loudly on an empty `runinfo` instead of caching it
- [ ] Rhodotorula repair job 28958145 -- record retained-pair fraction, then re-run the cell
- [ ] `collect_metrics.py` trace-merge fix before any runtime figure
