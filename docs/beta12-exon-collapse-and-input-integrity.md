---
topic: beta12-exon-collapse-and-input-integrity
description: RESOLVED (Thread 1). v1.9.0-beta.12/beta.13's exon-structure collapse was NOT TransDecoder (that attribution was wrong and is corrected below) -- it traced to a real duplicate-print-loop bug in PASA's rust_optimize branch (PASA_transcripts_and_assemblies_to_GFF3.dbi, introduced by commit bce776a), which doubled every validated-alignment GFF3/GTF line and fed a degraded, less-multiexon transcript pool into TransDecoder/getBestModel. Fixed at PASApipeline rust_optimize@4376a22 (tagged v2.6.1-rc.1); verified on a real pipeline run to restore PASA 2.5.3-matching structure. Three compounding Python-side bugs in funannotate-live (getBestModel tie-break, selectTrainingModels 200-gate, HiQ 0/0 conflation) were also found and fixed independently. Thread 2 (RNA-seq input integrity) is unchanged and still open.
created: 2026-09-20
last_updated: 2026-09-21
status: Thread 1 RESOLVED (fix verified); Thread 2 still OPEN
---

# v1.9.0-beta.12/beta.13 exon-structure collapse, and RNA-seq input integrity

> Opened 2026-09-20 from BFD/Funannotate_benchmarking. Two independent threads,
> recorded together because they were found in the same pass and both bear on
> whether the benchmark's accuracy numbers can be published.
>
> Thread 1 (**was blocking, now resolved**): beta.12/beta.13 predicted genes in
> the right places but with the wrong internal structure on a subset of
> genomes. Root cause found and fixed 2026-09-21.
>
> Thread 2 (**input integrity, still open**): the benchmark's RNA-seq
> provenance is weaker than assumed, and one genome's reads were wrong
> outright.

---

## Thread 1 — exon-structure collapse: root cause and fix (2026-09-21 update)

### Correction: the TransDecoder attribution below was wrong

The original version of this doc (2026-09-20) concluded the collapse traced to
TransDecoder 5.7.1 -> 6.0.0, based on comparing `pasa.step1.gff3` (believed
pre-TransDecoder) against `funannotate_train.pasa.gff3` (post-TransDecoder).
That reasoning was wrong on a factual point: **`pasa.step1.gff3` is itself
POST-TransDecoder output** (`funannotate/train.py` writes it after PASA's own
`--TRANSDECODER` pass, not before). Re-reading `train.py` caught this.

A controlled test settled it: TransDecoder 5.7.1 and 6.0.0 run on byte-identical
input Trinity transcripts produced **byte-identical output** (job 28958410).
TransDecoder is exonerated. The "TransDecoder retention ratio" table further
down in this doc measured something real (a structural drop between two
already-different stages) but attributed it to the wrong stage boundary.

### The real root cause: a duplicate-print-loop bug in PASA itself

Isolating PASA version as a variable (genome, Trinity assembly, and
`--aligners` all held identical via md5-verified inputs) showed PASA 2.6.0_rust
producing a genuinely less-multiexon transcript pool than PASA 2.5.3 at the
`pasa.step1.gff3` stage — this was real, just mis-attributed to TransDecoder
rather than PASA itself:

| | PASA 2.5.3 (1.8.17) | PASA 2.6.0_rust (broken beta.13) |
|---|---|---|
| `step1.gff3` mean CDS/model | 2.51 | 1.87 |
| `step1.gff3` % >=3 exon | 39.3% | 21.1% |

Tracing PASA's own intermediate files found every validated-alignment record
in `valid_gmap_alignments.gff3` / `valid_custom_alignments.gff3` printed
**exactly twice** (55,205 unique lines -> 110,410 total, verified via
exact-duplicate-line counting). This was not a database, threading, or
DB-backend effect (all ruled out empirically: forcing `-T 1` on the clustering
script did not fix it; SQLite backend reproduced the identical 2x duplication;
no retry/reconnect events fired in any log).

`git log`/`git blame` on the `rust_optimize` branch (hyphaltip/PASApipeline)
found the actual cause: commit `bce776a` ("parallelize per-asmbl_id scripts,
fix N+1 query in GFF3 output", 2026-06-29, self-authored) added a new
batch-query print loop to `scripts/PASA_transcripts_and_assemblies_to_GFF3.dbi`
to replace an old per-alignment N+1-query loop, but never deleted the old loop.
Both ran unconditionally for GFF3/GTF output, so every validated alignment
segment was printed twice. BED output was unaffected (both loops wrote into
the same `%ALIGNMENTS` hash key, an overwrite not an append), which is why the
bug was invisible until GFF3 structure was compared directly against PASA
2.5.3 under matched inputs. **This is not an upstream PASA bug** — the file is
byte-identical to `master` aside from this branch's own commits; nothing needed
reporting upstream.

A second, related defect was found in the same pass: `Launch_PASA_pipeline.pl`
had a copy-pasted duplicate `@cmds` block queuing 5 of 6 output-writing
commands a second time, and a checkpoint-filename typo (`chkpt =>
"...failed_${map_program}_alignments.gff3.ok"` on the failed-alignments `.bed`
writer, colliding with the `.gff3` writer's checkpoint) that silently skipped
the `.bed` writer on every run.

### Fix

Both bugs fixed in `PASApipeline` (`rust_optimize` branch), commit `4376a22`,
tagged `v2.6.1-rc.1`:

1. `PASA_transcripts_and_assemblies_to_GFF3.dbi`: deleted the legacy
   per-alignment loop, kept only the batch-query loop.
2. `Launch_PASA_pipeline.pl`: deleted the duplicate `@cmds` block; fixed the
   checkpoint typo.

Independently reviewed end-to-end by a second model pass (control flow, other
callers of the GFF3 script, checkpoint-name collisions across the whole file)
— PASS, no regressions.

A separate, unrelated build issue was found and fixed while rebuilding: the
`rust_optimize` branch's `pasa_rust/Cargo.toml` had `edition = "2026"` (not a
real Rust edition), breaking `cargo build --release` outright on current
toolchains. Fixed to `edition = "2021"`. `v2.6.1-rc.1` was retagged to include
this fix (`2a4aeae`).

### Verification: fix restores PASA 2.5.3-matching structure

Rebuilt the conda-rust PASA install from the fixed source (Perl scripts
verified byte-identical to the patched checkout; rust binaries freshly
compiled from the corrected `Cargo.toml`) and reran the isolation cell with the
same genome/Trinity/aligner inputs used throughout this investigation:

| | PASA 2.5.3<br>(1.8.17) | PASA 2.6.0_rust<br>(broken, beta.13) | PASA 2.6.1-rc.1<br>(fixed) |
|---|---|---|---|
| `step1.gff3` mean CDS/model | 2.51 | 1.87 | **2.51** |
| `step1.gff3` % >=2 exon | 66.0% | 52.1% | **66.0%** |
| `step1.gff3` % >=3 exon | 39.3% | 21.1% | **39.3%** |
| final training models | 8,633 | 7,600 | **8,729** |
| final mean CDS/model | 2.55 | 2.02 | **2.55** |
| final % >=2 exon | 66.6% | 60.5% | **66.7%** |
| final % >=3 exon | 41.2% | 25.2% | **41.2%** |

`valid_gmap_alignments.gff3` confirmed 55,205 total = 55,205 unique lines (no
duplication) on the fixed build, vs 110,410/55,205 (2x) on the broken build.

Also resolves the "second, separate defect" flagged in the original version of
this doc (beta.12 producing 30-40% fewer PASA assemblies than 1.8.17 from
equal/larger Trinity input, e.g. Botrytis 23,402 vs 15,266) — that was very
likely the same root cause (the corrupted, duplicate-inflated alignment pool
feeding the assembler), not a separate defect. Not re-verified across all
genomes yet (see Open items).

### Compounding, independently-found bugs in funannotate-live itself

Found and fixed in the same investigation, upstream of/parallel to the PASA
bug, each real on its own:

- `getBestModel()` in `train.py`: TPM tie-break blind to CDS structure — when
  TransDecoder emits multiple ORFs per PASA assembly sharing identical exon
  sets (same TPM), the winner was decided by Python's stable-sort iteration
  order, not model quality. Fixed to break ties on `(nCDS, cdsLen)`.
- `selectTrainingModels()` in `library.py`: the `countKeeperCDS >= min_models`
  gate, when the filtered set was large enough but had too few multi-CDS
  genes, dropped the multi-CDS requirement entirely rather than falling back
  to the unfiltered set. Fixed.
- HiQ Augustus model selection in `predict.py`: an intronless model
  (`# CDS introns: 0/0`) triggered a `ZeroDivisionError` caught and scored
  `support=0`, making HiQ structurally impossible for any intronless gene when
  `--rna_bam` was used. Fixed to fall back to hint-support percentage.

These are real, independent fixes (all committed to `funannotate-live`
`local/beta.13`, fast-forward-merged into `target_1.9/rust_EVM_trinity_PASA`)
and should NOT be un-fixed now that the PASA bug is also found — both classes
of defect were compounding on the same collapsed genomes.

### Ruled out during the investigation (do not re-litigate)

- **Aligner choice** (gmap vs blat vs minimap2): no effect on multi-exon
  structure once PASA version and Trinity input are held constant
  (md5-verified controlled test).
- **PASA clustering-script threading** (`-T` thread count on
  `assign_clusters_by_stringent_alignment_overlap.dbi`): forcing `-T 1` did
  not eliminate the duplication; this script runs *after* the file that
  showed duplication is already written, so it was structurally impossible
  for it to be the cause anyway.
- **DB backend** (MySQL vs SQLite): both reproduced the identical 2x
  duplication before the fix, confirming it is a pure output-loop bug, not a
  database/connection-retry effect.
- **`DB_connect.pm`'s SQLite-busy-retry / reconnect logic**: confirmed dead
  code for MySQL-backend runs (zero retry/reconnect events in any log).
- **The rust-vs-C++ PASA assembler question**: `PASA_alignment_assembler.pm`
  defaults to the C++ `pasa` binary (Rust only used if `PASA_ASSEMBLER` is set
  or no C++ binary exists), confirmed via `PerlLib/CDNA/PASA_alignment_assembler.pm`
  and `docs/Cpp_PASA_Assembler_Optimization.md` in the PASApipeline repo. Our
  pipeline never sets `PASA_ASSEMBLER`, so the core assembly algorithm has been
  C++ in every "_rust" cell all along — irrelevant to this bug (which is a
  pure Perl output-loop issue) but relevant context for any future "is rust
  faster/more-accurate" claim about the assembler specifically. `slclust`
  (clustering) genuinely defaults to Rust; `pasa` (assembly) does not. Both are
  part of the `rust_EVM_trinity_PASA` effort along with Trinity/EVM Rust work,
  so the "_rust" cell label is not invalidated by this — see the parent branch
  name for the intended scope.

### Bearing on the manuscript

Original conclusion ("1.9.0 loses gene structure, cause looks like a dependency
version bump") is superseded. The cause was a self-inflicted bug in this
project's own PASA fork, now fixed, now verified to restore PASA
2.5.3-matching structure on the one genome tested end-to-end (Botrytis). A full
re-run of the accuracy comparison against the fixed build (`v1.9.0-rc.1`,
tagged in `cells.tsv` as `v1.9.0-rc1_container` / `v1.9.0-rc1_container_rust`)
is needed before any "is 1.9.0 as accurate as 1.8.17" claim — see Open items.

---

## Thread 2 — RNA-seq input integrity

*(Unchanged from 2026-09-20 — still open. See prior content below.)*

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

### Rhodotorula toruloides: reads replaced, then read-length-mismatch fixed

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

**Update 2026-09-21**: the accepted accession (SRR29721378) then failed
`bbnorm` deterministically on any chunk of the data with
`AssertionError: List size mismatch`. Root-caused (after ruling out
multi-member gzip and R1/R2 stream desync) to genuine per-mate read-length
asymmetry (2.3% of pairs, e.g. R1=150bp/R2=89bp) — the accession was deposited
already asymmetrically pre-trimmed per mate, not raw uniform-150bp sequencer
output. Fixed with `fix_fastq_header_trinity` (strip embedded `length=NNN`
deflines) run *before* `enforce_seqpair_readlen` (truncate longer mate to
match shorter, drop pairs below `minlen=75`) — order matters, since the
truncation tool compares full FASTQ headers including the length-dependent
embedded description text. Verified on a 20,000-pair slice (0 mismatches),
then run on the full 24,996,630-pair accession: job 28962314, completed,
2,307,195 final normalized pairs. **Not yet wired back into any benchmark
cell's `rnaseq_data/` for a re-run of FUNANNOTATE_TRAIN/PREDICT** — see Open
items.

### bbnorm pairing failure during the rebuild (superseded, see above)

Both accessions downloaded clean and mate-balanced (SRR29721378: 24,996,630
pairs; SRR29721379: 24,123,142 pairs; mate1 == mate2 within each). bbnorm then
died:

```
java.lang.AssertionError: List size mismatch: 6 vs 5
  at stream.PairStreamer.nextList(PairStreamer.java:88)
  at jgi.BBMerge.makeInsertHistogram(BBMerge.java:1339)
```

This was the SAME symptom, root-caused differently at the time (blamed on
`parallel-fastq-dump`'s 15-way block reassembly). The read-length-asymmetry
root cause above (found 2026-09-21) is the correct explanation; the
`rename_headers`-blacklist / `repair.sh` workaround described in the original
version of this section did not actually fix it (confirmed: it failed again
under every prior hypothesis before the header-order fix worked).

---

## Open items

- [x] ~~Controlled TransDecoder 5.7.1 vs 6.0.0 test on one identical assembly set~~ — done; TransDecoder exonerated, real cause found and fixed (see Thread 1)
- [ ] Re-verify the "30-40% fewer PASA assemblies" defect is resolved by the same fix, across all affected genomes (only Botrytis confirmed so far)
- [ ] Re-run the full accuracy comparison (gffcompare Sn/Sp, CDS/mRNA) against `v1.9.0-rc1_container` / `v1.9.0-rc1_container_rust` before any "is 1.9.0 as accurate as 1.8.17" claim
- [ ] Decide Thread 2 provenance option 1 vs 2
- [ ] Fix `sra_query.nf` to fail loudly on an empty `runinfo` instead of caching it
- [ ] Wire the fixed Rhodotorula toruloides RNA-seq (2,307,195 pairs, job 28962314) into the actual benchmark cells and re-run FUNANNOTATE_TRAIN/PREDICT
- [ ] `collect_metrics.py` trace-merge fix before any runtime figure
- [ ] Fix `scripts/preseed_clean_genomes.py`'s unscoped whole-tree race condition (same `--cell` scoping already applied to `preseed_pasa_checkpoints.py`)
