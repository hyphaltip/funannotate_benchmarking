#!/usr/bin/env python3
"""Tag species with a 0-transcript (or missing) Trinity-GG assembly for
future rebuild/debug triage, distinguishing genuine undiagnosed failures from
already-explained or intentional cases.

Motivation: 5,822 of 7,499 rnaseq_data/*.trinity-GG.fasta entries are
0-transcript/empty project-wide. Most of those are legitimately "no RNA-seq
available" (SRA query found nothing) or intentionally excluded (documented in
rnaseq_skip.csv / rnaseq_force_no_rnaseq.csv) or composite-hybrid placeholders
(e.g. Cryptococcus_deneoformans -> a *.composite-parents.parent_prefix_map.tsv
pointer, not a real assembly attempt -- confirmed NOT a failure, working as
designed via BUILD_HYBRID_COMPOSITE_TRINITY). A small subset, though, are
genuine undiagnosed Trinity failures where real reads exist but assembly
produced nothing (the pattern already found and root-caused for Rhodotorula
toruloides this session, and for the two MAG species already documented in
rnaseq_force_no_rnaseq.csv) -- THOSE are what this script surfaces.

Categories assigned, in priority order:
  KNOWN_SKIP           in rnaseq_skip.csv (already diagnosed, ab-initio in use)
  KNOWN_FORCE_NO_RNASEQ in rnaseq_force_no_rnaseq.csv (already diagnosed)
  COMPOSITE_HYBRID     trinity-GG.fasta content matches a composite-parents
                       pointer file rather than real FASTA records (>= 1
                       '>' record required to NOT be flagged this way)
  NO_READS_FOUND       no non-trivial (>1KB) norm_R1/R2 or norm_SE read file
                       found in rnaseq_reads/ or rnaseq_normalized/ -- SRA
                       query legitimately found nothing, or reads were never
                       fetched; ab-initio-only is the correct, non-bug outcome
  UNDIAGNOSED_FAILURE  real, non-trivial read files exist, species is not in
                       either known-exclusion list, and is not a composite
                       placeholder, yet trinity-GG.fasta has zero transcripts
                       -- THIS is the priority list for rebuild/debug, same
                       failure shape as Rhodotorula toruloides (read-length
                       mismatch) and the two already-documented MAG cases.
  HAS_TRANSCRIPTS      not empty; not in scope for this script

Output: results/trinity_failure_tags.tsv, one row per species in
rnaseq_data/*.trinity-GG.fasta (all of them, not just empty ones, so
HAS_TRANSCRIPTS rows are there for completeness/cross-checking).

Usage:
    python3 scripts/tag_trinity_failures.py \
        --fungi-bfd-runs /bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs \
        --out results/trinity_failure_tags.tsv
"""
from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

MIN_READ_BYTES = 1024  # placeholders are 0 bytes; real reads are MB+


def load_skip_species(path: Path) -> list[str]:
    """Returns the raw 'out' tags (species_strain), used as PREFIXES against
    the bare species name from rnaseq_data/*.trinity-GG.fasta -- the skip
    list's tags carry a strain suffix the bare species name doesn't."""
    tags = []
    if not path.is_file():
        return tags
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh)
        if reader.fieldnames and "out" in reader.fieldnames:
            for row in reader:
                tag = (row.get("out") or "").strip()
                if tag:
                    tags.append(tag)
    return tags


def species_in_skip_list(species: str, skip_tags: list[str]) -> bool:
    return any(tag == species or tag.startswith(species + "_") for tag in skip_tags)


def load_force_no_rnaseq_species(path: Path) -> set[str]:
    species = set()
    if not path.is_file():
        return species
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            tag = (row.get("species_tag") or "").strip()
            if tag:
                species.add(tag)
    return species


def is_composite_placeholder(fasta_path: Path) -> bool:
    """True if the file has no real FASTA '>' records (composite-parent
    pointer files are a bare TSV line, not FASTA)."""
    try:
        with open(fasta_path, "rb") as fh:
            head = fh.read(4096)
    except OSError:
        return False
    return b">" not in head


def real_read_file_present(reads_dirs: list[Path], species: str) -> bool:
    for d in reads_dirs:
        for suffix in ("_norm_R1.fastq.gz", "_norm_R2.fastq.gz", "_norm_SE.fastq.gz"):
            p = d / f"{species}{suffix}"
            try:
                if p.is_file() and p.stat().st_size > MIN_READ_BYTES:
                    return True
            except OSError:
                continue
    return False


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--fungi-bfd-runs", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    args = ap.parse_args()

    root = args.fungi_bfd_runs
    rnaseq_data = root / "rnaseq_data"
    reads_dirs = [root / "rnaseq_reads", root / "rnaseq_normalized"]

    skip_tags = load_skip_species(root / "rnaseq_skip.csv")
    force_no_species = load_force_no_rnaseq_species(root / "rnaseq_force_no_rnaseq.csv")
    print(f"[INFO] {len(skip_tags)} tags in rnaseq_skip.csv, "
          f"{len(force_no_species)} species in rnaseq_force_no_rnaseq.csv", file=sys.stderr)

    fasta_files = sorted(rnaseq_data.glob("*.trinity-GG.fasta"))
    print(f"[INFO] {len(fasta_files)} trinity-GG.fasta entries to classify", file=sys.stderr)

    rows = []
    for i, fasta in enumerate(fasta_files, 1):
        if i % 1000 == 0:
            print(f"[INFO] ... {i}/{len(fasta_files)}", file=sys.stderr)
        species = fasta.name[: -len(".trinity-GG.fasta")]
        size = fasta.stat().st_size

        if size > 0 and not is_composite_placeholder(fasta):
            rows.append({"species": species, "size_bytes": size, "category": "HAS_TRANSCRIPTS", "note": ""})
            continue

        if species_in_skip_list(species, skip_tags):
            rows.append({"species": species, "size_bytes": size, "category": "KNOWN_SKIP", "note": "see rnaseq_skip.csv"})
            continue
        if species in force_no_species:
            rows.append({"species": species, "size_bytes": size, "category": "KNOWN_FORCE_NO_RNASEQ", "note": "see rnaseq_force_no_rnaseq.csv"})
            continue
        if size > 0 and is_composite_placeholder(fasta):
            rows.append({"species": species, "size_bytes": size, "category": "COMPOSITE_HYBRID", "note": "composite-parents pointer, not a real assembly attempt"})
            continue

        has_reads = real_read_file_present(reads_dirs, species)
        if not has_reads:
            rows.append({"species": species, "size_bytes": size, "category": "NO_READS_FOUND", "note": "no non-trivial read file in rnaseq_reads/ or rnaseq_normalized/"})
        else:
            rows.append({"species": species, "size_bytes": size, "category": "UNDIAGNOSED_FAILURE",
                         "note": "real reads present but 0-transcript assembly, not in either known-exclusion list -- rebuild/debug candidate"})

    args.out.parent.mkdir(parents=True, exist_ok=True)
    # Priority sort: undiagnosed failures first.
    priority = {"UNDIAGNOSED_FAILURE": 0, "COMPOSITE_HYBRID": 1, "KNOWN_SKIP": 2,
                "KNOWN_FORCE_NO_RNASEQ": 3, "NO_READS_FOUND": 4, "HAS_TRANSCRIPTS": 5}
    rows.sort(key=lambda r: (priority.get(r["category"], 9), r["species"]))

    with open(args.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["species", "size_bytes", "category", "note"], delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    from collections import Counter
    counts = Counter(r["category"] for r in rows)
    print(f"[INFO] wrote {len(rows)} rows -> {args.out}", file=sys.stderr)
    for cat in ("UNDIAGNOSED_FAILURE", "COMPOSITE_HYBRID", "KNOWN_SKIP", "KNOWN_FORCE_NO_RNASEQ", "NO_READS_FOUND", "HAS_TRANSCRIPTS"):
        print(f"[SUMMARY] {cat}: {counts.get(cat, 0)}", file=sys.stderr)


if __name__ == "__main__":
    main()
