#!/usr/bin/env python3
"""Pre-seed per-cell RNA-seq reads + SRA-query cache so cells never hit SRA.

nf_funannotate1's FETCH_RNASEQ subworkflow only skips the network entirely
when BOTH of these are true for a species_tag (Genus_species, spaces ->
underscores):

  1. `<launchDir>/rnaseq_reads/<tag>_norm_R1.fastq.gz` (+ _R2, _SE) already
     exist -- SRA_FETCH/SRA_FETCH_SE declare these as storeDir outputs, so
     Nextflow skips the process outright when they're already present
     (modules/local/sra_fetch.nf).
  2. `params.skip_sra_query=true` AND
     `<launchDir>/rnaseq_reads/sra_query/<tag>.sra_query.csv` already exists
     -- otherwise SRA_QUERY_BATCH still makes NCBI network calls every run
     (subworkflows/local/fetch_rnaseq.nf), even though its result is only
     used for routing (PE vs SE vs no-data), not the actual download.

This benchmark already has both sources of truth on disk, just not in the
cell-local layout the pipeline checks:

  * conf/benchmark.yaml pool.rnaseq -- the master SRA manifest
    (species_tag,taxonid,sra_accession,spots,platform,layout), already
    covering every species Fungi_BFD has queried.
  * fetch.rnaseq_reads_dir -- the shared pool of already-downloaded,
    normalized reads (<tag>_norm_R1/_R2/_SE.fastq.gz, including zero-byte
    placeholders for species with no RNA-seq at all).

This script symlinks/derives the cell-local layout from those two sources so
every cell starts from already-masked genomes (scripts/preseed_clean_genomes.py)
AND already-downloaded RNA-seq -- SRA_FETCH and SRA_QUERY_BATCH should never
run for a species this script has covered.

Usage:
    python3 scripts/preseed_rnaseq_reads.py [--config conf/benchmark.yaml] [--runs-dir runs] [--dry-run]

Rerunnable/idempotent: existing symlinks/files are left untouched; missing
ones are created. A species with no rows in pool.rnaseq and no reads in
fetch.rnaseq_reads_dir is reported as "uncovered" (falls through to whatever
--skip_sra_query does with a missing cache CSV: SRA_QUERY_BATCH is skipped
with a warning and that species gets no RNA-seq -- ab-initio training only).
"""
from __future__ import annotations

import argparse
import csv
import os
import sys
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lib  # noqa: E402

READ_SUFFIXES = ("_norm_R1.fastq.gz", "_norm_R2.fastq.gz", "_norm_SE.fastq.gz")
SRA_QUERY_HEADER = ["species_tag", "taxonid", "sra_accession", "spots", "platform", "layout"]


def parse_args():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--config", default="conf/benchmark.yaml")
    ap.add_argument("--runs-dir", default="runs")
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--verbose", "-v", action="count", default=0)
    return ap.parse_args()


def species_tag(species: str) -> str:
    return "_".join(species.strip().split())


def load_manifest_by_tag(manifest_path: str):
    rows_by_tag = defaultdict(list)
    if not os.path.exists(manifest_path):
        lib.LOG.warning("pool.rnaseq manifest not found: %s", manifest_path)
        return rows_by_tag
    for row in lib.read_table(manifest_path):
        tag = row.get("species_tag", "").strip()
        if tag:
            rows_by_tag[tag].append(row)
    return rows_by_tag


def symlink_reads(tag: str, src_dir: str, dest_dir: str, dry_run: bool) -> str:
    """Returns 'linked' | 'present' | 'relinked' | 'missing' for this
    species_tag's read trio.

    The master pool (src_dir) is always authoritative when it has the full
    trio: a dest that's already the correct symlink is left alone ('present'
    -- the common case, cheap no-op), but anything else at dest (a stale
    regular file from a cell's own pre-preseed SRA_FETCH, or a symlink to
    somewhere else) is replaced ('relinked'). This matters in practice: cells
    launched before this script existed can have truncated/corrupt
    leftovers (e.g. a handful of bytes from an aborted chunked download) that
    `os.path.exists()` alone can't tell apart from a real fetch -- trusting
    "already there" caused a real Trinity failure on Allomyces_macrogynus
    (2026-09-07, a 23-byte R1/R2 pair silently reused instead of the real
    120MB pool copy). Only fall back to whatever the cell already has when
    the master pool doesn't cover this species at all.
    """
    srcs = [os.path.join(src_dir, f"{tag}{suf}") for suf in READ_SUFFIXES]
    if not all(os.path.exists(s) for s in srcs):
        dests = [os.path.join(dest_dir, f"{tag}{suf}") for suf in READ_SUFFIXES]
        return "present" if all(os.path.exists(d) or os.path.islink(d) for d in dests) else "missing"

    dests = [os.path.join(dest_dir, f"{tag}{suf}") for suf in READ_SUFFIXES]
    status = "present"
    if not dry_run:
        os.makedirs(dest_dir, exist_ok=True)
    for s, d in zip(srcs, dests):
        if os.path.islink(d) and os.readlink(d) == s:
            continue  # already correct
        if os.path.exists(d) or os.path.islink(d):
            status = "relinked"
            if not dry_run:
                os.remove(d)
        else:
            status = "linked" if status == "present" else status
        if not dry_run:
            os.symlink(s, d)
    return status


def write_sra_query_csv(tag: str, rows, dest_dir: str, dry_run: bool) -> str:
    """Returns 'written' | 'present'. Always writes a header even with 0 rows
    (species with no RNA-seq at all) so skip_sra_query finds a cache hit
    instead of silently dropping the species."""
    dest = os.path.join(dest_dir, f"{tag}.sra_query.csv")
    if os.path.exists(dest):
        return "present"
    if not dry_run:
        os.makedirs(dest_dir, exist_ok=True)
        with open(dest, "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(SRA_QUERY_HEADER)
            for r in rows:
                w.writerow([r.get(c, "") for c in SRA_QUERY_HEADER])
    return "written"


def main():
    args = parse_args()
    lib.setup_logging(args.verbose)
    cfg = lib.load_yaml(args.config)
    fetch = cfg.get("fetch", {})
    pool = cfg.get("pool", {})

    reads_src_dir = fetch.get("rnaseq_reads_dir")
    manifest_path = pool.get("rnaseq")
    if not reads_src_dir:
        lib.LOG.error("conf/benchmark.yaml fetch.rnaseq_reads_dir is not set")
        sys.exit(1)
    manifest_by_tag = load_manifest_by_tag(manifest_path) if manifest_path else {}

    total_cells = total_linked = total_present = total_relinked = total_missing = 0
    total_csv_written = total_csv_present = 0
    for cell in sorted(os.listdir(args.runs_dir)):
        cell_dir = os.path.join(args.runs_dir, cell)
        samples_csv = os.path.join(cell_dir, "samples.csv")
        if not os.path.isdir(cell_dir) or not os.path.exists(samples_csv):
            continue
        rows = lib.read_table(samples_csv)
        tags = sorted({species_tag(r.get("SPECIES", "")) for r in rows if r.get("SPECIES", "").strip()})

        reads_dest = os.path.join(cell_dir, "rnaseq_reads")
        query_dest = os.path.join(reads_dest, "sra_query")

        n_linked = n_present = n_relinked = n_missing = n_csv_written = n_csv_present = 0
        for tag in tags:
            status = symlink_reads(tag, reads_src_dir, reads_dest, args.dry_run)
            if status == "linked":
                n_linked += 1
            elif status == "present":
                n_present += 1
            elif status == "relinked":
                n_relinked += 1
                lib.LOG.warning("%s: replaced stale/incorrect reads with master pool copy (%s)", cell, tag)
            else:
                n_missing += 1
                lib.LOG.warning("%s: no cached reads under %s (%s)", cell, reads_src_dir, tag)

            csv_status = write_sra_query_csv(tag, manifest_by_tag.get(tag, []), query_dest, args.dry_run)
            if csv_status == "written":
                n_csv_written += 1
            else:
                n_csv_present += 1

        total_cells += 1
        total_linked += n_linked
        total_relinked += n_relinked
        total_present += n_present
        total_missing += n_missing
        total_csv_written += n_csv_written
        total_csv_present += n_csv_present
        print(f"{cell:<32} reads: {n_linked} linked, {n_present} already present, "
              f"{n_relinked} relinked (stale/incorrect replaced), {n_missing} uncovered  |  "
              f"sra_query csv: {n_csv_written} written, {n_csv_present} already present")

    print(f"total: {total_cells} cells, {total_linked} read-trios linked, "
          f"{total_present} already present, {total_relinked} relinked, "
          f"{total_missing} uncovered species, "
          f"{total_csv_written} sra_query csvs written, {total_csv_present} already present")
    if total_missing:
        print(f"[WARN] {total_missing} species have no cached reads in {reads_src_dir} "
              "-- they will get ab-initio-only training (no cache CSV entry to route SRA_FETCH to)")


if __name__ == "__main__":
    main()
