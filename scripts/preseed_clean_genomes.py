#!/usr/bin/env python3
"""Pre-seed per-cell cleaned genomes so benchmark cells start after GENOME_CLEAN.

nf_funannotate1's CLEAN_GENOMES subworkflow skips the GENOME_CLEAN process for
any sample whose cleaned genome already exists as
`<launchDir>/input_clean_genomes/<ASMID>.fa.gz` (see
subworkflows/local/clean_genomes.nf). The benchmark consumes already-masked
pool genomes that ARE the post-GENOME_CLEAN product, so the 6 cells should
never re-run FCS-GX/length filtering -- they start at RNASEQ_PREP from the
pre-seeded cleaned genome.

This script symlinks each cell's `<runs-dir>/<cell>/input_clean_genomes/<ASMID>.fa.gz`
to the shared masked genome the cell's samples.csv GENOME column points at,
exactly matching the storeDir target the pipeline checks for. Symlinks (not
copies) keep storage flat across cells.

Usage:
    python3 scripts/preseed_clean_genomes.py [--runs-dir runs] [--dry-run]

Rerunnable/idempotent: existing symlinks/files are left untouched; missing
ones are created.
"""
from __future__ import annotations

import argparse
import csv
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lib  # noqa: E402


def parse_args():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--runs-dir", default="runs")
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--verbose", "-v", action="count", default=0)
    return ap.parse_args()


def main():
    args = parse_args()
    lib.setup_logging(args.verbose)

    total_cells = 0
    total_link = 0
    total_skip = 0
    total_missing = 0
    for cell in sorted(os.listdir(args.runs_dir)):
        cell_dir = os.path.join(args.runs_dir, cell)
        samples_csv = os.path.join(cell_dir, "samples.csv")
        if not os.path.isdir(cell_dir) or not os.path.exists(samples_csv):
            continue
        with open(samples_csv, newline="") as fh:
            rows = list(csv.DictReader(fh))
        clean_dir = os.path.join(cell_dir, "input_clean_genomes")
        if not args.dry_run:
            os.makedirs(clean_dir, exist_ok=True)
        n_link = n_skip = n_missing = 0
        for r in rows:
            asmid = r.get("ASMID", "").strip()
            genome = r.get("GENOME", "").strip()
            if not asmid or not genome:
                continue
            target = os.path.join(clean_dir, f"{asmid}.fa.gz")
            if os.path.exists(target) or os.path.islink(target):
                n_skip += 1
                continue
            if not os.path.exists(genome):
                lib.LOG.warning("%s: GENOME missing: %s", asmid, genome)
                n_missing += 1
                continue
            if not args.dry_run:
                os.symlink(genome, target)
            n_link += 1
        total_cells += 1
        total_link += n_link
        total_skip += n_skip
        total_missing += n_missing
        print(f"{cell:<32} {n_link} linked, {n_skip} already present"
              + (f", {n_missing} missing GENOME" if n_missing else ""))

    print(f"total: {total_cells} cells, {total_link} symlinks created, "
          f"{total_skip} already present, {total_missing} missing")


if __name__ == "__main__":
    main()
