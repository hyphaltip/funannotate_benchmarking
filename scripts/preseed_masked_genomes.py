#!/usr/bin/env python3
"""Pre-seed a cell's tantan-masked genomes from a sibling cell's completed output.

nf_funannotate1's MASKREPEAT_TANTAN_RUN uses storeDir and is skipped whenever
`<launchDir>/input_clean_genomes/<ASMID>.masked.fasta.gz` already exists (see
modules/local/maskrepeat_tantan_run.nf). Tantan masking depends only on the
funannotate/tantan version, not on conda-vs-container provisioning, so a cell
that already masked a genome under the same funannotate version can donate
that file to a sibling cell stuck re-running the same masking step (e.g. the
v1.9.0-beta10 container cells, which were failing MASKREPEAT_TANTAN_RUN on
squashfuse-less nodes -- see conf/site_ucr_hpcc_singularity.config).

This script symlinks `<runs-dir>/<target-cell>/input_clean_genomes/<ASMID>.masked.fasta.gz`
to the same file under `<runs-dir>/<source-cell>/input_clean_genomes/`, for
every ASMID in the target cell's samples.csv that's missing its masked genome
but present in the source cell.

Usage:
    python3 scripts/preseed_masked_genomes.py --source-cell v1.9.0-beta10_conda \\
        --target-cells v1.9.0-beta10_container v1.9.0-beta10_container_rust

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
    ap.add_argument("--source-cell", required=True,
                     help="cell whose input_clean_genomes/<ASMID>.masked.fasta.gz to reuse")
    ap.add_argument("--target-cells", nargs="+", required=True,
                     help="cells to seed with the source cell's masked genomes")
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--verbose", "-v", action="count", default=0)
    return ap.parse_args()


def main():
    args = parse_args()
    lib.setup_logging(args.verbose)

    source_masked_dir = os.path.join(args.runs_dir, args.source_cell, "input_clean_genomes")
    if not os.path.isdir(source_masked_dir):
        lib.LOG.error("source cell masked dir not found: %s", source_masked_dir)
        sys.exit(1)

    total_link = total_skip = total_missing = 0
    for cell in args.target_cells:
        cell_dir = os.path.join(args.runs_dir, cell)
        samples_csv = os.path.join(cell_dir, "samples.csv")
        if not os.path.exists(samples_csv):
            lib.LOG.warning("%s: no samples.csv, skipping", cell)
            continue
        with open(samples_csv, newline="") as fh:
            rows = list(csv.DictReader(fh))
        target_dir = os.path.join(cell_dir, "input_clean_genomes")
        if not args.dry_run:
            os.makedirs(target_dir, exist_ok=True)
        n_link = n_skip = n_missing = 0
        for r in rows:
            asmid = r.get("ASMID", "").strip()
            if not asmid:
                continue
            fname = f"{asmid}.masked.fasta.gz"
            target = os.path.join(target_dir, fname)
            if os.path.exists(target) or os.path.islink(target):
                n_skip += 1
                continue
            source = os.path.join(source_masked_dir, fname)
            if not os.path.exists(source):
                lib.LOG.warning("%s: %s not yet masked in source cell %s", cell, asmid, args.source_cell)
                n_missing += 1
                continue
            if not args.dry_run:
                os.symlink(os.path.abspath(source), target)
            n_link += 1
        total_link += n_link
        total_skip += n_skip
        total_missing += n_missing
        print(f"{cell:<32} {n_link} linked, {n_skip} already present"
              + (f", {n_missing} not available in source" if n_missing else ""))

    print(f"total: {total_link} symlinks created, {total_skip} already present, "
          f"{total_missing} missing from source")


if __name__ == "__main__":
    main()
