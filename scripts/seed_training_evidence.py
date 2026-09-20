#!/usr/bin/env python3
"""Copy per-genome funannotate training evidence from one cell into another.

Why this exists
---------------
The v1.9.0-beta.12 container pair (`-norust.sif` / `.sif`) is byte-identical
apart from FUNANNOTATE_EVM_ENGINE (perl vs rust) -- same Trinity 2.16.1_rust,
same salmon 2.7.0, same PASA, same TransDecoder (verified 2026-09-19). EVM runs
in `funannotate predict`, NOT in `funannotate train`, so training both cells
independently would:

  * burn ~16 redundant multi-hour Trinity assemblies, and
  * WEAKEN the comparison, by letting run-to-run nondeterminism differ between
    the two arms.

Instead: train once in the source cell, copy the finished evidence into the
target cell, and let the target's FUNANNOTATE_TRAIN short-circuit on the
published `funannotate_train.pasa.gff3` ("training already resolved"). Both arms
then run `predict` over provably identical input, so any difference in the
output is attributable to the EVM engine alone.

Only the `funannotate_train.*` set is copied (~142 MB/genome dereferenced), not
the whole training directory (~18 GB/genome). Symlinks are dereferenced so the
target is self-contained.

Usage
-----
    python3 scripts/seed_training_evidence.py \
        --from-cell v1.9.0-beta12_container \
        --to-cell   v1.9.0-beta12_container_rust [--dry-run]

Re-runnable: a genome already seeded with a matching checksum is skipped.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import os
import shutil
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lib  # noqa: E402

PREFIX = "funannotate_train."
# The file the train module's "already resolved" check keys on.
SENTINEL = "funannotate_train.pasa.gff3"


def md5(path, blocks=1 << 20):
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(blocks), b""):
            h.update(chunk)
    return h.hexdigest()


def genomes_of(cell_dir):
    sc = os.path.join(cell_dir, "samples.csv")
    if not os.path.exists(sc):
        return []
    out = []
    for r in csv.DictReader(open(sc)):
        out.append("{}_{}".format(r["SPECIES"].replace(" ", "_"),
                                  r["STRAIN"].replace(" ", "_")))
    return out


def parse_args():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--runs-dir", default="runs")
    ap.add_argument("--from-cell", required=True)
    ap.add_argument("--to-cell", required=True)
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--verbose", "-v", action="count", default=0)
    return ap.parse_args()


def main():
    args = parse_args()
    lib.setup_logging(args.verbose)

    src_cell = os.path.join(args.runs_dir, args.from_cell)
    dst_cell = os.path.join(args.runs_dir, args.to_cell)
    for d in (src_cell, dst_cell):
        if not os.path.isdir(d):
            lib.LOG.error("cell dir missing: %s", d)
            return 2

    genomes = genomes_of(dst_cell) or genomes_of(src_cell)
    seeded = skipped = missing = 0

    for g in genomes:
        sdir = os.path.join(src_cell, "genome_annotation_training", g, "training")
        ddir = os.path.join(dst_cell, "genome_annotation_training", g, "training")
        sentinel_src = os.path.join(sdir, SENTINEL)

        if not os.path.exists(sentinel_src):
            print(f"  {g:<44} source not trained yet -- skip")
            missing += 1
            continue

        files = sorted(f for f in os.listdir(sdir) if f.startswith(PREFIX))
        dst_sent = os.path.join(ddir, SENTINEL)
        if os.path.exists(dst_sent) and md5(dst_sent) == md5(sentinel_src):
            print(f"  {g:<44} already seeded (checksum match) -- skip")
            skipped += 1
            continue

        total = 0
        if not args.dry_run:
            os.makedirs(ddir, exist_ok=True)
        for f in files:
            s = os.path.join(sdir, f)
            if not os.path.isfile(s):      # follows symlinks
                continue
            total += os.path.getsize(os.path.realpath(s))
            if not args.dry_run:
                # copyfile() reads through the symlink -> self-contained target
                shutil.copyfile(s, os.path.join(ddir, f))

        verdict = "DRY-RUN"
        if not args.dry_run:
            verdict = "OK" if md5(dst_sent) == md5(sentinel_src) else "CHECKSUM MISMATCH"
        print(f"  {g:<44} {len(files)} files, {total/1048576:.0f} MB  {verdict}")
        if verdict == "CHECKSUM MISMATCH":
            lib.LOG.error("%s: seeded evidence does not match source", g)
            return 1
        seeded += 1

    print(f"\n  seeded={seeded} already_present={skipped} source_untrained={missing}")
    if missing:
        print("  NOTE: genomes whose source cell has not finished training were skipped;"
              "\n        re-run this script after the source cell completes.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
