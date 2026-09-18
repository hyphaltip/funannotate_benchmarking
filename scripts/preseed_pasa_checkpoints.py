#!/usr/bin/env python3
"""Clear stale PASA MySQL checkpoints before a cell (re)launch.

FUNANNOTATE_TRAIN spins up a fresh MariaDB sidecar for every task attempt
("Initializing fresh MariaDB system tables via sidecar image's
mysql_install_db" -- no data persists between attempts). PASA's own
checkpoint markers, though, live on disk under
`<training_dir>/pasa*/__pasa_<genome>_pasa_mysql_chkpts/*.ok` and DO persist
across relaunches/`-resume`. Once a genome's PASA run has gotten past
`create_db.ok` once, any later re-attempt sees that checkpoint, skips
re-creating the database against the new (empty) MariaDB instance, and dies
at `update_fli_status.dbi` with `Unknown database '<genome>_pasa'`.

Confirmed 2026-09-14 on v1.8.17_conda: every genome re-attempted after an
earlier successful-looking run hit this identical error (Malassezia_globosa,
Saccharomyces_cerevisiae, Umbelopsis_ramanniana, Aspergillus_fumigatus, ...) --
including genomes whose predict_results/*.gff3 already existed from a prior
success, since something upstream still invalidated nextflow's task cache and
forced FUNANNOTATE_TRAIN to re-execute. Existing predict output is therefore
NOT a reliable signal that a genome's TRAIN task won't re-run.

Fix: unconditionally delete every genome's PASA checkpoint dir(s) -- both the
live `pasa/` and any `pasa.stale_*/` funannotate itself renamed aside -- before
each launch, so PASA always starts clean against that task's fresh MySQL
sidecar. This is harmless for a genome nextflow ends up skipping via -resume
(a cached task never looks at this directory again).

Usage:
    python3 scripts/preseed_pasa_checkpoints.py [--runs-dir runs] [--dry-run]

Rerunnable/idempotent: a genome with no pasa* checkpoint dir is a no-op.
"""
from __future__ import annotations

import argparse
import glob
import os
import shutil
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

    total_cells = total_cleared = 0
    for cell in sorted(os.listdir(args.runs_dir)):
        cell_dir = os.path.join(args.runs_dir, cell)
        samples_csv = os.path.join(cell_dir, "samples.csv")
        training_root = os.path.join(cell_dir, "genome_annotation_training")
        if not os.path.isdir(cell_dir) or not os.path.exists(samples_csv):
            continue
        if not os.path.isdir(training_root):
            continue

        n_cleared = 0
        for genome in sorted(os.listdir(training_root)):
            pasa_dirs = sorted(glob.glob(
                os.path.join(training_root, genome, "training", "pasa*")))
            if not pasa_dirs:
                continue
            n_cleared += 1
            for d in pasa_dirs:
                lib.LOG.info("%s/%s: removing PASA checkpoint dir %s", cell, genome, d)
                if not args.dry_run:
                    shutil.rmtree(d)

        total_cells += 1
        total_cleared += n_cleared
        print(f"{cell:<32} {n_cleared} genomes' PASA checkpoints cleared")

    print(f"total: {total_cells} cells, {total_cleared} genomes' PASA checkpoints cleared")


if __name__ == "__main__":
    main()
