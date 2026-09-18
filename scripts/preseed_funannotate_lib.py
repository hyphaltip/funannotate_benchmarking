#!/usr/bin/env python3
"""Pre-seed per-cell lib/swissprot_fungi.faa + lib/template.sbt.

nf_funannotate1's params.proteins and params.sbt_template default to
`<launchDir>/lib/swissprot_fungi.faa` and `<launchDir>/lib/template.sbt`
(conf/profile_annotate.config). Unlike `lib/augustus/3.5/config` (auto-seeded
by the SETUP_AUGUSTUS_CONFIG process, storeDir-cached), nothing in the
pipeline populates these two for a new cell -- every prior working cell had
them hand-symlinked once and forgotten about.

This mostly went unnoticed under the conda profile, which reads the real
filesystem path directly, so a missing file there only matters if something
downstream actually opens it. The container/singularity profile bind-mounts
`<launchDir>/lib` (params.proteins' parent dir) into the container
(conf/provision_singularity.config funannotate_binds), and
FUNANNOTATE_PREDICT hard-fails with "<cell>/lib/swissprot_fungi.faa is not a
valid file, exiting" when the file isn't actually there -- confirmed on
v1.9.0-beta10_container(_rust) 2026-09-06 and v1.9.0-beta11_container_rust
2026-09-07 (a "fix" landed once already, 2026-09-04, as a symlink at the
project ROOT lib/ instead of any cell's own runs/<cell>/lib/ -- which is
never a Nextflow launchDir, so it never actually took effect; see the
_stale_pre_swissprot_fix_2026-09-04 directory left behind by that attempt).

Usage:
    python3 scripts/preseed_funannotate_lib.py [--config conf/benchmark.yaml] [--runs-dir runs] [--dry-run]

Rerunnable/idempotent: existing symlinks/files are left untouched.
"""
from __future__ import annotations

import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lib  # noqa: E402

FILES = ["swissprot_fungi.faa", "template.sbt"]


def parse_args():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--config", default="conf/benchmark.yaml")
    ap.add_argument("--runs-dir", default="runs")
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--verbose", "-v", action="count", default=0)
    return ap.parse_args()


def main():
    args = parse_args()
    lib.setup_logging(args.verbose)
    cfg = lib.load_yaml(args.config)
    src_dir = cfg.get("fetch", {}).get("funannotate_lib_dir")
    if not src_dir:
        lib.LOG.error("conf/benchmark.yaml fetch.funannotate_lib_dir is not set")
        sys.exit(1)

    total_cells = total_linked = total_present = total_missing = 0
    for cell in sorted(os.listdir(args.runs_dir)):
        cell_dir = os.path.join(args.runs_dir, cell)
        samples_csv = os.path.join(cell_dir, "samples.csv")
        if not os.path.isdir(cell_dir) or not os.path.exists(samples_csv):
            continue
        lib_dir = os.path.join(cell_dir, "lib")

        n_linked = n_present = n_missing = 0
        for fname in FILES:
            src = os.path.join(src_dir, fname)
            dest = os.path.join(lib_dir, fname)
            if not os.path.exists(src):
                n_missing += 1
                lib.LOG.warning("%s: source missing: %s", cell, src)
                continue
            if os.path.exists(dest) or os.path.islink(dest):
                n_present += 1
                continue
            n_linked += 1
            if not args.dry_run:
                os.makedirs(lib_dir, exist_ok=True)
                os.symlink(src, dest)

        total_cells += 1
        total_linked += n_linked
        total_present += n_present
        total_missing += n_missing
        print(f"{cell:<32} {n_linked} linked, {n_present} already present"
              + (f", {n_missing} missing source" if n_missing else ""))

    print(f"total: {total_cells} cells, {total_linked} symlinks created, "
          f"{total_present} already present, {total_missing} missing source")


if __name__ == "__main__":
    main()
