#!/usr/bin/env python3
"""Post-hoc BUSCO completeness for delivered predict_results proteins.

DESIGN.md "Gene-content comparison" calls for BUSCO completeness per cell.
funannotate's own logfiles/busco.log only captures INTERMEDIATE ab-initio
training checkpoints -- BUSCO in genome mode run BEFORE any gene prediction
(to assess assembly completeness / pick augustus training loci), and BUSCO in
protein mode run against augustus's raw training predictions specifically.
Neither reflects the completeness of the FINAL delivered gene set in
predict_results/, and (confirmed) Fungi_BFD_runs doesn't even keep that log
for older genomes -- so reusing busco.log values would silently compare two
different things across cells/trees, not the annotation quality DESIGN.md
actually wants.

This script runs BUSCO directly against each cell's
predict_results/<genome>.proteins.fa, so completeness is a genuine, uniform
measurement across every cell/version, independent of what each funannotate
release happens to log internally.

Output lands under
  runs/<cell>/genome_annotation/<genome>/busco_completeness/
so lib.find_busco_summary() (which walks the whole genome dir looking for
short_summary*.txt) picks it up automatically -- compare_predictions.py and
validate_harness.py need no changes to start using these results.

Idempotent/resumable: skips a genome as soon as scripts/lib.py's
find_busco_summary() already finds a short_summary*.txt anywhere under its
genome_annotation/<genome>/ dir (whether from a prior run of this script or,
in principle, from funannotate itself). Re-run any time more genomes finish.

Requires `busco` on PATH (UCR HPCC: `module load busco`).

Usage:
  module load busco
  python3 scripts/run_busco_completeness.py --cells-dir runs --samples samples.csv
  python3 scripts/run_busco_completeness.py --cells v1.8.17_conda --genomes Foo_bar --cpus 8
"""
from __future__ import annotations

import argparse
import glob
import os
import shutil
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lib  # noqa: E402


def find_predict_proteins(genome_dir: str) -> str | None:
    hits = glob.glob(os.path.join(genome_dir, "predict_results", "*.proteins.fa"))
    return hits[0] if hits else None


def run_busco(proteins_fa: str, lineage_path: str, out_dir: str, out_name: str, cpus: int) -> bool:
    os.makedirs(out_dir, exist_ok=True)
    cmd = [
        "busco", "-i", proteins_fa, "-l", lineage_path, "-o", out_name,
        "--out_path", out_dir, "-m", "proteins", "-c", str(cpus), "-f",
    ]
    try:
        subprocess.run(cmd, check=True, capture_output=True, text=True, timeout=3600)
    except (subprocess.CalledProcessError, subprocess.TimeoutExpired) as e:
        lib.LOG.warning("busco failed for %s: %s", proteins_fa, e)
        return False
    return True


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--cells-dir", default="runs")
    ap.add_argument("--cells", nargs="*", default=None, help="restrict to these cell names (default: all under --cells-dir with a genome_annotation/ dir)")
    ap.add_argument("--genomes", nargs="*", default=None, help="restrict to these genome tags (default: every completed genome in each cell)")
    ap.add_argument("--samples", default="samples.csv", help="pool metadata (for busco_lineage_used per genome)")
    ap.add_argument("--busco-lineages-dir", default="/srv/projects/db/BUSCO/v10/lineages")
    ap.add_argument("--cpus", type=int, default=8)
    ap.add_argument("-v", "--verbose", action="count", default=0)
    args = ap.parse_args()
    lib.setup_logging(args.verbose)

    if shutil.which("busco") is None:
        lib.LOG.error("busco not on PATH -- module load busco first")
        sys.exit(1)
    if not os.path.exists(args.samples):
        lib.LOG.error("samples.csv not found: %s", args.samples)
        sys.exit(1)
    samples_idx = lib.load_samples_index(args.samples)

    cells = args.cells or sorted(
        n for n in os.listdir(args.cells_dir)
        if os.path.isdir(os.path.join(args.cells_dir, n, "genome_annotation"))
    )
    if not cells:
        lib.LOG.error("no cells with a genome_annotation/ dir found under %s", args.cells_dir)
        sys.exit(1)
    lib.LOG.info("cells: %s", ", ".join(cells))

    n_ran = n_skipped_done = n_skipped_no_predict = n_skipped_no_lineage = n_failed = 0
    for cell in cells:
        ga_dir = os.path.join(args.cells_dir, cell, "genome_annotation")
        genomes = args.genomes or sorted(os.listdir(ga_dir))
        for genome in genomes:
            genome_dir = os.path.join(ga_dir, genome)
            if not os.path.isdir(genome_dir):
                continue
            if lib.find_busco_summary(genome_dir):
                n_skipped_done += 1
                continue
            proteins_fa = find_predict_proteins(genome_dir)
            if not proteins_fa:
                n_skipped_no_predict += 1
                continue
            srow = samples_idx.get(genome)
            lineage = srow.get("busco_lineage_used", "") if srow else ""
            if not lineage:
                lib.LOG.warning("%s/%s: no busco_lineage_used in samples.csv -- skipping", cell, genome)
                n_skipped_no_lineage += 1
                continue
            lineage_path = os.path.join(args.busco_lineages_dir, lineage)
            if not os.path.isdir(lineage_path):
                lib.LOG.warning("%s/%s: lineage dir not found: %s -- skipping", cell, genome, lineage_path)
                n_skipped_no_lineage += 1
                continue
            out_dir = os.path.join(genome_dir, "busco_completeness")
            out_name = genome
            lib.LOG.info("%s/%s: running busco (%s)", cell, genome, lineage)
            if run_busco(proteins_fa, lineage_path, out_dir, out_name, args.cpus):
                n_ran += 1
            else:
                n_failed += 1

    print(
        f"ran={n_ran} already_done={n_skipped_done} no_predict_output={n_skipped_no_predict} "
        f"no_lineage={n_skipped_no_lineage} failed={n_failed}",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
