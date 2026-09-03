#!/usr/bin/env python3
"""Validation gate: nf_funannotate1 output vs Fungi_BFD's known-good output.

DESIGN.md "Validation gate (must pass before trusting the N=65 run)": before
trusting a full benchmark cell, spot-check that this newer
nf_funannotate1-driven harness reproduces Fungi_BFD's established results
(same funannotate engine, different wrapper) on a handful of genomes
Fungi_BFD has already annotated well. This exists because nf_funannotate1 is
newer/less battle-tested than Fungi_BFD's production pipeline — the
benchmark must not end up measuring wrapper immaturity instead of
funannotate itself.

This script does NOT launch the pipeline (see run_annotate.sh / sbatch
launchers for that) — it diffs whichever benchmark-cell output already
exists against Fungi_BFD_runs/genome_annotation/ for the same genome tag:
  * gene count delta
  * BUSCO completeness delta (short_summary*.txt, if present on both sides)
  * gffcompare sensitivity/precision between the two GFF3s directly (not vs
    RefSeq — this is a wrapper-vs-wrapper check, RefSeq comparison is
    compare_predictions.py's job)

A genome only enters the gate if funannotate predict has finished for it in
BOTH trees. Pass/fail thresholds are a loose sanity check, not a strict
regression test, and are overridable.

Usage:
  module load gffcompare
  python3 scripts/validate_harness.py --cell v1.8.17_conda
  python3 scripts/validate_harness.py --cell v1.8.17_conda \\
      --genomes Aspergillus_fumigatus_Af293 Saccharomyces_cerevisiae_S288C
"""
from __future__ import annotations

import argparse
import os
import sys
import tempfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lib  # noqa: E402


def pick_genomes(cells_dir: str, cell: str, fungi_bfd_dir: str, n: int) -> list:
    ours = os.path.join(cells_dir, cell, "genome_annotation")
    theirs = os.path.join(fungi_bfd_dir, "genome_annotation")
    if not (os.path.isdir(ours) and os.path.isdir(theirs)):
        return []
    ours_done = {d for d in os.listdir(ours) if lib.find_predict_gff3(os.path.join(ours, d))}
    theirs_done = {d for d in os.listdir(theirs) if lib.find_predict_gff3(os.path.join(theirs, d))}
    return sorted(ours_done & theirs_done)[:n]


def compare_one(genome: str, our_dir: str, their_dir: str, tmp: str) -> dict:
    row = {"genome": genome}
    q_gff = lib.find_predict_gff3(our_dir)
    r_gff = lib.find_predict_gff3(their_dir)
    row["ours_gff3"] = q_gff or ""
    row["fungi_bfd_gff3"] = r_gff or ""
    if not q_gff or not r_gff:
        row["compare_status"] = "missing_output"
        return row

    row["ours_gene_count"] = lib.count_gff3_genes(q_gff)
    row["fungi_bfd_gene_count"] = lib.count_gff3_genes(r_gff)
    if row["fungi_bfd_gene_count"]:
        row["gene_count_delta_pct"] = round(
            100.0 * (row["ours_gene_count"] - row["fungi_bfd_gene_count"]) / row["fungi_bfd_gene_count"], 2)
    else:
        row["gene_count_delta_pct"] = ""

    ours_busco = lib.find_busco_summary(our_dir)
    theirs_busco = lib.find_busco_summary(their_dir)
    if ours_busco:
        row.update({f"ours_{k}": v for k, v in lib.parse_busco_summary(ours_busco).items()})
    if theirs_busco:
        row.update({f"fungi_bfd_{k}": v for k, v in lib.parse_busco_summary(theirs_busco).items()})
    try:
        delta = float(row["ours_busco_complete_pct"]) - float(row["fungi_bfd_busco_complete_pct"])
        row["busco_complete_pct_delta"] = round(delta, 2)
    except (TypeError, ValueError, KeyError):
        row["busco_complete_pct_delta"] = ""

    gffc = lib.run_gffcompare(q_gff, r_gff, os.path.join(tmp, genome))
    if gffc:
        row.update(gffc)

    row["compare_status"] = "compared"
    return row


def gate_verdict(row: dict, max_busco_delta: float, max_gene_pct: float) -> str:
    if row.get("compare_status") != "compared":
        return "SKIP"
    problems = []
    try:
        if abs(float(row["busco_complete_pct_delta"])) > max_busco_delta:
            problems.append("busco")
    except (TypeError, ValueError, KeyError):
        pass
    try:
        if abs(float(row["gene_count_delta_pct"])) > max_gene_pct:
            problems.append("gene_count")
    except (TypeError, ValueError, KeyError):
        pass
    return "FAIL" if problems else "PASS"


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--cell", required=True, help="benchmark cell to validate, e.g. v1.8.17_conda")
    ap.add_argument("--cells-dir", default="runs")
    ap.add_argument("--fungi-bfd-dir", default="/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs")
    ap.add_argument("--genomes", nargs="*", default=None, help="genome tags to check (default: auto-pick)")
    ap.add_argument("--n-genomes", type=int, default=3)
    ap.add_argument("--max-busco-delta", type=float, default=2.0, help="max abs BUSCO complete%% delta before FAIL")
    ap.add_argument("--max-gene-count-delta-pct", type=float, default=10.0)
    ap.add_argument("--out", default="results/validate_harness.tsv")
    ap.add_argument("-v", "--verbose", action="count", default=0)
    args = ap.parse_args()
    lib.setup_logging(args.verbose)

    genomes = args.genomes or pick_genomes(args.cells_dir, args.cell, args.fungi_bfd_dir, args.n_genomes)
    if not genomes:
        lib.LOG.error(
            "no genome has completed predict in both %s/%s/genome_annotation and "
            "%s/genome_annotation yet — nothing to validate", args.cells_dir, args.cell, args.fungi_bfd_dir)
        sys.exit(1)
    lib.LOG.info("validating %d genome(s): %s", len(genomes), ", ".join(genomes))

    rows = []
    with tempfile.TemporaryDirectory(prefix="validate_harness_") as tmp:
        for genome in genomes:
            our_dir = os.path.join(args.cells_dir, args.cell, "genome_annotation", genome)
            their_dir = os.path.join(args.fungi_bfd_dir, "genome_annotation", genome)
            row = compare_one(genome, our_dir, their_dir, tmp)
            row["verdict"] = gate_verdict(row, args.max_busco_delta, args.max_gene_count_delta_pct)
            rows.append(row)

    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    lib.write_table(args.out, rows, sep="\t")

    n_pass = sum(1 for r in rows if r["verdict"] == "PASS")
    n_fail = sum(1 for r in rows if r["verdict"] == "FAIL")
    n_skip = sum(1 for r in rows if r["verdict"] == "SKIP")
    print(f"validate_harness[{args.cell}]: {n_pass} PASS, {n_fail} FAIL, {n_skip} SKIP "
          f"(of {len(rows)}) -> {args.out}", file=sys.stderr)
    sys.exit(1 if n_fail else 0)


if __name__ == "__main__":
    main()
