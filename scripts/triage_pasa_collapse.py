#!/usr/bin/env python3
"""Triage production Fungi_BFD_runs genomes for the PASA duplicate-output-loop
collapse (see BFD/Funannotate_benchmarking docs/beta12-exon-collapse-and-input-integrity.md).

Population: every species in Fungi_BFD_runs/rnaseq_data with a non-empty
<species>.trinity-GG.fasta (i.e. it actually had RNA-seq evidence and would
have invoked PASA training) -- ab-initio-only genomes never ran PASA and are
not in scope for this specific bug.

For each such species, searches every do_annotation_* run directory for a
genome_annotation_training/<species>_<strain>/ dir (species can appear in more
than one run dir -- backfills, re-runs -- all matches are reported), and from
its training/ subdirectory:
  - parses funannotate_train.pasa.gff3 (final PASA/TransDecoder training set):
    n_models, mean CDS/model, %>=2exon, %>=3exon
  - flags presence/absence of stringtie.gtf, junctions.bed (the second,
    compounding gap: --stop_after_trinity exited before StringTie ran on
    every production build before today's PASA fix, per the same doc)
  - flags presence/absence of funannotate_train.coordSorted.bam (the hisat2
    BAM needed to regenerate junctions.bed via bam2juncbed() WITHOUT
    re-running Trinity -- see doc for why this makes remediation cheap)

Output: results/pasa_collapse_triage.tsv, one row per (species, run_dir) match,
sorted worst-first by mean CDS/model, so the most likely-collapsed genomes
sort to the top. This is a HEURISTIC SCREEN, not a definitive per-genome
verdict -- some genomes (e.g. yeasts) legitimately have low CDS/model with no
PASA damage at all; a "low" number here is a prioritization signal for manual
follow-up, not proof of collapse. There is no per-genome RefSeq baseline to
compare against for most of these 1,677 species, unlike the controlled
benchmark cells.

Usage:
    python3 scripts/triage_pasa_collapse.py \
        --rnaseq-dir /bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/rnaseq_data \
        --runs-root /bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs \
        --out results/pasa_collapse_triage.tsv
"""
from __future__ import annotations

import argparse
import csv
import os
import re
import sys
from collections import defaultdict
from pathlib import Path

RUN_DIR_PREFIX = "do_annotation"


def find_run_dirs(runs_root: Path) -> list[Path]:
    return sorted(
        p for p in runs_root.iterdir()
        if p.is_dir() and p.name.startswith(RUN_DIR_PREFIX)
    )


def real_trinity_species(rnaseq_dir: Path) -> list[str]:
    """Species tags with a non-empty <species>.trinity-GG.fasta."""
    species = []
    for f in rnaseq_dir.glob("*.trinity-GG.fasta"):
        try:
            if f.stat().st_size > 0:
                species.append(f.name[: -len(".trinity-GG.fasta")])
        except OSError:
            continue
    return sorted(species)


def find_training_dirs(run_dir: Path, species: str) -> list[Path]:
    """genome_annotation_training/<species>_<strain>/ -- species is a prefix,
    strain suffix varies. Matches training/ subdir must exist."""
    base = run_dir / "genome_annotation_training"
    if not base.is_dir():
        return []
    hits = []
    prefix = species + "_"
    try:
        entries = os.listdir(base)
    except OSError:
        return []
    for name in entries:
        if name == species or name.startswith(prefix):
            d = base / name / "training"
            if d.is_dir():
                hits.append(d)
    return hits


def parse_cds_per_model(gff3_path: Path) -> dict | None:
    """Parse CDS features by Parent= from a GFF3 file. Returns None if the
    file is missing or has zero CDS features."""
    if not gff3_path.is_file():
        return None
    models = defaultdict(int)
    parent_re = re.compile(r"Parent=([^;]+)")
    try:
        with open(gff3_path) as fh:
            for line in fh:
                if line.startswith("#") or not line.strip():
                    continue
                f = line.rstrip("\n").split("\t")
                if len(f) < 9 or f[2] != "CDS":
                    continue
                m = parent_re.search(f[8])
                if not m:
                    continue
                models[m.group(1).split(",")[0]] += 1
    except OSError:
        return None
    n = len(models)
    if n == 0:
        return None
    vals = list(models.values())
    return {
        "n_models": n,
        "mean_cds_per_model": sum(vals) / n,
        "pct_ge2exon": 100.0 * sum(1 for v in vals if v >= 2) / n,
        "pct_ge3exon": 100.0 * sum(1 for v in vals if v >= 3) / n,
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--rnaseq-dir", required=True, type=Path)
    ap.add_argument("--runs-root", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--flag-mean-cds-below", type=float, default=2.0,
                     help="heuristic screen threshold: mean CDS/model below this is flagged (default: 2.0)")
    args = ap.parse_args()

    run_dirs = find_run_dirs(args.runs_root)
    if not run_dirs:
        print(f"ERROR: no {RUN_DIR_PREFIX}* directories found under {args.runs_root}", file=sys.stderr)
        sys.exit(1)
    print(f"[INFO] scanning {len(run_dirs)} run directories: {', '.join(d.name for d in run_dirs)}", file=sys.stderr)

    species_list = real_trinity_species(args.rnaseq_dir)
    print(f"[INFO] {len(species_list)} species have a non-empty Trinity assembly", file=sys.stderr)

    rows = []
    n_no_training_dir = 0
    n_no_pasa_gff3 = 0
    for i, species in enumerate(species_list, 1):
        if i % 200 == 0:
            print(f"[INFO] ... {i}/{len(species_list)}", file=sys.stderr)
        matched_any = False
        for run_dir in run_dirs:
            for training_dir in find_training_dirs(run_dir, species):
                matched_any = True
                gff3 = training_dir / "funannotate_train.pasa.gff3"
                metrics = parse_cds_per_model(gff3)
                has_stringtie = (training_dir / "funannotate_train.stringtie.gtf").is_file()
                has_junctions = any(training_dir.glob("*junction*"))
                has_bam = (training_dir / "funannotate_train.coordSorted.bam").is_file()
                if metrics is None:
                    n_no_pasa_gff3 += 1
                    rows.append({
                        "species": species,
                        "run_dir": run_dir.name,
                        "training_dir": str(training_dir),
                        "n_models": "",
                        "mean_cds_per_model": "",
                        "pct_ge2exon": "",
                        "pct_ge3exon": "",
                        "flagged": "NO_PASA_GFF3",
                        "has_stringtie_gtf": has_stringtie,
                        "has_junctions_bed": has_junctions,
                        "has_coordsorted_bam": has_bam,
                    })
                    continue
                flagged = metrics["mean_cds_per_model"] < args.flag_mean_cds_below
                rows.append({
                    "species": species,
                    "run_dir": run_dir.name,
                    "training_dir": str(training_dir),
                    "n_models": metrics["n_models"],
                    "mean_cds_per_model": round(metrics["mean_cds_per_model"], 3),
                    "pct_ge2exon": round(metrics["pct_ge2exon"], 1),
                    "pct_ge3exon": round(metrics["pct_ge3exon"], 1),
                    "flagged": "LOW_CDS_PER_MODEL" if flagged else "",
                    "has_stringtie_gtf": has_stringtie,
                    "has_junctions_bed": has_junctions,
                    "has_coordsorted_bam": has_bam,
                })
        if not matched_any:
            n_no_training_dir += 1
            rows.append({
                "species": species,
                "run_dir": "",
                "training_dir": "",
                "n_models": "",
                "mean_cds_per_model": "",
                "pct_ge2exon": "",
                "pct_ge3exon": "",
                "flagged": "NO_TRAINING_DIR_FOUND",
                "has_stringtie_gtf": "",
                "has_junctions_bed": "",
                "has_coordsorted_bam": "",
            })

    # Sort: real numeric rows worst (lowest mean_cds_per_model) first, then
    # missing-data rows grouped at the end.
    def sort_key(r):
        if r["mean_cds_per_model"] == "":
            return (1, 0.0)
        return (0, r["mean_cds_per_model"])

    rows.sort(key=sort_key)

    args.out.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "species", "run_dir", "training_dir", "n_models", "mean_cds_per_model",
        "pct_ge2exon", "pct_ge3exon", "flagged",
        "has_stringtie_gtf", "has_junctions_bed", "has_coordsorted_bam",
    ]
    with open(args.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    n_flagged_low = sum(1 for r in rows if r["flagged"] == "LOW_CDS_PER_MODEL")
    n_missing_stringtie = sum(1 for r in rows if r["has_stringtie_gtf"] is False)
    print(f"[INFO] wrote {len(rows)} rows -> {args.out}", file=sys.stderr)
    print(f"[SUMMARY] species with no training dir found in any run: {n_no_training_dir}", file=sys.stderr)
    print(f"[SUMMARY] matches with a training dir but no funannotate_train.pasa.gff3: {n_no_pasa_gff3}", file=sys.stderr)
    print(f"[SUMMARY] flagged LOW_CDS_PER_MODEL (mean CDS/model < {args.flag_mean_cds_below}): {n_flagged_low}", file=sys.stderr)
    print(f"[SUMMARY] missing stringtie.gtf: {n_missing_stringtie}", file=sys.stderr)


if __name__ == "__main__":
    main()
