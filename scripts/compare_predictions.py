#!/usr/bin/env python3
"""Gene-content comparison vs RefSeq (DESIGN.md "Gene-content comparison").

For every (cell, genome) with both a funannotate predict_results/*.gff3
output and a fetched RefSeq GFF (results/reference_annotations/, from
scripts/fetch_genomes.py), computes:
  * gffcompare sensitivity/precision (transcript/exon/intron level) and
    missed/novel loci-exons-introns percentages
  * a bedtools-based coordinate-overlap tally at the gene level (catches
    same-locus evidence-source swaps that BUSCO/count deltas alone would
    hide — the funannotate-abinitio-reuse-validation.md failure mode)
  * BUSCO short_summary completeness, if a short_summary*.txt is found
    anywhere under the genome's output directory
  * raw gene counts (query vs reference)

Writes results/gene_content_comparison.tsv, one row per (cell, genome), each
annotated with refseq_annotation_release from samples.csv — RefSeq GFFs from
different NCBI pipeline eras (PGAP vs legacy) are not treated as one uniform
gold standard, so keep that column around for any downstream slicing.

Requires `gffcompare` and `bedtools` on PATH (UCR HPCC: `module load
gffcompare bedtools`). Either is skipped with a warning (not fatal) if
missing, and genomes with no predict_results/*.gff3 yet or no RefSeq GFF are
skipped and counted, not silently dropped.

Usage:
  module load gffcompare bedtools
  python3 scripts/compare_predictions.py --cells-dir runs --samples samples.csv
  python3 scripts/compare_predictions.py --cells v1.8.17_conda v1.8.17_container
"""
from __future__ import annotations

import argparse
import os
import sys
import tempfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lib  # noqa: E402


def iter_cell_genomes(cells_dir: str, cell: str):
    target = os.path.join(cells_dir, cell, "genome_annotation")
    if not os.path.isdir(target):
        return
    for name in sorted(os.listdir(target)):
        gdir = os.path.join(target, name)
        if os.path.isdir(gdir):
            yield name, gdir


def compare_one(cell: str, genome: str, gdir: str, srow: dict,
                 reference_dir: str, cache_dir: str, tmp: str) -> dict:
    row = {
        "cell": cell, "genome": genome,
        "accession": srow.get("accession", ""),
        "phylum": srow.get("phylum", ""),
        "refseq_annotation_release": srow.get("refseq_annotation_release", ""),
    }
    query_gff = lib.find_predict_gff3(gdir)
    row["query_gff3"] = query_gff or ""
    if not query_gff:
        row["compare_status"] = "no_predict_output"
        return row
    row["query_gene_count"] = lib.count_gff3_genes(query_gff)

    busco_path = lib.find_busco_summary(gdir)
    if busco_path:
        row.update(lib.parse_busco_summary(busco_path))

    asm_folder = lib.asm_folder_for(srow)
    ref_gz = lib.find_reference_gff(reference_dir, asm_folder)
    if not ref_gz:
        row["compare_status"] = "no_reference_gff"
        return row

    ref_gff = lib.gunzip_to_cache(ref_gz, cache_dir)
    row["ref_gff3"] = ref_gff
    row["ref_gene_count"] = lib.count_gff3_genes(ref_gff)

    # gffcompare's own `-o` handling silently drops the ".stats" suffix
    # whenever the prefix already contains a "." (confirmed: "v1.8.17_conda"
    # cell names produce a stats file named EXACTLY the given prefix, no
    # extension, while .loci/.tracking/.annotated.gtf still get suffixed
    # normally) -- sanitize dots out of the prefix so run_gffcompare's
    # `f"{out_prefix}.stats"` lookup actually matches what gets written.
    safe_prefix = f"{cell}__{genome}".replace(".", "_")
    out_prefix = os.path.join(tmp, safe_prefix)
    gffc = lib.run_gffcompare(query_gff, ref_gff, out_prefix)
    if gffc:
        row.update(gffc)

    q_bed = os.path.join(tmp, f"{cell}__{genome}.query.bed")
    r_bed = os.path.join(tmp, f"{cell}__{genome}.ref.bed")
    lib.gff3_gene_bed(query_gff, q_bed)
    lib.gff3_gene_bed(ref_gff, r_bed)
    overlap = lib.run_bedtools_overlap(q_bed, r_bed)
    if overlap:
        row.update({f"overlap_{k}": v for k, v in overlap.items()})

    row["compare_status"] = "compared"
    return row


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--cells-dir", default="runs")
    ap.add_argument("--cells", nargs="*", default=None, help="restrict to these cell names (default: all under --cells-dir)")
    ap.add_argument("--samples", default="samples.csv")
    ap.add_argument("--reference-dir", default="results/reference_annotations")
    ap.add_argument("--out", default="results/gene_content_comparison.tsv")
    ap.add_argument("--cache-dir", default=None, help="decompressed-GFF cache (default: <reference-dir>/.cache)")
    ap.add_argument("-v", "--verbose", action="count", default=0)
    args = ap.parse_args()
    lib.setup_logging(args.verbose)

    if not os.path.exists(args.samples):
        lib.LOG.error("samples.csv not found: %s", args.samples)
        sys.exit(1)
    samples_idx = lib.load_samples_index(args.samples)
    cache_dir = args.cache_dir or os.path.join(args.reference_dir, ".cache")

    cells = args.cells or sorted(
        n for n in os.listdir(args.cells_dir)
        if os.path.isdir(os.path.join(args.cells_dir, n, "genome_annotation"))
    )
    if not cells:
        lib.LOG.error("no cells with a genome_annotation/ dir found under %s", args.cells_dir)
        sys.exit(1)
    lib.LOG.info("comparing cells: %s", ", ".join(cells))

    rows = []
    counts = {"compared": 0, "no_predict_output": 0, "no_reference_gff": 0, "no_samples_row": 0}
    with tempfile.TemporaryDirectory(prefix="compare_predictions_") as tmp:
        for cell in cells:
            for genome, gdir in iter_cell_genomes(args.cells_dir, cell):
                srow = samples_idx.get(genome)
                if not srow:
                    lib.LOG.warning("no samples.csv row for genome tag %r (cell %s) — skipping", genome, cell)
                    counts["no_samples_row"] += 1
                    rows.append({"cell": cell, "genome": genome, "compare_status": "no_samples_row"})
                    continue
                row = compare_one(cell, genome, gdir, srow, args.reference_dir, cache_dir, tmp)
                counts[row["compare_status"]] = counts.get(row["compare_status"], 0) + 1
                rows.append(row)

    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    lib.write_table(args.out, rows, sep="\t")
    lib.LOG.info("wrote %d rows -> %s", len(rows), args.out)
    print(f"compared={counts['compared']} no_predict_output={counts['no_predict_output']} "
          f"no_reference_gff={counts['no_reference_gff']} no_samples_row={counts['no_samples_row']}",
          file=sys.stderr)


if __name__ == "__main__":
    main()
