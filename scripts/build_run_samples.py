#!/usr/bin/env python3
"""Generate per-cell nf_funannotate1 samplesheets (build_run_samples.py).

Reads the append-only benchmark samples.csv plus conf/cells.tsv and writes
one nf_funannotate1 `samples.csv` per cell under <runs-dir>/<cell>/:

  SPECIES,STRAIN,ASMID,LOCUSTAG,BUSCO_LINEAGE,TRANSL_TABLE,NCBI_TAXONID,GENOME

Mapping (nf_funannotate1 contract, see subworkflows/local/input_check.nf and
lib/SampleUtils.groovy):
  SPECIES          <- samples.csv species
  STRAIN           <- samples.csv strain
  ASMID            <- samples.csv accession (GCF_/GCA_; also the output key)
  LOCUSTAG         <- collision-free 6-char prefix: 3 chars genus + 3 chars
                      species, uppercased, alnum only; on collision the genus
                      part is extended until unique (then digit-suffixed)
  BUSCO_LINEAGE    <- samples.csv busco_lineage_used (e.g. sordariomycetes_odb10)
  TRANSL_TABLE     <- 1 (fungal nuclear genomes; adjust here if an exception
                      sample needs another table)
  NCBI_TAXONID     <- samples.csv taxid
  GENOME           <- absolute path to input_clean_genomes/<ASM_FOLDER>.masked.fasta.gz

The same accession/ASMID and species_strain id are used across all cells so a
genome's results are directly comparable cell-to-cell.
"""
from __future__ import annotations

import argparse
import csv
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lib  # noqa: E402

FIELDS = ["SPECIES", "STRAIN", "ASMID", "LOCUSTAG", "BUSCO_LINEAGE",
          "TRANSL_TABLE", "NCBI_TAXONID", "GENOME"]


def parse_args():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--config", default="conf/benchmark.yaml")
    ap.add_argument("--samples", default="samples.csv")
    ap.add_argument("--cells", default="conf/cells.tsv")
    ap.add_argument("--runs-dir", default="runs")
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--verbose", "-v", action="count", default=0)
    return ap.parse_args()


def locustag_candidates(species: str) -> str:
    """6-char base prefix from species: 3 chars genus + 3 chars epithet."""
    parts = re.split(r"\s+", species.strip())
    genus = re.sub(r"[^A-Za-z]", "", parts[0]) if parts else ""
    epithet = re.sub(r"[^A-Za-z]", "", parts[1]) if len(parts) > 1 else ""
    return (genus[:3] + epithet[:3]).upper()


def unique_locustags(rows):
    """Assign LOCUSTAG per row, resolving collisions with longer genus prefixes."""
    base = [locustag_candidates(r.get("species", "")) for r in rows]
    used = {}
    tags = [""] * len(rows)
    for i, b in enumerate(base):
        if not b:
            tags[i] = "NULL" + str(i).zfill(2)
            used[tags[i]] = tags[i]
            continue
        candidate = b
        genus = re.sub(r"[^A-Za-z]", "", re.split(r"\s+", rows[i].get("species", ""))[0]).upper()
        k = 6
        while candidate in used or candidate != used.get(candidate, candidate):
            k += 1
            candidate = (genus[:k] + base[i][3:]).upper()
            if k >= len(genus) + 3:
                candidate = (base[i] + str(i)).upper()[:8]
                if candidate in used:
                    candidate += str(i)
        tags[i] = candidate
        used[candidate] = candidate
    return tags


def main():
    args = parse_args()
    lib.setup_logging(args.verbose)
    cfg = lib.load_yaml(args.config)
    fetch = cfg.get("fetch", {})
    masked_out = fetch.get("masked_out_dir", "input_clean_genomes")

    rows = lib.read_table(args.samples)
    cells = [c for c in lib.read_table(args.cells) if c.get("cell", "").strip()]
    tags = unique_locustags(rows)

    written = []
    for cell in cells:
        cell_name = cell["cell"]
        out_dir = os.path.join(args.runs_dir, cell_name)
        out_csv = os.path.join(out_dir, "samples.csv")
        if not args.dry_run:
            os.makedirs(out_dir, exist_ok=True)
        out_rows = []
        missing = []
        for r, tag in zip(rows, tags):
            acc = r.get("accession", "")
            asm_folder = f"{acc}_{r.get('asm_name', '')}".rstrip("_")
            genome = os.path.abspath(os.path.join(masked_out, f"{asm_folder}.masked.fasta.gz"))
            if not os.path.exists(genome):
                missing.append(genome)
            out_rows.append({
                "SPECIES": r.get("species", "").strip(),
                "STRAIN": r.get("strain", "").strip(),
                "ASMID": acc,
                "LOCUSTAG": tag,
                "BUSCO_LINEAGE": r.get("busco_lineage_used", "").strip(),
                "TRANSL_TABLE": "1",
                "NCBI_TAXONID": r.get("taxid", "").strip(),
                "GENOME": genome,
            })
        if missing:
            for m in missing:
                lib.LOG.warning("missing masked genome for %s (under %s)",
                                os.path.basename(m), masked_out)
        if not args.dry_run:
            with open(out_csv, "w", newline="") as fh:
                w = csv.DictWriter(fh, fieldnames=FIELDS)
                w.writeheader()
                w.writerows(out_rows)
        written.append((cell_name, out_csv, len(out_rows), len(missing)))

    for cell_name, out_csv, n, n_missing in written:
        print(f"{cell_name:<32} {n} samples -> {out_csv}"
              + (f"  ({n_missing} missing genome)" if n_missing else ""))


if __name__ == "__main__":
    main()
