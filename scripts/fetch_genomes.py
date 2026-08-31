#!/usr/bin/env python3
"""Provision inputs for the funannotate benchmark (fetch_genomes.py).

For every row of samples.csv:
  1. Genomes: symlink the pre-masked genome from the Fungi_BFD pool into
     <fetch.masked_out_dir>/ if one exists (preferred path). If a masked
     genome is missing, record it as needing download + masking (the raw FNA
     can be fetched with the 1KFG helper scripts/sync_ncbi_assembly.sh and
     masked via nf_funannotate1 GENOME_CLEAN/MASKREPEAT; see conf/benchmark.yaml
     fetch.mask_via_nf_funannotate1).
  2. Reference annotations: NCBI RefSeq/GenBank GFF3 (+ protein set, cds,
     sequence report) -- used later by scripts/compare_predictions.py. Files
     land under <fetch.reference_dir>/<ASM_FOLDER>/<ASM_FOLDER>_*.

Sources for the reference annotations, in order of preference:
  a) the local 1KFG NCBI_ASM mirror
     <fetch.reference_mirror_dir>/<ASM_FOLDER>/<ASM_FOLDER>_* -- hardlinked,
     no network (mirror ships *_sequence_report.jsonl.gz).
  b) the NCBI FTP path (see /bigdata/stajichlab/shared/projects/1KFG/2026/
     NCBI_fungi/scripts/sync_ncbi_assembly.sh --method aria2c)
       genomes/all/{PRE}/{d1}/{d2}/{d3}/{ASM_FOLDER}/{ASM_FOLDER}_<suffix>
     where PRE is GCF|GCA and d1/d2/d3 are the 3-digit groups of the 9-digit
     numeric accession, e.g. GCF_000001985.1 -> GCF/000/001/985.

The sequence report is normalized to plain "<ASM>_sequence_report.jsonl"
(decompressed when pulled from the mirror, which stores it gzipped).
"""
from __future__ import annotations

import argparse
import gzip
import os
import shutil
import sys
import time
import urllib.error
import urllib.request

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lib  # noqa: E402

FTP_BASE = "https://ftp.ncbi.nih.gov/genomes/all"

# kind -> (mirror basename, ftp suffix). Mirror stores the seq report gzipped.
REFERENCE_SUFFIXES = {
    "gff": ("_genomic.gff.gz", "_genomic.gff.gz"),
    "protein": ("_protein.faa.gz", "_protein.faa.gz"),
    "cds": ("_cds_from_genomic.fna.gz", "_cds_from_genomic.fna.gz"),
    "seq-report": ("_sequence_report.jsonl.gz", "_sequence_report.jsonl"),
}


def parse_args():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--config", default="conf/benchmark.yaml")
    ap.add_argument("--samples", default="samples.csv")
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--link-only", action="store_true",
                    help="only (re-)link masked genomes; skip reference provisioning")
    ap.add_argument("--no-link", action="store_true",
                    help="skip the masked-genome symlinking")
    ap.add_argument("--verbose", "-v", action="count", default=0)
    return ap.parse_args()


def ftp_dir_for(accession: str) -> str:
    """genomes/all path components from an accession (GCF_000001985.1 -> GCF/000/001/985)."""
    pre, rest = accession.split("_", 1)
    num = rest.split(".")[0]
    return f"{pre}/{num[0:3]}/{num[3:6]}/{num[6:9]}"


def download(url: str, dest: str, retries: int = 3, sleep: float = 4.0,
             timeout: int = 120) -> bool:
    """Download url to dest (skips if exists). Returns True on success."""
    if os.path.exists(dest) and os.path.getsize(dest) > 0:
        return True
    tmp = dest + ".part"
    for attempt in range(1, retries + 1):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "funannotate-benchmark/0.1"})
            with urllib.request.urlopen(req, timeout=timeout) as resp, \
                    open(tmp, "wb") as fh:
                while True:
                    chunk = resp.read(1 << 20)
                    if not chunk:
                        break
                    fh.write(chunk)
            os.replace(tmp, dest)
            return True
        except (urllib.error.HTTPError, urllib.error.URLError, OSError) as e:
            lib.LOG.warning("download failed (%s) %s attempt %d/%d", e, url, attempt, retries)
            if attempt == retries:
                return False
            time.sleep(sleep * attempt)
    return False


def provide_reference(target_dir: str, asm_folder: str, kind: str,
                      mirror_dir: str, dest_name: str) -> str:
    """Populate <target_dir>/<dest_name> from mirror (hardlink) or FTP.

    Returns 'mirror', 'ftp', or 'missing'.
    """
    dest = os.path.join(target_dir, dest_name)
    if os.path.exists(dest) and os.path.getsize(dest) > 0:
        return "present"

    mirror_basename, ftp_suffix = REFERENCE_SUFFIXES[kind]
    mirror_src = os.path.join(mirror_dir, asm_folder, f"{asm_folder}{mirror_basename}") \
        if mirror_dir else ""
    if mirror_src and os.path.exists(mirror_src) and os.path.getsize(mirror_src) > 0:
        os.makedirs(target_dir, exist_ok=True)
        if kind == "seq-report":  # mirror stores it gzipped; normalize to plain
            with gzip.open(mirror_src, "rb") as fin, open(dest, "wb") as fout:
                shutil.copyfileobj(fin, fout)
        else:
            os.link(mirror_src, dest)
        return "mirror"

    base = f"{FTP_BASE}/{ftp_dir_for(asm_folder)}/{asm_folder}/{asm_folder}"
    if download(f"{base}{ftp_suffix}", dest):
        return "ftp"
    return "missing"


def main():
    args = parse_args()
    lib.setup_logging(args.verbose)
    cfg = lib.load_yaml(args.config)
    fetch = cfg.get("fetch", {})
    masked_out = fetch.get("masked_out_dir", "input_clean_genomes")
    reference_dir = fetch.get("reference_dir", "results/reference_annotations")
    mirror_dir = fetch.get("reference_mirror_dir", "")

    rows = lib.read_table(args.samples)
    pool = lib.load_pool(cfg["pool"]["accessions"], cfg["pool"]["taxonomy"])
    pool_masked = fetch.get("pool_masked_genome_dir") or cfg["pool"].get("masked_genome_dir", "")
    masked_accs = lib.masked_accessions(pool_masked)

    report = []
    for row in rows:
        acc = row.get("accession", "")
        asm_folder = (pool.get(acc, {}) or {}).get("ASM_FOLDER", "") or acc

        masked_name = f"{asm_folder}.masked.fasta.gz"
        pool_masked_path = os.path.join(pool_masked, masked_name) if pool_masked else ""
        has_pool = pool_masked and os.path.exists(pool_masked_path)

        notes = []
        if not args.no_link:
            if has_pool:
                os.makedirs(masked_out, exist_ok=True)
                link = os.path.join(masked_out, masked_name)
                if not os.path.exists(link):
                    os.symlink(pool_masked_path, link)
            else:
                notes.append("masked_missing")
                lib.LOG.warning("%s: no pre-masked genome in pool (%s); "
                                "needs fetch+mask", acc, masked_name)

        ref_sources = []
        if not args.link_only:
            target_dir = os.path.join(reference_dir, asm_folder)
            for kind in REFERENCE_SUFFIXES:
                dest_name = f"{asm_folder}{REFERENCE_SUFFIXES[kind][1]}"
                src = provide_reference(target_dir, asm_folder, kind, mirror_dir, dest_name)
                ref_sources.append(f"{kind}={src}")
                if src == "missing":
                    notes.append(f"{kind}_missing")
                    lib.LOG.warning("%s: reference %s unavailable", acc, kind)

        report.append({
            "accession": acc,
            "asm_folder": asm_folder,
            "masked_linked": "yes" if has_pool else "no",
            "reference_sources": ",".join(ref_sources) if ref_sources else "",
            "notes": ";".join(notes),
        })

    if args.dry_run:
        print("(dry run; nothing written)")
        return
    os.makedirs("results", exist_ok=True)
    lib.write_table("results/fetch_report.tsv", report, sep="\t")
    n_mask = sum(1 for r in report if r["masked_linked"] == "yes")
    n_mirror = sum(1 for r in report if "mirror" in r["reference_sources"])
    n_all_ref = sum(1 for r in report
                    if r["reference_sources"]
                    and all(s.endswith("=mirror") or s.endswith("=ftp") or s.endswith("=present")
                            for s in r["reference_sources"].split(",")))
    print(f"linked masked genomes: {n_mask}/{len(rows)}; "
          f"assemblies with all references: {n_all_ref}/{len(rows)} "
          f"(at least one from mirror: {n_mirror})")
    print("see results/fetch_report.tsv")


if __name__ == "__main__":
    main()
