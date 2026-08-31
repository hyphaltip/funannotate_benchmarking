#!/usr/bin/env python3
"""Shared helpers for the Funannotate benchmarking harness.

Conventions:
  * Every script is stdlib-first; requests/pandas/pyyaml are used when present
    and degrade gracefully otherwise.
  * Paths to the local Fungi_BFD candidate pool come from conf/benchmark.yaml
    (see scripts/select_genomes.py --config).
"""
from __future__ import annotations

import csv
import io
import json
import logging
import os
import re
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from typing import Any, Dict, List, Optional

LOG = logging.getLogger("funannotate_bench")

DATASETS_API = "https://api.ncbi.nlm.nih.gov/datasets/v2alpha"
EUTILS_API = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"


def setup_logging(verbose: int = 0, stream=None) -> None:
    level = logging.WARNING
    if verbose == 1:
        level = logging.INFO
    elif verbose >= 2:
        level = logging.DEBUG
    logging.basicConfig(level=level, format="%(levelname)s %(name)s: %(message)s",
                        stream=stream or sys.stderr)


# ── generic csv/tsv ────────────────────────────────────────────────────────
def read_table(path: str, sep: Optional[str] = None) -> List[Dict[str, str]]:
    """Read a CSV/TSV into a list of dicts (values stripped, '' -> '').

    Lines starting with '#' (comments) are ignored; the header is the first
    non-comment line.
    """
    if not os.path.exists(path):
        raise FileNotFoundError(f"table not found: {path}")
    with open(path, newline="") as fh:
        raw = [ln for ln in fh if not ln.lstrip().startswith("#")]
        if not raw:
            return []
        sample = "".join(raw[:8])
        if sep is None:
            dialect = csv.Sniffer().sniff(sample, delimiters=",\t;")
            sep = dialect.delimiter
        reader = csv.DictReader(raw, delimiter=sep)
        rows = []
        for row in reader:
            rows.append({k: (v or "").strip() for k, v in row.items()})
    return rows


def write_table(path: str, rows: List[Dict[str, str]], sep: str = ",") -> None:
    if not rows:
        LOG.warning("no rows to write for %s", path)
        return
    cols = list(rows[0].keys())
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=cols, delimiter=sep,
                                extrasaction="ignore")
        writer.writeheader()
        for r in rows:
            writer.writerow({k: r.get(k, "") for k in cols})


# ── yaml (stdlib fallback) ─────────────────────────────────────────────────
def load_yaml(path: str) -> Dict[str, Any]:
    try:
        import yaml  # type: ignore
        with open(path) as fh:
            return yaml.safe_load(fh) or {}
    except ImportError:
        # Minimal fallback for the flat subset of yaml we use.
        out: Dict[str, Any] = {}
        with open(path) as fh:
            for line in fh:
                line = line.split("#", 1)[0].rstrip()
                if not line.strip() or line[0].isspace():
                    continue
                if ":" in line:
                    k, _, v = line.partition(":")
                    out[k.strip()] = v.strip()
        return out


# ── candidate pool (Fungi_BFD accession + taxonomy tables) ─────────────────
def load_pool(accessions_csv: str, taxonomy_csv: str) -> Dict[str, Dict[str, str]]:
    """Join ncbi_accessions.csv with ncbi_accessions_taxonomy.csv on ASM_FOLDER.

    Returns {accession: row} where row carries both accession-level fields
    (strain, bioproject, ASM_LENGTH, N50, …) and taxonomy fields (phylum,
    class, order, genus, species, subphylum when present).
    """
    taxon = {}
    for row in read_table(taxonomy_csv):
        taxon[row.get("ASM_ACCESSION", "")] = row

    pool = {}
    for row in read_table(accessions_csv):
        acc = row.get("ACCESSION", "")
        asm_folder = row.get("ASM_FOLDER", "")
        t = taxon.get(asm_folder) or taxon.get(acc) or {}
        merged = dict(row)
        merged.update({k: t.get(k, "") for k in
                       ("PHYLUM", "SUBPHYLUM", "CLASS", "SUBCLASS", "ORDER",
                        "FAMILY", "GENUS", "SPECIES")})
        merged["taxid"] = merged.get("NCBI_TAXID", "") or merged.get("NCBI_TAXONID", "")
        pool[acc] = merged
    return pool


def load_subphylum_map(path: str) -> Dict[str, Dict[str, str]]:
    """class -> {subphylum, phylum} from conf/subphylum_map.tsv."""
    out = {}
    for row in read_table(path):
        out[row["class"]] = {"subphylum": row["subphylum"], "phylum": row["phylum"]}
    return out


def load_lineage_map(path: str) -> Dict[str, str]:
    out = {}
    for row in read_table(path):
        out[row["class"]] = row["lineage"]
    return out


def resolve_subphylum(row: Dict[str, str],
                      class_map: Optional[Dict[str, Dict[str, str]]] = None) -> str:
    sp = (row.get("SUBPHYLUM") or "").strip()
    if sp:
        return sp
    if class_map and row.get("CLASS"):
        return class_map.get(row["CLASS"], {}).get("subphylum", "")
    return ""


def resolve_lineage(phylum: str, subphylum: str, klass: str,
                    lineage_map: Dict[str, str]) -> str:
    for key in (klass, subphylum, phylum, "DEFAULT"):
        if key and key in lineage_map:
            return lineage_map[key]
    return "fungi_odb10"


def load_rnaseq(path: str) -> Dict[str, List[str]]:
    """species_tag (underscored species name) -> list of SRA accessions.

    From Fungi_BFD's samples.rnaseq_sra.csv (species_tag,taxonid,sra_accession,...).
    """
    out: Dict[str, List[str]] = {}
    for row in read_table(path):
        tag = (row.get("species_tag") or "").strip()
        sra = (row.get("sra_accession") or "").strip()
        if tag and sra:
            out.setdefault(tag, []).append(sra)
    return out


def masked_accessions(masked_dir: str) -> set:
    """Set of ASM_FOLDER base names that have a masked genome file.

    Uses the Fungi_BFD naming <ASM_FOLDER>.masked.fasta.gz where ASM_FOLDER is
    <ACCESSION>_<ASM_NAME> (e.g. GCF_010015735.1_Aaoar1.masked.fasta.gz).
    """
    if not os.path.isdir(masked_dir):
        return set()
    out = set()
    suffixes = (".masked.fasta.gz", ".fa.gz", ".fna.gz", ".fasta.gz")
    try:
        for n in os.listdir(masked_dir):
            for suf in suffixes:
                if n.endswith(suf):
                    out.add(n[: -len(suf)])
                    break
    except OSError:
        pass
    return out


def masked_genome_for(masked_dir: str, accession: str) -> Optional[str]:
    """Find an already-masked genome file for an accession in masked_dir.

    Looks for <acc>_<something>.masked.fasta.gz (Fungi_BFD naming), falling
    back to <acc>_<something>.fa.gz. Returns the absolute filename or None.
    """
    if not os.path.isdir(masked_dir):
        return None
    pref = f"{accession}_"
    try:
        names = os.listdir(masked_dir)
    except OSError:
        return None
    hits = sorted(n for n in names if n.startswith(pref))
    for n in hits:
        if n.endswith(".masked.fasta.gz"):
            return os.path.join(masked_dir, n)
    for n in hits:
        if n.endswith(".fa.gz") or n.endswith(".fna.gz") or n.endswith(".fasta.gz"):
            return os.path.join(masked_dir, n)
    return None


# ── NCBI datasets API ──────────────────────────────────────────────────────
def _http_get_json(url: str, retries: int = 3, sleep: float = 3.0) -> Optional[Any]:
    for attempt in range(1, retries + 1):
        try:
            req = urllib.request.Request(url, headers={"Accept": "application/json",
                                                       "User-Agent": "funannotate-benchmark/0.1"})
            with urllib.request.urlopen(req, timeout=60) as resp:
                data = json.loads(resp.read().decode())
                return data if data else None
        except (urllib.error.HTTPError, urllib.error.URLError, json.JSONDecodeError) as e:
            LOG.warning("http error %s on %s (attempt %d/%d)", e, url, attempt, retries)
            if attempt == retries:
                return None
            time.sleep(sleep * attempt)
    return None


def dataset_report(accession: str, ncbi_email: str = "", retries: int = 3) -> Optional[Dict[str, Any]]:
    """GET /genome/accession/{acc}/dataset_report. Returns the report dict or None."""
    url = f"{DATASETS_API}/genome/accession/{accession}/dataset_report"
    if ncbi_email:
        sep = "&" if "?" in url else "?"
        url += f"{sep}tool=funannotate-benchmark&email={urllib.parse.quote(ncbi_email)}"
    d = _http_get_json(url, retries=retries)
    if not d or "reports" not in d or not d["reports"]:
        return None
    return d["reports"][0]


def enrich_assembly_metadata(accession: str, ncbi_email: str = "",
                             retries: int = 2, sleep: float = 3.0) -> Dict[str, str]:
    """Extract DESIGN metadata fields for one accession via the datasets API.

    Returns a dict of extra columns; every field is a str and missing/unknown
    values are ''. refseq_annotation_release is inferred from the NCBI
    annotation naming convention: an 'RS_' release name marks PGAP's automated
    RefSeq re-annotation; a submitted/propagated annotation marks 'legacy'.

    The datasets API occasionally rate-limits with HTTP 200 + empty payload;
    empty results are retried (they are indistinguishable from a genuinely
    absent accession, but retrying with backoff makes the filled rate much
    higher). "Truly" absent accessions keep returning empty across retries.
    """
    out = {"assembly_level": "", "refseq_annotation_release": "",
           "ncbi_current_accession": "", "paired_accession": "",
           "assembly_status": "", "annotation_name": "", "annotation_date": "",
           "ncbi_species": "", "ncbi_genome_size_bp": ""}
    for _ in range(max(1, retries)):
        r = dataset_report(accession, ncbi_email=ncbi_email, retries=2)
        if r:
            ai = r.get("assembly_info", {})
            out["assembly_level"] = ai.get("assembly_level") or ""
            out["assembly_status"] = ai.get("assembly_status") or ""
            out["ncbi_current_accession"] = r.get("current_accession") or ""
            out["paired_accession"] = r.get("paired_accession") or ""
            out["ncbi_species"] = (r.get("organism") or {}).get("organism_name") or ""
            stats = r.get("assembly_stats", {}) or {}
            out["ncbi_genome_size_bp"] = str(stats.get("total_sequence_length") or "")

            ann = r.get("annotation_info") or {}
            name = ann.get("name") or ""
            provider = ann.get("provider") or ""
            out["annotation_name"] = name
            out["annotation_date"] = ann.get("release_date") or ""
            if name.startswith("Annotation submitted") or (name and "RS_" not in name):
                out["refseq_annotation_release"] = "legacy"
            elif name:
                out["refseq_annotation_release"] = "PGAP"
            elif provider == "NCBI":
                out["refseq_annotation_release"] = "PGAP"
            if out["assembly_level"] or out["refseq_annotation_release"]:
                break  # got a real payload
            if out["ncbi_current_accession"] or out["paired_accession"]:
                break
        time.sleep(sleep)
    return out


# ── eutils (rnaseq SRA presence is decided locally; this is for future use) ─
def esearch_sra_bioproject(bioproject: str, ncbi_email: str = "") -> int:
    q = urllib.parse.quote(f'"{bioproject}"[BioProject]')
    url = f"{EUTILS_API}/esearch.fcgi?db=sra&term={q}&retmode=json"
    if ncbi_email:
        url += f"&tool=funannotate-benchmark&email={urllib.parse.quote(ncbi_email)}"
    d = _http_get_json(url, retries=2)
    if not d:
        return 0
    return int((d.get("esearchresult") or {}).get("count") or 0)


def mb(size_bp_str: str) -> str:
    """Format a bp string as Mb (1 decimal). Returns '' for bad input."""
    try:
        return f"{int(size_bp_str) / 1e6:.1f}"
    except (ValueError, TypeError):
        return ""
