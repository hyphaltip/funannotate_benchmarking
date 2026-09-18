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
import gzip
import io
import json
import logging
import os
import re
import shutil
import subprocess
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
            try:
                dialect = csv.Sniffer().sniff(sample, delimiters=",\t;")
                sep = dialect.delimiter
            except csv.Error:
                # Sniffer chokes on TSVs with a ragged trailing column (rows
                # that omit trailing optional fields); fall back to the
                # extension, then to whichever delimiter the header uses.
                header = raw[0]
                if path.endswith(".tsv"):
                    sep = "\t"
                elif "\t" in header:
                    sep = "\t"
                elif ";" in header:
                    sep = ";"
                else:
                    sep = ","
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


# ── Nextflow trace parsing ──────────────────────────────────────────────────
# conf/profile_annotate.config pins trace.fields to:
#   task_id,hash,name,status,exit,realtime,%cpu,rss,tag
# i.e. no submit/complete timestamps, so per-genome "wall clock" below is
# really a sum of task realtime (tasks for one genome run concurrently across
# processes, so this overstates true wall time — it's still the right unit to
# compare processes/cells/genomes against each other).
_DURATION_RE = re.compile(r"(\d+(?:\.\d+)?)\s*(ms|s|m|h|d)\b")
_DURATION_UNIT_SECONDS = {"ms": 0.001, "s": 1.0, "m": 60.0, "h": 3600.0, "d": 86400.0}


def parse_nf_duration(s: str) -> Optional[float]:
    """Nextflow duration string ('1h 2m 3.4s', '59.9s', '3ms') -> seconds."""
    if not s or s in ("-", "NA"):
        return None
    total = 0.0
    matched = False
    for val, unit in _DURATION_RE.findall(s):
        total += float(val) * _DURATION_UNIT_SECONDS[unit]
        matched = True
    return total if matched else None


_MEM_RE = re.compile(r"([\d.]+)\s*([KMGT]?B)")
_MEM_UNIT_GB = {"B": 1e-9, "KB": 1e-6, "MB": 1e-3, "GB": 1.0, "TB": 1e3}


def parse_mem_to_gb(s: str) -> Optional[float]:
    if not s or s in ("-", "NA"):
        return None
    m = _MEM_RE.search(s)
    if not m:
        return None
    val, unit = float(m.group(1)), m.group(2)
    return val * _MEM_UNIT_GB.get(unit, 1.0)


def parse_pct(s: str) -> Optional[float]:
    if not s or s in ("-", "NA"):
        return None
    try:
        return float(s.replace("%", "").strip())
    except ValueError:
        return None


def read_trace(path: str) -> List[Dict[str, str]]:
    return read_table(path, sep="\t")


# ── genome id <-> samples.csv row ───────────────────────────────────────────
def species_strain_tag(species: str, strain: str) -> str:
    """SPECIES+STRAIN -> the genome_annotation/<tag> directory name nf_funannotate1
    uses (build_run_samples.py / funannotate.nf): whitespace runs -> '_',
    everything else (including '.', '-') passed through unchanged.
    """
    raw = f"{species or ''} {strain or ''}".strip()
    return re.sub(r"\s+", "_", raw)


def load_samples_index(samples_csv: str) -> Dict[str, Dict[str, str]]:
    """genome tag (as it appears under runs/<cell>/genome_annotation/) -> row.

    Also indexed by 'acc:<accession>' for accession-keyed lookups.
    """
    idx: Dict[str, Dict[str, str]] = {}
    for row in read_table(samples_csv):
        tag = species_strain_tag(row.get("species", ""), row.get("strain", ""))
        if tag:
            idx[tag] = row
        acc = row.get("accession", "")
        if acc:
            idx.setdefault(f"acc:{acc}", row)
    return idx


def asm_folder_for(row: Dict[str, str]) -> str:
    """'<accession>_<asm_name>' — matches results/reference_annotations/<asm_folder>/."""
    acc = row.get("accession", "")
    asm = row.get("asm_name", "")
    return f"{acc}_{asm}" if asm else acc


# ── funannotate / reference GFF3 helpers ────────────────────────────────────
def find_predict_gff3(genome_dir: str) -> Optional[str]:
    """The funannotate predict output GFF3 under <genome_dir>/predict_results/."""
    pr = os.path.join(genome_dir, "predict_results")
    if not os.path.isdir(pr):
        return None
    hits = sorted(f for f in os.listdir(pr) if f.endswith(".gff3"))
    return os.path.join(pr, hits[0]) if hits else None


def find_reference_gff(reference_dir: str, asm_folder: str) -> Optional[str]:
    for suffix in ("_genomic.gff.gz", "_genomic.gff"):
        p = os.path.join(reference_dir, asm_folder, f"{asm_folder}{suffix}")
        if os.path.exists(p):
            return p
    return None


def gunzip_to_cache(path: str, cache_dir: str) -> str:
    """Return an uncompressed copy of a .gz file cached under cache_dir."""
    if not path.endswith(".gz"):
        return path
    os.makedirs(cache_dir, exist_ok=True)
    out = os.path.join(cache_dir, os.path.basename(path)[:-3])
    if not os.path.exists(out) or os.path.getmtime(out) < os.path.getmtime(path):
        with gzip.open(path, "rb") as fi, open(out, "wb") as fo:
            shutil.copyfileobj(fi, fo)
    return out


def _gff3_open(path: str):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def count_gff3_genes(gff3_path: str) -> int:
    n = 0
    with _gff3_open(gff3_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) >= 3 and cols[2] == "gene":
                n += 1
    return n


_GFF3_ID_RE = re.compile(r"ID=([^;\n]+)")


def gff3_gene_bed(gff3_path: str, out_bed: str) -> int:
    """Write a 6-col BED (0-based half-open) of 'gene' features. Returns count."""
    n = 0
    with _gff3_open(gff3_path) as fh, open(out_bed, "w") as out:
        for line in fh:
            if line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9 or cols[2] != "gene":
                continue
            chrom, start, end, strand, attrs = cols[0], cols[3], cols[4], cols[6], cols[8]
            m = _GFF3_ID_RE.search(attrs)
            gid = m.group(1) if m else f"gene_{n}"
            out.write(f"{chrom}\t{int(start) - 1}\t{end}\t{gid}\t.\t{strand}\n")
            n += 1
    return n


# ── BUSCO short_summary ─────────────────────────────────────────────────────
_BUSCO_RE = re.compile(
    r"C:([\d.]+)%\[S:([\d.]+)%,D:([\d.]+)%\],F:([\d.]+)%,M:([\d.]+)%,n:(\d+)"
)


def find_busco_summary(genome_dir: str) -> Optional[str]:
    for root, _dirs, files in os.walk(genome_dir):
        for f in files:
            if f.startswith("short_summary") and f.endswith(".txt"):
                return os.path.join(root, f)
    return None


def parse_busco_summary(path: str) -> Dict[str, str]:
    with open(path) as fh:
        text = fh.read()
    m = _BUSCO_RE.search(text)
    if not m:
        return {}
    keys = ("busco_complete_pct", "busco_single_pct", "busco_duplicated_pct",
            "busco_fragmented_pct", "busco_missing_pct", "busco_n")
    return dict(zip(keys, m.groups()))


# ── gffcompare / bedtools (external tools; module load gffcompare bedtools) ─
_GFFC_LEVEL_RE = re.compile(r"^\s*([A-Za-z][\w /-]*level):\s*([\d.]+)\s*\|\s*([\d.]+)")
_GFFC_MISSED_NOVEL_RE = re.compile(r"^\s*(Missed|Novel) (exons|introns|loci).*?:\s*(\d+)/(\d+)\s*\(\s*([\d.]+)%\)")


def parse_gffcompare_stats(stats_path: str) -> Dict[str, str]:
    out: Dict[str, str] = {}
    with open(stats_path) as fh:
        for line in fh:
            m = _GFFC_LEVEL_RE.match(line)
            if m:
                level, sens, prec = m.groups()
                key = re.sub(r"[^a-z0-9]+", "_", level.strip().lower()).strip("_")
                out[f"gffcompare_{key}_sensitivity"] = sens
                out[f"gffcompare_{key}_precision"] = prec
                continue
            m2 = _GFFC_MISSED_NOVEL_RE.match(line)
            if m2:
                kind, feat, _num, _denom, pct = m2.groups()
                out[f"gffcompare_{kind.lower()}_{feat}_pct"] = pct
    return out


def run_gffcompare(query_gff: str, ref_gff: str, out_prefix: str) -> Optional[Dict[str, str]]:
    if shutil.which("gffcompare") is None:
        LOG.warning("gffcompare not on PATH (module load gffcompare) — skipping accuracy metrics")
        return None
    os.makedirs(os.path.dirname(out_prefix) or ".", exist_ok=True)
    cmd = ["gffcompare", "-r", ref_gff, "-o", out_prefix, query_gff]
    try:
        subprocess.run(cmd, check=True, capture_output=True, text=True, timeout=900)
    except (subprocess.CalledProcessError, subprocess.TimeoutExpired) as e:
        LOG.warning("gffcompare failed for %s vs %s: %s", query_gff, ref_gff, e)
        return None
    stats_path = f"{out_prefix}.stats"
    if not os.path.exists(stats_path):
        LOG.warning("gffcompare produced no stats file for %s", query_gff)
        return None
    return parse_gffcompare_stats(stats_path)


def run_bedtools_overlap(query_bed: str, ref_bed: str) -> Optional[Dict[str, int]]:
    """Coordinate-aware gene-overlap tally between two BEDs (see gff3_gene_bed).

    Catches same-locus evidence-source swaps that BUSCO/count deltas alone
    would hide (funannotate-abinitio-reuse-validation.md failure mode).
    """
    if shutil.which("bedtools") is None:
        LOG.warning("bedtools not on PATH (module load bedtools) — skipping overlap diff")
        return None
    try:
        q_sorted, r_sorted = query_bed + ".sorted", ref_bed + ".sorted"
        for src, dst in ((query_bed, q_sorted), (ref_bed, r_sorted)):
            with open(dst, "w") as out:
                subprocess.run(["sort", "-k1,1", "-k2,2n", src], check=True, stdout=out, timeout=300)
        proc = subprocess.run(
            ["bedtools", "intersect", "-a", q_sorted, "-b", r_sorted, "-wo"],
            check=True, capture_output=True, text=True, timeout=900,
        )
    except (subprocess.CalledProcessError, subprocess.TimeoutExpired, OSError) as e:
        LOG.warning("bedtools intersect failed: %s", e)
        return None

    query_ids, ref_ids, same_strand, opp_strand = set(), set(), set(), set()
    for line in proc.stdout.splitlines():
        cols = line.split("\t")
        if len(cols) < 13:
            continue
        q_id, q_strand, r_id, r_strand = cols[3], cols[5], cols[9], cols[11]
        query_ids.add(q_id)
        ref_ids.add(r_id)
        (same_strand if q_strand == r_strand else opp_strand).add(q_id)

    with open(query_bed) as fh:
        n_query = sum(1 for _ in fh)
    with open(ref_bed) as fh:
        n_ref = sum(1 for _ in fh)
    return {
        "query_genes": n_query,
        "ref_genes": n_ref,
        "query_genes_overlapping_ref": len(query_ids),
        "query_genes_overlap_same_strand": len(same_strand),
        "query_genes_overlap_opposite_strand": len(opp_strand),
        "query_genes_no_overlap": n_query - len(query_ids),
        "ref_genes_not_recovered": n_ref - len(ref_ids),
    }
