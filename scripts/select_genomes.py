#!/usr/bin/env python3
"""Subsample genomes for the funannotate benchmark.

Reads the candidate pool (Fungi_BFD accessions + taxonomy + rnaseq tables),
applies the DESIGN.md selection constraints (taxonomic quota floors, genome
size cap, genus/species caps, RefSeq preference, class-weighted diversity),
and appends the chosen samples to samples.csv (NEVER removing rows already
present).

Metadata per row follows DESIGN.md: accession, species, strain, genus, phylum,
subphylum, class, order, genome_size_mb, rnaseq_available, refseq_annotation_release,
assembly_level, busco_lineage_used. The NCBI annotation/assembly metadata is
filled best-effort from the NCBI datasets API (cached under results/); pass
--no-enrich to skip network entirely.
"""
from __future__ import annotations

import argparse
import json
import os
import random
import sys
import time
from collections import Counter, defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lib  # noqa: E402

TARGET_PHYLA = ("Ascomycota", "Basidiomycota", "Mucoromycota",
                "Chytridiomycota", "Blastocladiomycota")


def parse_args():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--config", default="conf/benchmark.yaml")
    ap.add_argument("--samples", default="samples.csv")
    ap.add_argument("--cache-dir", default="results/metadata_cache")
    ap.add_argument("--seed", type=int, default=None,
                    help="override the seed in the config")
    ap.add_argument("--no-enrich", action="store_true",
                    help="skip NCBI metadata enrichment (offline mode)")
    ap.add_argument("--dry-run", action="store_true",
                    help="print the selection without writing samples.csv")
    ap.add_argument("--verbose", "-v", action="count", default=0)
    return ap.parse_args()


def bucket_key(phylum: str, subphylum: str) -> tuple:
    return (phylum, subphylum)


def species_tag(species: str) -> str:
    return species.replace(" ", "_") if species else ""


def pick_from_bucket(rows, caps, used, rng, class_exp,
                     refseq_only: bool):
    """Pick one accession from a bucket honoring genus/species caps.

    With refseq_only=True only GCF_ rows are eligible; when the bucket holds no
    eligible GCF_ row the caller may retry with False. Class-weighted: a class
    is drawn with probability proportional to (number of its candidates **
    class_exp). Returns (row, class_name) or (None, "").
    """
    eligible = [r for r in rows if not used[(r["ACCESSION"],)]
                and (not refseq_only or r["IS_REFSEQ"])
                and caps["genus"][r["GENUS"]] < caps["max_per_genus"]
                and caps["species"][r["SPECIES"]] < caps["max_per_species"]]
    if not eligible:
        return None, ""
    by_class = defaultdict(list)
    for r in eligible:
        by_class[r["CLASS"] or "unknown"].append(r)
    classes = list(by_class)
    w = [len(by_class[c]) ** class_exp for c in classes]
    for _ in range(len(classes) * 3):
        klass = rng.choices(classes, weights=w, k=1)[0]
        # pick the first eligible member of the drawn class
        r = by_class[klass][0]
        by_class[klass].pop(0)
        classes.remove(klass)
        if classes:
            w = [len(by_class[c]) ** class_exp for c in classes]
        return r, klass
    return None, ""


def bucket_has(key, buckets, used, refseq_required):
    """First eligible row in bucket key, optionally restricted to RefSeq."""
    for r in buckets.get(key, []):
        if used[(r["ACCESSION"],)]:
            continue
        if refseq_required and not r["IS_REFSEQ"]:
            continue
        return r
    return None


def main():
    args = parse_args()
    lib.setup_logging(args.verbose)
    cfg = lib.load_yaml(args.config)
    sel = cfg["selection"]

    seed = args.seed if args.seed is not None else int(sel["seed"])
    rng = random.Random(seed)
    max_size = int(sel["max_genome_size_mb"]) * 1_000_000
    refseq_only = sel["reference_preference"] == "refseq_only"
    class_exp = float(sel.get("class_weight_exponent", 0.5))

    pool = lib.load_pool(cfg["pool"]["accessions"], cfg["pool"]["taxonomy"])
    subph_map = lib.load_subphylum_map("conf/subphylum_map.tsv")
    lineage_map = lib.load_lineage_map("conf/busco_lineage_map.tsv")
    rnaseq = lib.load_rnaseq(cfg["pool"]["rnaseq"])
    masked_dir = cfg["pool"].get("masked_genome_dir", "")
    masked_accs = lib.masked_accessions(masked_dir) if masked_dir else set()

    # ── existing samples.csv is authoritative: never drop, never re-pick ──
    existing = []
    if os.path.exists(args.samples) and os.path.getsize(args.samples) > 0:
        existing = lib.read_table(args.samples)

    def already_selected(acc):
        return any(str(r.get("accession", "")).startswith(acc.split(".")[0])
                   for r in existing)

    caps = {
        "max_per_genus": max(1, int(sel["max_per_genus"])),
        "max_per_species": max(1, int(sel["max_per_species"])),
        "genus": defaultdict(int),
        "species": defaultdict(int),
    }
    for r in existing:
        caps["genus"][r.get("genus", "")] += 1
        caps["species"][r.get("species", "")] += 1

    # ── build candidate buckets: masked genome available + size cap + phylum ──
    buckets = {}
    for acc, row in pool.items():
        if already_selected(acc):
            continue
        phylum = row.get("PHYLUM", "")
        if phylum not in TARGET_PHYLA:
            continue
        try:
            if int(row.get("ASM_LENGTH", 0) or 0) > max_size:
                continue
        except ValueError:
            continue
        subphylum = lib.resolve_subphylum(row, subph_map)
        key = bucket_key(phylum, subphylum)
        entry = {
            "ACCESSION": acc,
            "SPECIES": row.get("SPECIES", ""),
            "STRAIN": row.get("STRAIN", ""),
            "NCBI_TAXID": row.get("NCBI_TAXID", ""),
            "BIOPROJECT": row.get("BIOPROJECT", ""),
            "GENUS": row.get("GENUS", ""),
            "CLASS": row.get("CLASS", ""),
            "ORDER": row.get("ORDER", ""),
            "PHYLUM": phylum,
            "SUBPHYLUM": subphylum,
            "ASM_LENGTH": row.get("ASM_LENGTH", ""),
            "ASM_NAME": row.get("ASM_NAME", ""),
            "N50": row.get("N50", ""),
            "IS_REFSEQ": acc.startswith("GCF_"),
            "masked": row.get("ASM_FOLDER", "") in masked_accs,
        }
        buckets.setdefault(key, []).append(entry)

    # ── sorted key order: refseq-first primary, masked-available secondary ──
    def sortkey(e):
        return (0 if e["IS_REFSEQ"] else 1, 0 if e["masked"] else 1, e["ACCESSION"])

    for key in buckets:
        buckets[key].sort(key=sortkey)

    selected = []          # ordered list of picks
    used = defaultdict(bool)
    selected_set = set()

    def add_pick(row):
        selected.append(row)
        selected_set.add(row["ACCESSION"])
        used[(row["ACCESSION"],)] = True
        caps["genus"][row["GENUS"]] += 1
        caps["species"][row["SPECIES"]] += 1

    # ── PHASE 1: quota floors (rare groups first so they are guaranteed) ──
    floors = sel.get("quota_floors", {})
    # order: Blastocladiomycota, Chytridiomycota, Mucoromycota (rarest) first,
    # then Basidiomycota subphylum floors, Ascomycota floors last.
    def floor_order(phylum):
        return {"Blastocladiomycota": 0, "Chytridiomycota": 1, "Mucoromycota": 2,
                "Basidiomycota": 3, "Ascomycota": 4}.get(phylum, 9)

    for phylum in sorted(floors, key=floor_order):
        spec = floors[phylum]
        n_existing_phylum = sum(1 for r in existing if r.get("phylum") == phylum)
        total = max(0, int(spec.get("total", 0)) - n_existing_phylum)
        sub_floors = spec.get("subphylum", {})
        # subphylum-specific floors first (refseq preferred, GCA fallback)
        for sp, n in sorted(sub_floors.items(), key=lambda kv: -int(kv[1])):
            n_existing_sp = sum(1 for r in existing
                                if r.get("phylum") == phylum and r.get("subphylum") == sp)
            need = max(0, int(n) - n_existing_sp)
            while need > 0:
                for required in ([refseq_only, False] if refseq_only else [False]):
                    row = pick_from_bucket(buckets[bucket_key(phylum, sp)], caps,
                                           used, rng, class_exp, required)[0]
                    if row:
                        break
                if row is None:
                    break
                add_pick(row)
                need -= 1
            total -= int(n)
        # phylum-level floor: draw across ALL of the phylum's buckets so an
        # empty-subphylum bucket cannot hog the fills
        while total > 0:
            bucket_keys = [k for k in buckets if k[0] == phylum
                           and bucket_has(k, buckets, used, False)]
            if not bucket_keys:
                break
            row = None
            for required in ([refseq_only, False] if refseq_only else [False]):
                if row:
                    break
                keys = [k for k in bucket_keys if bucket_has(k, buckets, used, required)]
                if not keys:
                    continue
                key = rng.choice(keys)   # uniform across non-empty phylum buckets
                row = pick_from_bucket(buckets[key], caps, used, rng,
                                       class_exp, required)[0]
            if row is None:
                break
            add_pick(row)
            total -= 1

    # ── PHASE 2: fill remaining slots proportional to fill weights ──
    fill_weights = sel.get("fill_weights", {})
    subph_weights = sel.get("subphylum_weights", {})
    n_target = int(sel["n_target"])
    remaining = max(0, n_target - len(existing) - len(selected_set))
    guard = 0
    while remaining > 0 and guard < n_target * 20:
        guard += 1
        # weighted phylum draw over phyla that still have candidates
        phyla = []
        w = []
        for phylum, fw in fill_weights.items():
            if any(bucket_has(k, buckets, used, refseq_only)
                   or (refseq_only and bucket_has(k, buckets, used, False))
                   for k in buckets if k[0] == phylum):
                phyla.append(phylum)
                w.append(fw)
        if not phyla:
            break
        phylum = rng.choices(phyla, weights=w, k=1)[0]
        # weighted subphylum draw
        subph_options = [k[1] for k in buckets if k[0] == phylum
                         and (bucket_has(k, buckets, used, refseq_only)
                              or (refseq_only and bucket_has(k, buckets, used, False)))]
        spw = subph_weights.get(phylum, {})
        subph_choices, subph_w = subph_options, [spw.get(sp, 0.0) for sp in subph_options]
        if not subph_choices or sum(subph_w) <= 0:
            continue
        subphylum = rng.choices(subph_choices, weights=subph_w, k=1)[0]
        key = bucket_key(phylum, subphylum)
        row, _ = pick_from_bucket(buckets[key], caps, used, rng,
                                  class_exp, refseq_only)
        if row is None and refseq_only:
            # fall back to GCA_ within this bucket
            row, _ = pick_from_bucket(buckets[key], caps, used, rng,
                                      class_exp, False)
        if row is None:
            continue
        add_pick(row)
        remaining -= 1

    # ── verify floors met ──
    def verify_floors(all_rows):
        counts = defaultdict(int)
        for r in all_rows:
            counts[(r.get("phylum", ""), r.get("subphylum", ""))] += 1
        ok = True
        for phylum, spec in sorted(floors.items()):
            got = sum(v for (p, _), v in counts.items() if p == phylum)
            if got < int(spec["total"]):
                lib.LOG.warning("phylum floor unmet: %s got %d want %d",
                                phylum, got, int(spec["total"]))
                ok = False
            for sp, n in spec.get("subphylum", {}).items():
                if counts[(phylum, sp)] < int(n):
                    lib.LOG.warning("subphylum floor unmet: %s/%s got %d want %d",
                                    phylum, sp, counts[(phylum, sp)], int(n))
                    ok = False
        return ok

    # ── metadata enrichment (best-effort, cached; idempotent for old rows) ──
    CAP_MAP = {"species": "SPECIES", "strain": "STRAIN", "genus": "GENUS",
               "phylum": "PHYLUM", "subphylum": "SUBPHYLUM", "class": "CLASS",
               "order": "ORDER", "taxid": "NCBI_TAXID", "bioproject": "BIOPROJECT",
               "asm_name": "ASM_NAME", "n50": "N50"}

    def enrich_row(row):
        acc = row.get("ACCESSION") or row.get("accession", "")
        row["accession"] = acc
        if not row.get("genome_size_mb"):
            row["genome_size_mb"] = lib.mb(row.get("ASM_LENGTH", ""))
        p = pool.get(acc, {})
        for k, cap in CAP_MAP.items():
            if not row.get(k):
                row[k] = row.get(cap, "") or p.get(cap, "")
        tag = species_tag(row["species"])
        sra = rnaseq.get(tag, [])
        row["rnaseq_available"] = "yes" if sra else "no"
        row["rnaseq_sra_ids"] = ",".join(sra[:8])
        if not row.get("busco_lineage_used"):
            row["busco_lineage_used"] = lib.resolve_lineage(
                row.get("phylum", ""), row.get("subphylum", ""), row.get("class", ""),
                lineage_map)
        row.setdefault("selected_date", time.strftime("%Y-%m-%d"))
        row.setdefault("status", "selected")
        if args.no_enrich:
            return row
        if row.get("refseq_annotation_release") or row.get("assembly_level") \
                or row.get("ncbi_current_accession") or row.get("paired_accession"):
            return row  # already enriched
        cache = os.path.join(args.cache_dir, acc + ".json")
        meta = None
        if os.path.exists(cache):
            meta = json.load(open(cache))
            if not any(meta.values()):   # stale empty cache (API hiccup)
                os.remove(cache)
                meta = None
        if meta is None:
            meta = lib.enrich_assembly_metadata(
                acc, ncbi_email=cfg.get("fetch", {}).get("ncbi_email", ""),
                retries=3, sleep=4.0)
            if any(meta.values()):       # never cache an empty/failed payload
                os.makedirs(args.cache_dir, exist_ok=True)
                json.dump(meta, open(cache, "w"), indent=1)
            time.sleep(0.5)
        for k, v in meta.items():
            if v:
                row.setdefault(k, v)
        row["refseq_annotation_release"] = meta.get("refseq_annotation_release", "")
        row["assembly_level"] = meta.get("assembly_level", "")
        row["ncbi_current_accession"] = meta.get("ncbi_current_accession", "")
        row["paired_accession"] = meta.get("paired_accession", "")
        return row

    cols = ["accession", "species", "strain", "taxid", "bioproject", "genus",
            "phylum", "subphylum", "class", "order", "genome_size_mb",
            "assembly_level", "refseq_annotation_release",
            "busco_lineage_used", "rnaseq_available", "rnaseq_sra_ids",
            "ncbi_current_accession", "paired_accession", "asm_name",
            "n50", "selected_date", "status"]
    rows = []
    for r in existing:
        enrich_row(r)
        rows.append({c: r.get(c, "") for c in cols})
    for r in selected:
        enrich_row(r)
        rows.append({c: r.get(c, "") for c in cols})
    ok = verify_floors(rows)

    # ── report ──
    phylum_totals = defaultdict(int)
    refseq_n = sum(1 for r in rows if r.get("accession", "").startswith("GCF_"))
    for r in rows:
        phylum_totals[r["phylum"]] += 1
    print(f"selected {len(rows)} rows total "
          f"(new this run: {len(selected)}), RefSeq: {refseq_n}/{len(rows)}")
    for phylum, n in sorted(phylum_totals.items()):
        subs = Counter(str(r["subphylum"]) for r in rows if r["phylum"] == phylum)
        detail = ", ".join(f"{sp}:{c}" for sp, c in sorted(subs.items()))
        print(f"  {phylum:<16} {n:>3}   [{detail}]")
    classes = Counter(r["class"] for r in rows)
    print("classes:", ", ".join(f"{c}:{n}" for c, n in sorted(classes.items())))
    print("floors met:", "YES" if ok else "NO (see warnings)")

    if args.dry_run:
        print("(dry run: not writing samples.csv)")
        return
    lib.write_table(args.samples, rows, sep=",")
    print(f"wrote {len(rows)} rows to {args.samples}")


if __name__ == "__main__":
    main()
