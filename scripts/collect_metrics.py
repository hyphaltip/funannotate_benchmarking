#!/usr/bin/env python3
"""Parse each cell's Nextflow trace into results/metrics.tsv.

DESIGN.md "Metrics collection": wall-clock time, CPU-hours, peak memory per
process and pipeline-total, normalized by genome size (Mb).

Each benchmark cell (runs/<cell>/) writes its own trace at
logs/nextflow/annotate_trace.txt with a fixed field list (see
nf_funannotate1's conf/profile_annotate.config `trace.fields`):
    task_id,hash,name,status,exit,realtime,%cpu,rss,tag
`tag` is the genome id (matches runs/<cell>/genome_annotation/<tag>/).

Outputs:
  results/metrics_by_process.tsv   cell x genome x process rollup
  results/metrics.tsv              cell x genome pipeline totals

Caveat: the trace has no submit/complete timestamps, so "wall clock" here is
a SUM of task realtime, not true elapsed time — tasks for one genome run
concurrently across processes/SLURM jobs, so this overstates true wall time.
It's still the right unit to compare processes/cells/genomes against each
other. Container pull/build time is a one-time cost and is intentionally not
included (DESIGN.md tracks that separately).

Usage:
  python3 scripts/collect_metrics.py --cells-dir runs --samples samples.csv
  python3 scripts/collect_metrics.py --cells-dir runs --cells v1.8.17_conda v1.8.17_container
"""
from __future__ import annotations

import argparse
import os
import sys
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lib  # noqa: E402


def discover_cells(cells_dir: str) -> list:
    cells = []
    for name in sorted(os.listdir(cells_dir)):
        trace = os.path.join(cells_dir, name, "logs", "nextflow", "annotate_trace.txt")
        if os.path.isfile(trace):
            cells.append(name)
    return cells


def process_name(full_name: str) -> str:
    """'TRAIN_PREDICT:GENEMARK_RUN (Foo_bar)' -> 'TRAIN_PREDICT:GENEMARK_RUN'."""
    return full_name.split(" (", 1)[0].strip()


def load_cell_tasks(cell: str, cells_dir: str) -> list:
    trace_path = os.path.join(cells_dir, cell, "logs", "nextflow", "annotate_trace.txt")
    rows = []
    for r in lib.read_trace(trace_path):
        realtime_s = lib.parse_nf_duration(r.get("realtime", ""))
        cpu_pct = lib.parse_pct(r.get("%cpu", ""))
        rss_gb = lib.parse_mem_to_gb(r.get("rss", ""))
        cpu_hours = None
        if cpu_pct is not None and realtime_s is not None:
            cpu_hours = (cpu_pct / 100.0) * (realtime_s / 3600.0)
        rows.append({
            "cell": cell,
            "genome": r.get("tag", "") or "",
            "process": process_name(r.get("name", "")),
            "task_name": r.get("name", ""),
            "status": r.get("status", ""),
            "exit": r.get("exit", ""),
            "realtime_s": realtime_s,
            "cpu_pct": cpu_pct,
            "peak_rss_gb": rss_gb,
            "cpu_hours": cpu_hours,
        })
    return rows


def rollup(all_tasks: list, key_fields) -> dict:
    agg = defaultdict(lambda: {"n_tasks": 0, "n_completed": 0, "n_failed": 0, "n_aborted": 0,
                                "realtime_s_sum": 0.0, "cpu_hours_sum": 0.0, "peak_rss_gb_max": 0.0})
    for t in all_tasks:
        key = tuple(t[f] for f in key_fields)
        a = agg[key]
        a["n_tasks"] += 1
        if t["status"] == "COMPLETED":
            a["n_completed"] += 1
        elif t["status"] == "FAILED":
            a["n_failed"] += 1
        elif t["status"] == "ABORTED":
            a["n_aborted"] += 1
        if t["realtime_s"] is not None:
            a["realtime_s_sum"] += t["realtime_s"]
        if t["cpu_hours"] is not None:
            a["cpu_hours_sum"] += t["cpu_hours"]
        if t["peak_rss_gb"] is not None:
            a["peak_rss_gb_max"] = max(a["peak_rss_gb_max"], t["peak_rss_gb"])
    return agg


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--cells-dir", default="runs")
    ap.add_argument("--cells", nargs="*", default=None, help="restrict to these cell names (default: auto-discover)")
    ap.add_argument("--samples", default="samples.csv")
    ap.add_argument("--out", default="results/metrics.tsv")
    ap.add_argument("--by-process-out", default="results/metrics_by_process.tsv")
    ap.add_argument("--tasks-out", default=None, help="optional: dump every raw task row here")
    ap.add_argument("-v", "--verbose", action="count", default=0)
    args = ap.parse_args()
    lib.setup_logging(args.verbose)

    if not os.path.isdir(args.cells_dir):
        lib.LOG.error("cells dir not found: %s", args.cells_dir)
        sys.exit(1)

    samples_idx = {}
    if os.path.exists(args.samples):
        samples_idx = lib.load_samples_index(args.samples)
    else:
        lib.LOG.warning("samples.csv not found (%s) — genome_size normalization skipped", args.samples)

    cells = args.cells or discover_cells(args.cells_dir)
    if not cells:
        lib.LOG.error("no cells with a trace file found under %s", args.cells_dir)
        sys.exit(1)
    lib.LOG.info("found %d cell(s): %s", len(cells), ", ".join(cells))

    all_tasks = []
    for cell in cells:
        cell_tasks = load_cell_tasks(cell, args.cells_dir)
        lib.LOG.info("%s: %d task rows", cell, len(cell_tasks))
        all_tasks.extend(cell_tasks)

    if args.tasks_out:
        lib.write_table(args.tasks_out, all_tasks, sep="\t")

    # ── per (cell, genome, process) ─────────────────────────────────────────
    by_process = rollup(all_tasks, ("cell", "genome", "process"))
    process_rows = []
    for (cell, genome, process), a in sorted(by_process.items()):
        process_rows.append({
            "cell": cell, "genome": genome, "process": process,
            "n_tasks": a["n_tasks"], "n_completed": a["n_completed"],
            "n_failed": a["n_failed"], "n_aborted": a["n_aborted"],
            "realtime_hours_sum": round(a["realtime_s_sum"] / 3600.0, 4),
            "cpu_hours_sum": round(a["cpu_hours_sum"], 4),
            "peak_rss_gb_max": round(a["peak_rss_gb_max"], 3),
        })
    os.makedirs(os.path.dirname(args.by_process_out) or ".", exist_ok=True)
    lib.write_table(args.by_process_out, process_rows, sep="\t")
    lib.LOG.info("wrote %d rows -> %s", len(process_rows), args.by_process_out)

    # ── per (cell, genome) pipeline totals ──────────────────────────────────
    by_genome = rollup(all_tasks, ("cell", "genome"))
    genome_rows = []
    for (cell, genome), a in sorted(by_genome.items()):
        status = "success" if a["n_failed"] == 0 and a["n_aborted"] == 0 else \
                 ("aborted" if a["n_aborted"] else "failed")
        row = {
            "cell": cell, "genome": genome, "status": status,
            "n_tasks": a["n_tasks"], "n_completed": a["n_completed"],
            "n_failed": a["n_failed"], "n_aborted": a["n_aborted"],
            "task_time_hours_sum": round(a["realtime_s_sum"] / 3600.0, 4),
            "cpu_hours_sum": round(a["cpu_hours_sum"], 4),
            "peak_rss_gb_max": round(a["peak_rss_gb_max"], 3),
        }
        srow = samples_idx.get(genome, {})
        row["accession"] = srow.get("accession", "")
        row["phylum"] = srow.get("phylum", "")
        size_mb = srow.get("genome_size_mb", "")
        row["genome_size_mb"] = size_mb
        try:
            size = float(size_mb)
            row["cpu_hours_per_mb"] = round(a["cpu_hours_sum"] / size, 5) if size else ""
            row["task_time_hours_per_mb"] = round(a["realtime_s_sum"] / 3600.0 / size, 5) if size else ""
        except (TypeError, ValueError):
            row["cpu_hours_per_mb"] = ""
            row["task_time_hours_per_mb"] = ""
        genome_rows.append(row)

    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    lib.write_table(args.out, genome_rows, sep="\t")
    lib.LOG.info("wrote %d rows -> %s", len(genome_rows), args.out)

    n_non_success = sum(1 for r in genome_rows if r["status"] != "success")
    print(f"cells={len(cells)} genome_rows={len(genome_rows)} non_success={n_non_success}",
          file=sys.stderr)


if __name__ == "__main__":
    main()
