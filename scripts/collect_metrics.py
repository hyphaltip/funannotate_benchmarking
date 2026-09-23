#!/usr/bin/env python3
"""Parse each cell's Nextflow trace into results/metrics.tsv.

DESIGN.md "Metrics collection": wall-clock time, CPU-hours, peak memory per
process and pipeline-total, normalized by genome size (Mb).

Each benchmark cell (runs/<cell>/) writes its own trace at
logs/nextflow/annotate_trace.<timestamp>.txt with a fixed field list (see
nf_funannotate1's conf/profile_annotate.config `trace.fields`):
    task_id,hash,name,status,exit,realtime,%cpu,rss,tag
plus, in traces written from 2026-09-23 on:
    attempt,cpus,memory,queue,hostname,submit,start,complete,peak_rss
`tag` is the genome id (matches runs/<cell>/genome_annotation/<tag>/). Older
traces lack the extra columns; those output columns are then left blank.

The filename carries a per-invocation timestamp rather than being fixed:
every `-resume` relaunch is a fresh `nextflow run`, and Nextflow's trace file
is written with overwrite=true, so a FIXED filename would silently lose every
previous invocation's rows -- including for genomes already cache-hit/
completed in the new run (confirmed: a cell resumed many times had a trace
with 5 rows despite 57/65 genomes having real published output). This script
globs every annotate_trace.*.txt snapshot per cell and merges them.

Which rows count as the measurement (select_final_tasks()):
  - Per (cell, task name), only the LAST COMPLETED row is the measured run
    ("final"). Snapshots are read in filename (= launch time) order, rows in
    file order, so a task re-run after a fix replaces its earlier completion.
  - CACHED rows are dropped: Nextflow repeats the original run's realtime/%cpu
    on a cache hit, so counting them would count that run twice.
  - FAILED/ABORTED rows, and COMPLETED rows superseded by a later completion,
    are real compute spent but not the measurement. They go into separate
    overhead_* columns, never into the timing sums.
  A task name with no COMPLETED row at all is "unresolved"; it sets the
  genome's status to failed/aborted.

runs/genemark_sidecar/ (launch/run_genemark_sidecar.sbatch) is discovered the
same way, from its own logs/nextflow/genemark_sidecar_trace.<timestamp>.txt
(same field list — conf/profile_genemark_sidecar.config). It is a ONE-TIME cost shared by
all 6 cells (GeneMark runs once per genome, not once per cell — see
DESIGN.md "GeneMark sidecar"), so it shows up as its own "cell" named
genemark_sidecar in every output table here rather than being folded into
each real cell's numbers; sum it in exactly once, not once per cell, when
computing a benchmark-wide total cost.

Outputs:
  results/metrics_by_process.tsv   cell x genome x process rollup
  results/metrics.tsv              cell x genome pipeline totals

task_time_hours_sum is a SUM of task realtime, not elapsed time -- tasks for
one genome run concurrently, so it overstates true wall time. When every
final task of a genome has start/complete timestamps, wall_clock_hours gives
the true elapsed span (first start -> last complete); otherwise it is blank.
Neither includes SLURM queue wait (submit -> start). Container pull/build time is a one-time cost and is intentionally not
included (DESIGN.md tracks that separately).

Usage:
  python3 scripts/collect_metrics.py --cells-dir runs --samples samples.csv
  python3 scripts/collect_metrics.py --cells-dir runs --cells v1.8.17_conda v1.8.17_container
"""
from __future__ import annotations

import argparse
import glob
import os
import sys
from collections import defaultdict
from datetime import datetime

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lib  # noqa: E402


# cell dir name -> trace filename GLOB it writes (see nf_funannotate1's
# conf/profile_annotate.config vs. conf/profile_genemark_sidecar.config).
# Each invocation's trace carries its own timestamp, so a cell accumulates
# one snapshot per `nextflow run`/`-resume` -- all of them are merged below.
TRACE_GLOBS = {
    "genemark_sidecar": "genemark_sidecar_trace*.txt",
}
DEFAULT_TRACE_GLOB = "annotate_trace*.txt"


def cell_trace_paths(cell: str, cells_dir: str) -> list:
    pattern = os.path.join(cells_dir, cell, "logs", "nextflow", TRACE_GLOBS.get(cell, DEFAULT_TRACE_GLOB))
    return sorted(glob.glob(pattern))


def discover_cells(cells_dir: str) -> list:
    cells = []
    for name in sorted(os.listdir(cells_dir)):
        if cell_trace_paths(name, cells_dir):
            cells.append(name)
    return cells


def process_name(full_name: str) -> str:
    """'TRAIN_PREDICT:GENEMARK_RUN (Foo_bar)' -> 'TRAIN_PREDICT:GENEMARK_RUN'."""
    return full_name.split(" (", 1)[0].strip()


def _field(r: dict, key: str) -> str:
    v = (r.get(key, "") or "").strip()
    return "" if v in ("-", "NA") else v


def load_cell_tasks(cell: str, cells_dir: str) -> list:
    rows = []
    for trace_path in cell_trace_paths(cell, cells_dir):
        for r in lib.read_trace(trace_path):
            realtime_s = lib.parse_nf_duration(r.get("realtime", ""))
            cpu_pct = lib.parse_pct(r.get("%cpu", ""))
            # `rss` is resident memory, not the peak; prefer `peak_rss` when
            # the trace has it (traces from 2026-09-23 on).
            peak_rss_gb = lib.parse_mem_to_gb(r.get("peak_rss", ""))
            rss_field = "peak_rss"
            if peak_rss_gb is None:
                peak_rss_gb = lib.parse_mem_to_gb(r.get("rss", ""))
                rss_field = "rss" if peak_rss_gb is not None else ""
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
                "peak_rss_gb": peak_rss_gb,
                "rss_field": rss_field,
                "cpu_hours": cpu_hours,
                "attempt": _field(r, "attempt"),
                "cpus": _field(r, "cpus"),
                "memory_req_gb": lib.parse_mem_to_gb(_field(r, "memory")),
                "queue": _field(r, "queue"),
                "hostname": _field(r, "hostname"),
                "submit": _field(r, "submit"),
                "start": _field(r, "start"),
                "complete": _field(r, "complete"),
                "trace_file": os.path.basename(trace_path),
                "role": "",
            })
    return rows


def select_final_tasks(all_tasks: list) -> None:
    """Set each row's "role" in place: final / overhead / cached / unresolved.

    See the module docstring. Rows must be in chronological order (trace files
    sorted by timestamped filename, rows in file order).
    """
    last_completed = {}
    for i, t in enumerate(all_tasks):
        if t["status"] == "COMPLETED":
            last_completed[(t["cell"], t["task_name"])] = i
    for i, t in enumerate(all_tasks):
        key = (t["cell"], t["task_name"])
        if t["status"] == "CACHED":
            t["role"] = "cached"
        elif last_completed.get(key) == i:
            t["role"] = "final"
        elif key in last_completed:
            t["role"] = "overhead"
        else:
            # Never completed. Still compute spent, so it is overhead, but it
            # also marks the task as unresolved for the genome status.
            t["role"] = "unresolved"


def _join_unique(values) -> str:
    return ",".join(sorted({v for v in values if v}, key=lambda x: (len(x), x)))


def _parse_ts(s: str):
    # Nextflow trace default date format: 'yyyy-MM-dd HH:mm:ss.SSS'
    for fmt in ("%Y-%m-%d %H:%M:%S.%f", "%Y-%m-%d %H:%M:%S"):
        try:
            return datetime.strptime(s, fmt)
        except (ValueError, TypeError):
            continue
    return None


def rollup(all_tasks: list, key_fields) -> dict:
    agg = defaultdict(lambda: {"n_tasks": 0, "n_failed": 0, "n_aborted": 0, "n_superseded": 0,
                               "n_unresolved": 0, "realtime_s_sum": 0.0, "cpu_hours_sum": 0.0,
                               "peak_rss_gb_max": 0.0, "overhead_realtime_s_sum": 0.0,
                               "overhead_cpu_hours_sum": 0.0, "final": []})
    unresolved_names = defaultdict(set)
    for t in all_tasks:
        if t["role"] == "cached":
            continue
        key = tuple(t[f] for f in key_fields)
        a = agg[key]
        if t["status"] == "FAILED":
            a["n_failed"] += 1
        elif t["status"] == "ABORTED":
            a["n_aborted"] += 1
        if t["role"] == "final":
            a["n_tasks"] += 1
            a["final"].append(t)
            if t["realtime_s"] is not None:
                a["realtime_s_sum"] += t["realtime_s"]
            if t["cpu_hours"] is not None:
                a["cpu_hours_sum"] += t["cpu_hours"]
            if t["peak_rss_gb"] is not None:
                a["peak_rss_gb_max"] = max(a["peak_rss_gb_max"], t["peak_rss_gb"])
            continue
        if t["role"] == "overhead" and t["status"] == "COMPLETED":
            a["n_superseded"] += 1
        if t["role"] == "unresolved":
            unresolved_names[key].add(t["task_name"])
        if t["realtime_s"] is not None:
            a["overhead_realtime_s_sum"] += t["realtime_s"]
        if t["cpu_hours"] is not None:
            a["overhead_cpu_hours_sum"] += t["cpu_hours"]
    for key, names in unresolved_names.items():
        agg[key]["n_unresolved"] = len(names)
    for a in agg.values():
        final = a.pop("final")
        a["attempts"] = _join_unique(t["attempt"] for t in final)
        a["cpus"] = _join_unique(t["cpus"] for t in final)
        a["queues"] = _join_unique(t["queue"] for t in final)
        a["rss_field"] = _join_unique(t["rss_field"] for t in final)
        starts = [_parse_ts(t["start"]) for t in final]
        completes = [_parse_ts(t["complete"]) for t in final]
        if final and all(starts) and all(completes):
            a["wall_clock_hours"] = round((max(completes) - min(starts)).total_seconds() / 3600.0, 4)
        else:
            a["wall_clock_hours"] = ""
    return agg


def _common_cols(a: dict) -> dict:
    return {
        "n_tasks": a["n_tasks"], "n_failed": a["n_failed"], "n_aborted": a["n_aborted"],
        "n_superseded": a["n_superseded"], "n_unresolved": a["n_unresolved"],
    }


def _resource_cols(a: dict) -> dict:
    return {
        "cpu_hours_sum": round(a["cpu_hours_sum"], 4),
        "peak_rss_gb_max": round(a["peak_rss_gb_max"], 3),
        "rss_field": a["rss_field"],
        "overhead_task_time_hours_sum": round(a["overhead_realtime_s_sum"] / 3600.0, 4),
        "overhead_cpu_hours_sum": round(a["overhead_cpu_hours_sum"], 4),
        "attempts": a["attempts"], "cpus": a["cpus"], "queues": a["queues"],
    }


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
    select_final_tasks(all_tasks)

    if args.tasks_out:
        lib.write_table(args.tasks_out, all_tasks, sep="\t")
    n_role = defaultdict(int)
    for t in all_tasks:
        n_role[t["role"]] += 1
    lib.LOG.info("task rows by role: %s", dict(n_role))

    # ── per (cell, genome, process) ─────────────────────────────────────────
    by_process = rollup(all_tasks, ("cell", "genome", "process"))
    process_rows = []
    for (cell, genome, process), a in sorted(by_process.items()):
        process_rows.append({
            "cell": cell, "genome": genome, "process": process,
            **_common_cols(a),
            "realtime_hours_sum": round(a["realtime_s_sum"] / 3600.0, 4),
            "wall_clock_hours": a["wall_clock_hours"],
            **_resource_cols(a),
        })
    os.makedirs(os.path.dirname(args.by_process_out) or ".", exist_ok=True)
    lib.write_table(args.by_process_out, process_rows, sep="\t")
    lib.LOG.info("wrote %d rows -> %s", len(process_rows), args.by_process_out)

    # ── per (cell, genome) pipeline totals ──────────────────────────────────
    by_genome = rollup(all_tasks, ("cell", "genome"))
    genome_rows = []
    for (cell, genome), a in sorted(by_genome.items()):
        # Status from UNRESOLVED tasks only: a task that failed, was retried
        # and then completed is a success (its failures are overhead).
        unresolved = [t for t in all_tasks if t["role"] == "unresolved"
                      and t["cell"] == cell and t["genome"] == genome]
        if not unresolved:
            status = "success"
        else:
            status = "failed" if any(t["status"] == "FAILED" for t in unresolved) else "aborted"
        row = {
            "cell": cell, "genome": genome, "status": status,
            **_common_cols(a),
            "task_time_hours_sum": round(a["realtime_s_sum"] / 3600.0, 4),
            "wall_clock_hours": a["wall_clock_hours"],
            **_resource_cols(a),
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

    # ── Total benchmark cost: per-cell cost + the ONE-TIME sidecar cost ──────
    # (not per-cell x 6 -- see module docstring). genemark_sidecar is just
    # another row in genome_rows above if its trace was found, so this is a
    # plain sum split by whether the row IS that shared cost.
    sidecar_cpu_hours = sum(r["cpu_hours_sum"] for r in genome_rows if r["cell"] == "genemark_sidecar")
    percell_cpu_hours = sum(r["cpu_hours_sum"] for r in genome_rows if r["cell"] != "genemark_sidecar")
    total_cpu_hours = sidecar_cpu_hours + percell_cpu_hours
    overhead_cpu_hours = sum(r["overhead_cpu_hours_sum"] for r in genome_rows)
    print(
        f"total_cpu_hours={round(total_cpu_hours, 2)} "
        f"(per-cell={round(percell_cpu_hours, 2)}, genemark_sidecar one-time={round(sidecar_cpu_hours, 2)}) "
        f"overhead_cpu_hours={round(overhead_cpu_hours, 2)} (failed/aborted/superseded, not in total)",
        file=sys.stderr,
    )
    if sidecar_cpu_hours == 0.0 and "genemark_sidecar" not in cells:
        lib.LOG.warning(
            "genemark_sidecar not found under %s -- its cost is NOT included above "
            "(run launch/run_genemark_sidecar.sbatch, or pass --cells including it)",
            args.cells_dir,
        )


if __name__ == "__main__":
    main()
