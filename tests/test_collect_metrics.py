"""collect_metrics.py: which trace rows count as the timing measurement."""
import csv
import os
import subprocess
import sys

SCRIPT = os.path.join(os.path.dirname(__file__), "..", "scripts", "collect_metrics.py")

OLD_HDR = "task_id\thash\tname\tstatus\texit\trealtime\t%cpu\trss\ttag"
NEW_HDR = OLD_HDR + "\tattempt\tcpus\tmemory\tqueue\thostname\tsubmit\tstart\tcomplete\tpeak_rss"


def write(path, header, rows):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as fh:
        fh.write(header + "\n")
        for r in rows:
            fh.write("\t".join(r) + "\n")


def run(tmp_path):
    out, byp = tmp_path / "m.tsv", tmp_path / "p.tsv"
    subprocess.run([sys.executable, SCRIPT, "--cells-dir", str(tmp_path / "runs"),
                    "--samples", str(tmp_path / "none.csv"), "--out", str(out),
                    "--by-process-out", str(byp)], check=True, capture_output=True)
    rd = lambda p: list(csv.DictReader(open(p), delimiter="\t"))
    return rd(out), rd(byp)


def test_final_attempt_only(tmp_path):
    d = tmp_path / "runs" / "cellA" / "logs" / "nextflow"
    # Old-format trace: TRAIN fails once (1h), then completes (2h). PREDICT
    # completes (1h) but is re-run after a fix in the newer trace.
    write(str(d / "annotate_trace.20260901_000000.txt"), OLD_HDR, [
        ["1", "aa/1", "TP:TRAIN (G1)", "FAILED", "-", "1h", "-", "-", "G1"],
        ["2", "aa/2", "TP:TRAIN (G1)", "COMPLETED", "0", "2h", "800%", "10 GB", "G1"],
        ["3", "aa/3", "TP:PREDICT (G1)", "COMPLETED", "0", "1h", "400%", "5 GB", "G1"],
    ])
    # New-format trace: TRAIN is a cache hit, PREDICT re-runs (3h, 8 cpus).
    write(str(d / "annotate_trace.20260902_000000.txt"), NEW_HDR, [
        ["2", "aa/2", "TP:TRAIN (G1)", "CACHED", "0", "2h", "800%", "10 GB", "G1",
         "2", "24", "192 GB", "epyc", "n1", "-", "-", "-", "12 GB"],
        ["4", "bb/4", "TP:PREDICT (G1)", "COMPLETED", "0", "3h", "100%", "6 GB", "G1",
         "1", "8", "32 GB", "epyc", "n2", "2026-09-02 00:00:00.000",
         "2026-09-02 01:00:00.000", "2026-09-02 04:00:00.000", "7 GB"],
    ])
    m, p = run(tmp_path)
    g = m[0]
    assert g["status"] == "success"
    assert g["n_tasks"] == "2" and g["n_failed"] == "1" and g["n_superseded"] == "1"
    # final = TRAIN 2h + PREDICT (re-run) 3h; the cached row is not counted twice
    assert float(g["task_time_hours_sum"]) == 5.0
    assert float(g["cpu_hours_sum"]) == 8 * 2 + 1 * 3
    # overhead = failed TRAIN 1h + superseded PREDICT 1h
    assert float(g["overhead_task_time_hours_sum"]) == 2.0
    # TRAIN's final row is old-format -> no start/complete -> no wall clock
    assert g["wall_clock_hours"] == ""
    pred = next(r for r in p if r["process"] == "TP:PREDICT")
    assert pred["wall_clock_hours"] == "3.0"
    assert pred["cpus"] == "8" and pred["queues"] == "epyc" and pred["attempts"] == "1"
    assert pred["rss_field"] == "peak_rss" and float(pred["peak_rss_gb_max"]) == 7.0


def test_never_completed_is_failed(tmp_path):
    d = tmp_path / "runs" / "cellB" / "logs" / "nextflow"
    write(str(d / "annotate_trace.20260901_000000.txt"), OLD_HDR, [
        ["1", "aa/1", "TP:TRAIN (G2)", "FAILED", "1", "30m", "100%", "1 GB", "G2"],
    ])
    m, _ = run(tmp_path)
    assert m[0]["status"] == "failed" and m[0]["n_unresolved"] == "1"
    assert float(m[0]["task_time_hours_sum"]) == 0.0
    assert float(m[0]["overhead_task_time_hours_sum"]) == 0.5
