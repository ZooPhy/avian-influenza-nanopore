from __future__ import annotations

import csv
import subprocess
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
SNAKEFILE = REPO_ROOT / "tests" / "integration2b" / "Snakefile"
WORK = REPO_ROOT / "tests" / "integration2b" / "work"


def read_tsv(path: Path) -> dict[str, str]:
    with path.open(newline="", encoding="utf-8") as handle:
        return next(csv.DictReader(handle, delimiter="\t"))


def read_blast(path: Path) -> dict[str, dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return {row["segment"]: row for row in csv.DictReader(handle)}


def test_segment_qc_to_blast_summary_wiring():
    for sample in ("qc_pass", "length_warning", "low_breadth", "no_hit"):
        (WORK / sample / "irma" / "segments").mkdir(parents=True, exist_ok=True)

    subprocess.run(
        [
            "snakemake",
            "--snakefile",
            str(SNAKEFILE),
            "--cores",
            "2",
            "--sdm",
            "conda",
            "--rerun-incomplete",
        ],
        cwd=REPO_ROOT,
        check=True,
    )

    pass_stats = read_tsv(WORK / "qc_pass" / "coverage_stats" / "HA.tsv")
    assert pass_stats["status"] == "PASS"
    assert pass_stats["coverage_status"] == "PASS"
    assert pass_stats["length_status"] == "PASS"

    warning_stats = read_tsv(WORK / "length_warning" / "coverage_stats" / "HA.tsv")
    assert warning_stats["status"] == "PASS"
    assert warning_stats["length_status"] == "WARNING"
    assert warning_stats["selection_status"] == "MULTIPLE_CANDIDATES"
    assert warning_stats["n_candidates"] == "2"

    fail_stats = read_tsv(WORK / "low_breadth" / "coverage_stats" / "HA.tsv")
    assert fail_stats["status"] == "FAIL"
    assert fail_stats["coverage_status"] == "FAIL"
    assert float(fail_stats["breadth_covered"]) < 0.95

    qc_pass = read_blast(WORK / "qc_pass" / "summary" / "blast_top_hits.csv")
    assert qc_pass["HA"]["hit_status"] == "HIGH_CONFIDENCE"
    assert qc_pass["HA"]["query_coverage"] == "100.000"

    warning = read_blast(WORK / "length_warning" / "summary" / "blast_top_hits.csv")
    assert warning["HA"]["hit_status"] == "LOW_CONFIDENCE"
    assert warning["HA"]["query_coverage"] == "89.261"

    low_breadth = read_blast(WORK / "low_breadth" / "summary" / "blast_top_hits.csv")
    assert low_breadth["HA"]["hit_status"] == "SKIPPED_QC"

    no_hit = read_blast(WORK / "no_hit" / "summary" / "blast_top_hits.csv")
    assert no_hit["HA"]["hit_status"] == "NO_HIT"

    for sample in ("qc_pass", "length_warning", "low_breadth", "no_hit"):
        rows = read_blast(WORK / sample / "summary" / "blast_top_hits.csv")
        assert set(rows) == {"HA", "NA", "PB2", "PB1", "PA", "NP", "MP", "NS"}
