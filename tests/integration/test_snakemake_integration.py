from __future__ import annotations

import csv
import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
CONFIG = REPO_ROOT / "tests" / "integration" / "config.yaml"
PREPARE = REPO_ROOT / "tests" / "integration" / "prepare_fixtures.py"
RESULTS = REPO_ROOT / "tests" / "integration" / "work" / "results"

CASES = {
    "it_detected": "DETECTED",
    "it_not_detected": "NOT_DETECTED",
    "it_indeterminate": "INDETERMINATE",
}


def read_one_tsv(path: Path) -> dict[str, str]:
    with path.open(newline="", encoding="utf-8") as handle:
        return next(csv.DictReader(handle, delimiter="\t"))


def test_main_snakefile_h5n1_to_sample_summary_integration():
    subprocess.run([sys.executable, str(PREPARE)], cwd=REPO_ROOT, check=True)

    targets = [
        f"tests/integration/work/results/{sample}/summary/{sample}.sample_summary.tsv"
        for sample in CASES
    ]

    command = [
        "snakemake",
        "--snakefile",
        "Snakefile",
        *targets,
        "--configfile",
        str(CONFIG),
        "--cores",
        "1",
        "--sdm",
        "conda",
        "--allowed-rules",
        "detect_h5n1",
        "sample_summary",
        "--rerun-incomplete",
    ]
    subprocess.run(command, cwd=REPO_ROOT, check=True)

    for sample, expected_status in CASES.items():
        flag = RESULTS / sample / "genoflu" / "h5n1.flag"
        summary = RESULTS / sample / "summary" / f"{sample}.sample_summary.tsv"
        assert flag.read_text(encoding="utf-8").strip() == expected_status
        row = read_one_tsv(summary)
        assert row["h5n1_screen"] == expected_status

    detected = read_one_tsv(RESULTS / "it_detected" / "summary" / "it_detected.sample_summary.tsv")
    assert "h5n1_screen_detected" in detected["review_flags"]

    not_detected = read_one_tsv(RESULTS / "it_not_detected" / "summary" / "it_not_detected.sample_summary.tsv")
    assert "h5n1_screen_detected" not in not_detected["review_flags"]
    assert "h5n1_screen_indeterminate" not in not_detected["review_flags"]

    indeterminate = read_one_tsv(RESULTS / "it_indeterminate" / "summary" / "it_indeterminate.sample_summary.tsv")
    assert "h5n1_screen_indeterminate" in indeterminate["review_flags"]
    assert "coverage_failures" in indeterminate["review_flags"]
    assert "multiple_irma_candidates" in indeterminate["review_flags"]
    assert indeterminate["multiple_irma_candidate_segments"] == "HA"
    assert indeterminate["max_irma_candidate_count"] == "2"
