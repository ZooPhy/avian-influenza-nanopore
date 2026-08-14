from __future__ import annotations

import csv
import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
CONFIG = REPO_ROOT / "tests" / "integration3" / "config.yaml"
PREPARE = REPO_ROOT / "tests" / "integration3" / "prepare_fixtures.py"
RESULTS = REPO_ROOT / "tests" / "integration3" / "work" / "results"
SAMPLE = "smoke_mixed"


def read_one_tsv(path: Path) -> dict[str, str]:
    with path.open(newline="", encoding="utf-8") as handle:
        return next(csv.DictReader(handle, delimiter="\t"))


def test_downstream_production_smoke_with_metadata():
    subprocess.run([sys.executable, str(PREPARE)], cwd=REPO_ROOT, check=True)

    target = (
        f"tests/integration3/work/results/{SAMPLE}/summary/"
        f"{SAMPLE}.sample_summary.tsv"
    )

    command = [
        "snakemake",
        "--snakefile",
        "Snakefile",
        target,
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

    flag = RESULTS / SAMPLE / "genoflu" / "h5n1.flag"
    summary = RESULTS / SAMPLE / "summary" / f"{SAMPLE}.sample_summary.tsv"

    assert flag.read_text(encoding="utf-8").strip() == "INDETERMINATE"

    row = read_one_tsv(summary)

    assert row["h5n1_screen"] == "INDETERMINATE"
    assert "h5n1_screen_indeterminate" in row["review_flags"]
    assert "coverage_failures" in row["review_flags"]

    assert row["sample_id"] == SAMPLE
    assert row["host"] == "MALL"
    assert row["host_common_name"] == "Mallard"
    assert row["host_species"] == "Anas platyrhynchos"
    assert row["collection_date"] == "2026-01-15"
    assert row["state"] == "Arizona"
    assert row["country"] == "USA"
    assert row["flyway"] == "Pacific"
    assert row["metadata_status"] == "COMPLETE"

    assert row["consensus_segments"] == "7"
    assert row["segments_pass"] == "7"
    assert row["failed_segment_names"] == "NA"
    assert row["ha_top_blast_hit"] == "TEST_H5_HIGH_CONFIDENCE"
    assert row["na_top_blast_hit"] == "SKIPPED_QC"

    na_inference = (
        RESULTS / SAMPLE / "medaka" / "NA" / "inference.status.tsv"
    ).read_text(encoding="utf-8")
    assert "SKIPPED_QC" in na_inference
