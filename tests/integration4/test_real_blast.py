from __future__ import annotations

import csv
import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
HERE = REPO_ROOT / "tests" / "integration4"
WORK = HERE / "work"


def test_real_blast_smoke():
    subprocess.run(
        [sys.executable, str(HERE / "prepare_fixtures.py")],
        cwd=REPO_ROOT,
        check=True,
    )

    subprocess.run(
        [
            "snakemake",
            "--snakefile",
            str(HERE / "Snakefile"),
            "--cores",
            "1",
            "--sdm",
            "conda",
            "--rerun-incomplete",
        ],
        cwd=REPO_ROOT,
        check=True,
    )

    raw = WORK / "blast" / "HA.blast.txt"
    assert raw.is_file()
    assert raw.stat().st_size > 0

    first = raw.read_text(encoding="utf-8").splitlines()[0].split("\t")
    assert first[0] == "QUERY_HA"
    assert first[1] == "TEST_HA_REF"
    assert float(first[3]) == 100.0
    assert int(first[4]) == 1700
    assert int(first[5]) == 1700

    summary = WORK / "summary" / "blast_top_hits.csv"
    with summary.open(newline="", encoding="utf-8") as handle:
        rows = {row["segment"]: row for row in csv.DictReader(handle)}

    assert rows["HA"]["top_hit"] == "TEST_HA_REF"
    assert rows["HA"]["hit_status"] == "HIGH_CONFIDENCE"
    assert rows["HA"]["percent_identity"] == "100.000"
    assert rows["HA"]["query_coverage"] == "100.000"

    assert rows["NA"]["hit_status"] == "SKIPPED_QC"
    assert rows["NA"]["top_hit"] == "SKIPPED_QC"

    for segment in ("PB2", "PB1", "PA", "NP", "MP", "NS"):
        assert rows[segment]["hit_status"] == "SKIPPED_QC"
