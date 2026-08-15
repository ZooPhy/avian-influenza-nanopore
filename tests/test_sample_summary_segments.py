from __future__ import annotations

import csv
import runpy
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]

sample_summary = runpy.run_path(
    str(REPO_ROOT / "scripts" / "sample_summary.py")
)
write_summary = sample_summary["write_summary"]


def test_segments_detected_excludes_missing_segment(tmp_path):
    coverage_rows = []

    for segment in ("HA", "NA", "PB2", "PB1", "PA", "NP", "MP", "NS"):
        missing = segment == "MP"

        coverage_rows.append(
            {
                "segment": segment,
                "coverage_flag": "MISSING" if missing else "PASS",
                "contig": "NA" if missing else f"{segment}_contig",
                "median_depth": "NA" if missing else "60.00",
                "candidate_count": "0" if missing else "1",
            }
        )

    output = tmp_path / "sample_summary.tsv"

    write_summary(
        output_path=output,
        sample="missing_mp",
        fastplong={},
        coverage_rows=coverage_rows,
        blast_hits={},
        h5n1_status="NOT_DETECTED",
        genoflu_status="DISABLED_BY_CONFIG",
        consensus_segments=7,
        metadata={},
    )

    with output.open(newline="", encoding="utf-8") as handle:
        row = next(csv.DictReader(handle, delimiter="\t"))

    assert row["segments_detected"] == "7"
    assert row["segments_pass"] == "7"
    assert row["segments_missing"] == "1"
    assert row["missing_segment_names"] == "MP"
    assert row["failed_segment_names"] == "NONE"
    assert "coverage_failures" not in row["review_flags"].split(";")
    assert "fewer_than_8_pass_segments" in row["review_flags"].split(";")
