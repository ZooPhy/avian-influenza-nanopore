from __future__ import annotations

import csv
import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
HERE = REPO_ROOT / "tests" / "integration7"
WORK = HERE / "work"


def manifest_rows() -> dict[str, dict[str, str]]:
    path = WORK / "irma" / "manifest.tsv"
    with path.open(newline="", encoding="utf-8") as handle:
        return {
            row["segment"]: row
            for row in csv.DictReader(handle, delimiter="\t")
        }


def test_checkpoint_aware_irma_normalization():
    subprocess.run(
        [sys.executable, str(HERE / "prepare_fixtures.py")],
        cwd=REPO_ROOT,
        check=True,
    )

    target = str(WORK / "downstream" / "ready_segments.tsv")

    subprocess.run(
        [
            "snakemake",
            "--snakefile",
            str(HERE / "Snakefile"),
            target,
            "--cores",
            "1",
            "--sdm",
            "conda",
            "--rerun-incomplete",
        ],
        cwd=REPO_ROOT,
        check=True,
    )

    rows = manifest_rows()
    assert list(rows) == ["HA", "NA", "PB2", "PB1", "PA", "NP", "MP", "NS"]

    ha = rows["HA"]
    assert ha["status"] == "READY"
    assert ha["contig"] == "A_HA_H5_HIGH"
    assert ha["candidate_count"] == "2"
    assert ha["selection_status"] == "MULTIPLE_CANDIDATES"
    assert ha["coverage_median"] == "100"

    na = rows["NA"]
    assert na["status"] == "FASTA_ONLY"
    assert na["contig"] == "A_NA_N1"
    assert na["candidate_count"] == "1"
    assert na["selection_status"] == "UNIQUE"

    pb2 = rows["PB2"]
    assert pb2["status"] == "BAM_ONLY"
    assert pb2["contig"] == "A_PB2_TEST"

    for segment in ("PB1", "PA", "NP", "MP", "NS"):
        assert rows[segment]["status"] == "MISSING"
        assert rows[segment]["candidate_count"] == "0"
        assert rows[segment]["selection_status"] == "MISSING"

    ha_dir = WORK / "irma" / "segments" / "HA"
    assert (ha_dir / "consensus.fasta").is_file()
    assert (ha_dir / "alignment.bam").is_file()
    assert (ha_dir / "irma_coverage.tsv").is_file()

    # This output can only be scheduled after Snakemake reevaluates the DAG
    # through checkpoint.get() and discovers the READY segment files.
    downstream = (
        WORK / "downstream" / "ready_segments.tsv"
    ).read_text(encoding="utf-8").splitlines()
    assert downstream == ["segment", "HA"]

    log = (WORK / "irma" / "normalize.log").read_text(
        encoding="utf-8", errors="replace"
    )
    assert "ESCAPE_STATUS=IRMA_NORMALIZATION_COMPLETED" in log
    assert "ESCAPE_READY_SEGMENT_COUNT=1" in log
    assert "ESCAPE_READY_SEGMENTS=HA" in log
