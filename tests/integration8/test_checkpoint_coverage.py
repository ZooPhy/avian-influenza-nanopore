from __future__ import annotations

import csv
import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
HERE = REPO_ROOT / "tests" / "integration8"
WORK = HERE / "work"
RESULTS = WORK / "results"
SAMPLE = "qc_checkpoint"
ALL_SEGMENTS = ("HA", "NA", "PB2", "PB1", "PA", "NP", "MP", "NS")
MISSING_SEGMENT = "MP"
READY_SEGMENTS = tuple(segment for segment in ALL_SEGMENTS if segment != MISSING_SEGMENT)


def read_tsv_rows(path: Path) -> dict[str, dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return {
            row["segment"]: row
            for row in csv.DictReader(handle, delimiter="\t")
        }


def test_production_checkpoint_to_coverage_qc():
    subprocess.run(
        [sys.executable, str(HERE / "prepare_fixtures.py")],
        cwd=REPO_ROOT,
        check=True,
    )

    # Build valid IRMA-like BAM/FASTA artifacts with the production coverage
    # environment, but do not run IRMA itself. MP is intentionally absent.
    subprocess.run(
        [
            "snakemake",
            "--snakefile",
            str(HERE / "prepare_project.smk"),
            "--cores",
            "1",
            "--sdm",
            "conda",
            "--rerun-incomplete",
        ],
        cwd=REPO_ROOT,
        check=True,
    )

    target = (
        f"tests/integration8/work/results/{SAMPLE}/coverage/coverage.tsv"
    )

    # This is the real production Snakefile. IRMA is intentionally excluded:
    # its completed project directory is a pre-seeded fixture. The real
    # normalize_irma_outputs checkpoint, check_coverage, and coverage_table
    # rules are allowed to execute.
    subprocess.run(
        [
            "snakemake",
            "--snakefile",
            "Snakefile",
            target,
            "--configfile",
            str(HERE / "config.yaml"),
            "--cores",
            "1",
            "--sdm",
            "conda",
            "--allowed-rules",
            "normalize_irma_outputs",
            "check_coverage",
            "coverage_table",
            "--rerun-incomplete",
        ],
        cwd=REPO_ROOT,
        check=True,
    )

    manifest = read_tsv_rows(
        RESULTS / SAMPLE / "irma" / "manifest.tsv"
    )
    assert list(manifest) == list(ALL_SEGMENTS)

    for segment in READY_SEGMENTS:
        row = manifest[segment]
        assert row["status"] == "READY", segment
        assert row["candidate_count"] == "1", segment
        assert row["selection_status"] == "UNIQUE", segment

    missing_manifest = manifest[MISSING_SEGMENT]
    assert missing_manifest["status"] == "MISSING"
    assert missing_manifest["candidate_count"] == "0"
    assert missing_manifest["selection_status"] == "MISSING"

    # The normalizer must not invent placeholder assembly files for a segment
    # that was genuinely absent from the IRMA-like project.
    missing_segment_dir = RESULTS / SAMPLE / "irma" / "segments" / MISSING_SEGMENT
    assert not (missing_segment_dir / "consensus.fasta").exists()
    assert not (missing_segment_dir / "alignment.bam").exists()

    coverage = read_tsv_rows(
        RESULTS / SAMPLE / "coverage" / "coverage.tsv"
    )
    assert list(coverage) == list(ALL_SEGMENTS)
    assert len(coverage) == 8

    # Straight PASS.
    assert coverage["HA"]["coverage_status"] == "PASS"
    assert coverage["HA"]["length_status"] == "PASS"
    assert coverage["HA"]["n_content_status"] == "PASS"
    assert coverage["HA"]["overall_status"] == "PASS"
    assert coverage["HA"]["median_depth"] == "60.00"
    assert coverage["HA"]["breadth_covered"] == "1.000"

    # Coverage failure: full breadth exists, but only 20x, below the 50x gate.
    assert coverage["NA"]["coverage_status"] == "FAIL"
    assert coverage["NA"]["length_status"] == "PASS"
    assert coverage["NA"]["n_content_status"] == "PASS"
    assert coverage["NA"]["overall_status"] == "FAIL"
    assert coverage["NA"]["median_depth"] == "20.00"
    assert coverage["NA"]["breadth_covered"] == "0.000"

    # Upper length bound is warning-only and does not block the segment.
    assert coverage["PB2"]["coverage_status"] == "PASS"
    assert coverage["PB2"]["length_status"] == "WARNING"
    assert coverage["PB2"]["overall_status"] == "PASS"
    assert coverage["PB2"]["length"] == "2450"

    # Lower length bound is a hard failure.
    assert coverage["PB1"]["coverage_status"] == "PASS"
    assert coverage["PB1"]["length_status"] == "FAIL"
    assert coverage["PB1"]["overall_status"] == "FAIL"
    assert coverage["PB1"]["length"] == "2100"

    # Excess N content is a hard failure.
    assert coverage["PA"]["coverage_status"] == "PASS"
    assert coverage["PA"]["length_status"] == "PASS"
    assert coverage["PA"]["n_content_status"] == "FAIL"
    assert coverage["PA"]["overall_status"] == "FAIL"
    assert float(coverage["PA"]["n_fraction"]) > 0.01

    # NP and NS remain ordinary passing segments. NS is intentionally retained
    # to demonstrate that the missing-segment regression is not NS-specific.
    for segment in ("NP", "NS"):
        assert coverage[segment]["overall_status"] == "PASS", segment

    # The deliberately absent MP segment must survive the production QC path as
    # explicit MISSING data rather than causing Snakemake DAG construction to fail.
    missing_coverage = coverage[MISSING_SEGMENT]
    assert missing_coverage["assembly_status"] == "MISSING"
    assert missing_coverage["coverage_status"] == "MISSING"
    assert missing_coverage["length_status"] == "MISSING"
    assert missing_coverage["n_content_status"] == "MISSING"
    assert missing_coverage["overall_status"] == "MISSING"
    assert missing_coverage["coverage_flag"] == "MISSING"
    assert missing_coverage["candidate_count"] == "0"
    assert missing_coverage["selection_status"] == "MISSING"

    expected_flags = {
        "HA": "PASS",
        "NA": "FAIL",
        "PB2": "PASS",
        "PB1": "FAIL",
        "PA": "FAIL",
        "NP": "PASS",
        "MP": "MISSING",
        "NS": "PASS",
    }
    for segment, expected in expected_flags.items():
        observed = (
            RESULTS / SAMPLE / "coverage_flags" / f"{segment}.flag"
        ).read_text(encoding="utf-8").strip()
        assert observed == expected, segment

    normalize_log = (
        RESULTS / SAMPLE / "irma" / "normalize.log"
    ).read_text(encoding="utf-8", errors="replace")
    assert "ESCAPE_STATUS=IRMA_NORMALIZATION_COMPLETED" in normalize_log
    assert "ESCAPE_READY_SEGMENT_COUNT=7" in normalize_log
