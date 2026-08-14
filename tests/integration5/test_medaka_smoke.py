from __future__ import annotations

import csv
import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
HERE = REPO_ROOT / "tests" / "integration5"
WORK = HERE / "work"
MODEL = "dna_r10.4.1_e8.2_400bps_hac@v5.0.0"


def read_one_tsv(path: Path) -> dict[str, str]:
    with path.open(newline="", encoding="utf-8") as handle:
        return next(csv.DictReader(handle, delimiter="\t"))


def test_medaka_model_and_status_smoke():
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

    version = (WORK / "medaka" / "tool_version.txt").read_text(
        encoding="utf-8", errors="replace"
    )
    assert "medaka" in version.lower()
    assert "2.2.2" in version

    model = read_one_tsv(WORK / "medaka" / "model.tsv")
    assert model["basecaller_model"] == MODEL
    assert model["medaka_consensus_selector"] == f"{MODEL}:consensus"
    assert model["medaka_variant_selector"] == f"{MODEL}:variant"
    assert model["records_with_model"] == "3"

    ha_inference = read_one_tsv(
        WORK / "medaka" / "HA" / "inference.status.tsv"
    )
    assert ha_inference["status"] == "FAILED"
    assert ha_inference["model_source"] == "FASTQ metadata"
    assert ha_inference["model"] == f"{MODEL}:consensus"
    assert ha_inference["reason"] == "medaka_inference_failed"

    ha_consensus_status = read_one_tsv(
        WORK / "medaka" / "HA" / "consensus.status.tsv"
    )
    assert ha_consensus_status["status"] == "FAILED"
    assert ha_consensus_status["consensus_source"] == "IRMA_FALLBACK"
    assert ha_consensus_status["reason"] == "medaka_inference_failed"

    assert (
        WORK / "medaka" / "HA" / "consensus.fasta"
    ).read_text(encoding="utf-8") == (
        WORK / "irma" / "HA.fasta"
    ).read_text(encoding="utf-8")

    na_inference = read_one_tsv(
        WORK / "medaka" / "NA" / "inference.status.tsv"
    )
    assert na_inference["status"] == "SKIPPED_QC"
    assert na_inference["model"] == f"{MODEL}:consensus"
    assert na_inference["reason"] == "segment_qc_failed"

    na_consensus_status = read_one_tsv(
        WORK / "medaka" / "NA" / "consensus.status.tsv"
    )
    assert na_consensus_status["status"] == "SKIPPED_QC"
    assert na_consensus_status["consensus_source"] == "NONE"
    assert na_consensus_status["reason"] == "segment_qc_failed"

    assert (WORK / "medaka" / "NA" / "features.hdf").stat().st_size == 0
    assert (WORK / "medaka" / "NA" / "consensus.fasta").stat().st_size == 0
