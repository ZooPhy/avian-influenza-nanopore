from __future__ import annotations

import csv
import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
HERE = REPO_ROOT / "tests" / "integration6"
WORK = HERE / "work"
MODEL = "dna_r10.4.1_e8.2_400bps_hac@v5.0.0"


def read_one_tsv(path: Path) -> dict[str, str]:
    with path.open(newline="", encoding="utf-8") as handle:
        return next(csv.DictReader(handle, delimiter="\t"))


def fasta_sequence(path: Path) -> str:
    lines = path.read_text(encoding="utf-8").splitlines()
    return "".join(line.strip() for line in lines if not line.startswith(">"))


def test_real_medaka_inference_and_sequence():
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

    model = read_one_tsv(WORK / "medaka" / "model.tsv")
    assert model["basecaller_model"] == MODEL
    assert model["medaka_consensus_selector"] == f"{MODEL}:consensus"

    bam = WORK / "HA.bam"
    bai = WORK / "HA.bam.bai"
    assert bam.stat().st_size > 0
    assert bai.stat().st_size > 0

    features = WORK / "medaka" / "HA" / "features.hdf"
    assert features.stat().st_size > 0

    inference = read_one_tsv(
        WORK / "medaka" / "HA" / "inference.status.tsv"
    )
    assert inference["status"] == "SUCCESS"
    assert inference["model_source"] == "FASTQ metadata"
    assert inference["model"] == f"{MODEL}:consensus"
    assert inference["reason"] == "inference_completed"

    consensus_path = WORK / "medaka" / "HA" / "consensus.fasta"
    consensus = fasta_sequence(consensus_path)
    assert len(consensus) >= 1200
    assert set(consensus.upper()) <= set("ACGTN")

    consensus_status = read_one_tsv(
        WORK / "medaka" / "HA" / "consensus.status.tsv"
    )
    assert consensus_status["status"] == "SUCCESS"
    assert consensus_status["consensus_source"] == "MEDAKA"
    assert consensus_status["reason"] == "polishing_completed"

    inference_log = (
        WORK / "medaka" / "HA" / "medaka_inference.log"
    ).read_text(encoding="utf-8", errors="replace")
    assert "Medaka consensus selector:" in inference_log
