from __future__ import annotations

import csv
import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
HERE = REPO_ROOT / "tests" / "integration9"
WORK = HERE / "work"
RESULTS = WORK / "results"
SAMPLE = "phase9"
MODEL = "dna_r10.4.1_e8.2_400bps_hac@v5.0.0"


def read_one_tsv(path: Path) -> dict[str, str]:
    with path.open(newline="", encoding="utf-8") as handle:
        return next(csv.DictReader(handle, delimiter="\t"))


def test_production_checkpoint_medaka_blast_chain():
    subprocess.run(
        [sys.executable, str(HERE / "prepare_fixtures.py")],
        cwd=REPO_ROOT,
        check=True,
    )

    # Create the valid IRMA-like BAM/FASTA project and tiny real BLAST database.
    subprocess.run(
        [
            "snakemake",
            "--snakefile",
            str(HERE / "prepare_assets.smk"),
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
        f"tests/integration9/work/results/{SAMPLE}/blast/HA.blast.txt"
    )

    # Run the real production chain. Full IRMA remains pre-seeded, while the
    # production checkpoint, QC, model resolution, Medaka, and BLAST rules run.
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
            "resolve_medaka_model",
            "medaka_inference",
            "medaka_consensus",
            "blastn",
            "--rerun-incomplete",
        ],
        cwd=REPO_ROOT,
        check=True,
    )

    sample_dir = RESULTS / SAMPLE

    manifest = read_one_tsv(sample_dir / "irma" / "manifest.tsv")
    assert manifest["segment"] == "HA"
    assert manifest["status"] == "READY"
    assert manifest["contig"] == "A_HA_H5_PHASE9"
    assert manifest["selection_status"] == "UNIQUE"

    qc = read_one_tsv(sample_dir / "coverage_stats" / "HA.tsv")
    assert qc["coverage_status"] == "PASS"
    assert qc["length_status"] == "PASS"
    assert qc["n_content_status"] == "PASS"
    assert qc["status"] == "PASS"
    assert qc["median_depth"] == "60.00"
    assert qc["breadth_covered"] == "1.000"

    assert (
        sample_dir / "coverage_flags" / "HA.flag"
    ).read_text(encoding="utf-8").strip() == "PASS"

    model = read_one_tsv(sample_dir / "medaka" / "model.tsv")
    assert model["basecaller_model"] == MODEL
    assert model["medaka_consensus_selector"] == f"{MODEL}:consensus"

    inference = read_one_tsv(
        sample_dir / "medaka" / "HA" / "inference.status.tsv"
    )
    assert inference["status"] == "SUCCESS"
    assert inference["model_source"] == "FASTQ metadata"
    assert inference["model"] == f"{MODEL}:consensus"

    features = sample_dir / "medaka" / "HA" / "features.hdf"
    assert features.stat().st_size > 0

    consensus_status = read_one_tsv(
        sample_dir / "medaka" / "HA" / "consensus.status.tsv"
    )
    assert consensus_status["status"] == "SUCCESS"
    assert consensus_status["consensus_source"] == "MEDAKA"

    consensus = sample_dir / "medaka" / "HA" / "consensus.fasta"
    assert consensus.stat().st_size > 0

    blast_path = sample_dir / "blast" / "HA.blast.txt"
    assert blast_path.stat().st_size > 0

    fields = blast_path.read_text(encoding="utf-8").splitlines()[0].split("\t")
    assert len(fields) == 12
    assert fields[1] == "PHASE9_H5_REFERENCE"
    assert float(fields[3]) >= 99.0

    qlen = int(fields[5])
    alignment_length = int(fields[4])
    qcov = alignment_length / qlen * 100.0
    assert qcov >= 99.0
