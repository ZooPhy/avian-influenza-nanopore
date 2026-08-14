from __future__ import annotations

import csv
import json
import os
import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
HERE = REPO_ROOT / "tests" / "integration13"
WORK = HERE / "work"
RESULTS = WORK / "results"
SAMPLE = "release_smoke"
SAMPLE_DIR = RESULTS / SAMPLE
BUNDLE = RESULTS / "wings_report_bundle.wings"


def run(cmd: list[str], env: dict[str, str] | None = None) -> None:
    subprocess.run(cmd, cwd=REPO_ROOT, env=env, check=True)


def read_one_tsv(path: Path) -> dict[str, str]:
    with path.open(newline="", encoding="utf-8") as handle:
        return next(csv.DictReader(handle, delimiter="\t"))


def test_end_to_end_release_smoke():
    run([sys.executable, str(HERE / "prepare_fixtures.py")])

    # Stage A: real production raw-read preprocessing. NanoPlot is independent
    # of the main assembly chain, so target it explicitly alongside seqtk rename.
    renamed = (
        "tests/integration13/work/results/release_smoke/"
        "fastplong/filtered_renamed.fastq.gz"
    )
    nanoplot_done = (
        "tests/integration13/work/results/release_smoke/nanoplot/done.txt"
    )
    run(
        [
            "snakemake",
            "--snakefile",
            "Snakefile",
            renamed,
            nanoplot_done,
            "--configfile",
            str(HERE / "config.yaml"),
            "--cores",
            "1",
            "--sdm",
            "conda",
            "--resources",
            "kaleido=1",
            "--rerun-incomplete",
        ]
    )

    assert (SAMPLE_DIR / "porechop" / "trimmed.fastq").stat().st_size > 0
    assert (SAMPLE_DIR / "fastplong" / "filtered.fastq.gz").stat().st_size > 0
    assert (SAMPLE_DIR / "fastplong" / "report.json").stat().st_size > 0
    assert (SAMPLE_DIR / "fastplong" / "filtered_renamed.fastq.gz").stat().st_size > 0
    assert (SAMPLE_DIR / "nanoplot" / "done.txt").stat().st_size > 0

    # Stage B: intentional integration boundary. Full IRMA assembly is the one
    # expensive production stage replaced with a deterministic IRMA-like project.
    # A tiny real BLAST database is built at the same time.
    run(
        [
            "snakemake",
            "--snakefile",
            str(HERE / "prepare_assets.smk"),
            "--cores",
            "1",
            "--sdm",
            "conda",
            "--rerun-incomplete",
        ]
    )

    project = SAMPLE_DIR / "irma" / "project"
    assert project.is_dir()

    # Stage C: drive the highest-level production release target. VADR is an
    # external boundary and is represented by a tiny executable stub; the
    # production VADR rule itself still executes its gating/runtime logic.
    env = os.environ.copy()
    env["PATH"] = str(HERE / "bin") + os.pathsep + env.get("PATH", "")

    run(
        [
            "snakemake",
            "--snakefile",
            "Snakefile",
            "tests/integration13/work/results/wings_report_bundle.wings",
            "--configfile",
            str(HERE / "config.yaml"),
            "--cores",
            "1",
            "--sdm",
            "conda",
            "--resources",
            "kaleido=1",
            "--rerun-incomplete",
        ],
        env=env,
    )

    # Metadata validation/extraction ran from the original user-facing metadata.
    validated = RESULTS / "metadata" / "validated_metadata.tsv"
    sample_metadata = SAMPLE_DIR / "metadata" / f"{SAMPLE}.metadata.tsv"
    assert validated.stat().st_size > 0
    assert sample_metadata.stat().st_size > 0
    metadata = read_one_tsv(sample_metadata)
    assert metadata["project_code"] == "P13-E2E"
    assert metadata["host_common_name"] == "Mallard"

    # The preseeded IRMA project must have been reused, not replaced by a real
    # IRMA invocation.
    assert not (SAMPLE_DIR / "irma" / "irma.log").exists()

    manifest = SAMPLE_DIR / "irma" / "manifest.tsv"
    assert manifest.stat().st_size > 0

    coverage = SAMPLE_DIR / "coverage" / "coverage.tsv"
    assert coverage.stat().st_size > 0

    ha_flag = (SAMPLE_DIR / "coverage_flags" / "HA.flag").read_text().strip()
    na_flag = (SAMPLE_DIR / "coverage_flags" / "NA.flag").read_text().strip()
    assert ha_flag == "PASS"
    assert na_flag == "FAIL"

    # Only HA should receive real Medaka polishing; low-depth segments are
    # explicitly skipped by production QC gating.
    ha_inference = read_one_tsv(
        SAMPLE_DIR / "medaka" / "HA" / "inference.status.tsv"
    )
    na_inference = read_one_tsv(
        SAMPLE_DIR / "medaka" / "NA" / "inference.status.tsv"
    )
    ha_consensus = read_one_tsv(
        SAMPLE_DIR / "medaka" / "HA" / "consensus.status.tsv"
    )

    assert ha_inference["status"] == "SUCCESS"
    assert na_inference["status"] == "SKIPPED_QC"
    assert ha_consensus["status"] == "SUCCESS"
    assert ha_consensus["consensus_source"] == "MEDAKA"

    blast_summary = SAMPLE_DIR / "summary" / "blast_top_hits.csv"
    assert blast_summary.stat().st_size > 0
    blast_text = blast_summary.read_text(encoding="utf-8")
    assert "PHASE13_H5_REFERENCE" in blast_text
    assert "SKIPPED_QC" in blast_text

    h5n1 = (SAMPLE_DIR / "genoflu" / "h5n1.flag").read_text().strip()
    assert h5n1 == "INDETERMINATE"

    # GenoFLU's executable is not invoked for an indeterminate screen; the
    # production gate records the state instead.
    genoflu = read_one_tsv(SAMPLE_DIR / "genoflu" / "GenoFLU.tsv")
    assert genoflu["status"] == "H5N1_INDETERMINATE"

    vadr_log = SAMPLE_DIR / "vadr" / f"{SAMPLE}.vadr.log"
    assert "PHASE13_VADR_STUB" in vadr_log.read_text(
        encoding="utf-8", errors="replace"
    )

    sample_tsv = SAMPLE_DIR / "summary" / f"{SAMPLE}.sample_summary.tsv"
    sample_html = SAMPLE_DIR / "summary" / f"{SAMPLE}.sample_summary.html"
    run_html = RESULTS / "run_summary" / "run_summary.html"
    provenance = RESULTS / "run_summary" / "run_provenance.json"

    for path in (sample_tsv, sample_html, run_html, provenance, BUNDLE):
        assert path.stat().st_size > 0, path

    summary = read_one_tsv(sample_tsv)
    assert summary["sample"] == SAMPLE
    assert summary["project_code"] == "P13-E2E"
    assert summary["segments_pass"] == "1"
    assert summary["h5n1_screen"] == "INDETERMINATE"
    assert summary["genoflu_status"] == "H5N1_INDETERMINATE"
    assert summary["consensus_segments"] == "1"

    sample_html_text = sample_html.read_text(encoding="utf-8", errors="replace")
    run_html_text = run_html.read_text(encoding="utf-8", errors="replace")
    assert "WINGS" in sample_html_text
    assert "Mallard" in sample_html_text
    assert "INDETERMINATE" in sample_html_text
    assert "WINGS" in run_html_text
    assert "Run Summary Report" in run_html_text
    assert "Partial" in run_html_text

    bundle = json.loads(BUNDLE.read_text(encoding="utf-8"))
    assert bundle["format"] == "WINGS_REPORT_BUNDLE"
    assert bundle["version"] == 1
    assert bundle["sample_count"] == 1
    assert set(bundle["samples"]) == {SAMPLE}
    assert "WINGS" in bundle["run_summary"]["html"]
    assert "Mallard" in bundle["samples"][SAMPLE]["html"]

    embedded_provenance = bundle["provenance"]["data"]
    assert embedded_provenance["workflow"]["name"] == "WINGS"
    assert embedded_provenance["workflow"]["sample_count"] == 1
