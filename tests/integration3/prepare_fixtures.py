#!/usr/bin/env python3
from __future__ import annotations

import csv
import gzip
import json
import shutil
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
WORK = ROOT / "tests" / "integration3" / "work"
DATA = WORK / "data"
RESULTS = WORK / "results"
SAMPLE = "smoke_mixed"

SEGMENTS = ("HA", "NA", "PB2", "PB1", "PA", "NP", "MP", "NS")


def write_tsv(path: Path, fieldnames: list[str], rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def write_csv(path: Path, fieldnames: list[str], rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    if WORK.exists():
        shutil.rmtree(WORK)

    DATA.mkdir(parents=True, exist_ok=True)
    sample_dir = RESULTS / SAMPLE

    with gzip.open(DATA / f"{SAMPLE}.fastq.gz", "wt", encoding="utf-8") as handle:
        handle.write(
            f"@{SAMPLE}_read1 basecall_model_version_id=dna_r10.4.1_e8.2_400bps_hac@v5.0.0\n"
            "ACGTACGTACGT\n+\nFFFFFFFFFFFF\n"
        )

    # Run-level metadata file exists only so config points to a realistic location.
    write_tsv(
        WORK / "metadata.tsv",
        [
            "sample_id",
            "host",
            "specimen_type",
            "collection_date",
            "state",
            "country",
        ],
        [{
            "sample_id": SAMPLE,
            "host": "MALL",
            "specimen_type": "Oral+Cloacal",
            "collection_date": "2026-01-15",
            "state": "Arizona",
            "country": "USA",
        }],
    )

    # Production sample_summary consumes the per-sample normalized metadata row.
    write_tsv(
        sample_dir / "metadata" / f"{SAMPLE}.metadata.tsv",
        [
            "sample_id",
            "host",
            "specimen_type",
            "collection_date",
            "state",
            "country",
            "host_common_name",
            "host_species",
            "flyway",
            "metadata_status",
        ],
        [{
            "sample_id": SAMPLE,
            "host": "MALL",
            "specimen_type": "Oral+Cloacal",
            "collection_date": "2026-01-15",
            "state": "Arizona",
            "country": "USA",
            "host_common_name": "Mallard",
            "host_species": "Anas platyrhynchos",
            "flyway": "Pacific",
            "metadata_status": "COMPLETE",
        }],
    )

    fastplong = {
        "summary": {
            "before_filtering": {
                "total_reads": 100,
                "total_bases": 100000,
                "q20_rate": 0.50,
                "q30_rate": 0.20,
            },
            "after_filtering": {
                "total_reads": 90,
                "total_bases": 90000,
                "q20_rate": 0.55,
                "q30_rate": 0.25,
            },
        },
        "filtering_result": {
            "passed_filter_reads": 90,
            "low_quality_reads": 5,
            "too_short_reads": 5,
        },
        "adapter_cutting": {"adapter_trimmed_reads": 0},
    }
    path = sample_dir / "fastplong" / "report.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(fastplong), encoding="utf-8")

    flags = sample_dir / "coverage_flags"
    stats = sample_dir / "coverage_stats"
    flags.mkdir(parents=True, exist_ok=True)
    stats.mkdir(parents=True, exist_ok=True)

    (flags / "HA.flag").write_text("PASS\n", encoding="utf-8")
    (flags / "NA.flag").write_text("FAIL\n", encoding="utf-8")

    stat_fields = ["segment", "selected_contig", "chosen_table"]
    write_tsv(
        stats / "HA.tsv",
        stat_fields,
        [{"segment": "HA", "selected_contig": "A_HA_H5", "chosen_table": "A_HA_H5"}],
    )
    write_tsv(
        stats / "NA.tsv",
        stat_fields,
        [{"segment": "NA", "selected_contig": "A_NA_N1", "chosen_table": "A_NA_N1"}],
    )

    coverage_rows = []
    for segment in SEGMENTS:
        status = "FAIL" if segment == "NA" else "PASS"
        contig = (
            "A_HA_H5"
            if segment == "HA"
            else "A_NA_N1"
            if segment == "NA"
            else f"A_{segment}_TEST"
        )
        coverage_rows.append({
            "segment": segment,
            "contig": contig,
            "coverage_flag": status,
            "median_depth": "100.00" if status == "PASS" else "10.00",
            "candidate_count": "1",
        })

    write_tsv(
        sample_dir / "coverage" / "coverage.tsv",
        ["segment", "contig", "coverage_flag", "median_depth", "candidate_count"],
        coverage_rows,
    )

    write_csv(
        sample_dir / "summary" / "blast_top_hits.csv",
        ["segment", "top_hit"],
        [
            {"segment": "HA", "top_hit": "TEST_H5_HIGH_CONFIDENCE"},
            {"segment": "NA", "top_hit": "SKIPPED_QC"},
        ],
    )

    genoflu_dir = sample_dir / "genoflu"
    genoflu_dir.mkdir(parents=True, exist_ok=True)
    (genoflu_dir / "GenoFLU.tsv").write_text(
        f"sample\tstatus\n{SAMPLE}\tH5N1_INDETERMINATE\n",
        encoding="utf-8",
    )

    merged = sample_dir / "merged" / "consensus_all_segments.fasta"
    merged.parent.mkdir(parents=True, exist_ok=True)
    with merged.open("w", encoding="utf-8") as handle:
        for i in range(7):
            handle.write(f">{SAMPLE}_{i+1}\nACGTACGT\n")

    # Pre-seed realistic Medaka status files. They are not consumed by
    # sample_summary.py, but make the fixture tree resemble a WINGS result tree.
    for segment in SEGMENTS:
        medaka_dir = sample_dir / "medaka" / segment
        medaka_dir.mkdir(parents=True, exist_ok=True)
        if segment == "NA":
            inference = "SKIPPED_QC"
            consensus = "SKIPPED_QC\tNONE"
            variants = "SKIPPED_QC"
        else:
            inference = "SUCCESS"
            consensus = "SUCCESS\tMEDAKA"
            variants = "SUCCESS"

        (medaka_dir / "inference.status.tsv").write_text(
            f"status\treason\n{inference}\tfixture\n",
            encoding="utf-8",
        )
        (medaka_dir / "consensus.status.tsv").write_text(
            f"status\tconsensus_source\treason\n{consensus}\tfixture\n",
            encoding="utf-8",
        )
        (medaka_dir / "variants.status.tsv").write_text(
            f"status\treason\n{variants}\tfixture\n",
            encoding="utf-8",
        )

    print(f"Prepared Phase 3 fixture under {WORK}")


if __name__ == "__main__":
    main()
