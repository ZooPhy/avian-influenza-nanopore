#!/usr/bin/env python3
"""Create tiny deterministic inputs for WINGS Snakemake integration tests."""

from __future__ import annotations

import csv
import gzip
import json
import shutil
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
WORK = ROOT / "tests" / "integration" / "work"
DATA = WORK / "data"
RESULTS = WORK / "results"

SEGMENTS = ("HA", "NA", "PB2", "PB1", "PA", "NP", "MP", "NS")

CASES = {
    "it_detected": {
        "ha_flag": "PASS",
        "na_flag": "PASS",
        "ha_contig": "A_HA_H5",
        "na_contig": "A_NA_N1",
        "candidate_count": 1,
    },
    "it_not_detected": {
        "ha_flag": "PASS",
        "na_flag": "PASS",
        "ha_contig": "A_HA_H7",
        "na_contig": "A_NA_N1",
        "candidate_count": 1,
    },
    "it_indeterminate": {
        "ha_flag": "FAIL",
        "na_flag": "PASS",
        "ha_contig": "A_HA_H5",
        "na_contig": "A_NA_N1",
        "candidate_count": 2,
    },
}


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


def make_dummy_fastq(sample: str) -> None:
    DATA.mkdir(parents=True, exist_ok=True)
    path = DATA / f"{sample}.fastq.gz"
    with gzip.open(path, "wt", encoding="utf-8") as handle:
        handle.write(
            f"@{sample}_read1 basecall_model_version_id=dna_r10.4.1_e8.2_400bps_hac@v5.0.0\n"
            "ACGTACGTACGT\n+\nFFFFFFFFFFFF\n"
        )


def make_case(sample: str, case: dict) -> None:
    sample_dir = RESULTS / sample
    make_dummy_fastq(sample)

    flags_dir = sample_dir / "coverage_flags"
    stats_dir = sample_dir / "coverage_stats"
    flags_dir.mkdir(parents=True, exist_ok=True)
    stats_dir.mkdir(parents=True, exist_ok=True)
    (flags_dir / "HA.flag").write_text(case["ha_flag"] + "\n", encoding="utf-8")
    (flags_dir / "NA.flag").write_text(case["na_flag"] + "\n", encoding="utf-8")

    stat_fields = ["segment", "selected_contig", "chosen_table"]
    write_tsv(
        stats_dir / "HA.tsv",
        stat_fields,
        [{"segment": "HA", "selected_contig": case["ha_contig"], "chosen_table": case["ha_contig"]}],
    )
    write_tsv(
        stats_dir / "NA.tsv",
        stat_fields,
        [{"segment": "NA", "selected_contig": case["na_contig"], "chosen_table": case["na_contig"]}],
    )

    metadata_dir = sample_dir / "metadata"
    metadata_dir.mkdir(parents=True, exist_ok=True)
    write_tsv(
        metadata_dir / f"{sample}.metadata.tsv",
        ["sample_id", "host", "collection_date", "country"],
        [{"sample_id": sample, "host": "TEST", "collection_date": "2026-01-01", "country": "USA"}],
    )

    fastplong_dir = sample_dir / "fastplong"
    fastplong_dir.mkdir(parents=True, exist_ok=True)
    fastplong = {
        "summary": {
            "before_filtering": {"total_reads": 100, "total_bases": 100000, "q20_rate": 0.50, "q30_rate": 0.20},
            "after_filtering": {"total_reads": 90, "total_bases": 90000, "q20_rate": 0.55, "q30_rate": 0.25},
        },
        "filtering_result": {"passed_filter_reads": 90, "low_quality_reads": 5, "too_short_reads": 5},
        "adapter_cutting": {"adapter_trimmed_reads": 0},
    }
    (fastplong_dir / "report.json").write_text(json.dumps(fastplong), encoding="utf-8")

    coverage_rows = []
    for segment in SEGMENTS:
        flag = "FAIL" if sample == "it_indeterminate" and segment == "HA" else "PASS"
        contig = case["ha_contig"] if segment == "HA" else case["na_contig"] if segment == "NA" else f"A_{segment}_TEST"
        coverage_rows.append(
            {
                "segment": segment,
                "contig": contig,
                "coverage_flag": flag,
                "median_depth": "100.00",
                "candidate_count": str(case["candidate_count"] if segment == "HA" else 1),
            }
        )
    write_tsv(
        sample_dir / "coverage" / "coverage.tsv",
        ["segment", "contig", "coverage_flag", "median_depth", "candidate_count"],
        coverage_rows,
    )

    write_csv(
        sample_dir / "summary" / "blast_top_hits.csv",
        ["segment", "top_hit"],
        [{"segment": "HA", "top_hit": "TEST_HA"}, {"segment": "NA", "top_hit": "TEST_NA"}],
    )

    genoflu_dir = sample_dir / "genoflu"
    genoflu_dir.mkdir(parents=True, exist_ok=True)
    (genoflu_dir / "GenoFLU.tsv").write_text(f"sample\tstatus\n{sample}\tTEST_FIXTURE\n", encoding="utf-8")

    merged_dir = sample_dir / "merged"
    merged_dir.mkdir(parents=True, exist_ok=True)
    records = 7 if sample == "it_indeterminate" else 8
    with (merged_dir / "consensus_all_segments.fasta").open("w", encoding="utf-8") as handle:
        for i in range(records):
            handle.write(f">{sample}_{i+1}\nACGTACGT\n")


def main() -> None:
    if WORK.exists():
        shutil.rmtree(WORK)
    for sample, case in CASES.items():
        make_case(sample, case)
    print(f"Prepared WINGS integration fixtures under {WORK}")


if __name__ == "__main__":
    main()
