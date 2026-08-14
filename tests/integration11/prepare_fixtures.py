#!/usr/bin/env python3
from __future__ import annotations

import csv
import gzip
import shutil
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
HERE = ROOT / "tests" / "integration11"
WORK = HERE / "work"
DATA = WORK / "data"
RESULTS = WORK / "results"

SEGMENTS = ("HA", "NA", "PB2", "PB1", "PA", "NP", "MP", "NS")

SAMPLES = {
    "complete": {
        "segments_pass": 8,
        "host": "Mallard",
        "host_species": "Anas platyrhynchos",
        "sample_type": "oropharyngeal swab",
        "collection_date": "2026-01-10",
        "state": "Arizona",
        "flyway": "Pacific",
        "h5n1": "NOT_DETECTED",
        "ha_contig": "A_HA_H3_COMPLETE",
        "na_contig": "A_NA_N2_COMPLETE",
        "ha_depth": "120.00",
        "na_depth": "115.00",
        "review_flags": "NONE",
        "project_code": "RUN11-A",
    },
    "near_complete": {
        "segments_pass": 7,
        "host": "Northern Pintail",
        "host_species": "Anas acuta",
        "sample_type": "cloacal swab",
        "collection_date": "2026-01-11",
        "state": "Arizona",
        "flyway": "Pacific",
        "h5n1": "DETECTED",
        "ha_contig": "A_HA_H5_NEAR",
        "na_contig": "A_NA_N1_NEAR",
        "ha_depth": "98.00",
        "na_depth": "91.00",
        "review_flags": "fewer_than_8_pass_segments;h5n1_screen_detected;coverage_failures",
        "project_code": "RUN11-B",
    },
    "partial": {
        "segments_pass": 3,
        "host": "Green-winged Teal",
        "host_species": "Anas crecca",
        "sample_type": "environmental swab",
        "collection_date": "2026-01-12",
        "state": "Nevada",
        "flyway": "Pacific",
        "h5n1": "INDETERMINATE",
        "ha_contig": "A_HA_H7_PARTIAL",
        "na_contig": "A_NA_N3_PARTIAL",
        "ha_depth": "72.00",
        "na_depth": "18.00",
        "review_flags": "fewer_than_8_pass_segments;h5n1_screen_indeterminate;coverage_failures",
        "project_code": "RUN11-C",
    },
    "failed": {
        "segments_pass": 0,
        "host": "Unknown wild bird",
        "host_species": "Unknown",
        "sample_type": "swab",
        "collection_date": "2026-01-13",
        "state": "California",
        "flyway": "Pacific",
        "h5n1": "INDETERMINATE",
        "ha_contig": "NA",
        "na_contig": "NA",
        "ha_depth": "NA",
        "na_depth": "NA",
        "review_flags": "fewer_than_8_pass_segments;h5n1_screen_indeterminate;coverage_failures;no_consensus_segments",
        "project_code": "RUN11-D",
    },
}


SUMMARY_FIELDS = [
    "sample", "sample_id", "host", "host_common_name", "host_species",
    "sample_type", "collection_date", "state", "flyway", "project_code",
    "reads_before", "reads_after", "bases_before", "bases_after",
    "q20_rate_before", "q20_rate_after", "q30_rate_before", "q30_rate_after",
    "reads_passed", "reads_low_quality", "reads_too_short", "reads_with_adapters",
    "segments_detected", "segments_pass", "pass_segment_names",
    "failed_segment_names", "ha_contig", "na_contig", "ha_median_depth",
    "na_median_depth", "h5n1_screen", "genoflu_status", "consensus_segments",
    "multiple_irma_candidate_segments", "multiple_irma_candidate_count",
    "max_irma_candidate_count", "review_flags", "ha_top_blast_hit",
    "na_top_blast_hit",
]


def write_tsv(path: Path, fields: list[str], rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    if WORK.exists():
        shutil.rmtree(WORK)
    DATA.mkdir(parents=True)

    # FASTQs are used only by the production Snakefile to discover SAMPLES.
    for sample in SAMPLES:
        with gzip.open(DATA / f"{sample}.fastq.gz", "wt", encoding="utf-8") as handle:
            handle.write("@fixture\nACGT\n+\nIIII\n")

    metadata_fields = [
        "sample_id", "host_common_name", "host_species", "sample_type",
        "collection_date", "state", "country", "flyway", "project_code",
    ]
    metadata_rows = []
    for sample, spec in SAMPLES.items():
        metadata_rows.append({
            "sample_id": sample,
            "host_common_name": spec["host"],
            "host_species": spec["host_species"],
            "sample_type": spec["sample_type"],
            "collection_date": spec["collection_date"],
            "state": spec["state"],
            "country": "USA",
            "flyway": spec["flyway"],
            "project_code": spec["project_code"],
        })

    write_tsv(WORK / "metadata.tsv", metadata_fields, metadata_rows)
    write_tsv(
        RESULTS / "metadata" / "validated_metadata.tsv",
        metadata_fields,
        metadata_rows,
    )

    for sample, spec in SAMPLES.items():
        sample_dir = RESULTS / sample
        pass_count = spec["segments_pass"]
        pass_segments = list(SEGMENTS[:pass_count])
        failed_segments = list(SEGMENTS[pass_count:])

        summary_row = {
            "sample": sample,
            "sample_id": sample,
            "host": spec["host"],
            "host_common_name": spec["host"],
            "host_species": spec["host_species"],
            "sample_type": spec["sample_type"],
            "collection_date": spec["collection_date"],
            "state": spec["state"],
            "flyway": spec["flyway"],
            "project_code": spec["project_code"],
            "reads_before": 1000,
            "reads_after": 900,
            "bases_before": 1500000,
            "bases_after": 1400000,
            "q20_rate_before": 0.90,
            "q20_rate_after": 0.95,
            "q30_rate_before": 0.70,
            "q30_rate_after": 0.80,
            "reads_passed": 900,
            "reads_low_quality": 60,
            "reads_too_short": 40,
            "reads_with_adapters": 0,
            "segments_detected": 8,
            "segments_pass": pass_count,
            "pass_segment_names": ",".join(pass_segments) if pass_segments else "NONE",
            "failed_segment_names": ",".join(failed_segments) if failed_segments else "NONE",
            "ha_contig": spec["ha_contig"],
            "na_contig": spec["na_contig"],
            "ha_median_depth": spec["ha_depth"],
            "na_median_depth": spec["na_depth"],
            "h5n1_screen": spec["h5n1"],
            "genoflu_status": "NOT_RUN",
            "consensus_segments": pass_count,
            "multiple_irma_candidate_segments": "NONE",
            "multiple_irma_candidate_count": 0,
            "max_irma_candidate_count": 1,
            "review_flags": spec["review_flags"],
            "ha_top_blast_hit": "fixture HA hit" if pass_count else "SKIPPED_QC",
            "na_top_blast_hit": "fixture NA hit" if pass_count >= 2 else "SKIPPED_QC",
        }
        write_tsv(
            sample_dir / "summary" / f"{sample}.sample_summary.tsv",
            SUMMARY_FIELDS,
            [summary_row],
        )

        report = sample_dir / "summary" / f"{sample}.sample_summary.html"
        report.parent.mkdir(parents=True, exist_ok=True)
        report.write_text(
            f"<html><body><h1>WINGS Sample {sample}</h1></body></html>\n",
            encoding="utf-8",
        )

        coverage_fields = [
            "sample", "segment", "contig", "length", "coverage_positions",
            "median_depth", "mean_depth", "breadth_covered", "n_count",
            "n_fraction", "expected_length_min", "expected_length_max",
            "maximum_n_fraction", "coverage_status", "length_status",
            "n_content_status", "overall_status", "coverage_flag",
            "assembly_status", "candidate_count", "selection_status",
            "selection_reason", "qc_reason",
        ]
        coverage_rows = []
        for i, segment in enumerate(SEGMENTS):
            passed = i < pass_count
            coverage_rows.append({
                "sample": sample,
                "segment": segment,
                "contig": (
                    spec["ha_contig"] if segment == "HA"
                    else spec["na_contig"] if segment == "NA"
                    else f"A_{segment}_{sample.upper()}"
                ),
                "length": 1700 if segment == "HA" else 1400 if segment == "NA" else 1000,
                "coverage_positions": 1700 if segment == "HA" else 1400 if segment == "NA" else 1000,
                "median_depth": "100.00" if passed else "10.00",
                "mean_depth": "100.00" if passed else "10.00",
                "breadth_covered": "1.000" if passed else "0.000",
                "n_count": 0,
                "n_fraction": "0.000000",
                "expected_length_min": 1,
                "expected_length_max": 9999,
                "maximum_n_fraction": "0.0100",
                "coverage_status": "PASS" if passed else "FAIL",
                "length_status": "PASS",
                "n_content_status": "PASS",
                "overall_status": "PASS" if passed else "FAIL",
                "coverage_flag": "PASS" if passed else "FAIL",
                "assembly_status": "READY",
                "candidate_count": 1,
                "selection_status": "UNIQUE",
                "selection_reason": "fixture",
                "qc_reason": "phase11 fixture",
            })
        write_tsv(
            sample_dir / "coverage" / "coverage.tsv",
            coverage_fields,
            coverage_rows,
        )

        genoflu_dir = sample_dir / "genoflu"
        genoflu_dir.mkdir(parents=True, exist_ok=True)
        (genoflu_dir / "GenoFLU.tsv").write_text(
            f"sample\tstatus\n{sample}\tNOT_RUN\n",
            encoding="utf-8",
        )

        vadr_dir = sample_dir / "vadr"
        vadr_dir.mkdir(parents=True, exist_ok=True)
        (vadr_dir / f"{sample}.vadr.log").write_text(
            "Phase 11 run-summary fixture; VADR not executed.\n",
            encoding="utf-8",
        )

        for i, segment in enumerate(SEGMENTS):
            passed = i < pass_count
            medaka_dir = sample_dir / "medaka" / segment
            medaka_dir.mkdir(parents=True, exist_ok=True)

            if passed:
                inference_status = "SUCCESS"
                inference_reason = "fixture_success"
                consensus_status = "SUCCESS"
                consensus_source = "MEDAKA"
                consensus_reason = "fixture_polished"
                variant_status = "SUCCESS"
                variant_reason = "fixture_vcf"
            else:
                inference_status = "SKIPPED_QC"
                inference_reason = "segment_qc_failed"
                consensus_status = "SKIPPED_QC"
                consensus_source = "NONE"
                consensus_reason = "segment_qc_failed"
                variant_status = "SKIPPED_QC"
                variant_reason = "segment_qc_failed"

            (medaka_dir / "inference.status.tsv").write_text(
                "status\tmodel_source\tmodel\treason\n"
                f"{inference_status}\tFASTQ metadata\tfixture:consensus\t{inference_reason}\n",
                encoding="utf-8",
            )
            (medaka_dir / "consensus.status.tsv").write_text(
                "status\tconsensus_source\treason\n"
                f"{consensus_status}\t{consensus_source}\t{consensus_reason}\n",
                encoding="utf-8",
            )
            (medaka_dir / "variants.status.tsv").write_text(
                "status\treason\n"
                f"{variant_status}\t{variant_reason}\n",
                encoding="utf-8",
            )

    print(
        "Prepared Phase 11 run-summary fixtures for: "
        + ", ".join(SAMPLES)
    )


if __name__ == "__main__":
    main()
