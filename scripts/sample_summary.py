#!/usr/bin/env python3

import csv
import json
from pathlib import Path


def read_fastplong_json(path: Path) -> dict:
    with path.open() as handle:
        data = json.load(handle)

    summary = data.get("summary", {})
    before = summary.get("before_filtering", {})
    after = summary.get("after_filtering", {})

    filtering = data.get("filtering_result", {})
    adapter_cutting = data.get("adapter_cutting", {})

    return {
        "reads_before": before.get("total_reads", "NA"),
        "bases_before": before.get("total_bases", "NA"),
        "q20_rate_before": before.get("q20_rate", "NA"),
        "q30_rate_before": before.get("q30_rate", "NA"),
        "reads_after": after.get("total_reads", "NA"),
        "bases_after": after.get("total_bases", "NA"),
        "q20_rate_after": after.get("q20_rate", "NA"),
        "q30_rate_after": after.get("q30_rate", "NA"),
        "reads_passed": filtering.get("passed_filter_reads", "NA"),
        "reads_low_quality": filtering.get("low_quality_reads", "NA"),
        "reads_too_short": filtering.get("too_short_reads", "NA"),
        "reads_with_adapters": adapter_cutting.get(
            "adapter_trimmed_reads",
            "NA",
        ),
    }


def read_coverage(path: Path) -> list[dict]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def read_blast(path: Path) -> dict[str, str]:
    hits = {}
    with path.open(newline="") as handle:
        for row in csv.DictReader(handle):
            hits[row["segment"]] = row["top_hit"]
    return hits


def read_h5n1_flag(path: Path) -> str:
    if not path.exists():
        return "MISSING"
    return path.read_text().strip() or "MISSING"


def read_genoflu(path: Path) -> str:
    if not path.exists() or path.stat().st_size == 0:
        return "MISSING"

    lines = [line.strip() for line in path.read_text().splitlines() if line.strip()]
    if not lines:
        return "MISSING"

    if len(lines) == 2 and lines[0].startswith("sample\tstatus"):
        return lines[1].split("\t", 1)[-1]

    return "COMPLETED"


def count_fasta_records(path: Path) -> int:
    if not path.exists():
        return 0

    count = 0
    with path.open() as handle:
        for line in handle:
            if line.startswith(">"):
                count += 1
    return count


def write_summary(
    output_path: Path,
    sample: str,
    fastplong: dict,
    coverage_rows: list[dict],
    blast_hits: dict[str, str],
    h5n1_status: str,
    genoflu_status: str,
    consensus_segments: int,
) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)

    pass_segments = [
        row["segment"]
        for row in coverage_rows
        if row.get("coverage_flag") == "PASS"
    ]

    failed_segments = [
        row["segment"]
        for row in coverage_rows
        if row.get("coverage_flag") not in {"PASS", ""}
    ]

    ha_rows = [row for row in coverage_rows if row.get("segment") == "HA"]
    na_rows = [row for row in coverage_rows if row.get("segment") == "NA"]

    ha_contig = ha_rows[0].get("contig", "NA") if ha_rows else "NA"
    na_contig = na_rows[0].get("contig", "NA") if na_rows else "NA"

    ha_median_depth = (
        ha_rows[0].get("median_depth", "NA") if ha_rows else "NA"
    )
    na_median_depth = (
        na_rows[0].get("median_depth", "NA") if na_rows else "NA"
    )

    review_flags = []

    if len(pass_segments) < 8:
        review_flags.append("fewer_than_8_pass_segments")

    if h5n1_status != "PASS":
        review_flags.append("h5n1_screen_not_pass")

    if failed_segments:
        review_flags.append("coverage_failures")

    if consensus_segments == 0:
        review_flags.append("no_consensus_segments")

    fieldnames = [
        "sample",
        "reads_before",
        "reads_after",
        "bases_before",
        "bases_after",
        "q20_rate_before",
        "q20_rate_after",
        "q30_rate_before",
        "q30_rate_after",
        "reads_passed",
        "reads_low_quality",
        "reads_too_short",
        "reads_with_adapters",
        "segments_detected",
        "segments_pass",
        "pass_segment_names",
        "failed_segment_names",
        "ha_contig",
        "na_contig",
        "ha_median_depth",
        "na_median_depth",
        "h5n1_screen",
        "genoflu_status",
        "consensus_segments",
        "review_flags",
        "ha_top_blast_hit",
        "na_top_blast_hit",
    ]

    row = {
        "sample": sample,
        **fastplong,
        "segments_detected": len(coverage_rows),
        "segments_pass": len(pass_segments),
        "pass_segment_names": ",".join(pass_segments) or "NONE",
        "failed_segment_names": ",".join(failed_segments) or "NONE",
        "ha_contig": ha_contig,
        "na_contig": na_contig,
        "ha_median_depth": ha_median_depth,
        "na_median_depth": na_median_depth,
        "h5n1_screen": h5n1_status,
        "genoflu_status": genoflu_status,
        "consensus_segments": consensus_segments,
        "review_flags": ";".join(review_flags) or "NONE",
        "ha_top_blast_hit": blast_hits.get("HA", "NO_HIT"),
        "na_top_blast_hit": blast_hits.get("NA", "NO_HIT"),
    }

    with output_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerow(row)


def main() -> None:
    sample = snakemake.wildcards.sample

    fastplong = read_fastplong_json(
        Path(snakemake.input.fastplong)
    )
    coverage_rows = read_coverage(
        Path(snakemake.input.coverage)
    )
    blast_hits = read_blast(
        Path(snakemake.input.blast)
    )
    h5n1_status = read_h5n1_flag(
        Path(snakemake.input.h5n1)
    )
    genoflu_status = read_genoflu(
        Path(snakemake.input.genoflu)
    )
    consensus_segments = count_fasta_records(
        Path(snakemake.input.consensus)
    )

    write_summary(
        output_path=Path(snakemake.output.tsv),
        sample=sample,
        fastplong=fastplong,
        coverage_rows=coverage_rows,
        blast_hits=blast_hits,
        h5n1_status=h5n1_status,
        genoflu_status=genoflu_status,
        consensus_segments=consensus_segments,
    )


if __name__ == "__main__":
    main()