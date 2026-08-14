#!/usr/bin/env python3
"""Summarize per-segment BLAST evidence while retaining all eight segments."""

import csv
from pathlib import Path

SEGMENTS = ("HA", "NA", "PB2", "PB1", "PA", "NP", "MP", "NS")
BLAST_FIELDS = (
    "query_id",
    "subject_accession",
    "subject_title",
    "percent_identity",
    "alignment_length",
    "query_length",
    "query_start",
    "query_end",
    "subject_start",
    "subject_end",
    "evalue",
    "bit_score",
)
OUTPUT_FIELDS = (
    "sample",
    "segment",
    "top_hit",
    "query_id",
    "subject_accession",
    "subject_title",
    "percent_identity",
    "alignment_length",
    "query_length",
    "query_coverage",
    "evalue",
    "bit_score",
    "hit_status",
)


def segment_from_path(path: Path, suffix: str) -> str:
    name = path.name
    if name.endswith(suffix):
        return name[: -len(suffix)]
    return path.stem


def read_qc_status(path: Path) -> str:
    if not path.is_file():
        return "FAIL"
    value = path.read_text(encoding="utf-8", errors="replace").strip().upper()
    return "PASS" if value == "PASS" else "FAIL"


def read_top_hit(path: Path):
    if not path.is_file() or path.stat().st_size == 0:
        return None

    with path.open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < len(BLAST_FIELDS):
                raise ValueError(
                    f"Malformed BLAST row in {path}: expected at least "
                    f"{len(BLAST_FIELDS)} columns, found {len(fields)}"
                )
            return dict(zip(BLAST_FIELDS, fields[: len(BLAST_FIELDS)]))
    return None


def blank_row(sample: str, segment: str, status: str):
    row = {field: "" for field in OUTPUT_FIELDS}
    row.update(
        {
            "sample": sample,
            "segment": segment,
            "top_hit": "NO_HIT" if status == "NO_HIT" else status,
            "hit_status": status,
        }
    )
    return row


def summarize_hit(sample: str, segment: str, hit: dict, min_identity: float, min_qcov: float):
    try:
        percent_identity = float(hit["percent_identity"])
        alignment_length = int(hit["alignment_length"])
        query_length = int(hit["query_length"])
        bit_score = float(hit["bit_score"])
        evalue = float(hit["evalue"])
    except (TypeError, ValueError) as exc:
        raise ValueError(f"Invalid numeric BLAST value for {sample} {segment}: {hit}") from exc

    if query_length <= 0:
        raise ValueError(f"BLAST query length must be positive for {sample} {segment}")

    query_coverage = 100.0 * alignment_length / query_length
    status = (
        "HIGH_CONFIDENCE"
        if percent_identity >= min_identity and query_coverage >= min_qcov
        else "LOW_CONFIDENCE"
    )

    return {
        "sample": sample,
        "segment": segment,
        "top_hit": hit["subject_accession"],
        "query_id": hit["query_id"],
        "subject_accession": hit["subject_accession"],
        "subject_title": hit["subject_title"],
        "percent_identity": f"{percent_identity:.3f}",
        "alignment_length": alignment_length,
        "query_length": query_length,
        "query_coverage": f"{query_coverage:.3f}",
        "evalue": f"{evalue:.6g}",
        "bit_score": f"{bit_score:.3f}",
        "hit_status": status,
    }


sample = str(snakemake.wildcards.sample)
min_identity = float(snakemake.params.min_identity)
min_query_coverage = float(snakemake.params.min_query_coverage)

by_segment = {
    segment_from_path(Path(path), ".blast.txt"): Path(path)
    for path in snakemake.input.blast_files
}
flags_by_segment = {
    segment_from_path(Path(path), ".flag"): Path(path)
    for path in snakemake.input.flags
}

output_path = Path(snakemake.output.csv)
output_path.parent.mkdir(parents=True, exist_ok=True)

with output_path.open("w", newline="", encoding="utf-8") as handle:
    writer = csv.DictWriter(handle, fieldnames=OUTPUT_FIELDS)
    writer.writeheader()

    for segment in SEGMENTS:
        qc_status = read_qc_status(flags_by_segment.get(segment, Path("__missing__")))
        if qc_status != "PASS":
            writer.writerow(blank_row(sample, segment, "SKIPPED_QC"))
            continue

        hit = read_top_hit(by_segment[segment]) if segment in by_segment else None
        if hit is None:
            writer.writerow(blank_row(sample, segment, "NO_HIT"))
            continue

        writer.writerow(
            summarize_hit(
                sample,
                segment,
                hit,
                min_identity,
                min_query_coverage,
            )
        )
