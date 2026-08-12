#!/usr/bin/env python3
"""Combine normalized manifest and per-segment QC statistics."""

import csv
from pathlib import Path

SEGMENTS = ("HA", "NA", "PB2", "PB1", "PA", "NP", "MP", "NS")


def read_manifest(path: Path) -> dict[str, dict[str, str]]:
    if not path.is_file():
        return {}
    with path.open(encoding="utf-8", errors="replace") as handle:
        return {
            row["segment"]: row
            for row in csv.DictReader(handle, delimiter="\t")
            if row.get("segment")
        }


def read_first_row(path: Path) -> dict[str, str]:
    if not path.is_file():
        return {}
    with path.open(encoding="utf-8", errors="replace") as handle:
        return next(csv.DictReader(handle, delimiter="\t"), {})


def read_flag(path: Path) -> str:
    try:
        return path.read_text(encoding="utf-8", errors="replace").strip() or "MISSING"
    except OSError:
        return "MISSING"


def main() -> None:
    sample = str(snakemake.wildcards.sample)
    manifest = read_manifest(Path(snakemake.input.manifest))
    stats_paths = {Path(path).stem: Path(path) for path in snakemake.input.stats}
    flag_paths = {Path(path).stem: Path(path) for path in snakemake.input.flags}
    output_path = Path(snakemake.output.tsv)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    rows = []
    for segment in SEGMENTS:
        manifest_row = manifest.get(segment, {})
        stats_row = read_first_row(stats_paths.get(segment, Path("/nonexistent")))
        overall_flag = read_flag(flag_paths.get(segment, Path("/nonexistent")))

        consensus_length = (
            stats_row.get("consensus_length")
            or manifest_row.get("fasta_length")
            or "0"
        )
        coverage_status = stats_row.get("coverage_status") or overall_flag
        length_status = stats_row.get("length_status") or "NA"
        n_content_status = stats_row.get("n_content_status") or "NA"
        overall_status = stats_row.get("status") or overall_flag

        rows.append(
            {
                "sample": sample,
                "segment": segment,
                "contig": stats_row.get("selected_contig")
                or stats_row.get("chosen_table")
                or manifest_row.get("contig")
                or "NA",
                "length": consensus_length,
                "coverage_positions": stats_row.get("positions") or "0",
                "median_depth": stats_row.get("median_depth") or "NA",
                "mean_depth": stats_row.get("mean_depth") or "NA",
                "breadth_covered": stats_row.get("breadth_covered") or "NA",
                "n_count": stats_row.get("n_count") or "0",
                "n_fraction": stats_row.get("n_fraction") or "NA",
                "expected_length_min": stats_row.get("expected_length_min") or "NA",
                "expected_length_max": stats_row.get("expected_length_max") or "NA",
                "maximum_n_fraction": stats_row.get("maximum_n_fraction") or "NA",
                "coverage_status": coverage_status,
                "length_status": length_status,
                "n_content_status": n_content_status,
                "overall_status": overall_status,
                # Legacy column name retained because sample_summary.py and older
                # reports use it as the passing-segment flag. It now represents
                # the overall segment QC decision, not coverage alone.
                "coverage_flag": overall_status,
                "assembly_status": manifest_row.get("status") or "MISSING",
                "candidate_count": manifest_row.get("candidate_count") or "0",
                "selection_reason": manifest_row.get("selection_reason") or "",
                "qc_reason": stats_row.get("selection_reason") or "",
            }
        )

    fieldnames = [
        "sample",
        "segment",
        "contig",
        "length",
        "coverage_positions",
        "median_depth",
        "mean_depth",
        "breadth_covered",
        "n_count",
        "n_fraction",
        "expected_length_min",
        "expected_length_max",
        "maximum_n_fraction",
        "coverage_status",
        "length_status",
        "n_content_status",
        "overall_status",
        "coverage_flag",
        "assembly_status",
        "candidate_count",
        "selection_reason",
        "qc_reason",
    ]
    with output_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
