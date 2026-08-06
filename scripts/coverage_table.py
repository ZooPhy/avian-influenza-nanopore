#!/usr/bin/env python3
"""Combine normalized manifest and per-segment coverage statistics."""


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
        flag = read_flag(flag_paths.get(segment, Path("/nonexistent")))

        positions = stats_row.get("positions") or manifest_row.get("fasta_length") or "0"
        rows.append(
            {
                "sample": sample,
                "segment": segment,
                "contig": stats_row.get("selected_contig")
                or stats_row.get("chosen_table")
                or manifest_row.get("contig")
                or "NA",
                "length": positions,
                "median_depth": stats_row.get("median_depth") or "NA",
                "mean_depth": stats_row.get("mean_depth") or "NA",
                "breadth_covered": stats_row.get("breadth_covered") or "NA",
                "coverage_flag": flag,
                "assembly_status": manifest_row.get("status") or "MISSING",
                "candidate_count": manifest_row.get("candidate_count") or "0",
                "selection_reason": manifest_row.get("selection_reason") or "",
            }
        )

    fieldnames = [
        "sample",
        "segment",
        "contig",
        "length",
        "median_depth",
        "mean_depth",
        "breadth_covered",
        "coverage_flag",
        "assembly_status",
        "candidate_count",
        "selection_reason",
    ]
    with output_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
