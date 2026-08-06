#!/usr/bin/env python3
"""Summarize per-segment BLAST output while retaining all eight segments."""


import csv
from pathlib import Path

SEGMENTS = ("HA", "NA", "PB2", "PB1", "PA", "NP", "MP", "NS")


def top_hit(path: Path) -> str:
    if not path.is_file() or path.stat().st_size == 0:
        return "NO_HIT"
    with path.open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 2:
                return fields[1]
    return "NO_HIT"


sample = str(snakemake.wildcards.sample)
by_segment = {
    Path(path).name.split(".blast.txt")[0]: Path(path)
    for path in snakemake.input.blast_files
}
output_path = Path(snakemake.output.csv)
output_path.parent.mkdir(parents=True, exist_ok=True)

with output_path.open("w", newline="", encoding="utf-8") as handle:
    writer = csv.DictWriter(handle, fieldnames=["sample", "segment", "top_hit"])
    writer.writeheader()
    for segment in SEGMENTS:
        writer.writerow(
            {
                "sample": sample,
                "segment": segment,
                "top_hit": top_hit(by_segment[segment]) if segment in by_segment else "NO_HIT",
            }
        )
