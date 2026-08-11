#!/usr/bin/env python3
"""Extract one sample row from the validated WINGS metadata table."""

import csv
from pathlib import Path

source = Path(str(snakemake.input.metadata))
out = Path(str(snakemake.output.tsv))
sample = str(snakemake.wildcards.sample)

with source.open(encoding="utf-8", newline="") as handle:
    reader = csv.DictReader(handle, delimiter="\t")
    fieldnames = reader.fieldnames or []
    row = next((r for r in reader if r.get("sample_id") == sample), None)

if row is None:
    raise ValueError(f"Validated metadata has no row for sample {sample!r}")

out.parent.mkdir(parents=True, exist_ok=True)
with out.open("w", encoding="utf-8", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
    writer.writeheader()
    writer.writerow(row)
