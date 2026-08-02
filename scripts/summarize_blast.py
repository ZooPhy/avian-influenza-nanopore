#!/usr/bin/env python3
"""
Summarise per-segment BLASTn results into a single CSV.

Expected Snakemake variables
----------------------------
* snakemake.input  : list of BLAST 6-column files (one per segment)
* snakemake.output : length-1 list; where to write the CSV
"""

import pandas as pd
from pathlib import Path

rows = []

for blast_file in snakemake.input:
    blast_path = Path(blast_file)
    # barcode20/blast/PB2.blast.txt  → sample=barcode20, segment=PB2
    sample  = blast_path.parents[1].name
    segment = blast_path.stem.split(".")[0]

    # BLAST -outfmt 6 columns (qseqid sseqid pident etc.)
    # Keep the first (best-scoring) hit; BLAST sorts by bit-score by default.
    try:
        top_hit = pd.read_csv(blast_file, sep="\t", header=None).iloc[0, 1]
    except (pd.errors.EmptyDataError, IndexError):
        top_hit = "NO_HIT"

    rows.append({"sample": sample, "segment": segment, "top_hit": top_hit})

# Combine and write
out_csv = snakemake.output[0]
Path(out_csv).parent.mkdir(parents=True, exist_ok=True)
pd.DataFrame(rows).to_csv(out_csv, index=False)

print(f"★ Wrote BLAST summary {out_csv}")
