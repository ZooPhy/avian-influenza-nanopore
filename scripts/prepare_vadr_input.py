#!/usr/bin/env python3
"""Concatenate coverage-qualified Medaka consensus FASTAs for VADR."""

from pathlib import Path


def is_pass(path: Path) -> bool:
    try:
        return path.read_text(encoding="utf-8", errors="replace").strip() == "PASS"
    except OSError:
        return False


output_path = Path(snakemake.output.fasta)
output_path.parent.mkdir(parents=True, exist_ok=True)
records: list[str] = []

for consensus_name, flag_name in zip(snakemake.input.consensus, snakemake.input.flags):
    consensus = Path(consensus_name)
    flag = Path(flag_name)
    if not is_pass(flag) or not consensus.is_file() or consensus.stat().st_size == 0:
        continue
    text = consensus.read_text(encoding="utf-8", errors="replace").strip()
    if text:
        records.append(text)

output_path.write_text(("\n".join(records) + "\n") if records else "", encoding="utf-8")
