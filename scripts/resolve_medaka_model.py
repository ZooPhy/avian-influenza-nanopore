#!/usr/bin/env python3

"""
Resolve Medaka model selectors from Oxford Nanopore FASTQ metadata.

The script reads FASTQ record headers and extracts:

    basecall_model_version_id=<model>

For example:

    basecall_model_version_id=dna_r10.4.1_e8.2_400bps_hac@v5.0.0

It writes a TSV containing the detected basecaller model and Medaka-compatible
model selectors:

    dna_r10.4.1_e8.2_400bps_hac@v5.0.0:consensus
    dna_r10.4.1_e8.2_400bps_hac@v5.0.0:variant

Medaka performs the actual mapping from the basecaller model to the
corresponding Medaka model. This avoids maintaining a duplicate model mapping
inside WINGS.
"""

from __future__ import annotations

import argparse
import gzip
import re
import sys
from collections import Counter
from pathlib import Path


MODEL_PATTERN = re.compile(
    r"(?:^|\s)basecall_model_version_id=([^\s]+)"
)


def open_text(path: Path):
    """Open plain-text or gzip-compressed FASTQ."""
    if path.suffix.lower() == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return path.open("rt", encoding="utf-8", errors="replace")


def detect_basecaller_model(
    fastq: Path,
    max_records: int = 100,
) -> tuple[str, int]:
    """
    Extract the basecaller model from FASTQ headers.

    Returns
    -------
    tuple[str, int]
        The unique basecaller model and number of records containing it.

    Raises
    ------
    RuntimeError
        If no model is found or multiple models are detected.
    """
    models: list[str] = []

    with open_text(fastq) as handle:
        record_number = 0

        while record_number < max_records:
            header = handle.readline()

            if not header:
                break

            sequence = handle.readline()
            plus = handle.readline()
            quality = handle.readline()

            if not sequence or not plus or not quality:
                raise RuntimeError(
                    f"FASTQ appears truncated near record "
                    f"{record_number + 1}: {fastq}"
                )

            record_number += 1

            if not header.startswith("@"):
                raise RuntimeError(
                    f"Invalid FASTQ header at record {record_number}: "
                    f"{header.rstrip()!r}"
                )

            match = MODEL_PATTERN.search(header)
            if match:
                models.append(match.group(1))

    if not models:
        raise RuntimeError(
            "No basecall_model_version_id metadata was found in the first "
            f"{max_records} FASTQ records of {fastq}"
        )

    counts = Counter(models)

    if len(counts) != 1:
        details = ", ".join(
            f"{model} ({count} reads)"
            for model, count in sorted(counts.items())
        )
        raise RuntimeError(
            "Multiple basecaller models were detected in the FASTQ: "
            f"{details}. WINGS will not choose a Medaka model automatically "
            "for mixed-basecaller input."
        )

    model, count = next(iter(counts.items()))
    return model, count


def write_tsv(
    output: Path,
    fastq: Path,
    basecaller_model: str,
    records_with_model: int,
) -> None:
    """Write resolved model metadata in a Snakemake-friendly TSV."""
    output.parent.mkdir(parents=True, exist_ok=True)

    consensus_selector = f"{basecaller_model}:consensus"
    variant_selector = f"{basecaller_model}:variant"

    with output.open("w", encoding="utf-8") as handle:
        handle.write(
            "\t".join(
                [
                    "fastq",
                    "basecaller_model",
                    "medaka_consensus_selector",
                    "medaka_variant_selector",
                    "records_with_model",
                ]
            )
            + "\n"
        )

        handle.write(
            "\t".join(
                [
                    str(fastq),
                    basecaller_model,
                    consensus_selector,
                    variant_selector,
                    str(records_with_model),
                ]
            )
            + "\n"
        )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Extract Oxford Nanopore basecaller metadata from FASTQ and "
            "generate Medaka model selectors."
        )
    )

    parser.add_argument(
        "--fastq",
        required=True,
        type=Path,
        help="Input FASTQ or FASTQ.GZ file.",
    )

    parser.add_argument(
        "--output",
        required=True,
        type=Path,
        help="Output TSV file.",
    )

    parser.add_argument(
        "--max-records",
        type=int,
        default=100,
        help="Maximum FASTQ records to inspect (default: 100).",
    )

    return parser.parse_args()


def main() -> int:
    args = parse_args()

    if args.max_records < 1:
        print(
            "ERROR: --max-records must be at least 1.",
            file=sys.stderr,
        )
        return 2

    if not args.fastq.is_file():
        print(
            f"ERROR: FASTQ does not exist: {args.fastq}",
            file=sys.stderr,
        )
        return 2

    try:
        basecaller_model, records_with_model = detect_basecaller_model(
            args.fastq,
            args.max_records,
        )

        write_tsv(
            args.output,
            args.fastq,
            basecaller_model,
            records_with_model,
        )

    except (OSError, RuntimeError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1

    print(
        f"Detected basecaller model: {basecaller_model}",
        file=sys.stderr,
    )
    print(
        f"Medaka consensus selector: "
        f"{basecaller_model}:consensus",
        file=sys.stderr,
    )
    print(
        f"Medaka variant selector: "
        f"{basecaller_model}:variant",
        file=sys.stderr,
    )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())