#!/usr/bin/env python3
"""Validate WINGS sample metadata and emit a normalized table plus QC report."""

import csv
import re
from datetime import date
from pathlib import Path

REQUIRED_COLUMNS = ("sample_id",)
RECOMMENDED_COLUMNS = (
    "barcode",
    "sample_name",
    "host_species",
    "host_common_name",
    "sample_type",
    "collection_date",
    "collection_location",
    "state",
    "country",
    "latitude",
    "longitude",
    "flyway",
    "age",
    "sex",
    "wildlife_agency",
    "notes",
)


def clean(value):
    return "" if value is None else str(value).strip()


def parse_bool(value):
    return str(value).strip().lower() in {"1", "true", "yes", "y", "on"}


def valid_iso_date(value):
    if not value:
        return True
    try:
        date.fromisoformat(value)
        return True
    except ValueError:
        return False


def valid_coordinate(value, low, high):
    if not value:
        return True
    try:
        number = float(value)
    except ValueError:
        return False
    return low <= number <= high


def main():
    source = Path(str(snakemake.input.metadata))
    validated = Path(str(snakemake.output.validated))
    report = Path(str(snakemake.output.report))
    expected_samples = [s for s in str(snakemake.params.samples).split(",") if s]
    expected_set = set(expected_samples)
    require_all = parse_bool(snakemake.params.require_all_samples)

    if not source.is_file():
        raise ValueError(
            f"Metadata file not found: {source}. Create it from metadata.example.tsv "
            "or update metadata_file in config.yaml."
        )

    with source.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"Metadata file has no header: {source}")
        fieldnames = [clean(name) for name in reader.fieldnames]
        missing_required = [c for c in REQUIRED_COLUMNS if c not in fieldnames]
        if missing_required:
            raise ValueError(
                "Metadata file is missing required column(s): " + ", ".join(missing_required)
            )
        rows = []
        for raw in reader:
            row = {name: clean(raw.get(name, "")) for name in fieldnames}
            if not any(row.values()):
                continue
            rows.append(row)

    errors = []
    warnings = []
    seen = set()
    for i, row in enumerate(rows, start=2):
        sample = row.get("sample_id", "")
        if not sample:
            errors.append(f"line {i}: sample_id is empty")
            continue
        if sample in seen:
            errors.append(f"line {i}: duplicate sample_id {sample!r}")
        seen.add(sample)
        if sample not in expected_set:
            warnings.append(f"line {i}: metadata sample {sample!r} has no matching FASTQ")
        if not valid_iso_date(row.get("collection_date", "")):
            errors.append(
                f"line {i}: collection_date must be ISO YYYY-MM-DD: {row.get('collection_date')!r}"
            )
        if not valid_coordinate(row.get("latitude", ""), -90.0, 90.0):
            errors.append(f"line {i}: invalid latitude {row.get('latitude')!r}")
        if not valid_coordinate(row.get("longitude", ""), -180.0, 180.0):
            errors.append(f"line {i}: invalid longitude {row.get('longitude')!r}")
        barcode = row.get("barcode", "")
        if barcode and not re.fullmatch(r"barcode\d+", barcode, flags=re.IGNORECASE):
            warnings.append(f"line {i}: unusual barcode value {barcode!r}")

    missing_samples = [sample for sample in expected_samples if sample not in seen]
    if missing_samples:
        message = "FASTQ sample(s) missing from metadata: " + ", ".join(missing_samples)
        if require_all:
            errors.append(message)
        else:
            warnings.append(message)

    output_fields = list(fieldnames)
    for column in RECOMMENDED_COLUMNS:
        if column not in output_fields:
            output_fields.append(column)
    if "metadata_status" not in output_fields:
        output_fields.append("metadata_status")

    normalized_rows = []
    by_sample = {row.get("sample_id", ""): row for row in rows if row.get("sample_id")}
    for sample in expected_samples:
        source_row = by_sample.get(sample, {})
        row = {field: clean(source_row.get(field, "")) for field in output_fields}
        row["sample_id"] = sample
        row["metadata_status"] = "COMPLETE" if source_row else "MISSING"
        normalized_rows.append(row)

    validated.parent.mkdir(parents=True, exist_ok=True)
    with validated.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=output_fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(normalized_rows)

    report.parent.mkdir(parents=True, exist_ok=True)
    with report.open("w", encoding="utf-8") as handle:
        handle.write(f"source\t{source}\n")
        handle.write(f"expected_samples\t{len(expected_samples)}\n")
        handle.write(f"metadata_rows\t{len(rows)}\n")
        handle.write(f"matched_samples\t{len(expected_set & seen)}\n")
        handle.write(f"missing_samples\t{len(missing_samples)}\n")
        handle.write(f"warnings\t{len(warnings)}\n")
        handle.write(f"errors\t{len(errors)}\n")
        for warning in warnings:
            handle.write(f"WARNING\t{warning}\n")
        for error in errors:
            handle.write(f"ERROR\t{error}\n")

    if errors:
        raise ValueError(
            "Metadata validation failed:\n- " + "\n- ".join(errors) + f"\nSee {report}."
        )

    print(
        f"[metadata] validated {len(expected_samples)} FASTQ sample(s); "
        f"matched={len(expected_set & seen)}, missing={len(missing_samples)}, warnings={len(warnings)}"
    )


if __name__ == "__main__":
    main()
