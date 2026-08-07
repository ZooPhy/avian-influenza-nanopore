#!/usr/bin/env python3
"""Build a portable WINGS report bundle from rendered HTML reports.

The resulting .wings file is JSON and contains only the rendered run-summary
and sample-summary HTML files plus a small manifest. It is designed to be read
locally in a web browser without uploading report contents to a server.
"""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path

FORMAT_NAME = "WINGS_REPORT_BUNDLE"
FORMAT_VERSION = 1
SAMPLE_SUFFIX = ".sample_summary.html"


def read_text(path: Path) -> str:
    if not path.is_file():
        raise FileNotFoundError(f"Report file not found: {path}")
    text = path.read_text(encoding="utf-8")
    if not text.strip():
        raise ValueError(f"Report file is empty: {path}")
    return text


def sample_id_from_path(path: Path) -> str:
    name = path.name
    if not name.endswith(SAMPLE_SUFFIX):
        raise ValueError(
            f"Sample report filename must end with '{SAMPLE_SUFFIX}': {path}"
        )
    sample_id = name[: -len(SAMPLE_SUFFIX)]
    if not sample_id:
        raise ValueError(f"Could not determine sample ID from: {path}")
    return sample_id


def build_bundle(run_summary: Path, sample_reports: list[Path]) -> dict:
    samples = {}
    for report_path in sample_reports:
        sample_id = sample_id_from_path(report_path)
        if sample_id in samples:
            raise ValueError(f"Duplicate sample report for sample '{sample_id}'")
        samples[sample_id] = {
            "filename": report_path.name,
            "html": read_text(report_path),
        }

    if not samples:
        raise ValueError("At least one sample report is required")

    return {
        "format": FORMAT_NAME,
        "version": FORMAT_VERSION,
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "sample_count": len(samples),
        "run_summary": {
            "filename": run_summary.name,
            "html": read_text(run_summary),
        },
        "samples": samples,
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build a portable WINGS report bundle from rendered HTML reports."
    )
    parser.add_argument(
        "--run-summary",
        required=True,
        type=Path,
        help="Rendered WINGS run_summary.html file",
    )
    parser.add_argument(
        "--output",
        required=True,
        type=Path,
        help="Output .wings bundle path",
    )
    parser.add_argument(
        "sample_reports",
        nargs="+",
        type=Path,
        help="Rendered *.sample_summary.html files",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    bundle = build_bundle(args.run_summary, args.sample_reports)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    temp_output = args.output.with_suffix(args.output.suffix + ".tmp")
    with temp_output.open("w", encoding="utf-8") as handle:
        json.dump(bundle, handle, ensure_ascii=False, separators=(",", ":"))
        handle.write("\n")
    temp_output.replace(args.output)

    size_mb = args.output.stat().st_size / (1024 * 1024)
    print(
        f"WINGS report bundle created: {args.output} "
        f"({bundle['sample_count']} samples, {size_mb:.2f} MB)"
    )


if __name__ == "__main__":
    main()
