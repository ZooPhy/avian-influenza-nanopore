#!/usr/bin/env python3
"""Normalize IRMA outputs into a stable per-segment interface.

The script scans an IRMA project recursively, selects at most one contig for each
influenza A segment family, copies the selected FASTA/BAM/coverage artifacts into
fixed paths, and writes a manifest that downstream rules can consume without
knowing IRMA's internal directory layout.
"""


import argparse
import csv
import re
import shutil
import statistics
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable

try:
    import pysam
except ImportError:  # The normalizer still works; BAM indexing is then skipped.
    pysam = None

SEGMENTS = ("HA", "NA", "PB2", "PB1", "PA", "NP", "MP", "NS")
FASTA_SUFFIXES = {".fasta", ".fa", ".fas", ".fna"}


@dataclass
class Candidate:
    contig: str
    segment: str
    fasta: Path | None = None
    bam: Path | None = None
    bam_index: Path | None = None
    coverage: Path | None = None
    fasta_length: int = 0
    coverage_median: float | None = None
    alternates: list[str] = field(default_factory=list)

    @property
    def ready(self) -> bool:
        return self.fasta is not None and self.bam is not None


def segment_from_contig(contig: str) -> str | None:
    """Return the canonical segment family for an IRMA contig identifier."""
    if not contig.startswith("A_"):
        return None
    remainder = contig[2:]
    if remainder.startswith("HA"):
        return "HA"
    if remainder.startswith("NA"):
        return "NA"
    for segment in ("PB2", "PB1", "PA", "NP", "MP", "NS"):
        if remainder == segment or remainder.startswith(f"{segment}_"):
            return segment
    return None


def fasta_length(path: Path) -> int:
    length = 0
    try:
        with path.open(encoding="utf-8", errors="replace") as handle:
            for line in handle:
                if not line.startswith(">"):
                    length += len(line.strip())
    except OSError:
        return 0
    return length


def norm_col(name: str) -> str:
    return re.sub(r"[^a-z0-9]+", "_", name.lower()).strip("_")


def coverage_median(path: Path) -> float | None:
    """Read median depth from an IRMA coverage table when possible."""
    try:
        with path.open(encoding="utf-8", errors="replace") as handle:
            reader = csv.reader(handle, delimiter="\t")
            header = next(row for row in reader if row)
            normalized = [norm_col(value) for value in header]
            depth_index = None
            for key in ("coverage_depth", "depth", "read_depth", "readdepth"):
                if key in normalized:
                    depth_index = normalized.index(key)
                    break
            if depth_index is None:
                return None

            depths: list[float] = []
            for row in reader:
                if not row or depth_index >= len(row):
                    continue
                value = row[depth_index].strip()
                if not value or value.upper() == "NA":
                    continue
                try:
                    depths.append(float(value))
                except ValueError:
                    continue
    except (OSError, StopIteration):
        return None
    return statistics.median(depths) if depths else None


def relative_or_absolute(path: Path | None, base: Path) -> str:
    if path is None:
        return ""
    try:
        return str(path.resolve().relative_to(base.resolve()))
    except ValueError:
        return str(path.resolve())


def choose_existing(paths: Iterable[Path]) -> Path | None:
    existing = sorted((path for path in paths if path.is_file()), key=lambda p: p.as_posix())
    return existing[0] if existing else None


def collect_candidates(project: Path) -> dict[str, Candidate]:
    candidates: dict[str, Candidate] = {}

    for path in sorted(project.rglob("*"), key=lambda p: p.as_posix()):
        if not path.is_file():
            continue

        name = path.name
        lower_name = name.lower()
        contig: str | None = None
        kind: str | None = None

        if path.suffix.lower() in FASTA_SUFFIXES:
            contig = path.stem
            kind = "fasta"
        elif path.suffix.lower() == ".bam":
            contig = path.stem
            kind = "bam"
        elif lower_name.endswith("-coverage.txt"):
            contig = name[: -len("-coverage.txt")]
            kind = "coverage"
        else:
            continue

        segment = segment_from_contig(contig)
        if segment is None:
            continue

        candidate = candidates.setdefault(contig, Candidate(contig=contig, segment=segment))
        current = getattr(candidate, kind)
        if current is None:
            setattr(candidate, kind, path)
        else:
            candidate.alternates.append(str(path))

    for candidate in candidates.values():
        if candidate.fasta is not None:
            candidate.fasta_length = fasta_length(candidate.fasta)
        if candidate.coverage is not None:
            candidate.coverage_median = coverage_median(candidate.coverage)
        if candidate.bam is not None:
            candidate.bam_index = choose_existing(
                [
                    Path(f"{candidate.bam}.bai"),
                    candidate.bam.with_suffix(".bai"),
                ]
            )

    return candidates


def candidate_rank(candidate: Candidate) -> tuple[float, ...]:
    """Prefer complete pairs, then measured depth, BAM size, and FASTA length."""
    bam_size = candidate.bam.stat().st_size if candidate.bam and candidate.bam.exists() else 0
    measured_depth = candidate.coverage_median if candidate.coverage_median is not None else -1.0
    return (
        1.0 if candidate.ready else 0.0,
        1.0 if candidate.coverage_median is not None else 0.0,
        measured_depth,
        float(bam_size),
        float(candidate.fasta_length),
    )


def selection_reason(candidate: Candidate, count: int) -> str:
    parts = [f"selected from {count} candidate(s)"]
    if candidate.ready:
        parts.append("FASTA+BAM pair available")
    elif candidate.fasta is not None:
        parts.append("FASTA only")
    elif candidate.bam is not None:
        parts.append("BAM only")
    if candidate.coverage_median is not None:
        parts.append(f"IRMA median depth={candidate.coverage_median:g}")
    else:
        parts.append("no parseable IRMA coverage table")
    return "; ".join(parts)


def copy_candidate(candidate: Candidate, segment_dir: Path) -> dict[str, str]:
    segment_dir.mkdir(parents=True, exist_ok=True)
    copied = {"fasta": "", "bam": "", "bam_index": "", "coverage": ""}

    if candidate.fasta is not None:
        destination = segment_dir / "consensus.fasta"
        shutil.copy2(candidate.fasta, destination)
        copied["fasta"] = str(destination)

    if candidate.bam is not None:
        destination = segment_dir / "alignment.bam"
        shutil.copy2(candidate.bam, destination)
        copied["bam"] = str(destination)
        index_destination = segment_dir / "alignment.bam.bai"
        if candidate.bam_index is not None:
            shutil.copy2(candidate.bam_index, index_destination)
            copied["bam_index"] = str(index_destination)
        elif pysam is not None:
            try:
                pysam.index(str(destination))
                if index_destination.is_file():
                    copied["bam_index"] = str(index_destination)
            except Exception:
                # Coverage can still be streamed without an index, and Medaka's
                # fail-soft behavior will handle BAMs it cannot consume.
                pass

    if candidate.coverage is not None:
        destination = segment_dir / "irma_coverage.tsv"
        shutil.copy2(candidate.coverage, destination)
        copied["coverage"] = str(destination)

    return copied


def normalize(project: Path, segments_dir: Path, manifest: Path, sample: str) -> None:
    if segments_dir.exists():
        shutil.rmtree(segments_dir)
    segments_dir.mkdir(parents=True, exist_ok=True)
    manifest.parent.mkdir(parents=True, exist_ok=True)

    candidates = collect_candidates(project)
    by_segment: dict[str, list[Candidate]] = {segment: [] for segment in SEGMENTS}
    for candidate in candidates.values():
        by_segment[candidate.segment].append(candidate)

    rows: list[dict[str, str | int | float]] = []
    log_lines = [f"sample={sample}", f"project={project}", f"candidate_contigs={len(candidates)}"]

    for segment in SEGMENTS:
        segment_candidates = sorted(by_segment[segment], key=lambda c: c.contig)
        if not segment_candidates:
            rows.append(
                {
                    "sample": sample,
                    "segment": segment,
                    "status": "MISSING",
                    "contig": "",
                    "fasta": "",
                    "bam": "",
                    "bam_index": "",
                    "coverage_table": "",
                    "fasta_length": 0,
                    "coverage_median": "",
                    "candidate_count": 0,
                    "selection_reason": "no segment-specific IRMA FASTA or BAM found",
                    "source_fasta": "",
                    "source_bam": "",
                    "source_coverage": "",
                }
            )
            log_lines.append(f"{segment}: MISSING")
            continue

        selected = sorted(
            segment_candidates,
            key=lambda candidate: (
                -candidate_rank(candidate)[0],
                -candidate_rank(candidate)[1],
                -candidate_rank(candidate)[2],
                -candidate_rank(candidate)[3],
                -candidate_rank(candidate)[4],
                candidate.contig,
            ),
        )[0]
        copied = copy_candidate(selected, segments_dir / segment)

        if selected.ready:
            status = "READY"
        elif selected.fasta is not None:
            status = "FASTA_ONLY"
        elif selected.bam is not None:
            status = "BAM_ONLY"
        else:
            status = "MISSING"

        rows.append(
            {
                "sample": sample,
                "segment": segment,
                "status": status,
                "contig": selected.contig,
                "fasta": relative_or_absolute(Path(copied["fasta"]) if copied["fasta"] else None, manifest.parent),
                "bam": relative_or_absolute(Path(copied["bam"]) if copied["bam"] else None, manifest.parent),
                "bam_index": relative_or_absolute(Path(copied["bam_index"]) if copied["bam_index"] else None, manifest.parent),
                "coverage_table": relative_or_absolute(Path(copied["coverage"]) if copied["coverage"] else None, manifest.parent),
                "fasta_length": selected.fasta_length,
                "coverage_median": "" if selected.coverage_median is None else f"{selected.coverage_median:.6g}",
                "candidate_count": len(segment_candidates),
                "selection_reason": selection_reason(selected, len(segment_candidates)),
                "source_fasta": relative_or_absolute(selected.fasta, project),
                "source_bam": relative_or_absolute(selected.bam, project),
                "source_coverage": relative_or_absolute(selected.coverage, project),
            }
        )
        log_lines.append(
            f"{segment}: {status}; contig={selected.contig}; "
            f"candidates={len(segment_candidates)}"
        )

    fieldnames = [
        "sample",
        "segment",
        "status",
        "contig",
        "fasta",
        "bam",
        "bam_index",
        "coverage_table",
        "fasta_length",
        "coverage_median",
        "candidate_count",
        "selection_reason",
        "source_fasta",
        "source_bam",
        "source_coverage",
    ]
    with manifest.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    (segments_dir / "normalization.log").write_text("\n".join(log_lines) + "\n", encoding="utf-8")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project", required=True, type=Path)
    parser.add_argument("--segments", required=True, type=Path)
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--sample", required=True)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    normalize(args.project.resolve(), args.segments.resolve(), args.manifest.resolve(), args.sample)


if __name__ == "__main__":
    main()
