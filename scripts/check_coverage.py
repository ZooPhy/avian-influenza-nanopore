#!/usr/bin/env python3
"""Compute segment coverage directly from a normalized IRMA BAM."""


import csv
import statistics
from pathlib import Path

import pysam


def read_manifest_row(manifest_path: Path, segment: str) -> dict[str, str]:
    if not manifest_path.is_file():
        return {}
    with manifest_path.open(encoding="utf-8", errors="replace") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row.get("segment") == segment:
                return row
    return {}


def write_outputs(
    flag_out: Path,
    stats_out: Path,
    *,
    segment: str,
    threshold: float,
    status: str,
    contig: str = "NA",
    positions: int = 0,
    median: str = "NA",
    mean: str = "NA",
    minimum: str = "NA",
    maximum: str = "NA",
    breadth: str = "NA",
    candidate_count: int = 0,
    reason: str = "",
) -> None:
    flag_out.parent.mkdir(parents=True, exist_ok=True)
    stats_out.parent.mkdir(parents=True, exist_ok=True)
    flag_out.write_text(f"{status}\n", encoding="utf-8")

    fieldnames = [
        "segment",
        "chosen_table",
        "selected_contig",
        "positions",
        "median_depth",
        "mean_depth",
        "min_depth",
        "max_depth",
        "breadth_covered",
        "threshold",
        "status",
        "n_candidates",
        "selection_reason",
    ]
    row = {
        "segment": segment,
        # Retained for compatibility with the H5N1 screen and existing reports.
        "chosen_table": contig,
        "selected_contig": contig,
        "positions": positions,
        "median_depth": median,
        "mean_depth": mean,
        "min_depth": minimum,
        "max_depth": maximum,
        "breadth_covered": breadth,
        "threshold": f"{threshold:.2f}",
        "status": status,
        "n_candidates": candidate_count,
        "selection_reason": reason,
    }
    with stats_out.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerow(row)


def depths_without_index(bam: pysam.AlignmentFile) -> tuple[list[int], str]:
    """Compute per-position depth by streaming a BAM without requiring an index."""
    references = list(bam.references)
    lengths = list(bam.lengths)
    depth_by_reference = {
        reference: [0] * length for reference, length in zip(references, lengths)
    }

    for read in bam.fetch(until_eof=True):
        if read.is_unmapped or read.reference_id < 0:
            continue
        reference = bam.get_reference_name(read.reference_id)
        values = depth_by_reference.get(reference)
        if values is None:
            continue
        for position in read.get_reference_positions(full_length=False):
            if position is not None and 0 <= position < len(values):
                values[position] += 1

    all_depths: list[int] = []
    for reference in references:
        all_depths.extend(depth_by_reference[reference])
    contig = references[0] if len(references) == 1 else ",".join(references)
    return all_depths, contig


def depths_with_index(bam: pysam.AlignmentFile) -> tuple[list[int], str]:
    all_depths: list[int] = []
    for reference, length in zip(bam.references, bam.lengths):
        a, c, g, t = bam.count_coverage(reference, start=0, stop=length, quality_threshold=0)
        all_depths.extend(
            int(a[index] + c[index] + g[index] + t[index])
            for index in range(length)
        )
    contig = bam.references[0] if len(bam.references) == 1 else ",".join(bam.references)
    return all_depths, contig


def main() -> None:
    segments_dir = Path(snakemake.input.segments_dir)
    manifest_path = Path(snakemake.input.manifest)
    segment = str(snakemake.wildcards.segment)
    flag_out = Path(snakemake.output.flag)
    stats_out = Path(snakemake.output.stats)
    threshold = float(snakemake.params.min_median_depth)

    manifest_row = read_manifest_row(manifest_path, segment)
    manifest_status = manifest_row.get("status", "MISSING")
    contig = manifest_row.get("contig") or "NA"
    candidate_count = int(manifest_row.get("candidate_count") or 0)
    bam_path = segments_dir / segment / "alignment.bam"

    if manifest_status != "READY" or not bam_path.is_file() or bam_path.stat().st_size == 0:
        write_outputs(
            flag_out,
            stats_out,
            segment=segment,
            threshold=threshold,
            status="MISSING",
            contig=contig,
            candidate_count=candidate_count,
            reason=f"normalized segment status={manifest_status}; usable BAM unavailable",
        )
        print(f"[check_coverage] {segment}: no usable normalized BAM -> MISSING")
        return

    try:
        with pysam.AlignmentFile(str(bam_path), "rb") as bam:
            if not bam.references:
                raise ValueError("BAM has no reference sequences")
            try:
                has_index = bam.has_index()
            except (ValueError, OSError):
                has_index = False
            if has_index:
                depths, bam_contig = depths_with_index(bam)
                method = "pysam count_coverage using BAM index"
            else:
                depths, bam_contig = depths_without_index(bam)
                method = "streamed BAM alignments without index"
    except Exception as exc:
        write_outputs(
            flag_out,
            stats_out,
            segment=segment,
            threshold=threshold,
            status="ERROR",
            contig=contig,
            candidate_count=candidate_count,
            reason=f"could not read normalized BAM: {exc}",
        )
        print(f"[check_coverage] {segment}: BAM error -> ERROR: {exc}")
        return

    if not depths:
        write_outputs(
            flag_out,
            stats_out,
            segment=segment,
            threshold=threshold,
            status="MISSING",
            contig=contig,
            candidate_count=candidate_count,
            reason="BAM contained no reference positions",
        )
        print(f"[check_coverage] {segment}: no reference positions -> MISSING")
        return

    median = statistics.median(depths)
    mean = sum(depths) / len(depths)
    minimum = min(depths)
    maximum = max(depths)
    breadth = sum(depth > 0 for depth in depths) / len(depths)
    status = "PASS" if median >= threshold else "FAIL"
    selected_contig = contig if contig != "NA" else bam_contig

    write_outputs(
        flag_out,
        stats_out,
        segment=segment,
        threshold=threshold,
        status=status,
        contig=selected_contig,
        positions=len(depths),
        median=f"{median:.2f}",
        mean=f"{mean:.2f}",
        minimum=f"{minimum:.2f}",
        maximum=f"{maximum:.2f}",
        breadth=f"{breadth:.3f}",
        candidate_count=candidate_count,
        reason=f"coverage computed from normalized BAM; {method}",
    )
    print(
        f"[check_coverage] {segment}: contig={selected_contig}; positions={len(depths)}; "
        f"median={median:.2f}; mean={mean:.2f}; breadth={breadth:.3f} -> {status}"
    )


if __name__ == "__main__":
    main()
