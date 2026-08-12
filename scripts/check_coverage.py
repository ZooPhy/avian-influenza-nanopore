#!/usr/bin/env python3
"""Compute segment coverage from a normalized IRMA BAM using samtools depth."""

import csv
import statistics
import subprocess
from pathlib import Path


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
    min_breadth: float,
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
        "minimum_breadth",
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
        "minimum_breadth": f"{min_breadth:.3f}",
        "status": status,
        "n_candidates": candidate_count,
        "selection_reason": reason,
    }
    with stats_out.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerow(row)


def samtools_depths(bam_path: Path) -> tuple[list[int], str]:
    """Return per-reference-position depths from ``samtools depth -aa``.

    ``-aa`` is required so zero-depth positions are included in the denominator
    used for coverage breadth. Base-quality and mapping-quality thresholds are
    set explicitly to zero to make the intended counting behavior reproducible.
    """
    command = [
        "samtools",
        "depth",
        "-aa",
        "-q",
        "0",
        "-Q",
        "0",
        str(bam_path),
    ]

    depths: list[int] = []
    contigs: list[str] = []
    seen_contigs: set[str] = set()

    try:
        process = subprocess.Popen(
            command,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            encoding="utf-8",
            errors="replace",
        )
    except OSError as exc:
        raise RuntimeError(f"could not start samtools depth: {exc}") from exc

    assert process.stdout is not None
    for line_number, line in enumerate(process.stdout, start=1):
        line = line.rstrip("\n")
        if not line:
            continue
        fields = line.split("\t")
        if len(fields) < 3:
            process.kill()
            raise RuntimeError(
                f"unexpected samtools depth output at line {line_number}: {line!r}"
            )
        reference = fields[0]
        try:
            depth = int(fields[2])
        except ValueError as exc:
            process.kill()
            raise RuntimeError(
                f"non-integer depth from samtools at line {line_number}: {fields[2]!r}"
            ) from exc

        if reference not in seen_contigs:
            seen_contigs.add(reference)
            contigs.append(reference)
        depths.append(depth)

    stderr = process.stderr.read() if process.stderr is not None else ""
    return_code = process.wait()
    if return_code != 0:
        message = stderr.strip() or f"samtools depth exited with status {return_code}"
        raise RuntimeError(message)

    contig = "NA" if not contigs else (contigs[0] if len(contigs) == 1 else ",".join(contigs))
    return depths, contig


def main() -> None:
    segments_dir = Path(snakemake.input.segments_dir)
    manifest_path = Path(snakemake.input.manifest)
    segment = str(snakemake.wildcards.segment)
    flag_out = Path(snakemake.output.flag)
    stats_out = Path(snakemake.output.stats)
    threshold = float(snakemake.params.min_median_depth)
    min_breadth = float(snakemake.params.min_breadth)

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
            min_breadth=min_breadth,
            status="MISSING",
            contig=contig,
            candidate_count=candidate_count,
            reason=f"normalized segment status={manifest_status}; usable BAM unavailable",
        )
        print(f"[check_coverage] {segment}: no usable normalized BAM -> MISSING")
        return

    try:
        depths, bam_contig = samtools_depths(bam_path)
    except Exception as exc:
        write_outputs(
            flag_out,
            stats_out,
            segment=segment,
            threshold=threshold,
            min_breadth=min_breadth,
            status="ERROR",
            contig=contig,
            candidate_count=candidate_count,
            reason=f"samtools depth failed for normalized BAM: {exc}",
        )
        print(f"[check_coverage] {segment}: samtools depth error -> ERROR: {exc}")
        return

    if not depths:
        write_outputs(
            flag_out,
            stats_out,
            segment=segment,
            threshold=threshold,
            min_breadth=min_breadth,
            status="MISSING",
            contig=contig,
            candidate_count=candidate_count,
            reason="samtools depth returned no reference positions",
        )
        print(f"[check_coverage] {segment}: no reference positions -> MISSING")
        return

    median = statistics.median(depths)
    mean = sum(depths) / len(depths)
    minimum = min(depths)
    maximum = max(depths)
    breadth = sum(depth >= threshold for depth in depths) / len(depths)
    status = "PASS" if median >= threshold and breadth >= min_breadth else "FAIL"
    selected_contig = contig if contig != "NA" else bam_contig

    write_outputs(
        flag_out,
        stats_out,
        segment=segment,
        threshold=threshold,
        min_breadth=min_breadth,
        status=status,
        contig=selected_contig,
        positions=len(depths),
        median=f"{median:.2f}",
        mean=f"{mean:.2f}",
        minimum=f"{minimum:.2f}",
        maximum=f"{maximum:.2f}",
        breadth=f"{breadth:.3f}",
        candidate_count=candidate_count,
        reason=(
            "coverage computed from normalized BAM with samtools depth -aa -q 0 -Q 0; "
            f"PASS requires median depth >= {threshold:.2f}x and "
            f"breadth at >= {threshold:.2f}x >= {min_breadth:.3f}"
        ),
    )
    print(
        f"[check_coverage] {segment}: contig={selected_contig}; positions={len(depths)}; "
        f"median={median:.2f}; mean={mean:.2f}; "
        f"breadth_at_{threshold:g}x={breadth:.3f} (min={min_breadth:.3f}) -> {status}"
    )


if __name__ == "__main__":
    main()
