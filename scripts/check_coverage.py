#!/usr/bin/env python3
"""Compute segment QC from normalized IRMA BAM and consensus FASTA.

Coverage is calculated with ``samtools depth -aa -q 0 -Q 0``. A segment passes
the hard QC gate when median depth, breadth at the configured depth threshold,
minimum consensus length, and maximum N fraction pass. A consensus above the
configured upper length guide receives a warning but is not rejected.
"""

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


def read_consensus_sequence(fasta_path: Path) -> str:
    """Read exactly one non-empty FASTA record and return its sequence."""
    if not fasta_path.is_file() or fasta_path.stat().st_size == 0:
        raise ValueError("consensus FASTA is missing or empty")

    record_count = 0
    sequence_parts: list[str] = []
    with fasta_path.open(encoding="utf-8", errors="replace") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                record_count += 1
                if record_count > 1:
                    raise ValueError(
                        "consensus FASTA contains more than one record; "
                        f"second header encountered at line {line_number}"
                    )
                continue
            if record_count == 0:
                raise ValueError(
                    f"sequence data encountered before FASTA header at line {line_number}"
                )
            sequence_parts.append("".join(line.split()))

    if record_count != 1:
        raise ValueError(f"expected one FASTA record, found {record_count}")

    sequence = "".join(sequence_parts).upper()
    if not sequence:
        raise ValueError("consensus FASTA record has an empty sequence")
    return sequence


def write_outputs(
    flag_out: Path,
    stats_out: Path,
    *,
    segment: str,
    threshold: float,
    min_breadth: float,
    expected_length_min: int,
    expected_length_max: int,
    max_n_fraction: float,
    status: str,
    coverage_status: str = "NA",
    length_status: str = "NA",
    n_content_status: str = "NA",
    contig: str = "NA",
    positions: int = 0,
    consensus_length: int = 0,
    n_count: int = 0,
    n_fraction: str = "NA",
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
        "consensus_length",
        "n_count",
        "n_fraction",
        "median_depth",
        "mean_depth",
        "min_depth",
        "max_depth",
        "breadth_covered",
        "threshold",
        "minimum_breadth",
        "expected_length_min",
        "expected_length_max",
        "maximum_n_fraction",
        "coverage_status",
        "length_status",
        "n_content_status",
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
        "consensus_length": consensus_length,
        "n_count": n_count,
        "n_fraction": n_fraction,
        "median_depth": median,
        "mean_depth": mean,
        "min_depth": minimum,
        "max_depth": maximum,
        "breadth_covered": breadth,
        "threshold": f"{threshold:.2f}",
        "minimum_breadth": f"{min_breadth:.3f}",
        "expected_length_min": expected_length_min,
        "expected_length_max": expected_length_max,
        "maximum_n_fraction": f"{max_n_fraction:.4f}",
        "coverage_status": coverage_status,
        "length_status": length_status,
        "n_content_status": n_content_status,
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
    consensus_path = Path(snakemake.input.consensus)
    segment = str(snakemake.wildcards.segment)
    flag_out = Path(snakemake.output.flag)
    stats_out = Path(snakemake.output.stats)
    threshold = float(snakemake.params.min_median_depth)
    min_breadth = float(snakemake.params.min_breadth)
    expected_length_min = int(snakemake.params.expected_length_min)
    expected_length_max = int(snakemake.params.expected_length_max)
    max_n_fraction = float(snakemake.params.max_n_fraction)

    manifest_row = read_manifest_row(manifest_path, segment)
    manifest_status = manifest_row.get("status", "MISSING")
    contig = manifest_row.get("contig") or "NA"
    candidate_count = int(manifest_row.get("candidate_count") or 0)
    bam_path = segments_dir / segment / "alignment.bam"

    common = dict(
        segment=segment,
        threshold=threshold,
        min_breadth=min_breadth,
        expected_length_min=expected_length_min,
        expected_length_max=expected_length_max,
        max_n_fraction=max_n_fraction,
        contig=contig,
        candidate_count=candidate_count,
    )

    if (
        manifest_status != "READY"
        or not bam_path.is_file()
        or bam_path.stat().st_size == 0
        or not consensus_path.is_file()
        or consensus_path.stat().st_size == 0
    ):
        write_outputs(
            flag_out,
            stats_out,
            **common,
            status="MISSING",
            coverage_status="MISSING",
            length_status="MISSING",
            n_content_status="MISSING",
            reason=(
                f"normalized segment status={manifest_status}; usable BAM and/or "
                "consensus FASTA unavailable"
            ),
        )
        print(f"[check_coverage] {segment}: normalized inputs unavailable -> MISSING")
        return

    try:
        sequence = read_consensus_sequence(consensus_path)
    except Exception as exc:
        write_outputs(
            flag_out,
            stats_out,
            **common,
            status="ERROR",
            coverage_status="NA",
            length_status="ERROR",
            n_content_status="ERROR",
            reason=f"could not read normalized consensus FASTA: {exc}",
        )
        print(f"[check_coverage] {segment}: consensus FASTA error -> ERROR: {exc}")
        return

    consensus_length = len(sequence)
    n_count = sequence.count("N")
    n_fraction_value = n_count / consensus_length
    if consensus_length < expected_length_min:
        length_status = "FAIL"
    elif consensus_length > expected_length_max:
        # Sequences above the configured upper guide are retained for downstream
        # analysis. Extra terminal sequence can be biologically valid or assay-
        # specific, so this is a review warning rather than a hard failure.
        length_status = "WARNING"
    else:
        length_status = "PASS"
    n_content_status = "PASS" if n_fraction_value <= max_n_fraction else "FAIL"

    try:
        depths, bam_contig = samtools_depths(bam_path)
    except Exception as exc:
        write_outputs(
            flag_out,
            stats_out,
            **common,
            status="ERROR",
            coverage_status="ERROR",
            length_status=length_status,
            n_content_status=n_content_status,
            consensus_length=consensus_length,
            n_count=n_count,
            n_fraction=f"{n_fraction_value:.6f}",
            reason=f"samtools depth failed for normalized BAM: {exc}",
        )
        print(f"[check_coverage] {segment}: samtools depth error -> ERROR: {exc}")
        return

    if not depths:
        write_outputs(
            flag_out,
            stats_out,
            **common,
            status="MISSING",
            coverage_status="MISSING",
            length_status=length_status,
            n_content_status=n_content_status,
            consensus_length=consensus_length,
            n_count=n_count,
            n_fraction=f"{n_fraction_value:.6f}",
            reason="samtools depth returned no reference positions",
        )
        print(f"[check_coverage] {segment}: no reference positions -> MISSING")
        return

    median = statistics.median(depths)
    mean = sum(depths) / len(depths)
    minimum = min(depths)
    maximum = max(depths)
    breadth = sum(depth >= threshold for depth in depths) / len(depths)
    coverage_status = "PASS" if median >= threshold and breadth >= min_breadth else "FAIL"
    status = (
        "PASS"
        if coverage_status == "PASS"
        and length_status in {"PASS", "WARNING"}
        and n_content_status == "PASS"
        else "FAIL"
    )
    selected_contig = contig if contig != "NA" else bam_contig

    write_outputs(
        flag_out,
        stats_out,
        **{**common, "contig": selected_contig},
        status=status,
        coverage_status=coverage_status,
        length_status=length_status,
        n_content_status=n_content_status,
        positions=len(depths),
        consensus_length=consensus_length,
        n_count=n_count,
        n_fraction=f"{n_fraction_value:.6f}",
        median=f"{median:.2f}",
        mean=f"{mean:.2f}",
        minimum=f"{minimum:.2f}",
        maximum=f"{maximum:.2f}",
        breadth=f"{breadth:.3f}",
        reason=(
            "coverage computed from normalized BAM with samtools depth -aa -q 0 -Q 0; "
            f"segment PASS requires median depth >= {threshold:.2f}x, "
            f"breadth at >= {threshold:.2f}x >= {min_breadth:.3f}, "
            f"consensus length >= {expected_length_min} nt (>{expected_length_max} nt is warning-only), and "
            f"N fraction <= {max_n_fraction:.4f}"
        ),
    )
    print(
        f"[check_coverage] {segment}: contig={selected_contig}; positions={len(depths)}; "
        f"median={median:.2f}; mean={mean:.2f}; "
        f"breadth_at_{threshold:g}x={breadth:.3f}; length={consensus_length}; "
        f"N_fraction={n_fraction_value:.4f}; coverage={coverage_status}; "
        f"length_qc={length_status}; N_qc={n_content_status} -> {status}"
    )


if __name__ == "__main__":
    main()
