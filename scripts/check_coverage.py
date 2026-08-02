#!/usr/bin/env python3
from pathlib import Path
import csv
import re
import statistics as stats
from typing import Optional


def norm_col(name: str) -> str:
    return re.sub(r"[^a-z0-9]+", "_", name.lower()).strip("_")


def find_depth_column(header_row):
    """Find a sensible depth column in an IRMA coverage table header."""
    norm = [norm_col(c) for c in header_row]
    for key in ("coverage_depth", "depth", "read_depth", "readdepth"):
        if key in norm:
            return norm.index(key)
    raise ValueError(f"Could not find a depth column in header: {header_row}")


def read_depths(table_path: Path):
    """Read numeric depth values from a coverage table."""
    with open(table_path) as f:
        rdr = csv.reader(f, delimiter="\t")
        hdr = next(r for r in rdr if r)
        dcol = find_depth_column(hdr)

        depths = []
        for r in rdr:
            if not r or dcol >= len(r):
                continue
            v = r[dcol].strip()
            if not v or v.upper() == "NA":
                continue
            try:
                depths.append(float(v))
            except ValueError:
                continue

    if not depths:
        raise ValueError(f"No numeric depth values found in {table_path}")
    return depths


def derive_irma_key_from_segment_file(seg_file: Optional[Path], segment: str) -> Optional[str]:
    """
    Given a segment FASTA/BAM like:
      A_HA_H5.fasta  -> key = A_HA_H5
      A_NA_N1.bam    -> key = A_NA_N1
      A_PB2.fasta    -> key = A_PB2
    Returns None if it can't parse or seg_file is missing.
    """
    if seg_file is None:
        return None
    # IMPORTANT: optional() may hand us a placeholder that doesn't exist.
    if not seg_file.exists():
        return None

    # Match: A_<segment> or A_<segment>_<subtype>, up to the file extension
    m = re.search(rf"^(A_{re.escape(segment)}(?:_[A-Za-z0-9]+)?)\.", seg_file.name)
    return m.group(1) if m else None


def write_outputs(flag_out: Path, stats_out: Path, segment: str, min_med: float,
                  status: str, chosen: str = "NA", n: int = 0,
                  median="NA", mean="NA", dmin="NA", dmax="NA",
                  n_candidates: int = 0, reason: str = ""):
    flag_out.parent.mkdir(parents=True, exist_ok=True)
    stats_out.parent.mkdir(parents=True, exist_ok=True)

    with open(flag_out, "w") as f:
        f.write(f"{status}\n")

    with open(stats_out, "w") as f:
        f.write(
            "segment\tchosen_table\tpositions\tmedian_depth\tmean_depth\tmin_depth\tmax_depth\t"
            "threshold\tstatus\tn_candidates\tselection_reason\n"
        )
        f.write(
            f"{segment}\t{chosen}\t{n}\t{median}\t{mean}\t{dmin}\t{dmax}\t"
            f"{min_med:.2f}\t{status}\t{n_candidates}\t{reason}\n"
        )


def main():
    tables_dir = Path(snakemake.input.table_dir)
    segment = snakemake.wildcards.segment
    flag_out = Path(snakemake.output.flag)
    stats_out = Path(snakemake.output.stats)
    min_med = float(snakemake.params.min_median_depth)

    # segment_file is optional in Snakefile
    seg_file = None
    if hasattr(snakemake.input, "segment_file") and snakemake.input.segment_file:
        seg_file = Path(snakemake.input.segment_file)

    # Candidate tables (HA/NA families can have multiple subtypes)
    cands = sorted(tables_dir.glob(f"A_{segment}*-coverage.txt"))
    if not cands:
        cands = sorted(tables_dir.glob(f"*_{segment}*-coverage.txt"))

    if not cands:
        write_outputs(
            flag_out, stats_out, segment, min_med,
            status="MISSING", n_candidates=0,
            reason=f"no coverage tables in {tables_dir}"
        )
        print(f"[check_coverage] {segment}: no coverage tables in {tables_dir} -> MISSING")
        return

    chosen = None
    reason = ""

    # 1) Best: match subtype-specific table to the actual IRMA segment FASTA used downstream
    key = derive_irma_key_from_segment_file(seg_file, segment)
    if key:
        matches = [p for p in cands if p.name.startswith(key) and p.name.endswith("-coverage.txt")]
        if len(matches) == 1:
            chosen = matches[0]
            reason = f"matched segment_file ({seg_file.name}) -> {key}"
        elif len(matches) > 1:
            chosen = matches[0]
            reason = f"multiple tables matched {key}; picked first"
        else:
            reason = f"no table matched {key}; will fallback"

    # 2) If only one candidate overall, use it
    if chosen is None and len(cands) == 1:
        chosen = cands[0]
        reason = "single candidate"

    # 3) Fallback: pick table with highest median depth
    if chosen is None:
        best = None
        best_med = None
        for p in cands:
            try:
                depths = read_depths(p)
                med = stats.median(depths)
            except Exception:
                continue
            if best is None or med > best_med:
                best = p
                best_med = med

        if best is None:
            write_outputs(
                flag_out, stats_out, segment, min_med,
                status="AMBIGUOUS", n_candidates=len(cands),
                reason="candidates exist but none parseable"
            )
            print(f"[check_coverage] {segment}: {len(cands)} tables but none parseable -> AMBIGUOUS")
            return

        chosen = best
        reason = "fallback: highest median depth"

    # Compute stats & PASS/FAIL
    depths = read_depths(chosen)
    median = stats.median(depths)
    mean = sum(depths) / len(depths)
    n = len(depths)
    dmin = min(depths)
    dmax = max(depths)

    status = "PASS" if median >= min_med else "FAIL"

    write_outputs(
        flag_out, stats_out, segment, min_med,
        status=status, chosen=chosen.name, n=n,
        median=f"{median:.2f}", mean=f"{mean:.2f}",
        dmin=f"{dmin:.2f}", dmax=f"{dmax:.2f}",
        n_candidates=len(cands), reason=reason
    )

    print(
        f"[check_coverage] {segment}: table={chosen.name} ({reason}); "
        f"n={n} median={median:.2f} mean={mean:.2f} min={dmin:.2f} max={dmax:.2f} -> {status}"
    )


if __name__ == "__main__":
    main()
