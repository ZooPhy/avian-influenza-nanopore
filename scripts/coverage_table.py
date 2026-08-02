import csv, re
from pathlib import Path
import statistics as stats
import pysam

sample  = snakemake.wildcards.sample
project = Path(snakemake.input.project)
out_tsv = Path(snakemake.output.tsv)
out_tsv.parent.mkdir(parents=True, exist_ok=True)

# Read PASS/FAIL from provided flag files
flag_map = {}
for flag_path in snakemake.input.flags:
    p = Path(flag_path)
    seg = p.stem  # e.g., "HA" from ".../HA.flag"
    status = "FAIL"
    try:
        with open(p) as fh:
            if "PASS" in fh.read():
                status = "PASS"
    except FileNotFoundError:
        status = "MISSING"
    flag_map[seg] = status

SEGSET = ["PB2","PB1","PA","NP","MP","NS"]

def segment_from_name(name: str) -> str:
    # Works for BAM names ("A_HA_H5.bam") and coverage files ("A_HA_H5-coverage.txt")
    n = name
    if n.startswith("A_"):
        n = n[2:]
    n = re.sub(r"-coverage\.txt$", "", n)
    n = re.sub(r"\.bam$", "", n)
    if n.startswith("HA"):
        return "HA"
    if n.startswith("NA"):
        return "NA"
    for s in SEGSET:
        if n.startswith(s):
            return s
    return n.split("_")[0]

def norm_header(h):
    return re.sub(r'[^a-z0-9]+', '_', h.lower()).strip('_')

def stats_from_cov_table(cov_path: Path):
    # Parse IRMA tables: A_* -coverage.txt
    with cov_path.open() as f:
        # Robust header handling
        header = f.readline().rstrip("\n").split("\t")
        keys = [norm_header(c) for c in header]
        try:
            depth_idx = keys.index("coverage_depth")
        except ValueError:
            raise RuntimeError(f"'Coverage Depth' column not found in {cov_path}")
        # Reference_Name usually present; otherwise infer from filename
        ref_idx = keys.index("reference_name") if "reference_name" in keys else None

        depths = []
        ref_name = None
        for line in f:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if ref_idx is not None and ref_name is None:
                ref_name = parts[ref_idx]
            try:
                d = float(parts[depth_idx])
            except (ValueError, IndexError):
                continue
            depths.append(d)

    if not depths:
        return None  # let caller fall back to BAM

    length = len(depths)
    mean_depth = sum(depths) / length
    median_depth = stats.median(depths)
    breadth = sum(1 for d in depths if d > 0) / length
    contig_label = ref_name if ref_name else cov_path.stem.replace("-coverage", "")
    return contig_label, length, median_depth, mean_depth, breadth

def stats_from_bam(bam_path: Path):
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        total_len = 0
        depths_all = []
        contig_label = Path(bam_path).stem  # e.g., "A_HA_H5"
        for ridx, ref in enumerate(bam.references):
            length = bam.lengths[ridx]
            total_len += length
            A, C, G, T = bam.count_coverage(ref, quality_threshold=0)
            for i in range(length):
                depths_all.append(A[i] + C[i] + T[i] + G[i])
        if not depths_all or total_len == 0:
            return contig_label, 0, 0.0, 0.0, 0.0
        mean_depth = sum(depths_all) / total_len
        median_depth = stats.median(depths_all)
        breadth = sum(1 for d in depths_all if d > 0) / total_len
        return contig_label, total_len, median_depth, mean_depth, breadth

rows = []

# Iterate over BAMs to ensure we list exactly the assembled contigs,
# but prefer IRMA coverage tables for the stats.
for bam_path in sorted(project.glob("A_*.bam"), key=lambda p: p.name):
    seg = segment_from_name(bam_path.name)
    cov_path = project / "tables" / f"{bam_path.stem}-coverage.txt"

    if cov_path.exists():
        stats_tuple = stats_from_cov_table(cov_path)
        if stats_tuple is None:
            # empty or unparsable table -> fall back to BAM
            contig_label, length, med, mean, br = stats_from_bam(bam_path)
        else:
            contig_label, length, med, mean, br = stats_tuple
    else:
        contig_label, length, med, mean, br = stats_from_bam(bam_path)

    rows.append({
        "sample": sample,
        "segment": seg,
        "contig": contig_label,
        "length": length,
        "median_depth": f"{med:.2f}",
        "mean_depth": f"{mean:.2f}",
        "breadth_covered": f"{br:.3f}",
        "coverage_flag": flag_map.get(seg, "NA"),
    })

order = {"HA":0,"NA":1,"PB2":2,"PB1":3,"PA":4,"NP":5,"MP":6,"NS":7}
rows.sort(key=lambda r: (order.get(r["segment"], 99), r["contig"]))

with out_tsv.open("w", newline="") as f:
    w = csv.writer(f, delimiter="\t")
    w.writerow(["sample","segment","contig","length","median_depth","mean_depth","breadth_covered","coverage_flag"])
    for r in rows:
        w.writerow([r[k] for k in ["sample","segment","contig","length","median_depth","mean_depth","breadth_covered","coverage_flag"]])
