import csv
from collections import defaultdict
from pathlib import Path

SEGMENT_ORDER = ["PB2", "PB1", "PA", "HA", "NP", "NA", "MP", "NS"]


def segment_from_name(name: str) -> str:
    token = name.strip().split()[0]
    for segment in SEGMENT_ORDER:
        if token == segment or token.endswith("_" + segment):
            return segment
    return token


def read_names(path: Path):
    if not path.exists() or path.stat().st_size == 0:
        return []
    names = []
    for line in path.read_text(errors="replace").splitlines():
        line = line.strip()
        if line and not line.startswith("#"):
            names.append(line.split()[0])
    return names


outdir = Path(str(snakemake.input.outdir))
sample = str(snakemake.wildcards.sample)
output_path = Path(str(snakemake.output.tsv))
output_path.parent.mkdir(parents=True, exist_ok=True)

pass_names = read_names(outdir / f"{sample}.vadr.pass.list")
fail_names = read_names(outdir / f"{sample}.vadr.fail.list")
status = {segment_from_name(name): "PASS" for name in pass_names}
status.update({segment_from_name(name): "FAIL" for name in fail_names})

alert_counts = defaultdict(int)
fatal_counts = defaultdict(int)
alt_path = outdir / f"{sample}.vadr.alt"
if alt_path.exists() and alt_path.stat().st_size > 0:
    for raw in alt_path.read_text(errors="replace").splitlines():
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        fields = line.split()
        if not fields:
            continue
        segment = segment_from_name(fields[0])
        alert_counts[segment] += 1
        lower = line.lower()
        if "fatal" in lower or " fail" in lower:
            fatal_counts[segment] += 1

segments = sorted(set(status) | set(alert_counts), key=lambda x: SEGMENT_ORDER.index(x) if x in SEGMENT_ORDER else 99)
fieldnames = ["sample_id", "segment", "vadr_status", "alert_count", "fatal_alert_count"]
with output_path.open("w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
    writer.writeheader()
    for segment in segments:
        writer.writerow({
            "sample_id": sample,
            "segment": segment,
            "vadr_status": status.get(segment, "REVIEW" if alert_counts[segment] else "NOT_REPORTED"),
            "alert_count": alert_counts[segment],
            "fatal_alert_count": fatal_counts[segment],
        })
