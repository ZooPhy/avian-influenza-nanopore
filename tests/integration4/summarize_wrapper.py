from pathlib import Path
import runpy
from types import SimpleNamespace

repo_root = Path(snakemake.scriptdir).resolve().parents[1]

fake = SimpleNamespace(
    wildcards=SimpleNamespace(sample="blast_smoke"),
    params=SimpleNamespace(
        min_identity=float(snakemake.params.min_identity),
        min_query_coverage=float(snakemake.params.min_query_coverage),
    ),
    input=SimpleNamespace(
        blast_files=[str(p) for p in snakemake.input.blast_files],
        flags=[str(p) for p in snakemake.input.flags],
    ),
    output=SimpleNamespace(csv=str(snakemake.output.csv)),
)

runpy.run_path(
    str(repo_root / "scripts" / "summarize_blast.py"),
    init_globals={"snakemake": fake},
)
