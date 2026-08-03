################################################################################
# Avian influenza Nanopore workflow
#
# FASTQ -> NanoPlot -> Porechop -> fastplong -> IRMA -> coverage filtering
#       -> Medaka consensus/VCF -> BLAST -> summary -> H5N1 screen -> GenoFLU
#
# The workflow is architecture-neutral. Tool portability is controlled through
# Conda environments plus an IRMA runtime that can use Apptainer/Singularity on
# Linux clusters or Docker on laptops.
################################################################################

configfile: "config.yaml"

import csv
import os
import shlex
import shutil
import statistics as stats
import subprocess
from pathlib import Path

from snakemake.io import glob_wildcards

shell.executable("/bin/bash")

# -----------------------------------------------------------------------------
# Configuration helpers
# -----------------------------------------------------------------------------
def as_bool(value):
    """Parse YAML booleans and common string representations safely."""
    if isinstance(value, bool):
        return value
    if value is None:
        return False
    if isinstance(value, str):
        normalized = value.strip().lower()
        if normalized in {"true", "yes", "y", "1", "on"}:
            return True
        if normalized in {"false", "no", "n", "0", "off", ""}:
            return False
    return bool(value)


def config_path(key, default):
    """Return a normalized path string while preserving relative-path behavior."""
    return os.path.normpath(str(config.get(key, default)))


READS = config_path("reads_dir", "data")
RESULTS = config_path("results_dir", "results")
READ_PATTERN = str(config.get("reads_pattern", "{sample}.fastq.gz"))

if "{sample}" not in READ_PATTERN:
    raise ValueError("config key 'reads_pattern' must contain the wildcard {sample}")

READ_GLOB = os.path.join(READS, READ_PATTERN)

COVERAGE_MIN = float(config.get("coverage_min_depth", 50.0))

IRMA_IMAGE = str(config.get("irma_image", "docker://ghcr.io/cdcgov/irma:latest"))
IRMA_MODULE = str(config.get("irma_module", "FLU-minion"))
IRMA_RUNTIME = str(config.get("irma_runtime", "auto")).strip().lower()

BLAST_DB = config_path("blast_db", "data/flu_db/flu")
PORECHOP_CMD = config.get("porechop_command", "porechop_abi")

FASTPLONG_MEAN_QUAL = int(config.get("fastplong_mean_quality", 10))
FASTPLONG_MIN_LENGTH = int(config.get("fastplong_min_length", 500))

MEDAKA_MODEL = config.get("medaka_model")
MEDAKA_MODEL_ARG = (
    f"--model {shlex.quote(str(MEDAKA_MODEL))}" if MEDAKA_MODEL else ""
)
MEDAKA_FAIL_SOFT = as_bool(config.get("medaka_fail_soft", True))

NANOPLOT_INSTALL_CHROME = as_bool(
    config.get("nanoplot_install_chrome", True)
)
RUN_GENOFLU = as_bool(config.get("run_genoflu", True))

# Rule-level threads can be overridden in config or by a Snakemake profile.
NANOPLOT_THREADS = int(config.get("nanoplot_threads", 1))
FASTPLONG_THREADS = int(config.get("fastplong_threads", 4))
IRMA_THREADS = int(config.get("irma_threads", 4))
MEDAKA_THREADS = int(config.get("medaka_threads", 2))
BLAST_THREADS = int(config.get("blast_threads", 2))

# Canonical influenza A segment families and deterministic reporting order.
SEGSET = {"PB2", "PB1", "PA", "NP", "MP", "NS"}
SEGMENT_ORDER = {
    "HA": 0,
    "NA": 1,
    "PB2": 2,
    "PB1": 3,
    "PA": 4,
    "NP": 5,
    "MP": 6,
    "NS": 7,
}

# Discover samples from the configurable FASTQ pattern.
SAMPLES = sorted(set(glob_wildcards(READ_GLOB).sample))
if not SAMPLES:
    print(
        f"WARNING: no samples matched {READ_GLOB!r}. "
        "Check reads_dir and reads_pattern in config.yaml."
    )


# -----------------------------------------------------------------------------
# Deterministic path and subtype-selection utilities
# -----------------------------------------------------------------------------
def _sorted_rglob(base: Path, pattern: str):
    return sorted(base.rglob(pattern), key=lambda path: path.as_posix())


def _normalized_header(value: str):
    return "_".join(
        part for part in "".join(
            char.lower() if char.isalnum() else " " for char in value
        ).split()
        if part
    )


def _median_depth_from_coverage_table(table_path: Path):
    """Return median depth from an IRMA coverage table, or None if unparseable."""
    try:
        with table_path.open() as handle:
            reader = csv.reader(handle, delimiter="\t")
            header = next(row for row in reader if row)
            normalized = [_normalized_header(column) for column in header]

            depth_index = None
            for key in ("coverage_depth", "depth", "read_depth", "readdepth"):
                if key in normalized:
                    depth_index = normalized.index(key)
                    break
            if depth_index is None:
                return None

            depths = []
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

        return stats.median(depths) if depths else None
    except (OSError, StopIteration):
        return None


def segments_for_sample(wildcards):
    """Resolve segment families after the IRMA checkpoint completes."""
    project = Path(checkpoints.irma.get(sample=wildcards.sample).output.project)
    fastas = sorted(project.glob("A_*.fasta"), key=lambda path: path.name)

    segments = set()
    for fasta in fastas:
        name = fasta.stem[2:]  # remove leading A_
        if name in SEGSET:
            segments.add(name)
        elif name.startswith("HA"):
            segments.add("HA")
        elif name.startswith("NA"):
            segments.add("NA")

    return sorted(
        segments,
        key=lambda segment: (SEGMENT_ORDER.get(segment, 99), segment),
    )


def _selected_irma_stem(wildcards):
    """
    Select one IRMA contig consistently for a segment family.

    For HA/NA with multiple subtype-specific outputs, choose the candidate whose
    matching IRMA coverage table has the highest median depth. This keeps the BAM,
    FASTA, coverage decision, Medaka, and downstream subtype screen aligned.
    """
    project = Path(checkpoints.irma.get(sample=wildcards.sample).output.project)
    segment = wildcards.segment

    exact_fasta = project / f"A_{segment}.fasta"
    if exact_fasta.exists():
        return exact_fasta.stem

    candidates = sorted(project.glob(f"A_{segment}*.fasta"), key=lambda path: path.name)
    if not candidates:
        raise ValueError(
            f"No IRMA FASTA found for sample={wildcards.sample}, segment={segment}"
        )
    if len(candidates) == 1:
        return candidates[0].stem

    table_dir = project / "tables"
    scored = []
    for fasta in candidates:
        table = table_dir / f"{fasta.stem}-coverage.txt"
        median = _median_depth_from_coverage_table(table)
        if median is not None:
            scored.append((median, fasta.stem))

    if scored:
        # Highest median depth; name is a deterministic tie-breaker.
        scored.sort(key=lambda item: (-item[0], item[1]))
        return scored[0][1]

    # Deterministic fallback if coverage tables are absent or unparseable.
    return candidates[0].stem


def bam_path(wildcards):
    project = Path(checkpoints.irma.get(sample=wildcards.sample).output.project)
    stem = _selected_irma_stem(wildcards)

    exact = project / f"{stem}.bam"
    if exact.exists():
        return str(exact)

    matches = _sorted_rglob(project, f"{stem}.bam")
    if not matches:
        raise ValueError(
            f"No BAM found for sample={wildcards.sample}, segment={wildcards.segment}, "
            f"selected_contig={stem}"
        )
    return str(matches[0])


def fasta_path(wildcards):
    project = Path(checkpoints.irma.get(sample=wildcards.sample).output.project)
    stem = _selected_irma_stem(wildcards)

    exact = project / f"{stem}.fasta"
    if exact.exists():
        return str(exact)

    matches = _sorted_rglob(project, f"{stem}.fasta")
    if not matches:
        raise ValueError(
            f"No FASTA found for sample={wildcards.sample}, segment={wildcards.segment}, "
            f"selected_contig={stem}"
        )
    return str(matches[0])


def irma_table_dir(wildcards):
    project = Path(checkpoints.irma.get(sample=wildcards.sample).output.project)
    return str(project / "tables")


def irma_project_dir(wildcards):
    return str(checkpoints.irma.get(sample=wildcards.sample).output.project)


def blast_db_files(_wildcards):
    """Stage all files belonging to the configured BLAST database prefix."""
    prefix = Path(BLAST_DB)
    files = sorted(prefix.parent.glob(f"{prefix.name}.*"), key=lambda path: path.name)
    if not files:
        raise ValueError(
            f"No BLAST database files found for prefix {BLAST_DB!r}. "
            "Build the database or update blast_db in config.yaml."
        )
    return [str(path) for path in files]


# -----------------------------------------------------------------------------
# Final targets
# -----------------------------------------------------------------------------
FINAL_TARGETS = [
    expand(f"{RESULTS}/{{sample}}/irma/project", sample=SAMPLES),
    expand(f"{RESULTS}/{{sample}}/coverage/coverage.tsv", sample=SAMPLES),
    expand(f"{RESULTS}/{{sample}}/summary/blast_top_hits.csv", sample=SAMPLES),
    expand(
        f"{RESULTS}/{{sample}}/merged/consensus_all_segments.fasta",
        sample=SAMPLES,
    ),
    expand(
        f"{RESULTS}/{{sample}}/summary/{{sample}}.sample_summary.tsv",
        sample=SAMPLES,
    ),
    expand(
        f"{RESULTS}/{{sample}}/summary/{{sample}}.sample_summary.html",
        sample=SAMPLES,
    ),
]

if RUN_GENOFLU:
    FINAL_TARGETS.extend(
        expand(
            f"{RESULTS}/{{sample}}/genoflu/GenoFLU.tsv",
            sample=SAMPLES,
        )
    )


rule all:
    input:
        FINAL_TARGETS
	

# -----------------------------------------------------------------------------
# NanoPlot
# -----------------------------------------------------------------------------
rule nanoplot:
    input:
        fastq=READ_GLOB
    output:
        done=touch(f"{RESULTS}/{{sample}}/nanoplot/done.txt")
    log:
        f"{RESULTS}/{{sample}}/nanoplot/nanoplot.log"
    conda:
        "envs/nanoplot.yaml"
    threads:
        NANOPLOT_THREADS
    resources:
        kaleido=1
    params:
        install_chrome="true" if NANOPLOT_INSTALL_CHROME else "false"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {output.done:q})"
        exec > >(tee -a {log:q}) 2>&1

        if [[ {params.install_chrome:q} == "true" ]]; then
            python - <<'PY'
try:
    from choreographer.cli._cli_utils import get_chrome_sync
    get_chrome_sync()
except Exception as exc:
    print(f"Chrome setup skipped: {{exc}}")
PY
        fi

        NanoPlot \
            --fastq {input.fastq:q} \
            --threads {threads} \
            --tsv_stats \
            -o "$(dirname {output.done:q})" \
            -p {wildcards.sample:q}

        echo done > {output.done:q}
        """


# -----------------------------------------------------------------------------
# Porechop
# -----------------------------------------------------------------------------
rule porechop:
    input:
        fastq=READ_GLOB
    output:
        trimmed=f"{RESULTS}/{{sample}}/porechop/trimmed.fastq"
    log:
        f"{RESULTS}/{{sample}}/porechop/porechop.log"
    conda:
        "envs/porechop.yaml"
    threads:
        1
    resources:
        mem_mb=int(config.get("porechop_mem_mb", 64000)),
        time_min=int(config.get("porechop_time_min", 240))
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {output.trimmed:q})"
        exec > >(tee -a {log:q}) 2>&1
        {PORECHOP_CMD} -i {input.fastq} -o {output.trimmed}
        """


# -----------------------------------------------------------------------------
# fastplong
# -----------------------------------------------------------------------------
rule fastplong:
    input:
        trimmed=f"{RESULTS}/{{sample}}/porechop/trimmed.fastq"
    output:
        filtered=f"{RESULTS}/{{sample}}/fastplong/filtered.fastq.gz",
        html=f"{RESULTS}/{{sample}}/fastplong/report.html",
        json=f"{RESULTS}/{{sample}}/fastplong/report.json"
    log:
        f"{RESULTS}/{{sample}}/fastplong/fastplong.log"
    conda:
        "envs/fastplong.yaml"
    threads:
        FASTPLONG_THREADS
    params:
        mean_quality=FASTPLONG_MEAN_QUAL,
        minimum_length=FASTPLONG_MIN_LENGTH
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {output.filtered:q})"
        exec > >(tee -a {log:q}) 2>&1

        fastplong \
            -i {input.trimmed:q} \
            -o {output.filtered:q} \
            --mean_qual {params.mean_quality} \
            --length_required {params.minimum_length} \
            -h {output.html:q} \
            -j {output.json:q}
        """


# -----------------------------------------------------------------------------
# Rename FASTQ reads with seqtk
# -----------------------------------------------------------------------------
rule seqtk_rename:
    input:
        filtered=f"{RESULTS}/{{sample}}/fastplong/filtered.fastq.gz"
    output:
        renamed=f"{RESULTS}/{{sample}}/fastplong/filtered_renamed.fastq.gz"
    log:
        f"{RESULTS}/{{sample}}/fastplong/seqtk_rename.log"
    conda:
        "envs/seqtk.yaml"
    threads:
        1
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {output.renamed:q})"
        exec > >(tee -a {log:q}) 2>&1
        seqtk rename {input.filtered:q} | gzip -c > {output.renamed:q}
        """


# -----------------------------------------------------------------------------
# IRMA
#
# Runtime selection:
#   auto        -> Apptainer, Singularity, Docker, then local IRMA
#   apptainer   -> require apptainer
#   singularity -> require singularity
#   docker      -> require Docker
#   local       -> require IRMA on PATH
# -----------------------------------------------------------------------------
checkpoint irma:
    input:
        renamed=f"{RESULTS}/{{sample}}/fastplong/filtered_renamed.fastq.gz"
    output:
        project=directory(f"{RESULTS}/{{sample}}/irma/project"),
        segments=directory(f"{RESULTS}/{{sample}}/irma/segments")
    log:
        f"{RESULTS}/{{sample}}/irma/irma.log"
    threads:
        IRMA_THREADS
    params:
        image=IRMA_IMAGE,
        module=IRMA_MODULE,
        runtime=IRMA_RUNTIME
    run:
        input_path = Path(str(input.renamed)).resolve()
        project_path = Path(str(output.project)).resolve()
        segments_path = Path(str(output.segments)).resolve()
        log_path = Path(str(log[0])).resolve()

        project_path.parent.mkdir(parents=True, exist_ok=True)
        segments_path.parent.mkdir(parents=True, exist_ok=True)
        log_path.parent.mkdir(parents=True, exist_ok=True)

        if project_path.exists():
            shutil.rmtree(project_path)
        if segments_path.exists():
            shutil.rmtree(segments_path)

        runtime = str(params.runtime).lower()
        allowed = {"auto", "apptainer", "singularity", "docker", "local"}
        if runtime not in allowed:
            raise ValueError(
                f"Unsupported irma_runtime={runtime!r}; choose one of {sorted(allowed)}"
            )

        if runtime == "auto":
            if shutil.which("apptainer"):
                runtime = "apptainer"
            elif shutil.which("singularity"):
                runtime = "singularity"
            elif shutil.which("docker"):
                runtime = "docker"
            elif shutil.which("IRMA"):
                runtime = "local"
            else:
                raise RuntimeError(
                    "No IRMA runtime found. Install Apptainer/Singularity, Docker, "
                    "or IRMA locally; alternatively set irma_runtime in config.yaml."
                )

        mount_root = Path(
            os.path.commonpath([str(input_path.parent), str(project_path.parent)])
        )

        if runtime in {"apptainer", "singularity"}:
            executable = shutil.which(runtime)
            if executable is None:
                raise RuntimeError(f"Requested {runtime}, but it is not on PATH")
            command = [
                executable,
                "exec",
                "--bind",
                f"{mount_root}:{mount_root}",
                str(params.image),
                "IRMA",
                str(params.module),
                str(input_path),
                str(project_path),
            ]
        elif runtime == "docker":
            executable = shutil.which("docker")
            if executable is None:
                raise RuntimeError("Requested Docker, but docker is not on PATH")
            docker_image = str(params.image)
            if docker_image.startswith("docker://"):
                docker_image = docker_image[len("docker://"):]
            command = [
                executable,
                "run",
                "--rm",
                "-v",
                f"{mount_root}:{mount_root}",
                "-w",
                str(mount_root),
                docker_image,
                "IRMA",
                str(params.module),
                str(input_path),
                str(project_path),
            ]
        else:
            executable = shutil.which("IRMA")
            if executable is None:
                raise RuntimeError("Requested local IRMA, but IRMA is not on PATH")
            command = [
                executable,
                str(params.module),
                str(input_path),
                str(project_path),
            ]

        with log_path.open("w") as log_handle:
            log_handle.write(f"IRMA runtime: {runtime}\n")
            log_handle.write("Command: " + shlex.join(command) + "\n\n")
            log_handle.flush()
            subprocess.run(
                command,
                check=True,
                stdout=log_handle,
                stderr=subprocess.STDOUT,
            )

        segments_path.mkdir(parents=True, exist_ok=True)
        for fasta in sorted(project_path.glob("*.fasta"), key=lambda path: path.name):
            shutil.copy2(fasta, segments_path / fasta.name)


# -----------------------------------------------------------------------------
# Coverage decision per segment family
# -----------------------------------------------------------------------------
rule check_coverage:
    input:
        table_dir=irma_table_dir,
        project_dir=irma_project_dir
    output:
        flag=f"{RESULTS}/{{sample}}/coverage_flags/{{segment}}.flag",
        stats=f"{RESULTS}/{{sample}}/coverage_stats/{{segment}}.tsv"
    params:
        min_median_depth=COVERAGE_MIN
    conda:
        "envs/pysam.yaml"
    script:
        "scripts/check_coverage.py"


# -----------------------------------------------------------------------------
# Per-sample coverage summary
# -----------------------------------------------------------------------------
rule coverage_table:
    input:
        project=irma_project_dir,
        flags=lambda wildcards: [
            f"{RESULTS}/{wildcards.sample}/coverage_flags/{segment}.flag"
            for segment in segments_for_sample(wildcards)
        ]
    output:
        tsv=f"{RESULTS}/{{sample}}/coverage/coverage.tsv"
    conda:
        "envs/pysam.yaml"
    script:
        "scripts/coverage_table.py"


# -----------------------------------------------------------------------------
# Medaka inference
# -----------------------------------------------------------------------------
rule medaka_inference:
    input:
        bam=bam_path,
        fasta=fasta_path,
        flag=f"{RESULTS}/{{sample}}/coverage_flags/{{segment}}.flag"
    output:
        features=f"{RESULTS}/{{sample}}/medaka/{{segment}}/features.hdf"
    log:
        f"{RESULTS}/{{sample}}/medaka/{{segment}}/medaka_inference.log"
    conda:
        "envs/medaka.yaml"
    threads:
        MEDAKA_THREADS
    resources:
        mem_mb=int(config.get("medaka_inference_mem_mb", 16000)),
        time_min=int(config.get("medaka_inference_time_min", 240))
    params:
        model_arg=MEDAKA_MODEL_ARG,
        fail_soft="true" if MEDAKA_FAIL_SOFT else "false"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {output.features:q})"

        if grep -q '^PASS$' {input.flag:q}; then
            if medaka inference {input.bam:q} {output.features:q} \
                --threads {threads} {params.model_arg} \
                > {log:q} 2>&1; then
                [[ -e {output.features:q} ]] || : > {output.features:q}
            elif [[ {params.fail_soft:q} == "true" ]]; then
                echo "Medaka inference failed; creating an empty features file because medaka_fail_soft=true." >> {log:q}
                : > {output.features:q}
            else
                echo "Medaka inference failed and medaka_fail_soft=false." >> {log:q}
                exit 1
            fi
        else
            echo "Segment did not pass coverage; Medaka inference skipped." > {log:q}
            : > {output.features:q}
        fi
        """


# -----------------------------------------------------------------------------
# Medaka consensus
# -----------------------------------------------------------------------------
rule medaka_consensus:
    input:
        features=f"{RESULTS}/{{sample}}/medaka/{{segment}}/features.hdf",
        fasta=fasta_path,
        flag=f"{RESULTS}/{{sample}}/coverage_flags/{{segment}}.flag"
    output:
        consensus=f"{RESULTS}/{{sample}}/medaka/{{segment}}/consensus.fasta"
    log:
        f"{RESULTS}/{{sample}}/medaka/{{segment}}/medaka_sequence.log"
    conda:
        "envs/medaka.yaml"
    threads:
        MEDAKA_THREADS
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {output.consensus:q})"

        if grep -q '^PASS$' {input.flag:q} && [[ -s {input.features:q} ]]; then
            medaka sequence \
                --threads {threads} \
                {input.features:q} \
                {input.fasta:q} \
                {output.consensus:q} \
                > {log:q} 2>&1
        else
            echo "Segment did not pass coverage or features are empty; consensus skipped." > {log:q}
            : > {output.consensus:q}
        fi
        """


# -----------------------------------------------------------------------------
# Medaka VCF
# -----------------------------------------------------------------------------
rule medaka_vcf:
    input:
        features=f"{RESULTS}/{{sample}}/medaka/{{segment}}/features.hdf",
        fasta=fasta_path,
        flag=f"{RESULTS}/{{sample}}/coverage_flags/{{segment}}.flag"
    output:
        vcf=f"{RESULTS}/{{sample}}/medaka/{{segment}}/variants.vcf"
    log:
        f"{RESULTS}/{{sample}}/medaka/{{segment}}/medaka_variant.log"
    conda:
        "envs/medaka.yaml"
    threads:
        1
    resources:
        mem_mb=int(config.get("medaka_variant_mem_mb", 8000)),
        time_min=int(config.get("medaka_variant_time_min", 60))
    params:
        fail_soft="true" if MEDAKA_FAIL_SOFT else "false"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {output.vcf:q})"

        medaka --version >> {log:q} 2>&1 || true
        python -c "import medaka; print('python medaka __version__ =', getattr(medaka, '__version__', '<no __version__>'))" >> {log:q} 2>&1 || true

        if grep -q '^PASS$' {input.flag:q} && [[ -s {input.features:q} ]]; then
            if medaka variant \
                --features {input.features:q} \
                --reference {input.fasta:q} \
                --output {output.vcf:q} \
                >> {log:q} 2>&1; then
                if [[ ! -s {output.vcf:q} ]]; then
                    found_vcf="$(find "$(dirname {output.vcf:q})" -type f -name '*.vcf' -size +0c | LC_ALL=C sort | head -n 1 || true)"
                    if [[ -n "$found_vcf" && "$found_vcf" != {output.vcf:q} ]]; then
                        cp "$found_vcf" {output.vcf:q}
                    fi
                fi
                [[ -e {output.vcf:q} ]] || : > {output.vcf:q}
            elif [[ {params.fail_soft:q} == "true" ]]; then
                echo "Medaka variant failed; creating an empty VCF because medaka_fail_soft=true." >> {log:q}
                : > {output.vcf:q}
            else
                echo "Medaka variant failed and medaka_fail_soft=false." >> {log:q}
                exit 1
            fi
        else
            echo "Segment did not pass coverage or features are empty; VCF generation skipped." >> {log:q}
            : > {output.vcf:q}
        fi
        """


# -----------------------------------------------------------------------------
# BLAST
# -----------------------------------------------------------------------------
rule blastn:
    input:
        fasta=f"{RESULTS}/{{sample}}/medaka/{{segment}}/consensus.fasta",
        flag=f"{RESULTS}/{{sample}}/coverage_flags/{{segment}}.flag",
        database=blast_db_files
    output:
        txt=f"{RESULTS}/{{sample}}/blast/{{segment}}.blast.txt"
    log:
        f"{RESULTS}/{{sample}}/blast/{{segment}}.log"
    conda:
        "envs/blast.yaml"
    threads:
        BLAST_THREADS
    params:
        database_prefix=BLAST_DB
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {output.txt:q})"

        if [[ ! -s {input.fasta:q} ]] || ! grep -q '^PASS$' {input.flag:q}; then
            echo "Segment did not pass coverage or consensus is empty; BLAST skipped." > {log:q}
            : > {output.txt:q}
            exit 0
        fi

        blastn \
            -query {input.fasta:q} \
            -db {params.database_prefix:q} \
            -outfmt 6 \
            -num_threads {threads} \
            2> {log:q} \
          | LC_ALL=C sort -t $'\t' -k12,12gr -k11,11g -k3,3gr \
          > {output.txt:q}
        """


# -----------------------------------------------------------------------------
# BLAST summary
# -----------------------------------------------------------------------------
rule summarize_blast:
    input:
        blast_files=lambda wildcards: [
            f"{RESULTS}/{wildcards.sample}/blast/{segment}.blast.txt"
            for segment in segments_for_sample(wildcards)
        ]
    output:
        csv=f"{RESULTS}/{{sample}}/summary/blast_top_hits.csv"
    conda:
        "envs/py-tools.yaml"
    script:
        "scripts/summarize_blast.py"


# -----------------------------------------------------------------------------
# Concatenate polished segment consensuses
#
# The VCF list is an explicit dependency so variant files are generated during
# a normal rule-all run, even though only consensus FASTAs are concatenated.
# -----------------------------------------------------------------------------
rule concat_consensus:
    input:
        consensus=lambda wildcards: [
            f"{RESULTS}/{wildcards.sample}/medaka/{segment}/consensus.fasta"
            for segment in segments_for_sample(wildcards)
        ],
        vcfs=lambda wildcards: [
            f"{RESULTS}/{wildcards.sample}/medaka/{segment}/variants.vcf"
            for segment in segments_for_sample(wildcards)
        ]
    output:
        merged=f"{RESULTS}/{{sample}}/merged/consensus_all_segments.fasta"
    log:
        f"{RESULTS}/{{sample}}/merged/concat_consensus.log"
    run:
        output_path = Path(str(output.merged))
        log_path = Path(str(log[0]))
        output_path.parent.mkdir(parents=True, exist_ok=True)
        log_path.parent.mkdir(parents=True, exist_ok=True)

        consensus_paths = [Path(str(path)) for path in input.consensus]
        with output_path.open("wb") as destination:
            for consensus_path in consensus_paths:
                if not consensus_path.exists() or consensus_path.stat().st_size == 0:
                    continue
                with consensus_path.open("rb") as source:
                    shutil.copyfileobj(source, destination)

        log_path.write_text(
            "Concatenated consensus files:\n"
            + "".join(f"- {path}\n" for path in consensus_paths)
        )


# -----------------------------------------------------------------------------
# H5N1 screen
#
# This is an IRMA-supported H5/N1 screening criterion, not an independent
# definitive subtype call. It requires the selected HA and NA coverage tables
# to correspond to H5 and N1 and both to pass the configured depth threshold.
# -----------------------------------------------------------------------------
rule detect_h5n1:
    input:
        ha_flag=f"{RESULTS}/{{sample}}/coverage_flags/HA.flag",
        na_flag=f"{RESULTS}/{{sample}}/coverage_flags/NA.flag",
        ha_stats=f"{RESULTS}/{{sample}}/coverage_stats/HA.tsv",
        na_stats=f"{RESULTS}/{{sample}}/coverage_stats/NA.tsv"
    output:
        flag=f"{RESULTS}/{{sample}}/genoflu/h5n1.flag"
    log:
        f"{RESULTS}/{{sample}}/genoflu/detect_h5n1.log"
    params:
        threshold=COVERAGE_MIN
    run:
        def read_flag(path):
            try:
                return Path(str(path)).read_text().strip()
            except OSError:
                return "MISSING"

        def chosen_table(path):
            try:
                with Path(str(path)).open() as handle:
                    reader = csv.DictReader(handle, delimiter="\t")
                    row = next(reader, None)
                    return row.get("chosen_table", "NA") if row else "NA"
            except OSError:
                return "NA"

        ha_status = read_flag(input.ha_flag)
        na_status = read_flag(input.na_flag)
        ha_table = chosen_table(input.ha_stats)
        na_table = chosen_table(input.na_stats)

        is_h5 = ha_table.startswith("A_HA_H5")
        is_n1 = na_table.startswith("A_NA_N1")
        passed = ha_status == "PASS" and na_status == "PASS" and is_h5 and is_n1

        output_path = Path(str(output.flag))
        log_path = Path(str(log[0]))
        output_path.parent.mkdir(parents=True, exist_ok=True)
        log_path.parent.mkdir(parents=True, exist_ok=True)

        output_path.write_text("PASS\n" if passed else "FAIL\n")
        log_path.write_text(
            f"sample={wildcards.sample}\n"
            f"coverage_threshold={params.threshold}\n"
            f"HA_status={ha_status}\n"
            f"HA_chosen_table={ha_table}\n"
            f"NA_status={na_status}\n"
            f"NA_chosen_table={na_table}\n"
            f"H5N1_screen={'PASS' if passed else 'FAIL'}\n"
        )


# -----------------------------------------------------------------------------
# Produce the summary of the individual samples
# -----------------------------------------------------------------------------
rule sample_summary:
    input:
        fastplong=f"{RESULTS}/{{sample}}/fastplong/report.json",
        coverage=f"{RESULTS}/{{sample}}/coverage/coverage.tsv",
        blast=f"{RESULTS}/{{sample}}/summary/blast_top_hits.csv",
        h5n1=f"{RESULTS}/{{sample}}/genoflu/h5n1.flag",
        genoflu=f"{RESULTS}/{{sample}}/genoflu/GenoFLU.tsv",
        consensus=f"{RESULTS}/{{sample}}/merged/consensus_all_segments.fasta"
    output:
        tsv=f"{RESULTS}/{{sample}}/summary/{{sample}}.sample_summary.tsv"
    conda:
        "envs/reporting.yaml"
    script:
        "scripts/sample_summary.py"
        

# -----------------------------------------------------------------------------
# Produce the summary of the individual samples in html
# -----------------------------------------------------------------------------
rule sample_summary_html:
    input:
        summary=f"{RESULTS}/{{sample}}/summary/{{sample}}.sample_summary.tsv",
        coverage=f"{RESULTS}/{{sample}}/coverage/coverage.tsv",
        blast=f"{RESULTS}/{{sample}}/summary/blast_top_hits.csv"
    output:
        html=f"{RESULTS}/{{sample}}/summary/{{sample}}.sample_summary.html"
    params:
        template="scripts/sample_summary.qmd",
        coverage_threshold=COVERAGE_MIN
    conda:
        "envs/reporting.yaml"
    shell:
        r"""
        set -euo pipefail

        output_dir="$(cd "$(dirname {output.html:q})" && pwd)"
        summary_abs="$(cd "$(dirname {input.summary:q})" && pwd)/$(basename {input.summary:q})"
        coverage_abs="$(cd "$(dirname {input.coverage:q})" && pwd)/$(basename {input.coverage:q})"
        blast_abs="$(cd "$(dirname {input.blast:q})" && pwd)/$(basename {input.blast:q})"
        template_abs="$(cd "$(dirname {params.template:q})" && pwd)/$(basename {params.template:q})"

        temp_qmd="$output_dir/.sample_summary.qmd"
        cp "$template_abs" "$temp_qmd"

        (
            cd "$output_dir"

            quarto render ".sample_summary.qmd" \
              --to html \
              --output "{wildcards.sample}.sample_summary.html" \
              -P "sample_id:{wildcards.sample}" \
              -P "summary_file:${{summary_abs}}" \
              -P "coverage_file:${{coverage_abs}}" \
              -P "blast_file:${{blast_abs}}" \
              -P "coverage_threshold:{params.coverage_threshold}"
        )

        rm -f "$temp_qmd"
        rm -rf "$output_dir/.sample_summary_files"
        """

# -----------------------------------------------------------------------------
# GenoFLU, gated by the H5N1 screen
# -----------------------------------------------------------------------------
rule genoflu:
    input:
        flag=f"{RESULTS}/{{sample}}/genoflu/h5n1.flag",
        fasta=f"{RESULTS}/{{sample}}/merged/consensus_all_segments.fasta",
        summary=f"{RESULTS}/{{sample}}/summary/blast_top_hits.csv"
    output:
        tsv=f"{RESULTS}/{{sample}}/genoflu/GenoFLU.tsv"
    log:
        f"{RESULTS}/{{sample}}/genoflu/genoflu.log"
    conda:
        "envs/genoflu.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {output.tsv:q})"

        if grep -q '^PASS$' {input.flag:q}; then
            fasta_dir="$(dirname {input.fasta:q})"
            fasta_base="$(basename {input.fasta:q})"
            (
                cd "$fasta_dir"
                genoflu.py -f "$fasta_base"
            ) | tee {log:q} > {output.tsv:q}
        else
            printf "sample\tstatus\n%s\tnot H5N1\n" {wildcards.sample:q} \
              | tee {log:q} > {output.tsv:q}
        fi
        """
