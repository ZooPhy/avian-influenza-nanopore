################################################################################
#  Influenza Nanopore → IRMA → Coverage (50×) → Medaka → BLAST → Summary → GenoFLU
################################################################################

configfile: "config.yaml"

import os
from pathlib import Path
from snakemake.io import glob_wildcards

# ─────────────────────────── Config & constants ───────────────────────────
READS   = config["reads_dir"]
RESULTS = config["results_dir"]
IRMA_IMAGE = config.get("irma_image", "docker://ghcr.io/cdcgov/irma:latest")
COVERAGE_MIN = float(config.get("coverage_min_depth", 50.0))
BASEDIR = Path(__file__).resolve().parent

# Discover sample names from fastq
SAMPLES, = glob_wildcards(os.path.join(READS, "{sample}.fastq.gz"))

# Accept canonical internal segments and HA/NA families
SEGSET = {"PB2", "PB1", "PA", "NP", "MP", "NS"}

# ─────────────────────── Deterministic path utilities ─────────────────────
def _sorted_rglob(base: Path, pat: str):
    return sorted(base.rglob(pat), key=lambda p: p.as_posix())

# Resolve segments from IRMA outputs deterministically
# NOTE: use FASTAs as the "segment exists" signal (avoids missing-FASTA crashes)
def segments_for_sample(wc):
    proj = checkpoints.irma.get(sample=wc.sample).output.project
    fastas = sorted(Path(proj).glob("A_*.fasta"), key=lambda p: p.name)

    segs = set()
    for fa in fastas:
        name = fa.stem[2:]  # strip leading "A_"
        if name in SEGSET:
            segs.add(name)
        elif name.startswith("HA"):
            segs.add("HA")
        elif name.startswith("NA"):
            segs.add("NA")
    return sorted(segs)

def bam_path(wc):
    proj = Path(checkpoints.irma.get(sample=wc.sample).output.project)
    exact = _sorted_rglob(proj, f"A_{wc.segment}.bam")
    if exact:
        return str(exact[0])
    wild = _sorted_rglob(proj, f"A_{wc.segment}*.bam")
    assert wild, f"No BAM found for {wc.sample}:{wc.segment}"
    return str(wild[0])

def fasta_path(wc):
    proj = Path(checkpoints.irma.get(sample=wc.sample).output.project)
    exact = _sorted_rglob(proj, f"A_{wc.segment}.fasta")
    if exact:
        return str(exact[0])
    wild = _sorted_rglob(proj, f"A_{wc.segment}*.fasta")
    assert wild, f"No FASTA found for {wc.sample}:{wc.segment}"
    return str(wild[0])

# ─────────────────────────────── Final Targets ───────────────────────────────
rule all:
    input:
        expand(f"{RESULTS}/{{sample}}/irma/project", sample=SAMPLES),
        expand(f"{RESULTS}/{{sample}}/coverage/coverage.tsv", sample=SAMPLES),
        expand(f"{RESULTS}/{{sample}}/summary/blast_top_hits.csv", sample=SAMPLES),
        expand(f"{RESULTS}/{{sample}}/merged/consensus_all_segments.fasta", sample=SAMPLES),
        expand(f"{RESULTS}/{{sample}}/genoflu/GenoFLU.tsv", sample=SAMPLES)

# ─────────────────────────────── NanoPlot ──────────────────────────────────
rule nanoplot:
    input: fastq=f"{READS}/{{sample}}.fastq.gz"
    output: done=touch(f"{RESULTS}/{{sample}}/nanoplot/done.txt")
    conda: "envs/nanoplot.yaml"
    resources: kaleido=1
    shell:
        r"""
        set -euo pipefail
        mkdir -p {RESULTS}/{wildcards.sample}/nanoplot
        python - <<'PY'
try:
    from choreographer.cli._cli_utils import get_chrome_sync
    get_chrome_sync()
except Exception:
    pass
PY
        NanoPlot --fastq {input.fastq} --only-report --tsv_stats \
                 -o {RESULTS}/{wildcards.sample}/nanoplot -p {wildcards.sample}
        echo done > {output.done}
        """

# ────────────────────────────── Porechop ─────────────────────────────────
rule porechop:
    input:
        fastq=f"{READS}/{{sample}}.fastq.gz"
    output:
        trimmed=f"{RESULTS}/{{sample}}/porechop/trimmed.fastq"
    conda: "envs/porechop.yaml"
    threads: 1
    resources:
        mem_mb = int(config.get("porechop_mem_mb", 64000)),
        time_min = int(config.get("porechop_time_min", 240))
    log:
        f"{RESULTS}/{{sample}}/porechop/porechop.log"
    shell:
        r"""
        set -euo pipefail
        mkdir -p {RESULTS}/{wildcards.sample}/porechop
        exec > >(tee -a {log}) 2>&1
        porechop -i {input.fastq} -o {output.trimmed}
        """

# ─────────────────────────────── fastplong ───────────────────────────────────
rule fastplong:
    input:
        trimmed = f"{RESULTS}/{{sample}}/porechop/trimmed.fastq"
    output:
        filtered = f"{RESULTS}/{{sample}}/fastplong/filtered.fastq.gz",
        html     = f"{RESULTS}/{{sample}}/fastplong/report.html",
        json     = f"{RESULTS}/{{sample}}/fastplong/report.json"
    conda: "envs/fastplong.yaml"
    threads: 4
    shell:
        r"""
        mkdir -p {RESULTS}/{wildcards.sample}/fastplong
        fastplong \
            -i {input.trimmed} \
            -o {output.filtered} \
            --mean_qual 10 \
            --length_required 500 \
            -h {output.html} \
            -j {output.json}
        """

# ─────────────────────────────── Rename FASTQ (seqtk) ───────────────────────
rule seqtk_rename:
    input:
        filtered = f"{RESULTS}/{{sample}}/fastplong/filtered.fastq.gz"
    output:
        renamed = f"{RESULTS}/{{sample}}/fastplong/filtered_renamed.fastq.gz"
    conda: "envs/seqtk.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p {RESULTS}/{wildcards.sample}/fastplong
        seqtk rename {input.filtered} | gzip -c > {output.renamed}
        """

# ─────────────────────────────── IRMA ────────────────────────────────────
checkpoint irma:
    input: renamed=f"{RESULTS}/{{sample}}/fastplong/filtered_renamed.fastq.gz"
    output:
        project = directory(f"{RESULTS}/{{sample}}/irma/project"),
        tables  = directory(f"{RESULTS}/{{sample}}/irma/project/tables"),
        segments= directory(f"{RESULTS}/{{sample}}/irma/segments")
    singularity: IRMA_IMAGE
    threads: 4
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.project})
        rm -rf {output.project}
        IRMA FLU-minion {input.renamed} {output.project}
        mkdir -p {output.segments}
        cp {output.project}/*.fasta {output.segments}/ || true
        """

# ──────────────────────────── Check Coverage ──────────────────────────────
# IMPORTANT: do NOT pass a possibly-missing FASTA as an input file.
# Instead pass project_dir (always exists), and the script will look up A_{segment}*.fasta if present.
rule check_coverage:
    input:
        table_dir   = lambda wc: checkpoints.irma.get(sample=wc.sample).output.tables,
        project_dir = lambda wc: checkpoints.irma.get(sample=wc.sample).output.project
    output:
        flag  = f"{RESULTS}/{{sample}}/coverage_flags/{{segment}}.flag",
        stats = f"{RESULTS}/{{sample}}/coverage_stats/{{segment}}.tsv"
    params:
        min_median_depth = float(config.get("coverage_min_depth", 50.0))
    conda: "envs/pysam.yaml"
    script: "scripts/check_coverage.py"

# ─────────────────────────── Coverage summary TSV ───────────────────────────
rule coverage_table:
    input:
        project = f"{RESULTS}/{{sample}}/irma/project",
        flags   = lambda wc: [
            f"{RESULTS}/{wc.sample}/coverage_flags/{seg}.flag"
            for seg in segments_for_sample(wc)
        ]
    output:
        tsv = f"{RESULTS}/{{sample}}/coverage/coverage.tsv"
    conda: "envs/pysam.yaml"
    script: "scripts/coverage_table.py"

# ─────────────────────────── Medaka Inference ─────────────────────────────
rule medaka_variant:
    input:
        bam=bam_path,
        fasta=fasta_path,
        flag=f"{RESULTS}/{{sample}}/coverage_flags/{{segment}}.flag"
    output:
        features=f"{RESULTS}/{{sample}}/medaka/{{segment}}/features.hdf"
    log:
        f"{RESULTS}/{{sample}}/medaka/{{segment}}/medaka_inference.log"
    conda: "envs/medaka.yaml"
    threads: 2
    resources:
        mem_mb=16000, time_min=240
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {output.features})" "{RESULTS}/{wildcards.sample}/benchmarks"

        if grep -q PASS {input.flag}; then
            if medaka inference {input.bam} {output.features} 2> {log}; then
                [[ -e {output.features} ]] || : > {output.features}
            else
                : > {output.features}
            fi
        else
            : > {output.features}
        fi
        """

# ───────────────────────────── Medaka Consensus ───────────────────────────
rule medaka_consensus:
    input:
        features=f"{RESULTS}/{{sample}}/medaka/{{segment}}/features.hdf",
        fasta=fasta_path,
        flag=f"{RESULTS}/{{sample}}/coverage_flags/{{segment}}.flag"
    output:
        consensus=f"{RESULTS}/{{sample}}/medaka/{{segment}}/consensus.fasta"
    conda: "envs/medaka.yaml"
    threads: 2
    shell:
        r"""
        set -euo pipefail
        if grep -q PASS {input.flag} && [ -s {input.features} ]; then
            medaka sequence --threads {threads} {input.features} {input.fasta} {output.consensus}
        else
            : > {output.consensus}
        fi
        """

# ───────────────────────────────── BLAST ─────────────────────────────────
rule blastn:
    input:
        fasta=f"{RESULTS}/{{sample}}/medaka/{{segment}}/consensus.fasta",
        flag=f"{RESULTS}/{{sample}}/coverage_flags/{{segment}}.flag"
    output:
        txt=f"{RESULTS}/{{sample}}/blast/{{segment}}.blast.txt"
    log:
        f"{RESULTS}/{{sample}}/blast/{{segment}}.log"
    conda: "envs/blast.yaml"
    threads: 2
    shell:
        r"""
        set -euo pipefail
        mkdir -p {RESULTS}/{wildcards.sample}/blast
        if [ ! -s {input.fasta} ] || ! grep -q PASS {input.flag}; then
            : > {output.txt}
            exit 0
        fi
        blastn -query {input.fasta} -db data/flu_db/flu -outfmt 6 -num_threads {threads} 2> {log} \
          | LC_ALL=C sort -t $'\t' -k12,12gr -k11,11g -k3,3gr > {output.txt}
        """

# ───────────────────────────── Medaka Variant (VCF) ───────────────────────────
rule medaka_vcf:
    input:
        features=f"{RESULTS}/{{sample}}/medaka/{{segment}}/features.hdf",
        fasta=fasta_path,
        flag=f"{RESULTS}/{{sample}}/coverage_flags/{{segment}}.flag"
    output:
        vcf=f"{RESULTS}/{{sample}}/medaka/{{segment}}/variants.vcf"
    log:
        f"{RESULTS}/{{sample}}/medaka/{{segment}}/medaka_variant.log"
    conda: "envs/medaka.yaml"
    threads: 1
    resources:
        mem_mb=8000, time_min=60
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.vcf}) "$(dirname {log})"

        medaka --version 2>> {log} || true
        python -c "import medaka,sys; print('python medaka __version__ =', getattr(medaka,'__version__','<no __version__>'))" 2>> {log} || true

        if grep -q PASS {input.flag} && [ -s {input.features} ]; then
            if medaka variant --features {input.features} --reference {input.fasta} --output {output.vcf} 2>> {log}; then
                if [ ! -s {output.vcf} ]; then
                    for f in $(dirname {output.vcf})/round_*/*.vcf $(dirname {output.vcf})/*.vcf; do
                        if [ -s "$f" ]; then
                            cp "$f" {output.vcf}
                            break
                        fi
                    done
                fi
            else
                : > {output.vcf}
            fi
        else
            : > {output.vcf}
        fi
        """

# ───────────────────────────── BLAST Summary ──────────────────────────────
rule summarize_blast:
    input:
        blast_files=lambda wc: [
            f"{RESULTS}/{wc.sample}/blast/{seg}.blast.txt"
            for seg in segments_for_sample(wc)
        ]
    output:
        csv=f"{RESULTS}/{{sample}}/summary/blast_top_hits.csv"
    conda: "envs/py-tools.yaml"
    script: "scripts/summarize_blast.py"

# ────────────────────────────── Concat Consensus ───────────────────────────
rule concat_consensus:
    input:
        cons=lambda wc: [
            f"{RESULTS}/{wc.sample}/medaka/{seg}/consensus.fasta"
            for seg in segments_for_sample(wc)
        ]
    output:
        merged=f"{RESULTS}/{{sample}}/merged/consensus_all_segments.fasta"
    conda: "envs/genoflu.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.merged})
        cat {input.cons} > {output.merged}
        """

# ─────────────────────────── H5N1 detector (requires 50× on HA & NA) ───────────────────────────
rule detect_h5n1:
    input:
        segments_dir=lambda wc: checkpoints.irma.get(sample=wc.sample).output.segments,
        ha_flag=f"{RESULTS}/{{sample}}/coverage_flags/HA.flag",
        na_flag=f"{RESULTS}/{{sample}}/coverage_flags/NA.flag",
    output:
        flag=f"{RESULTS}/{{sample}}/genoflu/h5n1.flag"
    log:
        f"{RESULTS}/{{sample}}/genoflu/detect_h5n1.log"
    conda: "envs/py-tools.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p {RESULTS}/{wildcards.sample}/genoflu
        exec > >(tee -a {log}) 2>&1

        echo "[detect_h5n1] segments_dir={input.segments_dir}"
        echo "[detect_h5n1] HA flag=$(cat {input.ha_flag})"
        echo "[detect_h5n1] NA flag=$(cat {input.na_flag})"

        if grep -q PASS {input.ha_flag} && grep -q PASS {input.na_flag} \
           && ls {input.segments_dir}/A_HA_H5*.fasta >/dev/null 2>&1 \
           && ls {input.segments_dir}/A_NA_N1*.fasta >/dev/null 2>&1; then
            echo PASS > {output.flag}
            echo "[detect_h5n1] PASS: found H5 and N1 with ≥50× coverage"
        else
            echo FAIL > {output.flag}
            echo "[detect_h5n1] FAIL: missing tags or coverage"
            echo "[detect_h5n1] present:"; ls -1 {input.segments_dir} || true
        fi
        """

# ──────────────────────────── GenoFLU (gated by detect_h5n1) ────────────────────────────
rule genoflu:
    input:
        flag   = f"{RESULTS}/{{sample}}/genoflu/h5n1.flag",
        fasta  = f"{RESULTS}/{{sample}}/merged/consensus_all_segments.fasta",
        summary= f"{RESULTS}/{{sample}}/summary/blast_top_hits.csv"
    output:
        tsv    = f"{RESULTS}/{{sample}}/genoflu/GenoFLU.tsv"
    log:
        f"{RESULTS}/{{sample}}/genoflu/genoflu.log"
    conda: "envs/genoflu.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p {RESULTS}/{wildcards.sample}/genoflu
        if grep -q PASS {input.flag}; then
            fasta_dir=$(dirname {input.fasta})
            fasta_base=$(basename {input.fasta})
            ( cd "$fasta_dir" && genoflu.py -f "$fasta_base" ) | tee {log} > {output.tsv}
        else
            printf "sample\tstatus\n%s\tnot H5N1\n" "{wildcards.sample}" | tee {log} > {output.tsv}
        fi
        """
