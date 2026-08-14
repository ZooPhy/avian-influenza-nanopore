from pathlib import Path

HERE = Path(workflow.basedir).resolve()
WORK = HERE / "work"
SOURCE = WORK / "source"
PROJECT = WORK / "results" / "phase9" / "irma" / "project"
DB = WORK / "blastdb" / "ha_fixture"


rule all:
    input:
        PROJECT,
        str(DB) + ".nhr",
        str(DB) + ".nin",
        str(DB) + ".nsq"


rule build_irma_like_project:
    input:
        sam=SOURCE / "HA.sam",
        fasta=SOURCE / "HA.fasta"
    output:
        project=directory(PROJECT)
    conda:
        "../../envs/coverage.yaml"
    shell:
        r"""
        set -euo pipefail

        rm -rf {output.project:q}
        mkdir -p \
          {output.project:q}/assemblies \
          {output.project:q}/alignments

        contig="$(head -n 1 {input.fasta:q} | sed 's/^>//')"

        cp {input.fasta:q} \
          {output.project:q}/assemblies/"$contig".fasta

        samtools view -bS {input.sam:q} \
          > {output.project:q}/alignments/"$contig".unsorted.bam

        samtools sort \
          -o {output.project:q}/alignments/"$contig".bam \
          {output.project:q}/alignments/"$contig".unsorted.bam

        rm {output.project:q}/alignments/"$contig".unsorted.bam

        samtools index \
          {output.project:q}/alignments/"$contig".bam

        samtools quickcheck -v \
          {output.project:q}/alignments/"$contig".bam
        """


rule build_blast_db:
    input:
        SOURCE / "blast_reference.fasta"
    output:
        nhr=str(DB) + ".nhr",
        nin=str(DB) + ".nin",
        nsq=str(DB) + ".nsq"
    conda:
        "../../envs/blast.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {DB:q})"
        makeblastdb \
          -in {input:q} \
          -dbtype nucl \
          -parse_seqids \
          -out {DB:q}
        """
