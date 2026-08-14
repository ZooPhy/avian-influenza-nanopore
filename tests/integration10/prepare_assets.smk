from pathlib import Path

HERE = Path(workflow.basedir).resolve()
WORK = HERE / "work"
SOURCE = WORK / "source"
PROJECT = WORK / "results" / "phase10" / "irma" / "project"
DB = WORK / "blastdb" / "ha_fixture"

SEGMENTS = ("HA", "NA", "PB2", "PB1", "PA", "NP", "MP", "NS")


rule all:
    input:
        PROJECT,
        str(DB) + ".nhr",
        str(DB) + ".nin",
        str(DB) + ".nsq"


rule build_irma_like_project:
    input:
        sams=[SOURCE / f"{segment}.sam" for segment in SEGMENTS],
        fastas=[SOURCE / f"{segment}.fasta" for segment in SEGMENTS]
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

        for segment in HA NA PB2 PB1 PA NP MP NS; do
            fasta={SOURCE}/"$segment".fasta
            sam={SOURCE}/"$segment".sam
            contig="$(head -n 1 "$fasta" | sed 's/^>//')"

            cp "$fasta" \
              {output.project:q}/assemblies/"$contig".fasta

            samtools view -bS "$sam" \
              > {output.project:q}/alignments/"$contig".unsorted.bam

            samtools sort \
              -o {output.project:q}/alignments/"$contig".bam \
              {output.project:q}/alignments/"$contig".unsorted.bam

            rm {output.project:q}/alignments/"$contig".unsorted.bam
            samtools index {output.project:q}/alignments/"$contig".bam
            samtools quickcheck -v {output.project:q}/alignments/"$contig".bam
        done
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
