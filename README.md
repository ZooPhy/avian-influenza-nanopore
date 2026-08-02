# Avian Influenza Nanopore Pipeline

A Snakemake workflow for Oxford Nanopore sequencing analysis of avian influenza viruses.

## Current workflow

FASTQ → NanoPlot → Porechop → fastplong → IRMA → coverage filtering →
Medaka → BLAST → H5N1 screening → GenoFLU

## Supported environments

- Linux ARM64 SLURM cluster
- Laptop support under development
- Future Linux AMD64 and Apple Silicon support

## Quick start on the ARM cluster

```bash
snakemake \
  --profile profiles/slurm-arm \
  --cores 4
