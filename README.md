# ESCAPE

**Enhanced Sequencing and Characterization of Avian Pathogens Engine**

ESCAPE is a portable Snakemake workflow for genomic analysis of avian influenza A virus from Oxford Nanopore sequencing reads. It performs read preprocessing, influenza assembly, segment-level quality assessment, consensus polishing, variant calling, subtype screening, genotype assignment, annotation, and generation of interactive HTML reports.

The workflow has been validated on:

- Apple Silicon macOS using Snakemake, Conda, and Docker Desktop
- Linux ARM64 on a SLURM cluster using Snakemake, Conda, and Apptainer or Singularity

Most tools run in rule-specific Conda environments. IRMA runs in a container selected for the host environment.

## Features

- Oxford Nanopore influenza A analysis
- Porechop ABI adapter trimming
- `fastplong` read-quality and length filtering
- IRMA `FLU-minion` assembly
- Segment-level coverage assessment
- Medaka consensus polishing and variant calling
- BLAST-based segment identification
- H5N1 analytical screening
- GenoFLU genotype assignment
- VADR sequence annotation and validation
- Interactive sample-level HTML reports
- Interactive sequencing-run summary report
- Apple Silicon and Linux ARM64 support
- Docker, Apptainer, Singularity, or local IRMA execution

## Workflow

```text
Nanopore FASTQ
    |
    v
Porechop ABI
(adapter trimming)
    |
    v
fastplong
(long-read quality and length filtering)
    |
    v
seqtk rename
(read identifier normalization)
    |
    v
IRMA FLU-minion
(influenza assembly)
    |
    v
Segment-level coverage assessment
(default median-depth threshold: 50x)
    |
    +--> Medaka consensus polishing
    +--> Medaka variant calling
    |
    v
BLAST segment identification
    |
    v
H5N1 screening
    |
    v
GenoFLU genotype assignment
    |
    v
VADR annotation and validation
    |
    v
Interactive HTML reports
    +--> Sample report
    +--> Run summary report
```

NanoPlot is also available as an optional raw-read quality-control target.

## Repository structure

```text
.
├── Snakefile
├── config.yaml                       # local configuration; not committed
├── config/
│   └── config.example.yaml
├── envs/
│   ├── blast.yaml
│   ├── fastplong.yaml
│   ├── genoflu.yaml
│   ├── medaka.yaml
│   ├── nanoplot.yaml
│   ├── porechop.yaml
│   ├── py-tools.yaml
│   ├── pysam.yaml
│   └── seqtk.yaml
├── scripts/
│   ├── build_blast_db.sh
│   ├── check_coverage.py
│   ├── coverage_table.py
│   ├── normalize_irma_outputs.py
│   ├── summarize_blast.py
│   ├── sample_summary.qmd
│   ├── run_summary.qmd
│   └── report/
│       ├── escape-report.html
│       ├── escape-report.js
│       └── sample-report.css
├── profiles/
│   └── slurm-arm/
├── data/                             # input FASTQ files; not committed
├── resources/
│   ├── fluA_reference.fasta.zip      # downloaded resource; not committed
│   └── flu_db/                       # generated BLAST database; not committed
├── results/                          # generated outputs; not committed
├── README.md
└── .gitignore
```

## Requirements

### All environments

- Git
- Bash
- Conda or Mamba
- Snakemake
- A supported IRMA execution method:
  - Docker
  - Apptainer
  - Singularity
  - a local IRMA installation

The BLAST database setup script additionally requires `curl`, `unzip`, and either local `makeblastdb` or one of the supported container runtimes.

### Apple Silicon laptop

The validated macOS configuration uses:

- Miniforge or another native `osx-arm64` Conda distribution
- Snakemake 9
- Docker Desktop
- native `osx-arm64` Conda environments for non-IRMA rules

Docker Desktop must be installed and running before IRMA executes.

#### Docker Desktop memory

Large or highly multiplexed influenza datasets can require substantially more memory than the Docker Desktop default. A low Docker memory limit can cause IRMA preprocessing to be killed even when the host Mac has ample RAM.

For a high-memory Apple Silicon workstation, a validated configuration is approximately:

- Memory: 96 GB
- Swap: 4 GB or more

Choose values appropriate for the physical RAM available on the host. Leave sufficient memory for macOS and other applications. Verify the memory visible inside Docker with:

```bash
docker run --rm alpine sh -c 'free -h'
```

ESCAPE now checks the IRMA log and expected outputs after execution. Internal failures such as a killed process, an out-of-memory condition, or `found no QC'd data` cause the workflow to stop instead of continuing to downstream reports.

### Linux ARM64 cluster

The validated cluster configuration uses:

- Mamba or Conda
- Snakemake
- Apptainer or Singularity
- SLURM

## Installation

### 1. Clone the repository

```bash
git clone https://github.com/ZooPhy/avian-influenza-nanopore.git
cd avian-influenza-nanopore
```

### 2. Install Miniforge on Apple Silicon macOS

Users who already have a working Conda or Mamba installation may skip this step.

```bash
curl -fsSLo Miniforge3.sh \
  "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-MacOSX-$(uname -m).sh"

bash Miniforge3.sh -b -p "$HOME/miniforge3"
"$HOME/miniforge3/bin/conda" init zsh
exec zsh
```

Configure strict channel priority:

```bash
conda config --set channel_priority strict
```

### 3. Create the Snakemake environment

```bash
mamba create -n snakemake_env \
  -c conda-forge \
  -c bioconda \
  snakemake
```

Activate it:

```bash
conda activate snakemake_env
```

Confirm the installation:

```bash
snakemake --version
```

### 4. Install and start Docker Desktop on macOS

IRMA runs in Docker on macOS. Docker may also be used by the BLAST database setup script when local BLAST+ is not installed.

After starting Docker Desktop, verify that the daemon is available:

```bash
docker info
```

Installing the Docker command-line program alone is not sufficient; the Docker daemon must be running.

## First-time workflow setup

### 1. Create a local configuration file

```bash
cp config/config.example.yaml config.yaml
```

`config.yaml` is intentionally excluded from Git so that paths and machine-specific settings can be changed locally.

### 2. Build the influenza BLAST database

Run the included setup script from the repository root:

```bash
./scripts/build_blast_db.sh
```

If needed, make the script executable first:

```bash
chmod +x scripts/build_blast_db.sh
```

The script:

1. Looks for `resources/fluA_reference.fasta.zip`.
2. Downloads the reference archive when it is absent.
3. Extracts the influenza A reference FASTA.
4. Removes an older database with the same prefix.
5. Builds a nucleotide BLAST database.
6. Writes the database under `resources/flu_db/`.

By default, the reference archive is obtained from the `v0.1.0` release of `ZooPhy/apgap-influenza-pipeline`. The source can be changed with environment variables when needed.

Expected output:

```text
resources/
├── fluA_reference.fasta.zip
└── flu_db/
    ├── fluA_reference.fasta
    ├── fluA_db.nhr
    ├── fluA_db.nin
    ├── fluA_db.nsq
    └── ...
```

The script tries the following database-building methods in order:

1. local `makeblastdb`
2. Docker
3. Apptainer
4. Singularity

Set the resulting database prefix in `config.yaml`:

```yaml
blast_db: "resources/flu_db/fluA_db"
```

The value must be the database prefix, not the directory and not an individual database file.

#### Override the reference source

Use another release tag:

```bash
RELEASE_TAG="vX.Y.Z" ./scripts/build_blast_db.sh
```

Use another GitHub repository and release:

```bash
GITHUB_REPO="owner/repository" \
RELEASE_TAG="vX.Y.Z" \
./scripts/build_blast_db.sh
```

Use a direct archive URL:

```bash
ZIP_URL="https://example.org/fluA_reference.fasta.zip" \
./scripts/build_blast_db.sh
```

For a private GitHub release, authenticate the GitHub CLI before running the script:

```bash
gh auth login
```

The archive may also be downloaded manually and placed at:

```text
resources/fluA_reference.fasta.zip
```

## Input data

Place one compressed Nanopore FASTQ file per sample in the configured input directory.

Example:

```text
data/
├── 24-0514.fastq.gz
├── sample02.fastq.gz
└── sample03.fastq.gz
```

The default pattern is:

```text
{sample}.fastq.gz
```

The text matched by `{sample}` becomes the sample identifier used in output paths.

## Configuration

### Apple Silicon laptop example

```yaml
reads_dir: "data"
reads_pattern: "{sample}.fastq.gz"
results_dir: "results"

coverage_min_depth: 50

porechop_command: "porechop_abi"
porechop_mem_mb: 12000
porechop_time_min: 240

fastplong_mean_quality: 10
fastplong_min_length: 500

irma_image: "docker://ghcr.io/cdcgov/irma:latest"
irma_module: "FLU-minion"
irma_runtime: "docker"

blast_db: "resources/flu_db/fluA_db"

medaka_model: null
medaka_fail_soft: true

run_genoflu: true
```

### Linux ARM64 cluster example

```yaml
reads_dir: "data"
reads_pattern: "{sample}.fastq.gz"
results_dir: "results"

coverage_min_depth: 50

porechop_command: "porechop_abi"
porechop_mem_mb: 80000
porechop_time_min: 720

fastplong_mean_quality: 10
fastplong_min_length: 500

irma_image: "docker://ghcr.io/cdcgov/irma:latest"
irma_module: "FLU-minion"
irma_runtime: "apptainer"

blast_db: "resources/flu_db/fluA_db"

medaka_model: null
medaka_fail_soft: true

run_genoflu: true
```

Use `irma_runtime: "singularity"` instead when Singularity is installed rather than Apptainer. `irma_runtime: "auto"` selects Apptainer, Singularity, Docker, or local IRMA in that order based on what is available.

### Medaka model

The Medaka model should match the Nanopore chemistry and basecalling model used to generate the reads. Leaving `medaka_model` unset preserves the current pipeline behavior, but specifying the correct model is recommended for reproducibility.

## Running on an Apple Silicon laptop

Activate the Snakemake environment and confirm that Docker Desktop is running:

```bash
conda activate snakemake_env
docker info
```

### Dry run

```bash
snakemake \
  --configfile config.yaml \
  --sdm conda \
  --cores 4 \
  --dry-run \
  --printshellcmds
```

`<TBD>` inputs are normal during a dry run of this workflow because IRMA is a checkpoint. Segment-specific jobs are added after IRMA completes and the available influenza segments are known.

### Pre-create Conda environments

```bash
snakemake \
  --configfile config.yaml \
  --sdm conda \
  --cores 4 \
  --conda-create-envs-only
```

Because some segment-specific jobs are created after the IRMA checkpoint, additional Conda environments may be created later during the first complete execution.

### Run the complete workflow

```bash
snakemake \
  --configfile config.yaml \
  --sdm conda \
  --cores 4 \
  --resources mem_mb=90000 kaleido=1 \
  --printshellcmds \
  --rerun-incomplete
```

Snakemake reuses completed outputs automatically when the command is rerun after a failure or interruption. The `mem_mb` value is a Snakemake scheduling resource; it does not configure Docker Desktop memory. Docker resources must be configured separately in Docker Desktop.

### Run NanoPlot optionally

NanoPlot is available as an optional target and is not required by the default final targets.

For one sample:

```bash
snakemake \
  --configfile config.yaml \
  --sdm conda \
  --cores 4 \
  results/24-0514/nanoplot/done.txt
```

## Running on a Linux ARM64 SLURM cluster

Example submission script:

```bash
#!/bin/bash

#SBATCH --job-name=avian-flu
#SBATCH --mem=200G
#SBATCH --partition=arm
#SBATCH --cpus-per-task=4
#SBATCH --qos=grp_mscotch
#SBATCH --error=snakemake.err
#SBATCH --output=snakemake.out
#SBATCH --time=7-00:00:00
#SBATCH --export=NONE
#SBATCH --nodelist=scgh003

set -euo pipefail

module purge
module load mamba/latest

eval "$(conda shell.bash hook)"
conda activate snakemake_env

cd /data/pipp2/Scotch/Snakemake/flu

snakemake \
  --configfile config.yaml \
  --sdm conda \
  --cores "$SLURM_CPUS_PER_TASK" \
  --resources mem_mb=200000 kaleido=1 \
  --conda-cleanup-pkgs cache \
  --latency-wait 300 \
  --printshellcmds \
  --rerun-incomplete
```

Set `irma_runtime` in the cluster `config.yaml` to `apptainer` or `singularity` as appropriate. The current Snakefile invokes the selected IRMA runtime directly, so a separate Snakemake container deployment flag is not required for IRMA.

Older Snakemake installations may use `--use-conda` instead of `--sdm conda`.

## Output structure

Outputs are organized by sample:

```text
results/<sample>/
├── nanoplot/                    # optional
├── porechop/
├── fastplong/
├── irma/
│   ├── project/
│   ├── segments/
│   ├── manifest.tsv
│   └── irma.log
├── coverage/
├── coverage_flags/
├── coverage_stats/
├── medaka/
├── blast/
├── merged/
├── genoflu/
├── vadr/
└── summary/
```

Important sample outputs include:

```text
results/<sample>/fastplong/report.html
results/<sample>/irma/manifest.tsv
results/<sample>/irma/segments/<segment>/consensus.fasta
results/<sample>/coverage/coverage.tsv
results/<sample>/coverage_stats/<segment>.tsv
results/<sample>/medaka/<segment>/consensus.fasta
results/<sample>/medaka/<segment>/variants.vcf
results/<sample>/blast/<segment>.blast.txt
results/<sample>/summary/blast_top_hits.csv
results/<sample>/merged/consensus_all_segments.fasta
results/<sample>/genoflu/h5n1.flag
results/<sample>/genoflu/GenoFLU.tsv
results/<sample>/vadr/<sample>.vadr.log
results/<sample>/summary/<sample>.sample_summary.tsv
results/<sample>/summary/<sample>.sample_summary.html
```

The sequencing-run report is written to:

```text
results/run_summary/run_summary.html
```

## Reports

### Sample report

Each sample report summarizes read filtering, segment recovery, coverage, BLAST assignments, H5N1 screening, GenoFLU results, VADR status, and review flags. The report is self-contained HTML and can be opened locally in a web browser.

Build one sample report with:

```bash
snakemake \
  --configfile config.yaml \
  --sdm conda \
  --cores 4 \
  --resources mem_mb=90000 kaleido=1 \
  --rerun-incomplete \
  results/<sample>/summary/<sample>.sample_summary.html
```

### Run summary report

The run summary combines all configured samples into a single HTML dashboard. Build it with:

```bash
snakemake \
  --configfile config.yaml \
  --sdm conda \
  --cores 4 \
  --resources mem_mb=90000 kaleido=1 \
  --rerun-incomplete \
  results/run_summary/run_summary.html
```

Because the run summary depends on all sample reports, target an individual sample report when rerunning only one barcode.

## Coverage criterion

Each influenza segment is evaluated independently. A segment passes when its median depth is greater than or equal to the configured threshold.

Default:

```yaml
coverage_min_depth: 50
```

Only passing segments proceed through Medaka and BLAST analysis.

For subtype-specific HA and NA outputs, the workflow selects the relevant IRMA coverage table deterministically and records the selected table and reason in the segment coverage statistics.

## H5N1 screening

The H5N1 rule is a screening criterion based on IRMA-supported HA and NA assignments and segment coverage. It requires:

- an H5-associated HA assignment
- an N1-associated NA assignment
- passing HA coverage
- passing NA coverage

A `PASS` result should be interpreted as an analytical screening flag rather than an independent confirmatory subtype test.

When the H5N1 screen does not pass, the GenoFLU output contains a status indicating that the sample was not classified as H5N1 by this screening criterion.

## Troubleshooting

### No samples are detected

Message:

```text
WARNING: no samples matched 'data/{sample}.fastq.gz'
```

Confirm that:

- FASTQ files are present under `reads_dir`
- filenames match `reads_pattern`
- `reads_pattern` contains `{sample}`

### Docker is installed but IRMA cannot start

Message:

```text
Cannot connect to the Docker daemon
```

Start Docker Desktop and verify:

```bash
docker info
```

Then rerun the complete Snakemake command. Completed upstream files will be reused.

### IRMA is killed or reports no QC'd data

Messages may include:

```text
Killed
found no QC'd data
```

This commonly indicates that the container runtime did not have enough memory. On macOS, increase Docker Desktop memory, restart Docker Desktop, remove the affected sample's incomplete IRMA and downstream outputs, and rerun only that sample report target.

Check Docker memory with:

```bash
docker run --rm alpine sh -c 'free -h'
```

Monitor the affected sample with:

```bash
tail -f results/<sample>/irma/irma.log
```

The workflow should stop on these failures rather than interpreting them as a biological negative result.

### Rerun one sample after an IRMA failure

Remove only that sample's IRMA and downstream outputs:

```bash
rm -rf \
  results/<sample>/irma \
  results/<sample>/coverage \
  results/<sample>/coverage_flags \
  results/<sample>/coverage_stats \
  results/<sample>/medaka \
  results/<sample>/blast \
  results/<sample>/merged \
  results/<sample>/genoflu \
  results/<sample>/vadr \
  results/<sample>/summary
```

Then target only that sample's HTML report:

```bash
snakemake \
  --configfile config.yaml \
  --sdm conda \
  --cores 4 \
  --resources mem_mb=90000 kaleido=1 \
  --rerun-incomplete \
  results/<sample>/summary/<sample>.sample_summary.html
```

### BLAST database files are not found

Message:

```text
No BLAST database files found for prefix
```

Build the database:

```bash
./scripts/build_blast_db.sh
```

Then confirm that `config.yaml` contains:

```yaml
blast_db: "resources/flu_db/fluA_db"
```

Verify the files:

```bash
ls -lh resources/flu_db/fluA_db.*
```

### Porechop command is not found on Apple Silicon

The workflow uses the maintained `porechop_abi` package on Apple Silicon and Linux ARM64. Confirm that the configuration contains:

```yaml
porechop_command: "porechop_abi"
```

The corresponding Conda environment should provide an executable named `porechop_abi`.

### A Conda environment fails to solve

First enable strict channel priority:

```bash
conda config --set channel_priority strict
```

Then remove only Snakemake's generated environments and allow them to be rebuilt:

```bash
rm -rf .snakemake/conda
```

Rerun the environment creation or complete workflow command.

## Reproducibility and data management

- `porechop_abi` is used instead of the original Porechop package for Apple Silicon and Linux ARM64 portability.
- Rule-specific Conda environments are stored under `.snakemake/conda/`.
- IRMA runs in Docker on macOS and in Apptainer or Singularity on the ARM64 cluster.
- The BLAST reference archive, generated database, input reads, results, local configuration, and Snakemake working files should not be committed to Git.
- Container tags and package versions should be pinned for a formal release after the validated versions are finalized.

Recommended `.gitignore` entries:

```text
config.yaml
.snakemake/
results/
data/*.fastq
data/*.fastq.gz
data/*.fq
data/*.fq.gz
resources/fluA_reference.fasta.zip
resources/flu_db/
*.log
.DS_Store
```

## Acknowledgements

ESCAPE integrates or builds on the following projects:

- Snakemake
- CDC IRMA
- Porechop ABI
- fastplong
- Medaka
- NCBI BLAST+
- GenoFLU
- VADR
- NanoPlot
- Oxford Nanopore Technologies sequencing software and file formats

Please cite the underlying tools used in an analysis according to their respective documentation and publications.

## Roadmap

Planned or under-development enhancements include:

- Expanded interactive run-level visualizations
- Segment-recovery and coverage heatmaps
- Improved genotype visualizations
- Automated public-health narrative summaries
- Additional export formats, including PDF

## License

Add the project license selected for this repository.

## Citation

Add a `CITATION.cff` file after the repository version, authorship, and preferred citation are finalized.
