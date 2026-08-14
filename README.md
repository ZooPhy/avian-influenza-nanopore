# Wild-bird Influenza Genomics and Surveillance (WINGS)

<p align="center">
  <img src="wings_logo.jpg" alt="WINGS logo" width="360">
</p>

**Overview**

WINGS is a portable Snakemake workflow for genomic analysis of avian influenza A virus from Oxford Nanopore sequencing reads. It performs read preprocessing, influenza assembly, segment-level quality assessment, consensus polishing, variant calling, subtype screening, genotype assignment, annotation, and generation of interactive HTML reports.

WINGS was developed in support of the [**Pandemic ESCAPE Center**](https://escape.engr.uky.edu/), with a focus on genomic epidemiology, bioinformatics, and surveillance of avian influenza viruses in wild birds.

The workflow has been validated on:

- Apple Silicon macOS using Snakemake, Conda, and Docker Desktop
- Linux ARM64 on a SLURM cluster using Snakemake, Conda, and Apptainer or Singularity

Most tools run in rule-specific Conda environments. IRMA runs in a container selected for the host environment.

## Features

- Oxford Nanopore influenza A analysis
- Porechop ABI adapter trimming
- `fastplong` read-quality and length filtering without a second adapter-trimming pass
- Sample metadata validation and integration
- IRMA `FLU-minion` assembly
- Segment-level QC using depth, breadth-at-depth, expected-length, and N-content criteria
- Medaka consensus polishing and variant calling
- BLAST-based segment identification with identity/query-coverage evidence and confidence classification
- H5N1 analytical screening
- GenoFLU genotype assignment
- VADR sequence annotation and validation
- Interactive sample-level HTML reports
- Interactive sequencing-run summary report
- Portable `.wings` report bundles containing the run summary, all sample reports, and embedded run-level provenance
- Run-level provenance capturing workflow state, configuration hashes, environment hashes, runtime details, and BLAST database provenance
- Browser-based local report viewing at `wings.scotchlab.org` with no sequencing-data upload
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
(long-read quality and length filtering; adapter trimming disabled)
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
Segment-level QC
(default: median depth >=50x; >=95% of positions at >=50x;
 minimum segment length; <=1% Ns; long segments warned)
    |
    +--> FASTQ basecaller-model detection
    |    (sample-specific Medaka model resolution)
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
    +--> Validated sample metadata
    +--> Sample report
    +--> Run summary report
    +--> Run-level provenance (TSV + JSON)
    +--> Portable WINGS report bundle (.wings; provenance embedded)
```

NanoPlot is also available as an optional raw-read quality-control target.

## Repository structure

```text
.
├── Snakefile
├── view_reports.sh                    # convenience launcher for local HTML reports
├── config.yaml                       # local configuration; not committed
├── config/
│   └── config.example.yaml
├── envs/
│   ├── blast.yaml
│   ├── coverage.yaml
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
│   ├── build_report_bundle.py
│   ├── check_coverage.py
│   ├── coverage_table.py
│   ├── normalize_irma_outputs.py
│   ├── resolve_medaka_model.py
│   ├── validate_metadata.py
│   ├── extract_sample_metadata.py
│   ├── serve_reports.py
│   ├── summarize_blast.py
│   ├── sample_summary.qmd
│   ├── run_summary.qmd
│   ├── write_run_provenance.py
│   └── report/
│       ├── escape-report.html
│       ├── escape-report.js
│       └── sample-report.css
├── profiles/
│   └── slurm-arm/
├── demo/
│   └── wings_demo.wings              # public demonstration bundle for the website
├── data/                             # input FASTQ files; not committed
├── metadata.tsv                      # sample metadata table
├── resources/
│   ├── fluA_reference.fasta.zip      # downloaded resource; not committed
│   └── flu_db/                       # generated BLAST database; not committed
│       └── database_manifest.tsv     # BLAST database provenance manifest
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

WINGS now checks the IRMA log and expected outputs after execution. Internal failures such as a killed process, an out-of-memory condition, or `found no QC'd data` cause the workflow to stop instead of continuing to downstream reports.

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
7. Writes `resources/flu_db/database_manifest.tsv` with the source release, archive and FASTA SHA-256 hashes, `makeblastdb` version, BLAST database format version, build method, container image, and database prefix.

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
    ├── database_manifest.tsv
    └── ...
```

The default container image for database construction is pinned to `ncbi/blast-static:2.17.0`. The script tries the following database-building methods in order:

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

### Sample metadata

WINGS supports a tab-delimited `metadata.tsv` file that is validated against the detected FASTQ samples before sample-level reporting. Metadata are matched to sequencing inputs by `sample_id`. By default, every detected FASTQ sample must have exactly one corresponding metadata record.

The supported schema is:

| Field | Requirement | Description |
|---|---|---|
| `sample_id` | Required | Must exactly match the `{sample}` identifier derived from the FASTQ filename |
| `host` | Required | Host species or host code; use `environmental` when there is no animal host |
| `collection_date` | Required | Collection date in ISO `YYYY-MM-DD` format |
| `country` | Required | Country of collection |
| `specimen_type` | Optional | Specimen or swab type |
| `state` | Optional | State, province, or equivalent first-level administrative area |
| `latitude` | Optional | Decimal latitude from -90 to 90 |
| `longitude` | Optional | Decimal longitude from -180 to 180 |

Example:

```tsv
sample_id	host	specimen_type	collection_date	state	country	latitude	longitude
12-11-2025_barcode03	BLVU	Dry (Oral+Cloacal)	2025-10-20	Kentucky	USA		
21-07-2026_barcode13	RTHA	Dry (Oral+Cloacal)	2026-04-08	Kentucky	USA	38.069739	-84.746138
```

Host values are preserved as supplied; WINGS does not silently expand or normalize species codes. Latitude and longitude may both be blank, but when supplied they must be valid numeric coordinates.

Metadata validation checks required columns and fields, unique `sample_id` values, ISO-formatted collection dates, coordinate ranges, and agreement between metadata records and detected FASTQ samples. The normalized table is written to `results/metadata/validated_metadata.tsv`, and each sample receives a sample-specific metadata table under `results/<sample>/metadata/`.

Configure metadata with:

```yaml
metadata_file: "metadata.tsv"
metadata_require_all_samples: true
```

Set `metadata_require_all_samples: false` only when intentionally allowing FASTQ samples without metadata.

## Configuration

### Apple Silicon laptop example

```yaml
reads_dir: "data"
reads_pattern: "{sample}.fastq.gz"
results_dir: "results"

coverage_min_depth: 50
coverage_min_breadth: 0.95
segment_max_n_fraction: 0.01

porechop_command: "porechop_abi"
porechop_mem_mb: 12000
porechop_time_min: 240

fastplong_mean_quality: 10
fastplong_min_length: 500
fastplong_threads: 4

metadata_file: "metadata.tsv"
metadata_require_all_samples: true

irma_image: "docker://ghcr.io/cdcgov/irma:latest"
irma_module: "FLU-minion"
irma_runtime: "docker"

blast_db: "resources/flu_db/fluA_db"
blast_min_identity: 95.0
blast_min_query_coverage: 90.0
blast_max_target_seqs: 10
blast_max_hsps: 1

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
coverage_min_breadth: 0.95
segment_max_n_fraction: 0.01

porechop_command: "porechop_abi"
porechop_mem_mb: 80000
porechop_time_min: 720

fastplong_mean_quality: 10
fastplong_min_length: 500
fastplong_threads: 4

metadata_file: "metadata.tsv"
metadata_require_all_samples: true

irma_image: "docker://ghcr.io/cdcgov/irma:latest"
irma_module: "FLU-minion"
irma_runtime: "apptainer"

blast_db: "resources/flu_db/fluA_db"
blast_min_identity: 95.0
blast_min_query_coverage: 90.0
blast_max_target_seqs: 10
blast_max_hsps: 1

medaka_model: null
medaka_fail_soft: true

run_genoflu: true
```

Use `irma_runtime: "singularity"` instead when Singularity is installed rather than Apptainer. `irma_runtime: "auto"` selects Apptainer, Singularity, Docker, or local IRMA in that order based on what is available.

### Medaka model

WINGS automatically determines the Medaka consensus model for each sample from the `basecall_model_version_id` metadata embedded in the original Oxford Nanopore FASTQ headers. The workflow inspects multiple FASTQ records, verifies that the detected basecaller model is consistent, and writes the result to `results/<sample>/medaka/model.tsv`.

The detected basecaller model is passed to Medaka using its `:consensus` model selector, allowing the installed Medaka version to resolve the corresponding supported consensus model. This avoids relying on automatic model detection from the normalized IRMA BAM, which may not retain the original Nanopore read-group metadata.

For example:

```text
FASTQ basecaller model:
dna_r10.4.1_e8.2_400bps_hac@v5.0.0

Medaka selector:
dna_r10.4.1_e8.2_400bps_hac@v5.0.0:consensus

Resolved Medaka model:
r1041_e82_400bps_hac_v5.0.0
```

`medaka_model` remains available as an expert configuration override. With `medaka_model: null`, WINGS uses automatic FASTQ-based model resolution. The number of FASTQ records inspected can be configured with `medaka_model_records` (default: 100).

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
├── metadata/
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
results/<sample>/metadata/<sample>.metadata.tsv
results/<sample>/irma/manifest.tsv
results/<sample>/irma/segments/<segment>/consensus.fasta
results/<sample>/coverage/coverage.tsv
results/<sample>/coverage_stats/<segment>.tsv
results/<sample>/medaka/model.tsv
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

Validated run-level metadata and provenance outputs include:

```text
results/metadata/validated_metadata.tsv
results/metadata/metadata_validation.tsv
results/run_summary/run_provenance.tsv
results/run_summary/run_provenance.json
```

The sequencing-run report is written to:

```text
results/run_summary/run_summary.html
```

A portable WINGS report bundle containing the run summary and all sample reports is written to:

```text
results/wings_report_bundle.wings
```

The `.wings` bundle is a self-contained JSON report package intended for browser-based viewing. It contains rendered HTML reports plus the embedded `run_provenance.json` record, but not the raw sequencing reads or intermediate analysis files.

## Run-level provenance

WINGS writes machine-readable provenance for each completed run to:

```text
results/run_summary/run_provenance.tsv
results/run_summary/run_provenance.json
```

The record captures the WINGS Git commit, branch and dirty/clean state; SHA-256 hashes of the `Snakefile` and active `config.yaml`; Snakemake, Python, operating-system and architecture information; the configured IRMA, Medaka, VADR, QC and BLAST settings; the BLAST database manifest and its checksum; and SHA-256 hashes of the Conda environment YAML files. The JSON provenance record is embedded directly in `results/wings_report_bundle.wings`, so the portable bundle carries its computational provenance with the rendered reports.

For a formal reproducible analysis, generate provenance from a clean committed repository state so `workflow.git_dirty` is `false`.

## Reports

### Sample report

Each sample report summarizes validated sample metadata, read filtering, segment recovery, coverage, BLAST assignments, H5N1 screening, GenoFLU results, VADR status, and review flags. Reports are generated as self-contained HTML and are also packaged into the portable `.wings` bundle for browser-based navigation.

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

### Build the portable WINGS report bundle

Build the portable report bundle with:

```bash
snakemake results/wings_report_bundle.wings \
  --configfile config.yaml \
  --sdm conda \
  --cores 4 \
  --resources mem_mb=90000 kaleido=1
```

The resulting `results/wings_report_bundle.wings` file contains the rendered run summary, all rendered sample reports, and the run-level provenance JSON in a single portable package. It can be opened at `wings.scotchlab.org` by selecting or dragging the `.wings` file into the report viewer. The browser reads the bundle locally; the sequencing results are not uploaded to the WINGS website.

For a public demonstration, a deliberately selected example bundle can be placed at:

```text
demo/wings_demo.wings
```

Only use non-sensitive data that are appropriate for public distribution in the demo bundle.

### View reports locally

The preferred way to review a completed analysis is to open `results/wings_report_bundle.wings` at `wings.scotchlab.org`. The site can display the run summary and navigate to individual sample reports directly from the local bundle without uploading the report contents.

The included local report server remains available as an alternative for development, offline use, or direct browsing of the generated HTML files. Some browsers, including Safari, restrict navigation between local `file://` HTML documents, so use the local server rather than opening the HTML files directly.

Launch the local report server from the repository root with:

```bash
./view_reports.sh
```

The wrapper starts `scripts/serve_reports.py`, which by default binds only to the local loopback interface and opens:

```text
http://127.0.0.1:4174/run_summary/run_summary.html
```

The report server serves the existing static HTML reports from the local `results/` directory. By default it binds only to the loopback interface, so it does not upload sequencing data or publish the reports to the network.

Press `Ctrl+C` to stop the server.

The Python launcher can also be run directly:

```bash
python scripts/serve_reports.py
```

To use another port:

```bash
python scripts/serve_reports.py --port 8080
```

To start the server without automatically opening a browser:

```bash
python scripts/serve_reports.py --no-browser
```

## Segment QC criterion

Each influenza segment is evaluated independently before Medaka polishing. Per-position depth is calculated from the normalized IRMA BAM with `samtools depth -aa -q 0 -Q 0`. The `-aa` option ensures that zero-depth reference positions are included in the denominator. Coverage breadth is the fraction of reference positions whose depth is greater than or equal to `coverage_min_depth`. WINGS also evaluates the normalized IRMA consensus FASTA for segment length and the fraction of ambiguous `N` bases. The configured lower length bound is a hard minimum; the upper bound is a review guide rather than a hard failure threshold.

Defaults:

```yaml
coverage_min_depth: 50
coverage_min_breadth: 0.95
segment_max_n_fraction: 0.01
segment_expected_lengths:
  PB2: [2200, 2400]
  PB1: [2200, 2400]
  PA: [2100, 2300]
  HA: [1600, 1800]
  NP: [1450, 1600]
  NA: [1300, 1500]
  MP: [950, 1050]
  NS: [800, 950]
```

With these defaults, a segment passes the hard QC gate when **all** of the following are true: median depth is at least **50x**; at least **95% of reference positions are at 50x or greater**; consensus length is at least the configured segment-specific minimum; and no more than **1% of consensus bases are `N`**. A consensus longer than the configured upper length guide remains eligible for downstream analysis but receives a length `WARNING` for review. The reported `breadth_covered` value therefore represents breadth at the configured depth threshold, not merely the fraction of positions with any coverage.

The segment-length bounds are configurable QC guardrails rather than subtype-confirmation criteria. The lower bound protects against truncated assemblies. The upper bound highlights unexpectedly long consensuses without automatically rejecting sequences that may contain valid terminal or assay-specific sequence. Bounds should be changed only when the assay design or validated biological targets justify different values.

The per-segment statistics record the individual `coverage_status`, `length_status`, and `n_content_status` values as well as the final `overall_status`. The existing `results/<sample>/coverage_flags/<segment>.flag` path is retained for workflow compatibility, but its `PASS` now means that the segment passed the complete segment-QC criterion rather than coverage alone.

Only segments passing the hard segment-QC criteria proceed through Medaka and BLAST analysis; an upper-length warning alone does not block downstream analysis.

For each segment, WINGS selects a normalized IRMA candidate deterministically and records the candidate count, selection status, selected contig, and selection reason in the manifest and segment QC outputs. When IRMA produces more than one candidate for a segment, the selected candidate remains eligible for the normal hard QC gate, but the sample receives a `multiple_irma_candidates` review flag and the ambiguity is surfaced in sample- and run-level reports.

## BLAST evidence and confidence

For each QC-passing segment, WINGS runs `blastn` against the configured influenza A nucleotide database and retains up to `blast_max_target_seqs` subject sequences with at most `blast_max_hsps` HSPs per subject. The default settings are:

```yaml
blast_min_identity: 95.0
blast_min_query_coverage: 90.0
blast_max_target_seqs: 10
blast_max_hsps: 1
```

The raw `results/<sample>/blast/<segment>.blast.txt` files retain the evidence rows. `results/<sample>/summary/blast_top_hits.csv` summarizes the top HSP for each segment with subject accession/title, percent identity, alignment length, query length, query coverage, E-value, bit score, and a confidence state. Query coverage is calculated as alignment length divided by query length for the selected top HSP; HSPs are not merged.

BLAST summary states are:

- `HIGH_CONFIDENCE`: a hit exists and meets both the identity and query-coverage thresholds.
- `LOW_CONFIDENCE`: a hit exists but fails one or both thresholds.
- `NO_HIT`: BLAST ran for a QC-passing segment but returned no hit.
- `SKIPPED_QC`: BLAST was not run because the segment failed the upstream hard QC gate.

The BLAST database provenance is recorded in `resources/flu_db/database_manifest.tsv` and is incorporated into the run-level provenance record.

## H5N1 screening

The H5N1 rule is a screening criterion based on IRMA-supported HA and NA assignments and segment coverage. It requires:

- an H5-associated HA assignment
- an N1-associated NA assignment
- passing HA segment QC
- passing NA segment QC

The H5N1 screen uses three states: `DETECTED` when QC-qualified HA and NA evidence supports H5 and N1; `NOT_DETECTED` when QC-qualified informative HA/NA evidence does not meet the H5+N1 criterion; and `INDETERMINATE` when HA or NA evidence is missing or does not pass segment QC. `DETECTED` is an analytical screening flag rather than an independent confirmatory subtype test. GenoFLU is run only for `DETECTED` samples.

## Troubleshooting

### Metadata validation fails

WINGS validates `metadata.tsv` before sample-level metadata are propagated into reports. Common causes of failure include missing required columns, duplicate `sample_id` values, sample identifiers that do not match FASTQ filenames, non-ISO collection dates, or invalid coordinates.

Check the configured metadata path:

```yaml
metadata_file: "metadata.tsv"
metadata_require_all_samples: true
```

Then inspect the identifiers derived from the input FASTQs and compare them with the first column of `metadata.tsv`:

```bash
printf "FASTQ samples:\n"
find data -maxdepth 1 -type f -name '*.fastq.gz' -print \
  | sed 's#^.*/##; s/\.fastq\.gz$//' \
  | sort

printf "\nMetadata sample_id values:\n"
cut -f1 metadata.tsv | tail -n +2 | sort
```

With `metadata_require_all_samples: true`, every detected FASTQ sample must have a matching metadata record.

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

- `porechop_abi` is used instead of the original Porechop package for Apple Silicon and Linux ARM64 portability, with ab-initio adapter inference enabled.
- `fastplong` performs long-read quality and length filtering with its adapter-trimming step disabled to avoid a second adapter-trimming pass after Porechop ABI.
- Segment depth is calculated with `samtools depth -aa -q 0 -Q 0`; breadth is the fraction of all reference positions meeting the configured depth threshold.
- Sample metadata are validated before report generation and propagated into sample-specific metadata outputs.
- The Oxford Nanopore basecaller model is detected from the original FASTQ metadata for each sample, and the resulting Medaka selector is recorded in `results/<sample>/medaka/model.tsv`.
- Segment QC requires the configured depth, breadth-at-depth, minimum-length, and N-content criteria; by default this is median depth >=50x, >=95% of positions at >=50x, the configured segment-specific minimum length, and <=1% Ns. Consensus lengths above the configured upper guide generate a warning rather than a hard failure.
- Rule-specific Conda environments are stored under `.snakemake/conda/`.
- IRMA runs in Docker on macOS and in Apptainer or Singularity on the ARM64 cluster.
- The BLAST database build records source and build provenance in `resources/flu_db/database_manifest.tsv`; the default BLAST build image is pinned to `ncbi/blast-static:2.17.0`.
- Run-level provenance is written to `results/run_summary/run_provenance.tsv` and `.json`, and the JSON record is embedded in the portable `.wings` bundle.
- The BLAST reference archive, generated database, input reads, results, local configuration, and Snakemake working files should not be committed to Git. Commit the BLAST provenance manifest only when intentionally maintaining a fixed reference build record in the repository.
- `results/wings_report_bundle.wings` is generated from local reports and should be treated as analysis output; do not publish it unless its contents are appropriate for public release.
- Remaining reproducibility gaps should be addressed before a formal release: IRMA and VADR currently use mutable `latest` container tags by default, and the GenoFLU Conda environment does not yet pin a specific GenoFLU version.

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

WINGS integrates or builds on the following projects:

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

- Additional interactive run-level visualizations and comparative views
- Improved genotype visualizations
- Automated public-health narrative summaries
- Additional export formats, including PDF

## Use of large language models and ChatGPT

Large language models, including OpenAI ChatGPT, were used during development of WINGS as a software-development and documentation assistant. Uses included brainstorming workflow design, reviewing and refining code, troubleshooting Snakemake and reporting behavior, and drafting or editing documentation.

All workflow logic, code changes, configuration decisions, and scientific interpretations remain the responsibility of the WINGS developers and should be independently reviewed and validated. ChatGPT is not used by the workflow to generate sequencing results, assemble influenza genomes, assign subtypes or genotypes, call variants, or replace the underlying bioinformatics tools described above.

Users adapting WINGS should apply the same standard to any LLM-assisted changes: review the generated code, verify tool parameters and dependencies, test changes on appropriate data, and document substantive LLM assistance when required by institutional, journal, or funding-agency policies.
