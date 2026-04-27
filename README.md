# Phage Assembly & Annotation Pipeline

Nextflow pipeline for end-to-end processing of Oxford Nanopore phage isolate reads.

## Steps

| # | Process | Tool(s) | Output |
|---|---------|---------|--------|
| 1 | QC & Trimming | fastplong | Trimmed reads + QC report |
| 2 | Assembly | Flye | Draft assembly |
| 3 | Polishing | Medaka | Refined contigs |
| 4 | Assembly QC | CheckV + viralVerify | Completeness / contamination |
| 5 | Annotation | Pharokka → Phold → Phynteny | GFF, PHROGs functions, structure-based calls, synteny context |
| 6 | Report | MultiQC | HTML summary report |

## Requirements

- Nextflow ≥ 23.04
- Micromamba (or conda)

## Setup

```bash
# Create the environment once (Nextflow will do this automatically,
# but you can pre-build it to avoid waiting per-run)
micromamba env create -f environment.yml

# Verify
micromamba run -n phage-pipeline pharokka --version
```

## Usage

```bash
# Test run (no tools required — uses stub outputs)
nextflow run main.nf -c nextflow.config --inputFile work/dummy.fastq -stub

# Real run
nextflow run main.nf -c nextflow.config --inputFile reads.fastq

## Output structure

results/
├── 01_qc/               # fastplong trimmed reads + HTML/JSON report
├── 02_assembly/         # Flye assembly.fasta + assembly_info.txt
├── 02install_databases.py -o <path/to/databse_dir>_medaka/        # viralFlye refined contigs + host prediction
├── 04_assembly_qc/
│   ├── checkv/          # CheckV quality summary
│   └── viralVerify/           # BUSCO completeness report
├── 05_annotation/
│   ├── pharokka/        # PHROGs-based GFF + summary
│   ├── phold/           # Structure-based GFF + summary
│   └── phynteny/        # Synteny-annotated GenBank + summary
└── 06_report/           # MultiQC HTML report
```