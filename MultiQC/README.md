# FastQC + FastP + MultiQC Nextflow Pipeline (DSL2)

## Overview

This Nextflow pipeline performs quality control and preprocessing of paired-end FASTQ files using:

* **FastQC**: For raw quality assessment
* **FastP**: For read trimming and filtering
* **MultiQC**: For aggregating FastQC and FastP reports into a unified summary

The pipeline is written in DSL2 and leverages official Docker containers for reproducibility.

---

## Features

* Fully compliant with **Nextflow DSL2** syntax
* Uses **official containers** from BioContainers and Quay.io
* Processes **paired-end FASTQ** files
* Generates:

  * FastQC HTML/ZIP reports
  * FastP HTML/JSON trimming reports
  * Aggregated MultiQC report
* Modular and easy to extend

---

## Usage

### Input

FASTQ files should follow this pattern:

```
data/SAMPLE_1.fastq.gz
 data/SAMPLE_2.fastq.gz
```

### Running the Pipeline

```bash
nextflow run main.nf \
    --reads "data/*{1,2}.fastq.gz" \
    --outdir "results"
```

### Output

* `results/fastqc/`: FastQC reports (`.html`, `.zip`)
* `results/fastp/`: FastP filtered files + QC reports (`.html`, `.json`)
* `results/multiqc/`: Combined `multiqc_report.html` and data folder

---

## Requirements

* Nextflow ≥ 21.10.0
* Docker (or Singularity/Podman, with minor changes)

---

## Containers Used

| Tool    | Container Image                                  |
| ------- | ------------------------------------------------ |
| FastQC  | `biocontainers/fastqc:v0.11.9_cv8`               |
| FastP   | `quay.io/biocontainers/fastp:0.23.2--h5f740d0_0` |
| MultiQC | `biocontainers/multiqc:v1.13--pyhdfd78af_0`      |

---

## Acknowledgements

* [Nextflow](https://www.nextflow.io/)
* [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
* [FastP](https://github.com/OpenGene/fastp)
* [MultiQC](https://multiqc.info/)

