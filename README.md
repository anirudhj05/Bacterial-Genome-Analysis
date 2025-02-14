# Bacterial-Genome-Analysis

## Project Overview
This project involved identifying an unknown species/strain of bacteria that caused the outbreak, pinpointing the origin and source of the outbreak, characterizing the functional profile of the virulent isolates (i.e. virulence factors and antimicrobial profile), and making specific recommendations about outbreak response and treatment.

# Sequence Cleaning and Genome Assembly
This pipeline is designed for the processing and quality assessment of genomic data. It involves read trimming, contaminant filtering, genome assembly, and assembly quality evaluation. The pipeline uses multiple environments to handle different parts of the workflow, ensuring the use of appropriate tools for each task.

## Prerequisites
This pipeline is a bash script that uses a Conda/Mamba environment with the appropriate packages and dependencies. It is recommended to  use Mamba over Conda due to its overall better performance in package installation and dependency resolution.

## Setup
1. Clone the repository:
```bash
git clone https://github.com/anirudhj05/Bacterial-Genome-Analysis.git
cd Bacterial-Genome-Analysis
```
2. Environment Setup

The required environments are set up as part of the main pipeline script (pipeline.sh). Here's what each environment is for and how they are configured:

- teamf environment:
This environment is specifically set up to handle the initial stages of genomic data processing which include read trimming and contaminant filtering.

Tools:
  - fastp: Used for high-performance read trimming, quality filtering, and quality control checks. It quickly processes large volumes of data and provides reports in HTML and JSON formats.
  - bbmap: Known for its powerful contaminant filtering capabilities. It can identify and eliminate unfavorable sequences using reference datasets such as phix, which is commonly utilized for identifying sequencing contamination.
  - skesa (Strategic K-mer Extension for Scrupulous Assemblies): It is used for assembling microbial genomes directly from read data.

```bash
conda create -n teamf -y
conda activate teamf
conda install -c bioconda fastp bbmap skesa -y
conda deactivate
```
- Pythonold environment:
This environment is intended to run older Python scripts that work exclusively with Python 2.7. It is critical to assure that older codebases can still be executed without modification.

Tools:
- 
