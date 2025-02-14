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

- 
