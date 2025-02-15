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

The required environments are set up as part of the main pipeline script (`sequence_cleaning_and_assembly_pipeline.sh`). Here's what each environment is for and how they are configured:

- ** teamf environment:**

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
- **Pythonold environment:**

This environment is intended to run older Python scripts that work exclusively with Python 2.7. It is critical to assure that older codebases can still be executed without modification.

Tools:
- Biopython: A set of freely available tools for biological computation. It's used here for manipulating biological data, within the custom script `filter.contigs.py`.
The `filter.contigs.py` script processes genomic data by filtering contigs based on length, GC content, coverage, and compositional complexity. Here's how it works:

### **Argument Parsing:**  
The script begins by defining the CLI, allowing the user to specify the input file, output file, and various filtering criteria such as minimum contig length and minimum coverage.

### **File Handling:**  
It handles both plain and gzip-compressed FASTA files as input.

### **Biopython Usage:**  
The script uses the Biopython library to read genomic sequence data from FASTA files. Biopython facilitates the manipulation of genomic data, notably the extraction and analysis of sequence information.

### **Filter Implementation:**
- **Length Filter:** Discards sequences shorter than a specified threshold.
- **Coverage Filter:** Uses regular expressions to extract coverage information from sequence headers and discards sequences with coverage below a specified threshold.
- **GC Content Filter:** Calculates the GC content of sequences using Biopython's utilities and filters out sequences with extreme GC content values.
- **Complexity Filter:** Ensures sequences contain a minimum number of different nucleotides, filtering out low-complexity sequences.

### **Output Generation:**  
Sequences that pass all filters are written to an output file in FASTA format. Optionally, it can also generate a file containing sequences that were discarded during the filtering process, including the reasons for their exclusion in the sequence headers.

```bash
conda create -n pythonold python=2.7 -y
conda activate pythonold
conda install biopython -y
conda deactivate
```

-**qual_eval Environment:**

This environment is used for the quality evaluation of the assembled genomes.  It assesses the quality and integrity of genome assemblies.

Tools:
- Quast: Analyses genome assemblies by computing metrics such as contig count, contig length, and N50 values. It aids in determining assembly quality.
- Matplotlib: A plotting library for Python. It's used within quast or other scripts to generate visualizations of the assembly quality metrics.

```bash
conda create -n qual_eval python=3.7.12 -y
conda activate qual_eval
pip install quast matplotlib
conda deactivate
```
The `sequence_cleaning_and_assembly_pipeline.sh` script runs the workflow by activating each environment only before the necessary tools are needed and deactivates them immediately after use. This ensures that the dependencies of one part of the workflow do not conflict with another, resulting in clean and well-organized process management.

## **Input Data:**  
The pipeline accepts paired-end FastQ files. Ensure your input directory contains files named in the format *R1* and *R2* for paired reads.

## **Running the Script:**  
To run the pipeline, use the following command:

```bash
sh pipeline.sh /path/to/input_dir /path/to/output_dir
```
- **input_dir:** Directory containing input FastQ files.
- **output_dir:** Directory where the pipeline will write its output.

## **Output Files and Structure:**
The pipeline organizes its outputs into several directories, each tailored to a specific processing stage, ensuring that the results are easy to manage and analyze:

- **fastp/**: Contains trimmed sequence files. Each subdirectory within `fastp/` is named after the isolate ID and includes detailed reports on the read trimming process:
  - `*_fastp_report.html`: Provides a visual summary of the trimming statistics, including quality scores and the effect of trimming on read lengths.
  - `*.json`: Offers structured data about the trimming process for potential automated analysis or detailed reviews.

- **bbduk/**: Stores the results of the contaminant filtering process. It includes filtered sequence files alongside statistical summaries:
  - `*_stats.txt`: Statistics files that detail the number and proportion of reads filtered out as contaminants, providing insights into the effectiveness of the filtering process.

- **skesa/**: Contains assembled contigs and their corresponding assembly logs. These are crucial for downstream genomic analyses:
  - `*.fasta`: The assembled genome sequences in FASTA format, representing the primary output from the genome assembly process.
  - `filtered_*.fna`: Filtered assembled sequences that have passed additional quality checks, ensuring that only high-quality genomic data is forwarded for further analysis.

- **quast/**: Holds the quality assessment reports generated by QUAST, which evaluate the quality of the assembled sequences:
  - `report.txt`: A comprehensive text summary of the assembly quality, including metrics such as N50, L50, and total length, which are essential for assessing the completeness and integrity of the assemblies.

# Gene Prediction and Annotation
