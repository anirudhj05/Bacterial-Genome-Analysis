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

- **teamf environment:**

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

1. **Argument Parsing:**  
The script begins by defining the CLI, allowing the user to specify the input file, output file, and various filtering criteria such as minimum contig length and minimum coverage.

2. **File Handling:**  
It handles both plain and gzip-compressed FASTA files as input.

3. **Biopython Usage:**  
The script uses the Biopython library to read genomic sequence data from FASTA files. Biopython facilitates the manipulation of genomic data, notably the extraction and analysis of sequence information.

4. **Filter Implementation:**
- **Length Filter:** Discards sequences shorter than a specified threshold.
- **Coverage Filter:** Uses regular expressions to extract coverage information from sequence headers and discards sequences with coverage below a specified threshold.
- **GC Content Filter:** Calculates the GC content of sequences using Biopython's utilities and filters out sequences with extreme GC content values.
- **Complexity Filter:** Ensures sequences contain a minimum number of different nucleotides, filtering out low-complexity sequences.

5. **Output Generation:**  
Sequences that pass all filters are written to an output file in FASTA format. Optionally, it can also generate a file containing sequences that were discarded during the filtering process, including the reasons for their exclusion in the sequence headers.

```bash
conda create -n pythonold python=2.7 -y
conda activate pythonold
conda install biopython -y
conda deactivate
```

- **qual_eval Environment:**

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
sh sequence_cleaning_and_assembly_pipeline.sh /path/to/input_dir /path/to/output_dir
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
This pipeline is designed for gene prediction and annotation in microbial genomes. It predicts genes, annotates genomic features, identifies antimicrobial resistance genes, and analyzes pathogenic potentials.

## Environment Setup

Before running the pipeline, ensure that the necessary environments and tools are set up:

- Prodigal for gene prediction: Fast and reliable microbial gene identification.
- Bakta for rapid annotation: Accurate and up-to-date microbial genome annotation using a curated database.
- AMRFinderPlus for AMR gene detection: Identifies AMR genes using the latest NCBI AMR gene database.
- PathoFact for pathogenomic analysis: Determines the pathogenic potential by integrating various analyses including virulence factors and resistance genes.

Environment for Prodigal:
```bash
#Prodigal
conda create -n gene_pred
conda activate gene_pred
conda install -c bioconda prodigal
conda deactivate
```

Environment for Bakta
```bash
#Bakta
conda create -n bakta
conda activate bakta
conda install -c conda-forge -c bioconda bakta
bakta_db download --output <output-path> --type full
conda deactivate
```

Environment for AMRfinder:
```bash
# env for amrfinder
conda create -n amrfinder
conda activate amrfinder
mamba install -c bioconda ncbi-amrfinderplus
conda deactivate
```

Environment for Pathofact:
```bash
#env for pathofact
git lfs install
git clone -b master --recursive https://git-r3lab.uni.lu/laura.denies/PathoFact.git
cd PathoFact
conda env create -f=envs/PathoFact.yaml
cd ..
```

### Note
Bakta requires assembly files in zipped fasta format. You can use the option of gunzip [Option] [archive name/file name] to zip the files.

## Running the Pipeline

This pipeline takes the assembled fasta files as input from the sequence cleaning and genome pipeline and those are used by Prodigal, Bakta and PathOFact. The prodigal generates the faa and gff files that are used by AMR for determining the genes responsible for pathogenicity.

To run the annotation pipeline, use the following command:
```bash
sh Gene_Prediction_and_Annotation.sh [input_dir] [output_dir]
```

- **`input_dir`:** Directory containing the input genome sequences.
- **`output_dir`:** Directory where the pipeline will write its outputs.

## Output Files

`prodigal/`- This directory contains the following output content:
- **GFF**: The file ([isolate]_gene.coords.gff) contains gene predictions in GFF format, including start and end positions, strand direction, and other information about each predicted gene.
- **FAA**: The file ([isolate]_protein.translations.faa) contains the amino acid sequences of the proteins encoded by the predicted genes, which can be used to further predict and analyze protein function.
- **FASTA**: The file ([isolate]_gene.predictions.fasta) contains the predicted genes' nucleotide sequences, which can be used for subsequent genomic study such as gene function studies or comparative genomics.
- **Log File**: A detailed log ([isolate]_logfile.log) containing the execution data as well as any potential warnings or errors discovered by Prodigal during the gene prediction process.

`bakta/` - This directory contains the following output content:
- Bakta provides detailed annotation files containing information on identified genes, their functions, and classifications. Annotations are frequently stored in formats such as **GenBank** and **JSON**, which work well with a variety of bioinformatics tools for further analysis.
- `bakta_output.log`: This log file contains detailed information about the annotation process, including any errors or warnings that occurred and statistics about the genomic features annotated.
- `bakta_error.log`: Specifically captures any errors encountered during the annotation process, helping in troubleshooting and ensuring the quality of the annotation.

`amrfinder/` - This directory contains the following output content:
- **AMR Gene List**: This file ([isolate]_amr) contains a list of all antimicrobial resistance genes found in the genome. The output comprises the gene's name, location on the genome, and the resistance mechanism or family it belongs to. This information is critical for understanding the organism's resistance profile and can help with  epidemiological studies. This output is typically in a plain text format

`PathoFact/` -  This directory contains the following output content:
PathoFact analyzes genomic data for virulence factors, antimicrobial resistance genes, and other pathogenic potentials. 
- **Virulence Factors**: Lists genes and elements identified as virulence factors, including their mechanisms and potential effects on host-pathogen interactions.
- **Resistance Genes**: Similar to AMRFinderPlus but integrated within the broader context of pathogenic potential.
- **Metabolic Capabilities**: Analyzes and reports on the metabolic pathways present in the organism, which can influence pathogenicity.
These outputs include files in TSV format.
- **Statistical Summary Files**: Typically in CSV format, providing an overview with numerical data about the detected pathogenic features.
- **Classification Summarie**s: Also in CSV format, these files break down the number of genes in different pathogenic categories, which aids in the quick assessment and comparison across samples.

# Genotyping and Taxonomic and Quality Assessments

This pipeline combines several bioinformatics tools to perform genotyping, taxonomic categorization, and quality assessments on microbial genomes. It uses high-throughput sequencing data to predict and identify microbial genes, as well as analyze antimicrobial resistance and genome assembly quality.

## Environment Setup
The required environments are set up as part of the main pipeline script (`Genotyping_taxonomic_quality_assessment.sh`). Here's what each environment is for and how they are configured:
## Genotyping

- **MLST Environment:**
This environment uses MLST (Multi-Locus Sequence Typing) for characterizing bacterial species based on sequences of internal housekeeping genes. It assigns allelic profiles to these genes and combines the profiles to define an isolate's sequence type. The tool extracts specific DNA fragments from the genome, sequences them, and compares the sequences to known alleles in a database.

- **MAAST Environment:**
MAAST (Metagenomic Analysis Software Tool) is intended for microbial genome analysis and typing. It integrates functions such as sequence alignment, phylogenetic analysis, and variant calling to provide information about microbial populations' genetic links and functional capacities.  The tool typically requires compilation from source code and is executed to analyze genome sequences to identify genetic types based on predefined criteria.

## Taxonomic Analysis
- **LisSero Environment:**
LisSero is specifically designed for rapid serotyping of Listeria species based on whole genome sequencing data. It uses known serotype-specific gene sequences and performs rapid comparisons to assign a serotype to Listeria isolates.

## Genome Quality Assessment
- **CheckM Environment:**
CheckM includes tools for evaluating the quality of microbial genomes retrieved from isolates, metagenomes, and single cells. It aids in determining the completeness and contamination of genome assemblies. It uses lineage-specific marker genes to evaluate the completeness and assess potential contamination of genome assemblies.

- **QUAST Environment:**
QUAST (Quality Assessment Tool for Genome Assemblies) evaluates genome assemblies by analyzing various metrics that reflect the quality of the assembly. It compares genome assemblies to reference sequences if available, and calculates statistical metrics such as N50, L50, and total length, among others.

## Running the Pipeline
- **Input Data:** Directory containing FASTA files (.fa or .fna)

- **Running the Pipeline:**
```bash
sh Genotyping_taxonomic_quality_assessment.sh [input_directory] [output_directory]
```
### Note: Ensure to provide the absolute path for input_directory. 

## Output Files:
- **MLST/**: Summary files in `.tsv` format detailing the MLST profiles.
- **MAAST/**: Detailed genotype data and phylogenetic trees in respective directories.
- **LisSero/**: Aggregated serotyping results are in `.txt` format.
- **CheckM/**: Outputs include detailed text and `.tsv` reports on genome quality, highlighting potential issues and overall assembly quality.
- **QUAST/**: Provides comprehensive assembly evaluation reports in HTML and text formats.

# Comparative Genomics

This pipeline is designed for comprehensive comparative genomic analysis including genomic similarity assessment, phylogenetic tree construction, pangenome analysis, and antimicrobial resistance (AMR) gene detection. Comparative genomics is an useful scientific approach that includes comparing the genetic information from different organisms to better comprehend their evolutionary links, find variations in genes, and reveal functional capabilities. In the context of pathogen outbreak research, comparative genomics can disclose how distinct strains differ from each other and from known reference genomes, which is critical for pinpointing the outbreak source and understanding the disease's evolution.

## Environment Setup
The required environments are set up as part of the main pipeline script (`comparative_genomics.sh`). Here's what each environment is for and how they are configured:

- **gsearch**: This tool is used to perform advanced genomic searches and comparisons using MinHash-like signatures and Hnsw (Hierarchical Navigable Small World) graphs.
- **parsnp**: This tool is used for phylogenetic reconstruction and identifying core genome single nucleotide polymorphisms (SNPs).
- **roary**: This tool is used for pangenome analysis, identifying the core and accessory genes across multiple genomes.
- **amrfinder**: This tool is used for identifying AMR genes, providing crucial data for understanding resistance mechanisms and guiding treatment options.
## Running the Pipeline
Before running the pipeline, update the base_dir variable in the script to reflect your local directory path.
```bash
base_dir="/home/USER/biol7210/group4_comp_genomics"
mkdir -p $base_dir/{gsearch_results,prokka,amrfinder_output,amrfinder_output_filtered,assemblies,parsnp_outdir,filtered_skesa_asm/outbreak_samples/{gff,faa,fasta}}
```
Replace /home/USER/biol7210/group4_comp_genomics with your desired base directory path.

To run the pipeline, use the following command:
```bash
./comparative_genomics.sh
```
- **Input data:**
  - FASTA files (.fa or .fna) stored in a  designated directory:  `$filtered_skesa_asm_dir`.
  - FASTA format for amino acids (.faa) used by AMRfinder to identify AMR genes Usually found within a directory structured, like `$base_dir/outbreak_samples/faa/`.
 
## Output files
- `gsearch_results/results.txt` - Contains the results of genomic database searches conducted using MinHash algorithms, which efficiently estimate the similarity across genomic datasets. This output file includes similarity scores and identities for genomes with substantial matches, which might help identify closely related strains or possible outbreaks.
- `$parsnp_outdir` -  Includes multiple files that together describe the core genome alignment and the resultant phylogenetic tree:
    - Phylogenetic Tree File: Typically in Newick format, which can be viewed using tools like FigTree. This tree helps visualize the evolutionary relationships between the analyzed genomes.
    - Alignment Files: Contain the aligned sequences of the core genomes, which are critical for identifying evolutionary conserved regions.
    - Summary Files: Provide a concise overview of the alignment metrics, such as coverage of the core genome and any notable genomic variations.
- `$prokka_dir/roary_results` - The results from Roary's pangenome analysis include:
    - Gene Presence and Absence Matrix: A comprehensive table that lists which genes are present in which genomes, providing insights into the core and accessory genes across the studied strains.
    - Graphs and Charts: Visual representations of the pangenome dynamics, such as the number of new genes discovered with each additional genome sequenced.
    - Detailed Statistics: Text files and spreadsheets detailing the number of core, soft-core, shell, and cloud genes, which help in understanding the genetic makeup of the microbial community or species under study.
- `$amrfinder_output` - Lists the genes associated with antimicrobial resistance, their possible functions, and the mechanisms by which they confer resistance and includes metadata about the conditions under which each gene confers resistance, which is crucial for clinical settings and resistance management.   
   
