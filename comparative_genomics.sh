#!/bin/bash

# Assuming Anaconda is already installed.
# Version used: 2022.05

# Base directory definition:
base_dir="/home/USER/biol7210/group4_comp_genomics"
gsearch_dir="$base_dir/gsearch_results"
gsearch_db_dir="$base_dir/GTDB/nucl/"
prokka_dir="$base_dir/prokka"
amrfinder_output="$base_dir/amrfinder_output"
amrfinder_output_filtered="$base_dir/amrfinder_output_filtered"
assemblies_dir="$base_dir/assemblies"
parsnp_outdir="$base_dir/parsnp_outdir"
filtered_skesa_asm_dir="$base_dir/filtered_skesa_asm"


# Create the necessary directories to store the outputs and our intermediate files as we process them.
mkdir -p $prokka_dir $amrfinder_output $amrfinder_output_filtered $assemblies_dir $parsnp_outdir $filtered_skesa_asm_dir/outbreak_samples/{gff,faa,fasta}

##---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------##

# Step 1: Setup the necessary Conda environments and install the required packages for different tools used in this pipeline. 

##---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------##

	## GSearch setup:
		### GSearch is used for genomic db searching using MinHash techniques to find genomic similarities more efficiently.
conda create -n gsearch -y
conda activate gsearch
conda install python=3.8
conda install -c bioconda gsearch
conda deactivate

	## ParSNP setup:
		### ParSNP is used for rapid core genome multi-alignment & making phylogenetic trees. In this case, it is being used to understand the evolutionary relationships between strains.
conda create -n harvestsuite -c bioconda parsnp harvesttools figtree -y
conda deactivate

	## Roary setup:
		### Roary is used for pangenome analysis. It allows us to compare whole genomes and identify shared / unique genes among bacterial isolates. 
conda create -n roary -y
conda activate roary
conda install -c bioconda roary
conda deactivate

	## AMRFinder setup:
		### AMRFinder is generally used to identify antimicrobial resistance genes in protein sequences. This is key in understanding the resistance profiles of bacterial isolates.
conda create -n amrfinder -y
conda activate amrfinder
conda install -c bioconda amrfinderplus
conda deactivate

	## Additional Python environment:
conda create python38 python=3.8 -y
conda activate python38
conda install -c bioconda gsearch pandas matplotlib biopython
conda deactivate

##---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------##

	# Intermediate step --> Prepare for ParSNP execution: in order to manage the directories / file preparations.

##---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------##

mkdir -p $assemblies_dir/new  # Create a new directory for assemblies if it does not exist
# Loop through each folder to prepare assemblies
for folder in $filtered_skesa_asm_dir/*/; do
    cp "${folder}filtered_contigs.fa" "$assemblies_dir/new/$(basename "$folder")_contigs.fa"
done

##---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------##

# Step 2: Run respective tools for comparative genome analysis.

##---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------##

	## GSearch: Rapid genomic comparisons to support further study on genomic epidemiology / evolutionary bio (e.g. tracking genetic variations across populations as pathogens). 
conda activate gsearch
gsearch --db $gsearch_db_dir --query $assemblies_dir/sample_genome.fna --output $gsearch_dir/results.txt
conda deactivate

	## ParSNP: Analyze the core genome / construct a phylogenetic tree to better visualize the relationships & evolutionary history of the given isolates. 
conda activate harvestsuite
parsnp -d $assemblies_dir/new -r $assemblies_dir/new/reference_genome.fna -o $parsnp_outdir -p 4
figtree $parsnp_outdir/parsnp.tree
conda deactivate

	## Roary: Perform pangenome analysis to find common / unique genes among the genomes.
conda activate roary
roary -f $prokka_dir/roary_results $prokka_dir/*.gff
roary -e --mafft -p 8 $prokka_dir/*.gff
query_pan_genome -a union $prokka_dir/*.gff
query_pan_genome -a intersection $prokka_dir/*.gff
query_pan_genome -a complement $prokka_dir/*.gff
cd $prokka_dir
query_pan_genome -a difference --input_set_one F1738908.gff,F1739691.gff,F2045925.gff,F0582884.gff --input_set_two F2049583.gff,F2049584.gff,F2052228.gff,F0784744.gff
conda deactivate

	## AMRFinder: Look through protein sequences to identify resistance genes so that we can better understand how to select effective treatments & potential resistance mechanisms. 
conda activate amrfinder
for predicted_genes in $base_dir/outbreak_samples/faa/*.faa; do
    base_name=$(basename "$predicted_genes" .faa)
    amrfinder -p "$predicted_genes" -o "${amrfinder_output}/${base_name}_amrfinder_output.txt"
done
conda deactivate

##---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------##