#!/bin/bash
source ~/.bash_profile

#env for prodigal
# conda create -n gene_pred
# conda activate gene_pred
# conda install -c bioconda prodigal
# conda deactivate

# env for bakta
# conda create -n bakta
# conda activate bakta
# conda install -c conda-forge -c bioconda bakta
# bakta_db download --output <output-path> --type full
# conda deactivate 

# env for amrfinder
# conda create -n amrfinder
# conda activate amrfinder
# mamba install -c bioconda ncbi-amrfinderplus
# conda deactivate

#env for pathofact
# git lfs install
# git clone -b master --recursive https://git-r3lab.uni.lu/laura.denies/PathoFact.git
# cd PathoFact
# conda env create -f=envs/PathoFact.yaml
# cd ..

# commandline usage: sh pipeline.sh [input_dir] [output_dir]
input_dir="$1"
output_dir="$2"

for file in "$input_dir"/*; do
	# isolate name saved as $isolate
	isolate=$(basename "$file")

	# create output dirs for fastp, bbduk, skesa and quast
	mkdir -p "$output_dir"/{prodigal, bakta, amrfinder}
	
	#gene prediction with prodigal
	conda activate gen_pred
	prodigal \
	-i "$isolate/filtered_contigs.fa" \
	-o "$output_dir/prodigal/${isolate}_gene.coords.gff" \
	-f gff \
	-a "$output_dir/prodigal/${isolate}_protein.translations.faa" \
	-d "$output_dir/prodigal/${isolate}_gene.predictions.fasta" \
	> "$output_dir/prodigal/${isolate}_logfile.log" 2>&1
	conda deactivate

	#annotation with bakta
	conda activate bakta
	bakta --db ~/bakta/db/ \
	"$isolate/filtered_contigs.fa.gz" \
	--output "$output_dir/bakta/bakta_output" \
	--threads 8 \
	> "$output_dir/bakta/bakta_output.log" \
	2> "$output_dir/bakta/bakta_error.log"
	conda deactivate

	#amr
	conda activate amrfinder
	amrfinder \
	-p "$output_dir/prodigal/${isolate}_gene.predictions.fasta" \
	> "$output_dir/amrfinder/${isolate}_amr"
	conda deactivate
done

#pathofact
conda activate PathoFact
sh ./pathofact.sh
conda deactivate
