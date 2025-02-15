#!/bin/bash
# usage: ./pathofact.sh

# activate git lfs
git lfs install
# clone branch incl. sub-modules
git clone -b master --recursive https://git-r3lab.uni.lu/laura.denies/PathoFact.git
# ONCE DOWNLOADED, MUST INSTALL SIGNALP (
# Adjust config.yaml with proper path names. Pipeline designed to work with these settings
# Sample = <file names>
# Project: PathoFact_results
# datadir:  project_data
# workflow: "Vir"
# signalp: <path to signalp>

cd PathoFact
conda env create -f=envs/PathoFact.yaml

# activate the env
conda activate PathoFact
# run the pipeline
# set --cores to the number of cores to use, e.g. 10
snakemake -s Snakefile --use-conda --reason --cores 4 -p 

# Processes summary of output data in TSV file
process_tsv_file() {
    file="$1"
    sample_id=$(basename "$file" | cut -f1 -d '.')
    total_entries=$(($(wc -l < "$file") - 1))
    pathogen_hmm=$(awk -F'\t' '$3 == "pathogenic" {count++} END {print count}' "$file")
    pathogen_class=$(awk -F'\t' '$4 == "pathogenic" {count++} END {print count}' "$file")
    uniq_hmm=$(awk -F'\t' '$3 == "pathogenic" && $4 == "negative" {count++} END {print count}' "$file")
    uniq_class=$(awk -F'\t' '$4 == "pathogenic" && $3 == "negative" {count++} END {print count}' "$file")
    final_pred_unclass=$(awk -F'\t' '$5 == "unclassified" {count++} END {print count}' "$file")
    final_pred_nonpath=$(awk -F'\t' '$5 == "non_pathogenic" {count++} END {print count}' "$file")
    final_pred_virulent=$(awk -F'\t' '$5 == "pathogenic" {count++} END {print count}' "$file")
	
    echo "$sample_id,$total_entries,$pathogen_hmm,$pathogen_class,$uniq_hmm,$uniq_class,$final_pred_unclass,$final_pred_nonpath,$final_pred_virulent"
}
# Processes classification count summary
count_categories() {
    file="$1"
    sample_id=$(basename "$file" | cut -f1 -d '.')
    count_negative=$(awk -F'\t' '$7 == "-" {count++} END {print count}' "$file")
    count_1=$(awk -F'\t' '$7 ~ /1: / {count++} END {print count}' "$file")
    count_2=$(awk -F'\t' '$7 ~ /2: / {count++} END {print count}' "$file")
    count_3=$(awk -F'\t' '$7 ~ /3: / {count++} END {print count}' "$file")
    count_4=$(awk -F'\t' '$7 ~ /4: / {count++} END {print count}' "$file")
    
    echo "$sample_id,$count_negative,$count_1,$count_2,$count_3,$count_4"
}
# Main function
main() {
    input_dir="./PathoFact/project_data/PathoFact_results/PathoFact_report"
    output_file="PathoFactSummary.csv"
    output_file2="category_counts.csv"
	echo "SAMPLE ID,TotalEntries,PathogenicHMM,PathogenicClassifier,UniqPathogenic_HMM,UniqPathogenic_Classifier,FinalPredUnclass,FinalPredNonPath,FinalPredVirulent" > "$output_file"
    for file in "$input_dir"/*.tsv; do
		if [ -f "$file" ]; then
            process_tsv_file "$file" >> "$output_file"
        fi
    done
	
	echo "SampleID,-,1,2,3,4" > "$output_file2"
    for file in "$input_dir"/*.tsv; do
        if [ -f "$file" ]; then
            count_categories "$file" >> "$output_file"
        fi
    done
}
main