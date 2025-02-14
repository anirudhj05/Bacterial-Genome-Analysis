#!/bin/bash
source ~/.bash_profile

# environment 1
conda create -n teamf -y
conda activate teamf
conda install -c bioconda fastp bbmap skesa -y
conda deactivate

# environment 2
conda create -n pythonold python=2.7 -y
conda activate pythonold
conda install biopython -y
conda deactivate

# environment 3
#create conda env 'qual_eval' with python=3.7.12
conda create -n qual_eval python=3.7.12 -y
conda activate qual_eval
pip install quast matplotlib
conda deactivate

# commandline usage: sh pipeline.sh [input_dir] [output_dir]
input_dir="$1"
output_dir="$2"

for file in "$input_dir"/*R1*; do
	R1=$(basename -- "$file")
	R2="${R1/R1/R2}"
	# isolate name saved as $isolate
	isolate="${R1:0:8}"

	# create output dirs for fastp, bbduk, skesa and quast
	mkdir -p "$output_dir"/{fastp,reports,bbduk,skesa,quast}/$isolate
	
	conda activate teamf
	#trimming with fastp
	fastp -w 6 \
	--qualified_quality_phred 30 \
	--unqualified_percent_limit 20 \
	--dont_eval_duplication \
	-i "$input_dir/${R1}" -I "$input_dir/${R2}"\
	-o "$output_dir/fastp/$isolate/trimmed_${R1}" \
	-O "$output_dir/fastp/$isolate/trimmed_${R2}" \
	-h "$output_dir/reports/${isolate}/${R1}_fastp_report.html" \
	-j "$output_dir/reports/${isolate}/${R2}_fastp.json" -z 6

	#bbduk
	bbduk.sh in="$output_dir/fastp/$isolate/trimmed_${R1}" out="$output_dir/bbduk/$isolate/unmatched_${R1}" outm="$output_dir/bbduk/$isolate/matched_${R1}" ref=phix k=31 hdist=1 overwrite=t stats="$output_dir/bbduk/$isolate/${R1}_stats.txt" 
    	bbduk.sh in="$output_dir/fastp/$isolate/trimmed_${R2}" out="$output_dir/bbduk/$isolate/unmatched_${R2}" outm="$output_dir/bbduk/$isolate/matched_${R2}" ref=phix k=31 hdist=1 overwrite=t stats="$output_dir/bbduk/$isolate/${R2}_stats.txt" 
	 
	#assembly with skesa
	skesa \
	--fastq "$output_dir/bbduk/$isolate/unmatched_${R1}" \
	"$output_dir/bbduk/$isolate/unmatched_${R2}" \
	--contigs_out "$output_dir/skesa/$isolate/$isolate.fasta" \
	1> "$output_dir/skesa/$isolate"/skesa.stdout.txt \
	2> "$output_dir/skesa/$isolate"/skesa.stderr.txt
	
	conda deactivate

	# filter
	
	conda activate pythonold
	python filter.contigs.py \
		-i "$output_dir/skesa/$isolate/$isolate.fasta" \
		-o "$output_dir/skesa/$isolate/filtered_$isolate.fna"
	conda deactivate
	
	#quality evaluation with quast
	conda activate qual_eval
	
	# Run Quast filtered_assemblies
	quast.py -o "$output_dir/quast/$isolate" "$output_dir/skesa/$isolate/filtered_$isolate.fna"

	# Deactivate the 'qual_eval' environment
	conda deactivate
done
