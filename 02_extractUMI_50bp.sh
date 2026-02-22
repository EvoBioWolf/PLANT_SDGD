#!/bin/bash

sp=$SPECIES
bc_file="/gpfs1/data/sp11/barcodes/$sp.txt"
individuals=($(awk '{print $2}' "$bc_file"))

# Get the current individual based on SLURM_ARRAY_TASK_ID
sample="${individuals[$SLURM_ARRAY_TASK_ID - 1]}"

input_dir="/data/sp11/01_demultiplex/species/$sp"
output_dir="/data/sp11/02_extractUMI/species/$sp"

# Ensure the output directory exists
mkdir -p "$output_dir"

fq1_in="${input_dir}/${sample}.1.fq.gz"
fq2_in="${input_dir}/${sample}.2.fq.gz"
fq1_out="${output_dir}/${sample}.1.fq.gz"
fq2_out="${output_dir}/${sample}.2.fq.gz"

cd "$output_dir" || exit 1

# Run the Python script for UMI extraction
python 02_extractUMI_50bp.py "$fq1_in" "$fq2_in" "$fq1_out" "$fq2_out"


# sbatch -a 1-$(wc -l <01_demultiplex.sp.list) 02_extractUMI_50bp.sh
