#!/bin/bash

sp=$SPECIES
input_dir="/data/sp11/04_cleanR2/species/$sp"
output_dir="/data/sp11/05_mapping/$sp"
mkdir -p "$output_dir"

bc_file="/gpfs1/data/sp11/barcodes/$sp.txt"
individuals=($(awk '{print $2}' "$bc_file"))
sample="${individuals[$SLURM_ARRAY_TASK_ID - 1]}"

ref="/data/sp11/reference_replaced/${sp}/${sp}.fasta"
fq2_in="${input_dir}/${sample}.2.cutadapted.fq.gz"
sorted="${output_dir}/${sample}.sorted.bam"
bwa mem -t 20 "$ref" "$fq2_in" | samtools view -b | samtools sort -o "$sorted"
samtools index $sorted
