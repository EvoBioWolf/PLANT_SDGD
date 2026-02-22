#!/bin/bash

sp=$SPECIES
input_dir="/data/sp11/05_mapping/${sp}"
output_dir="/data/sp11/05_mapping/isoformsID/${sp}"
mkdir -p "$output_dir"

bc_file="/gpfs1/data/sp11/barcodes/$sp.txt"
individuals=($(awk '{print $2}' "$bc_file"))
sample="${individuals[$SLURM_ARRAY_TASK_ID - 1]}"

cd $input_dir || exit
sortedBAM="${input_dir}/${sample}.sorted.bam"
ids="${output_dir}/${sample}.txt"
depth="${output_dir}/${sample}.pos_dep.txt"
depth10="${output_dir}/${sample}.pos_dep10.txt"
depth20="${output_dir}/${sample}.pos_dep20.txt"
depth30="${output_dir}/${sample}.pos_dep30.txt"

# samtools view $sortedBAM | awk '{print $3}' | sort | uniq >$ids
samtools depth $sortedBAM > $depth
# awk '{if($3>10){print $1","$2}}' $depth >$depth10
# awk '{if($3>20){print $1","$2}}' $depth >$depth20
# awk '{if($3>30){print $1","$2}}' $depth >$depth30
