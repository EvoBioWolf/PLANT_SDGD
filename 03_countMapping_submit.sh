#!/bin/bash

############################################################
#### count how many reads are mapped to the reference #####
############################################################

task="countMapping"
script="03_countMapping.sh"
runningTime="2:00:00"
memoryCPU="10G"
CPUs="1" 

while read -r species; do

    output_folder="/work/jiangc/SP11/2025/${task}/$species"
    mkdir -p $output_folder

    bc_file="/gpfs1/data/sp11/barcodes/$species.txt"
    num_samples=$(wc -l < "$bc_file")

    # Submit an array job for individuals within the species
    sbatch --job-name=${task} \
           --output="${output_folder}/%x-%a-%a.log" \
           --export=SPECIES="$species" \
           --time=$runningTime \
           --mail-user="chongyi.jiang@uni-jena.de" \
           --mail-type="BEGIN,END,FAIL" \
           --mem-per-cpu=$memoryCPU \
           --cpus-per-task=$CPUs \
           --array=1-"$num_samples" "$script"

done < 01_demultiplex.sp.list


