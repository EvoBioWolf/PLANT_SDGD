#!/bin/bash

##########################################################################
### Using cd-hit-dup to remove duplication; 16 (UMI) + 50 (R2) = 66bp ####
### all 66bp identical then consider it as duplication                ####
##########################################################################

task="cleanR2"
script="02_getR2.sh"
runningTime="4:00:00"
memoryCPU="4G"
CPUs="8" 

while read -r species; do

    bc_file="/gpfs1/data/sp11/barcodes/$species.txt"
    output_folder="/work/jiangc/SP11/2025/${task}/$species"
    mkdir -p $output_folder

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


