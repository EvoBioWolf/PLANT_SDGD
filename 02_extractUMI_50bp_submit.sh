#!/bin/bash
##########################################################################
###### Using UMI (first 16bp) in R1 and first 50bp in R2 to remove   #####
###### duplication via cdhit;                                        #####
##########################################################################

task="extractUMI"
runningTime="24:00:00"
memoryCPU="20G"
CPUs="1" 

script="02_extractUMI_50bp.sh"

while read -r species; do
    output_folder="/work/jiangc/SP11/2025/${task}/$species"
    mkdir -p $output_folder

    bc_file="/gpfs1/data/sp11/barcodes/$species.txt"
    num_samples=$(wc -l < "$bc_file")

    # Submit an array job for individuals within the species
    sbatch --job-name=${task} \
           --output="${output_folder}/%x-%a-%a.log" \
           --export=SPECIES="$species" \
           --mail-user="chongyi.jiang@uni-jena.de" \
           --mail-type="BEGIN,END,FAIL" \
           --time=$runningTime \
           --mem-per-cpu=$memoryCPU \
           --cpus-per-task=$CPUs \
           --array=1-"$num_samples" "$script"

done < 01_demultiplex.sp.list


