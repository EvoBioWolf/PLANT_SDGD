#!/bin/bash
########################################################
########## check genotype of replicates   ##############
########################################################

task="replicate_comparison"
runningTime="48:00:00"
memoryCPU="20G"
CPUs="1"

script="05_replicatesTest.sh"

while read -r species ploidy; do
    output_folder="/work/jiangc/SP11/2025/${task}/$species"
    mkdir -p "$output_folder"

    sbatch --job-name=${task} \
           --output="${output_folder}/%x-%A-%a.log" \
           --export=SPECIES="$species",PLOIDY="$ploidy" \
           --mail-user="chongyi.jiang@uni-jena.de" \
           --mail-type="BEGIN,END,FAIL" \
           --time=$runningTime \
           --mem-per-cpu=$memoryCPU \
           --cpus-per-task=$CPUs \
           "$script"

done < 20250905_replicates.list
