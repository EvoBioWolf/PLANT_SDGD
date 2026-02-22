#!/bin/bash

task="overlapAmongIndividuals"
script="03_countOverlapped.sh"
runningTime="4:00:00"
memoryCPU="100G"
CPUs="1" 

while read -r species; do

    output_folder="/work/jiangc/SP11/2025/${task}/$species"
    mkdir -p $output_folder

    sbatch --job-name=${task} \
           --output="${output_folder}/%x-%a-%a.log" \
           --export=SPECIES="$species" \
           --time=$runningTime \
           --mail-user="chongyi.jiang@uni-jena.de" \
           --mail-type="BEGIN,END,FAIL" \
           --mem-per-cpu=$memoryCPU \
           --cpus-per-task=$CPUs \
           "$script"

done < 01_demultiplex.sp.list



