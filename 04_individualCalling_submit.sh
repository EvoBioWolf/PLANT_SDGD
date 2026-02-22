#!/bin/bash
########################################################
###### calling variants with GTAK for individual   #####
########################################################

task="gatk_individual"
runningTime="24:00:00"
memoryCPU="10G"
CPUs="1"

script="04_individualCalling.sh"

while read -r species ploidy; do
    output_folder="/work/jiangc/SP11/2025/${task}/$species"
    mkdir -p "$output_folder"

    toCallInd="/data/sp11/script/toCallInds/$species.txt"
    num_samples=$(wc -l < "$toCallInd")

    sbatch --job-name=${task} \
           --output="${output_folder}/%x-%a-%a.log" \
           --export=SPECIES="$species",PLOIDY="$ploidy" \
           --mail-user="chongyi.jiang@uni-jena.de" \
           --mail-type="BEGIN,END,FAIL" \
           --time=$runningTime \
           --mem-per-cpu=$memoryCPU \
           --cpus-per-task=$CPUs \
           --array=1-"$num_samples" "$script"

done < 04_calling.list

