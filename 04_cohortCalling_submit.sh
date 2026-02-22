#!/bin/bash
##########################################
###### do cohort calling with GATK   #####
##########################################

task="gatk_cohort"
runningTime="480:00:00"
memoryCPU="620G"
CPUs="1"

script="04_cohortCalling.sh"

while read -r species ploidy; do
    output_folder="/work/jiangc/SP11/2025/${task}/$species"
    mkdir -p "$output_folder"

    sbatch --job-name=${task} \
           --output="${output_folder}/%x-%j-%a.log" \
           --export=SPECIES="$species" \
           --mail-user="chongyi.jiang@uni-jena.de" \
           --mail-type="BEGIN,END,FAIL" \
           --time=$runningTime \
           --mem-per-cpu=$memoryCPU \
           --cpus-per-task=$CPUs \
           "$script"

done < 13_calling.list
