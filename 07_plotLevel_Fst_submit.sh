#!/bin/bash
########################################################
          ###### calculating Fst   #####
########################################################

task="Fst"
runningTime="36:00:00"
memoryCPU="30G"
CPUs="1"

script="07_plotLevel_Fst.sh"

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

done < 07_plotLevel_Fst.list
