#!/bin/bash
########################################################################
###### prepare tables to calculate within species heterozygosity   #####
########################################################################

task="vcfToTable"
runningTime="24:00:00"
memoryCPU="40G"
CPUs="1"

script="06_withinSp_heterozygosity.sh"

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

done < 20250905_VCF_to_table1.list
