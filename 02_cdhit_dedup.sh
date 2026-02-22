#!/bin/bash

sp=$SPECIES
bc_file="/gpfs1/data/sp11/barcodes/$sp.txt"
individuals=($(awk '{print $2}' "$bc_file"))

# Get the current individual based on SLURM_ARRAY_TASK_ID
sample="${individuals[$SLURM_ARRAY_TASK_ID - 1]}"
input_dir="/data/sp11/02_extractUMI/species/$sp"
output_dir="/data/sp11/03_cdhit_dedup/species/$sp"
mkdir -p "$output_dir"

fq1_in="${input_dir}/${sample}.1.fq.gz"
fq2_in="${input_dir}/${sample}.2.fq.gz"
fq1_="${output_dir}/${sample}.1.fq"
fq2_="${output_dir}/${sample}.2.fq"

fq1_out="${output_dir}/${sample}.1.deDup.fq"
fq2_out="${output_dir}/${sample}.2.deDup.fq"
# fa1="${input_dir}/${sample}.1.fa"
# fa1_out="${input_dir}/${sample}.1.removeUMI.fa"

gunzip -c $fq1_in > $fq1_
gunzip -c $fq2_in > $fq2_

cd "$output_dir" || exit 1
# sed -n '1~4s/^@/>/p;2~4p' $fq1_ > $fa1
# cd-hit -i $fa1 -o $fa1_out -c 1
/home/jiangc/cd-hit-v4.8.1-2019-0228/cd-hit-auxtools/cd-hit-dup -i $fq1_ -i2 $fq2_ -o $fq1_out -o2 $fq2_out
rm -f $fq1_ $fq2_



# sbatch -a 1-$(wc -l <01_demultiplex.sp.list) 02_cdhit_dedup.sh
