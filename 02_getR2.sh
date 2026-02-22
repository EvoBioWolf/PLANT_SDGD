#!/bin/bash

sp=$SPECIES
input_dir1="/data/sp11/03_cdhit_dedup/species/$sp"
input_dir2="/data/sp11/01_demultiplex/species/$sp"
output_dir="/data/sp11/04_cleanR2/species/$sp"
mkdir -p "$output_dir"

bc_file="/gpfs1/data/sp11/barcodes/$sp.txt"
individuals=($(awk '{print $2}' "$bc_file"))
sample="${individuals[$SLURM_ARRAY_TASK_ID - 1]}"

fq2_in="${input_dir1}/${sample}.2.deDup.fq"
fq2_name="${input_dir1}/${sample}.2.name.txt"
grep '^@' $fq2_in | sed 's/^@//' >$fq2_name

fq2="${input_dir2}/${sample}.2.fq.gz"
fq2_out="${output_dir}/${sample}.2.fq.gz"
seqkit grep -f $fq2_name $fq2 -o $fq2_out && rm -f $fq2_name

fq2_clean1="${output_dir}/${sample}.2.clean1.fq.gz"
fq2_clean2="${output_dir}/${sample}.2.clean2.fq.gz"
fq2_clean3="${output_dir}/${sample}.2.cutadapted.fq.gz"

cd "$output_dir" || exit 1
cutadapt -q 15,20 -a "A{30};anywhere;e=0.1;min_overlap=10" -u 1 --minimum-length 60 -j 8 -o $fq2_clean1 $fq2_out && rm -f $fq2_out
cutadapt -q 15,20 -a "A{10};e=0.1" --minimum-length 60 -j 8 -o $fq2_clean2 $fq2_clean1 && rm -f $fq2_clean1
cutadapt -q 15,20 -b "G{20};e=0.1" --minimum-length 60 -j 8 -o $fq2_clean3 $fq2_clean2 && rm -f $fq2_clean2


# sbatch -a 1-$(wc -l <01_demultiplex.sp.list) 02_getR2.sh
