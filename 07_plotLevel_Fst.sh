#! /bin/bash

str=$SPECIES
# GalMol_2x
if [[ $str == *_* ]]; then
  sp="${str%%_*}"   # before _
  ploidy="${str#*_}"   # after _
else
  sp="$str"
  ploidy=""
fi

input_vcf="/data/sp11/13_cohortResult/vcf/${str}.cohort.vcf_all_sites.gz"
diversity_file="plot_div.txt"
filter_file="/data/sp11/12_replicates/${sp}/replicate.pair.tsv"
output_dir="/data/sp11/15_FstResult/"

mkdir -p $output_dir
output_table="${output_dir}/${str}.Fst.csv"

module load GCC/12.2.0 Pysam/0.21.0
python 07_plotLevel_Fst.py $input_vcf $filter_file $diversity_file $output_table
