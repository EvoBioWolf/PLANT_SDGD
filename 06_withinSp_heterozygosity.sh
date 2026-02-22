#! /bin/bash

sp=$SPECIES
input_vcf="/data/sp11/13_cohortResult/vcf/${sp}.cohort.vcf_all_sites.gz"
diversity_file="plot_div.txt"
expr_file="/data/sp11/10_expression_PacBio/${sp}.count"
filter_file="/data/sp11/12_replicates/${sp}/replicate.pair.tsv"
output_dir="/data/sp11/13_cohortResult/table"

mkdir -p $output_dir
output_table="${output_dir}/${sp}.all.tsv"
python 06_withinSp_heterozygosity.py $input_vcf $diversity_file $expr_file $filter_file $sp $output_table

# vcf_file, diversity_file, expr_file, filter_file, species_name, output_file = sys.argv[1:]
