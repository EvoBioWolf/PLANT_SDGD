#! /bin/bash

str=$SPECIES

if [[ $str == *_* ]]; then
  sp="${str%%_*}"   # before _
  ploidy="${str#*_}"   # after _
else
  sp="$str"
  ploidy=""
fi

vcf_file="/data/sp11/08_variantsCalling/gvcf/${str}.cohort.vcf_all_sites.gz"

module load GCC/12.2.0 Pysam/0.21.0
pairs_file="/data/sp11/12_replicates/${sp}/replicate.pair.tsv"
out_file="/data/sp11/12_replicates/${sp}/${str}_replicate.pair.result.tsv"
python 05_replicatesTest.py $pairs_file $vcf_file $out_file

pairs_file="/data/sp11/12_replicates/${sp}/random.pair.tsv"
out_file="/data/sp11/12_replicates/${sp}/${str}_random.pair.result.tsv"
python 05_replicatesTest.py $pairs_file $vcf_file $out_file
