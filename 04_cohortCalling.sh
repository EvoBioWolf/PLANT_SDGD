#!/bin/bash

sp=$SPECIES
module load GCCcore/12.2.0 GATK/4.3.0.0-Java-11
gvcf_dir="/data/sp11/08_variantsCalling/gvcf/${sp}"
gvcf_list="/data/sp11/08_variantsCalling/gvcf/${sp}.list"
cohort_vcf="${sp}.cohort.g.vcf.gz"
final_out="${sp}.cohort.vcf_all_sites.gz"
ref="/data/sp11/reference_replaced/${sp}/${sp}.fasta"

cd $gvcf_dir || exit

gatk --java-options "-Xmx610g" CombineGVCFs \
  -R "$ref" \
  $(awk '{print "-V", $1}' "$gvcf_list") \
  -O "$cohort_vcf"

# Genotype GVCFs
gatk --java-options "-Xmx610g" GenotypeGVCFs \
  -R "$ref" \
  -V "$cohort_vcf" \
  --include-non-variant-sites true \
  -O "$final_out"

