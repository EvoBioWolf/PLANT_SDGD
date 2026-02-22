#! /bin/bash

module load picard/2.25.7-Java-11
sp=$SPECIES
ploidy=$PLOIDY
toCallInd_file="/data/sp11/toCallInds/$sp.txt"
individual=$(awk "NR==$SLURM_ARRAY_TASK_ID" "$toCallInd_file")

input_bam="/data/sp11/05_mapping/${sp}/${individual}.sorted.bam"
output_dir="/data/sp11/08_variantsCalling/inputBAM/${sp}"
mkdir -p $output_dir

RG_bam="${output_dir}/${individual}.AddRG.bam"
sort_bam="${output_dir}/${individual}.sorted.bam"

java -Xmx10G -jar $EBROOTPICARD/picard.jar AddOrReplaceReadGroups -I $input_bam -O $RG_bam -RGID "$sp" -RGLB "$sp" -RGPL "Illumina" -RGPU "NovaSeq 6000" -RGSM "$individual"
java -Xmx10G -jar $EBROOTPICARD/picard.jar SortSam -I $RG_bam -O $sort_bam -SORT_ORDER coordinate && rm -f $RG_bam
samtools index $sort_bam

module load GCCcore/12.2.0 GATK/4.3.0.0-Java-11
ref="/data/sp11/reference_replaced/${sp}/${sp}.fasta"
output_dir="/data/sp11/08_variantsCalling/gvcf/${sp}"
mkdir -p $output_dir
gvcf="${output_dir}/${individual}.gvcf.gz"
gatk --java-options "-Xmx10G" HaplotypeCaller -R "$ref" -I "$sort_bam" -O "$gvcf" -ploidy "$ploidy" -ERC GVCF && rm -f "$sort_bam"


# sbatch -a 1-$(wc -l <13_calling.list) 04_variantCalling.sh
