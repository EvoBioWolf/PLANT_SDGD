#! /bin/bash

#SBATCH --job-name=fitModel
#SBATCH --output=/work/jiangc/modelFit/%x-%A-%a.log
#SBATCH --mail-user=chongyi.jiang@uni-jena.de
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --time=12:00:00
#SBATCH --mem-per-cpu=40G
#SBATCH --ntasks-per-node=1

module load GCC/12.2.0 OpenMPI/4.1.4 R/4.2.2
str=$(awk "NR==$SLURM_ARRAY_TASK_ID" 20251208_fitModel.list)

cd "/data/sp11/13_cohortResult/table/"
input_file="${str}.all.tsv"
output_summary_file="${str}.model_BiomassResid_withinDiv.summary"
output_RData="${str}.model_BiomassResid_withinDiv.RData"
Rscript 06_withinSp_heterozygosity_fitModel.R $input_file $output_summary_file $output_RData


# sbatch -a 1-$(wc -l <20251208_fitModel.list) 06_withinSp_heterozygosity_fitModel.sh
