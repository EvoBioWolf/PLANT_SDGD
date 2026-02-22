#!/bin/bash
#SBATCH --mail-user=chongyi.jiang@uni-jena.de
#SBATCH --mail-type=BEGIN,END,FAIL

sp=$(awk "NR==$SLURM_ARRAY_TASK_ID" 01_demultiplex.sp.list)
input_dir="/gpfs1/data/sp11/00_rawData/$sp"
output_dir="/gpfs1/data/sp11/01_demultiplex/$sp"
bc="/gpfs1/data/sp11/barcodes/$sp.txt"

module load GCC/10.2.0
mkdir -p $output_dir
cd $output_dir || exit

/home/jiangc/stack_2.64/bin/process_radtags --paired --inline-null --disable_rad_check -p $input_dir -o $output_dir -b $bc --threads 16 && rm -f *.rem.*

# sbatch -a 1-$(wc -l <01_demultiplex.sp.list) 01_demultiplex.sh
