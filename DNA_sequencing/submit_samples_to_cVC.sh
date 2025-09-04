#!/bin/bash
#SBATCH --time=00-10:00
#SBATCH --account=def-shapiro
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=100G
#SBATCH --cpus-per-task=10
chmod u+x /home/p1211536/scratch/S_Pombe/*.py /home/p1211536/scratch/S_Pombe/*.r /home/p1211536/scratch/S_Pombe/*.sh
module load StdEnv/2020 samtools/1.17
ls bam/*preprocessed_sorted_noduplicates.bam | sed 's/bam\///g' > lst_samples_bamfiles.txt
#depth report
samtools depth bam/*_preprocessed_sorted_noduplicates.bam -o coverage_analysis/common_depth_report.csv
