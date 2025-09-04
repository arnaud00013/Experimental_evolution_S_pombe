#!/bin/bash
#SBATCH --time=0-10:00
#SBATCH --account=def-shapiro
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=80G
#SBATCH --cpus-per-task=2
module load r/4.2.2
Rscript analyze_coverage_profile.r
