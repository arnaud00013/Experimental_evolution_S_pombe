#!/bin/bash
#SBATCH --time=00-10:00
#SBATCH --account=def-shapiro
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem=120G
#SBATCH --cpus-per-task=2
module load StdEnv/2020 gcc/9.3.0 fastp/0.23.1 bbmap/38.86 blast+/2.13.0 bwa/0.7.17 samtools/1.17 picard/2.26.3 varscan/2.4.2
chmod u+x /home/p1211536/scratch/S_Pombe/*.py /home/p1211536/scratch/S_Pombe/*.sh
/home/p1211536/Tensorflow_mod/bin/python3.9 run_iPMVC_in_parallel.py /home/p1211536/scratch/S_Pombe/ lst_samples.txt 2
#run the depth analysis
sbatch submit_samples_to_cVC.sh
