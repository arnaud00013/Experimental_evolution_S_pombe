# Experimental evolution in *S. pombe*
DNA and RNA-sequencing of 10 populations of the fission yeast evolved for 10,000 generations in the same conditions as a pre-existing budding yeast dataset, i.e. high-sugar media and hypoxic conditions, allowing us to observe repeatable evolutionary outcomes within species but diverse molecular mechanisms and targets of adaptation across species.

## DNA sequencing
We sequenced the DNA of each population every 1,000 generations. We successfully recovered 93 out of 100 possible timepoint DNA sequences, which were stored in NCBI SRA and analyzed with the scripts in the *DNA_sequencing* repository. The pipeline for the DNA analysis is customized for our Slurm environment (usernames, configuration, file paths, dependencies and resources need to be customized to your own computing environment):
- sbatch submit_samples_to_iPMVC.sh 
  - Dependencies: StdEnv/2020, gcc/9.3.0, fastp/0.23.1, bbmap/38.86, blast+/2.13.0, bwa/0.7.17, samtools/1.17, picard/2.26.3, varscan/2.4.2, r/4.2.2

-Rscript S_Pombe_evolution_analysis.R
  - Dependencies: r/4.2.2 and packages "ggplot2", "seqinr", "RColorBrewer", "randomcoloR", "FD", "vegan", "gplots", "lmPerm", "ggpubr", "gridExtra", and "tidyr"

-The reference fasta sequences for S. pombe and possible sources of contamination (other yeast species sequenced in parallel)

-The list of samples

## RNA sequencing
We also sequenced the RNA of the WT (time 0) and the evolved populations after 10,000 generations of evolution, which we stored in NCBI SRA and analyzed using the scripts in the *RNA_sequencing* repository. The pipeline for the RNA-seq analysis is customized for our Slurm environment (usernames, configuration, file paths, dependencies and resources need to be customized to your own computing environment):
- sbatch submit_RNA_seq_counts.sh 
  - Dependencies: StdEnv/2023, fastqc/0.12.1, fastp/0.24.0, bbmap/39.06, subread/2.0.6, star/2.7.11b, samtools/1.20 and r/4.2.2

-Rscript H2O2_binary_sensitivity_clusters_Differential_expression_analysis.r
  - Dependencies: r/4.2.2 and packages "DESeq2", "ggplot2", "dplyr" and "ggrepel"

-The reference fasta sequences and gtf files

-The list of samples
