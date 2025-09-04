#!/bin/bash
#SBATCH --time=01-20:00
#SBATCH --account=def-anguyen
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem=300G
#SBATCH --cpus-per-task=1

#loading modules and libraries
source ~/Tensorflow_mod/bin/activate
module load StdEnv/2023 fastqc/0.12.1 fastp/0.24.0 bbmap/39.06
#/home/p1211536/Tensorflow_mod/bin/python3.9

#Quality control on raw reads
fastqc --noextract --nogroup -o ./fastqc_results -t 8 raw_RNA_seq/*.fastq.gz

#Fastp trimmming and other pre-processing
cat lst_samples.txt | parallel -j4 \ 'fastp -i raw_RNA_seq/{}_R1_001.fastq.gz -I raw_RNA_seq/{}_R2_001.fastq.gz -o fastp_preprocessed/{}_trimmed_1.fastq.gz -O fastp_preprocessed/{}_trimmed_2.fastq.gz -h fastp_preprocessed/fastp_{}.html -j fastp_preprocessed/fastp_{}.json -l 151 --trim_poly_g --trim_poly_a --cut_front --cut_front_mean_quality 20 --cut_tail --cut_tail_mean_quality 20 --detect_adapter_for_pe --overrepresentation_analysis --thread 2'

#Quality control on preprocessed reads
fastqc --noextract --nogroup -o ./after_fastp_fastqc_results -t 8 fastp_preprocessed/*.fastq.gz

#import new modules for the next steps
module load StdEnv/2023 gcc/12.3 subread/2.0.6 star/2.7.11b samtools/1.20

#Build STAR index 
STAR --runMode genomeGenerate --runThreadN 8 --genomeDir S_pombe_star_index --genomeFastaFiles ${REFERENCE_FASTA_DIR}/${REFERENCE_FASTA_FILENAME}.fasta --sjdbGTFfile ${REFERENCE_GTF_DIR}/${REFERENCE_GTF_FILENAME}.gtf --sjdbOverhang 150

#Align reads with STAR
cat lst_samples.txt | parallel -j2 \ ' STAR --genomeDir S_pombe_star_index --runThreadN 4 --readFilesIn ${FASTP_files_path}/${}_trimmed_1.fastq.gz ${FASTP_files_path}/${}_trimmed_2.fastq.gz --readFilesCommand zcat --outFilterMultimapNmax 999 --outFileNamePrefix STAR_BAM/${}_ --outSAMtype BAM SortedByCoordinate --outSAMattributes NH HI NM MD AS --outReadsUnmapped Fastx'

#index output bam files
cat lst_samples.txt | parallel -j8 \ ' samtools index STAR_BAM/${}_Aligned.sortedByCoord.out.bam'

#Count fragments from paired-end mapping
cat lst_samples.txt | parallel -j2 \ 'featureCounts -T 4 -p --countReadPairs -B -C -a ${REFERENCE_GTF_DIR}/${REFERENCE_GTF_FILENAME}.gtf -o read_counts_summary/${}_counts.txt STAR_BAM/${}_Aligned.sortedByCoord.out.bam'