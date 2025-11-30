#import packages
library( "DESeq2" )
library("ggplot2")
library("ggpubr")
library("dplyr")
library("ggrepel")

#import arguments
output_workspace <- "C:/Users/arnau/Documents/S_Pombe/"
lst_sample_replicates <- read.csv2(file = paste0(output_workspace,"lst_replicate_samples.txt"),sep = "\t",header = FALSE,stringsAsFactors = FALSE)$V1
df_fragment_pairs_counts_ALL_replicates <- read.csv(paste0(output_workspace,"counts_ALL_SAMPLE_REPLICATES.csv"), header = TRUE, sep = "\t",skip = 1,check.names = F)
df_metadata_sensitivity_binary_clusters <- read.csv(paste0(output_workspace,"Metadata_sensitivity_binary_clusters.tsv"), header = TRUE, sep = "\t",check.names = F)
#make sure the columns representing the samples are the replicate names
colnames(df_fragment_pairs_counts_ALL_replicates)[7:ncol(df_fragment_pairs_counts_ALL_replicates)] <- (sapply(X = colnames(df_fragment_pairs_counts_ALL_replicates)[7:ncol(df_fragment_pairs_counts_ALL_replicates)],FUN=function(x) gsub(pattern = "alignments/preprocessed_mapped_",replacement = "",x,fixed = TRUE)))
colnames(df_fragment_pairs_counts_ALL_replicates)[7:ncol(df_fragment_pairs_counts_ALL_replicates)] <- (sapply(X = colnames(df_fragment_pairs_counts_ALL_replicates)[7:ncol(df_fragment_pairs_counts_ALL_replicates)],FUN=function(x) gsub(pattern = "alignments/decontaminated_preprocessed_mapped_",replacement = "",x,fixed = TRUE)))
colnames(df_fragment_pairs_counts_ALL_replicates)[7:ncol(df_fragment_pairs_counts_ALL_replicates)] <- (sapply(X = colnames(df_fragment_pairs_counts_ALL_replicates)[7:ncol(df_fragment_pairs_counts_ALL_replicates)],FUN=function(x) gsub(pattern = "_Aligned.sortedByCoord.out.bam",replacement = "",x,fixed = TRUE)))
print(all(names(df_fragment_pairs_counts_ALL_replicates)[7:ncol(df_fragment_pairs_counts_ALL_replicates)]==lst_sample_replicates))

#Get FPKM from fragment pairs counts
gene_lengths <- df_fragment_pairs_counts_ALL_replicates$Length
counts_data <- df_fragment_pairs_counts_ALL_replicates[, -c(1:6)]
#rownames(counts_data) <- df_fragment_pairs_counts_ALL_replicates$Geneid
library_sizes <- colSums(counts_data)
fpkm_matrix <- t(t(counts_data) / library_sizes) * 1e6  # Normalize by library size (in millions)
fpkm_matrix <- fpkm_matrix / (gene_lengths / 1e3) # Normalize by gene length (in kilobases)
log2p1_fpkm_matrix <- log2(fpkm_matrix+1)
fpkm_matrix <- as.data.frame(fpkm_matrix)
fpkm_matrix <- cbind(df_fragment_pairs_counts_ALL_replicates$Geneid,fpkm_matrix)
names(fpkm_matrix)[1] <- "Geneid"
v_fpkm_per_pop <- fpkm_matrix["SPAC4H3.10c",]
v_fpkm_per_pop
log2p1_fpkm_matrix <- as.data.frame(log2p1_fpkm_matrix)
log2p1_fpkm_matrix <- cbind(df_fragment_pairs_counts_ALL_replicates$Geneid,log2p1_fpkm_matrix)
names(log2p1_fpkm_matrix)[1] <- "Geneid"
#save & Remove genes without variation across samples (ALL 0 count?)
fpkm_nonvariable_genes <- fpkm_matrix[which(rowSds(as.matrix(fpkm_matrix[,-1]),na.rm = T)==0),]
fpkm_matrix_without_nonvariable_genes <- fpkm_matrix[which(rowSds(as.matrix(fpkm_matrix[,-1]),na.rm = T)!=0),]
log2_fpkm_matrix_without_nonvariable_genes <- log2p1_fpkm_matrix[which(rowSds(as.matrix(log2p1_fpkm_matrix[,-1]),na.rm = T)!=0),]

#new count_data with gene names as first column
new_counts_data <- df_fragment_pairs_counts_ALL_replicates[, -c(2:6)]

#log2(FPKM+1) PCA
mtx_log2_fpkm_matrix_without_nonvariable_genes <- as.matrix(log2_fpkm_matrix_without_nonvariable_genes[,-1])
rownames(mtx_log2_fpkm_matrix_without_nonvariable_genes) <- log2_fpkm_matrix_without_nonvariable_genes$Geneid
pca <- prcomp(t(mtx_log2_fpkm_matrix_without_nonvariable_genes), scale. = TRUE)
explained_variance <- pca$sdev^2 / sum(pca$sdev^2)
variance_df <- data.frame(
  PC = 1:length(explained_variance),
  Variance = explained_variance
)

#visualize first 2 PCs
pca_df <- as.data.frame(pca$x)
pca_df$Sample <- rownames(pca_df)
pca_df <- merge(pca_df, df_metadata_sensitivity_binary_clusters, by.x = "Sample", by.y = "id")

#PCA with both replicates
ggplot(pca_df, aes(x = PC1, y = PC2, color = H2O2_bin_cluster, label = Sample)) +
  geom_point(size = 4, alpha = 0.8) +
  geom_text_repel(size = 3) +
  labs(title = "PCA Scatter Plot", x = paste0("PC1 (", round(explained_variance[1] * 100, 1), "%)"),
       y = paste0("PC2 (", round(explained_variance[2] * 100, 1), "%)")) +
  scale_color_brewer(palette = "Set2") +
  theme_minimal()
ggsave(filename = "PCA_all_sample_replicates.svg", path=output_workspace, width = 20, height = 15, units = "cm",device = svg)


##############REPLICATE 2 SHOWS A BATCH EFFECT SO ONLY FOCUS ON 1ST REPLICATE##############################
#log2(FPKM+1) PCA on replicates #1
rep1_mtx_log2_fpkm_matrix_without_nonvariable_genes <- mtx_log2_fpkm_matrix_without_nonvariable_genes[,!grepl(pattern = "_RNA_2",x = colnames(mtx_log2_fpkm_matrix_without_nonvariable_genes),fixed = T)]
#   #focus on the more sensitive vs as sensitive
# rep1_mtx_log2_fpkm_matrix_without_nonvariable_genes <- rep1_mtx_log2_fpkm_matrix_without_nonvariable_genes[,!grepl(pattern = "H11_",x = colnames(rep1_mtx_log2_fpkm_matrix_without_nonvariable_genes),fixed = T)]
# rep1_mtx_log2_fpkm_matrix_without_nonvariable_genes <- rep1_mtx_log2_fpkm_matrix_without_nonvariable_genes[,!grepl(pattern = "G11_",x = colnames(rep1_mtx_log2_fpkm_matrix_without_nonvariable_genes),fixed = T)]
# rep1_mtx_log2_fpkm_matrix_without_nonvariable_genes <- rep1_mtx_log2_fpkm_matrix_without_nonvariable_genes[,!grepl(pattern = "G5_",x = colnames(rep1_mtx_log2_fpkm_matrix_without_nonvariable_genes),fixed = T)]

rep1_mtx_log2_fpkm_matrix_without_nonvariable_genes <- rep1_mtx_log2_fpkm_matrix_without_nonvariable_genes[which(rowSds(as.matrix(rep1_mtx_log2_fpkm_matrix_without_nonvariable_genes[,-1]),na.rm = T)!=0),]
rep1_pca <- prcomp(t(rep1_mtx_log2_fpkm_matrix_without_nonvariable_genes), scale. = T,center = T)
rep1_explained_variance <- rep1_pca$sdev^2 / sum(rep1_pca$sdev^2)
rep1_variance_df <- data.frame(
  PC = 1:length(rep1_explained_variance),
  Variance = rep1_explained_variance
)
#Scree Plot
ggplot(rep1_variance_df, aes(x = PC, y = Variance)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_line(group = 1, color = "darkblue") +
  geom_point(size = 3, color = "darkblue") +
  labs(title = "Scree Plot", x = "Principal Component", y = "Proportion of Variance") +
  theme_minimal()
#dataframe for the first 2 PCs
rep1_pca_df <- as.data.frame(rep1_pca$x)
rep1_pca_df$Sample <- rownames(rep1_pca_df)
rep1_pca_df <- merge(rep1_pca_df, df_metadata_sensitivity_binary_clusters, by.x = "Sample", by.y = "id")

#Plot the samples in the first 2PCs
ggplot(rep1_pca_df, aes(x = PC1, y = PC2, color = H2O2_bin_cluster, label = Sample)) +
  geom_point(size = 4, alpha = 0.8) +
  geom_text_repel(size = 3) +
  labs(title = "PCA Scatter Plot", x = paste0("PC1 (", round(rep1_explained_variance[1] * 100, 1), "%)"),
       y = paste0("PC2 (", round(rep1_explained_variance[2] * 100, 1), "%)")) +
  scale_color_brewer(palette = "Set2") +
  theme_minimal()
ggsave(filename = "PCA_replicate1.svg", path=output_workspace, width = 20, height = 15, units = "cm",device = svg)
#View genes loadings
#View(as.data.frame(rep1_pca$rotation))

#new count_data with gene names as first column
rep1_new_counts_data <- new_counts_data[,!grepl(pattern = "_RNA_2",x = colnames(new_counts_data),fixed = T)]

#transcriptomic clustering of samples' first replicates (rep1)
the_hclust_rep1<-hclust(vegan::vegdist(x = t(mtx_log2_fpkm_matrix_without_nonvariable_genes)[!grepl(pattern = "_RNA_2",x = colnames(mtx_log2_fpkm_matrix_without_nonvariable_genes),fixed = T),],method = "bray"),method = "ward.D2")
plot(the_hclust_rep1,xlab="Samples",main="Ward D2 clusters of expression-based Bray-Curtis distance")
rect.hclust(the_hclust_rep1, k = 2, border = "red") 

#Now we know the identity of WT-like populations (similar to WT in term of transcriptomic profile AND H2O2 sensitivity): G9 and H9
#AND the identity of populations that are more sensitive than WT to H2O2 and transcriptomically dissimilar to WT: G1 to G4.
#Therefore, we can perform a differential expression analysis to reveal the transcriptomic changes that are NOT IN THE WT 
#AND contribute to increasing the sensitivity to OS after adapting to HS.
#For the populations that are more sensitive than WT to H2O2 but transcriptomically similar to WT (G11 and H11), 
#the OS sensitivity is probably explained by the genomic hits. As for the single population that has the same sensitivity 
#to H2O2 than WT but is transcriptomically dissimilar to WT (G5), its adaptive route to HS did not lead to sensitivity to OS.
rep1_new_counts_data_clearly_as_sensitive_vs_more_sensitive <- rep1_new_counts_data[,colnames(rep1_new_counts_data)%in% c("Geneid",subset(df_metadata_sensitivity_binary_clusters,(!(grepl(pattern = "_RNA_2",x = id,fixed = T)|grepl(pattern = "H11_",x = id,fixed = T)|grepl(pattern = "G11_",x = id,fixed = T)|grepl(pattern = "G5_",x = id,fixed = T))) )$id )]
rep1_dds_clearly_as_sensitive_vs_more_sensitive <- DESeqDataSetFromMatrix(countData=rep1_new_counts_data_clearly_as_sensitive_vs_more_sensitive,colData=subset(df_metadata_sensitivity_binary_clusters,!(grepl(pattern = "_RNA_2",x = id,fixed = T)|grepl(pattern = "H11_",x = id,fixed = T)|grepl(pattern = "G11_",x = id,fixed = T)|grepl(pattern = "G5_",x = id,fixed = T)) ), design=~H2O2_bin_cluster, tidy = TRUE)
rep1_dds_clearly_as_sensitive_vs_more_sensitive
rep1_dds_clearly_as_sensitive_vs_more_sensitive <- DESeq(rep1_dds_clearly_as_sensitive_vs_more_sensitive)
rep1_res_clearly_as_sensitive_vs_more_sensitive <- results(rep1_dds_clearly_as_sensitive_vs_more_sensitive,alpha = 0.05)
rep1_res_clearly_as_sensitive_vs_more_sensitive <- rep1_res_clearly_as_sensitive_vs_more_sensitive[order(-abs(rep1_res_clearly_as_sensitive_vs_more_sensitive$log2FoldChange)),]
df_rep1_res_clearly_as_sensitive_vs_more_sensitive <- as.data.frame(rep1_res_clearly_as_sensitive_vs_more_sensitive@listData)
df_rep1_res_clearly_as_sensitive_vs_more_sensitive$Gene <- rep1_res_clearly_as_sensitive_vs_more_sensitive@rownames
df_rep1_res_clearly_as_sensitive_vs_more_sensitive$abs_log2FoldChange <- abs(df_rep1_res_clearly_as_sensitive_vs_more_sensitive$log2FoldChange)
summary(rep1_res_clearly_as_sensitive_vs_more_sensitive)
nrow(subset(df_rep1_res_clearly_as_sensitive_vs_more_sensitive,padj<0.05&log2FoldChange < -0.585))
nrow(subset(df_rep1_res_clearly_as_sensitive_vs_more_sensitive,padj<0.05&log2FoldChange > 0.585))
View(subset(df_rep1_res_clearly_as_sensitive_vs_more_sensitive,padj<0.05&log2FoldChange < -0.585))
View(subset(df_rep1_res_clearly_as_sensitive_vs_more_sensitive,padj<0.05&log2FoldChange > 0.585))

#save downregulated differential genes
df_downregulated_gene_in_more_sensitive_compared_to_as_sensitive <- subset(df_rep1_res_clearly_as_sensitive_vs_more_sensitive,padj<0.05&log2FoldChange < -0.585)[,c("Gene","log2FoldChange","padj")]
colnames(df_downregulated_gene_in_more_sensitive_compared_to_as_sensitive) <- c("Systematic_ID","log2FoldChange","padj")
write.table(x=df_downregulated_gene_in_more_sensitive_compared_to_as_sensitive,file = paste0(output_workspace,"Table_Downregulated_genes_in_more_sensitive_compared_to_as_sensitive.csv"),sep = "\t",na = "NA",row.names = FALSE,col.names = TRUE)
#save upregulated differential genes
write.table(x=subset(df_rep1_res_clearly_as_sensitive_vs_more_sensitive,padj<0.05&log2FoldChange>0.585)$Gene,file = paste0(output_workspace,"upregulated_differential_genes_clearly_as_sensitive_vs_more_sensitive.csv"),sep = "\t",na = "NA",row.names = FALSE,col.names = FALSE)
df_upregulated_gene_in_more_sensitive_compared_to_as_sensitive <- subset(df_rep1_res_clearly_as_sensitive_vs_more_sensitive,padj<0.05&log2FoldChange>0.585)[,c("Gene","log2FoldChange","padj")]
colnames(df_upregulated_gene_in_more_sensitive_compared_to_as_sensitive) <- c("Systematic_ID","log2FoldChange","padj")
Upregulated_genes_in_More_sensitive_pops_Compared_to_as_sensitive_pops <- read.delim(paste0(output_workspace,"Upregulated_genes_in_More_sensitive_pops_Compared_to_as_sensitive_pops.tsv"))
df_upregulated_gene_in_more_sensitive_compared_to_as_sensitive <- merge(x = df_upregulated_gene_in_more_sensitive_compared_to_as_sensitive, y = Upregulated_genes_in_More_sensitive_pops_Compared_to_as_sensitive_pops, by = "Systematic_ID")
write.table(x=df_upregulated_gene_in_more_sensitive_compared_to_as_sensitive,file = paste0(output_workspace,"Table_S1_upregulated_genes_in_more_sensitive_compared_to_as_sensitive.csv"),sep = "\t",na = "NA",row.names = FALSE,col.names = TRUE)

#Multi-hit genes
#bmc1 (SPBC2A9.10); #bop1 (SPAP32A8.03c); #clg1 (SPBC1D7.03); 
#cmr2 (SPAC56F8.02); #ctf18 (SPBC902.02c); #dep1 (SPBC21C3.02c);
#fet4 (SPBP26C9.03c); #fhn1 (SPBC1685.13); #iss9 (SPBC2A9.11c);
#pfl2 (SPAPB15E9.01c); #pfl5 (SPBC1289.15); #pof15 (SPAPB1A10.14);
#pqr1 (SPAC6B12.07c); #prr1 (SPAC8C9.14); #pyk1 (SPAC4H3.10c);
#pzh1 (SPAC57A7.08); #rpp1 (SPAC3A12.04c); #seb1 (SPAC222.09);
#ubp9 (SPBC1703.12)
View(subset(df_rep1_res_clearly_as_sensitive_vs_more_sensitive,padj<0.05&abs_log2FoldChange > 0.585))
#Multi-hit genes differential expression in populations that are more sensitive than WT to H2O2 (G1, G2, G3 and G4) vs WT-like populations (G9 and H9)
#bmc1 (downregulated); fhn1 (upregulated); iss9 (downregulated); 
#pof15 (upregulated); pqr1 (upregulated); pyk1 (downregulated);
#pzh1 (upregulated); rpp1 (downregulated); ubp9 (upregulated);

#Volcano plot
df_rep1_res_clearly_as_sensitive_vs_more_sensitive$Significance <- ifelse(df_rep1_res_clearly_as_sensitive_vs_more_sensitive$padj<0.05 , yes = "Significant (FDR<0.05)",no = "Not Significant")
df_rep1_res_clearly_as_sensitive_vs_more_sensitive$Differences <- ifelse(df_rep1_res_clearly_as_sensitive_vs_more_sensitive$padj<0.05&df_rep1_res_clearly_as_sensitive_vs_more_sensitive$log2FoldChange < -0.585 , yes = "Down-regulated",no =ifelse(df_rep1_res_clearly_as_sensitive_vs_more_sensitive$padj<0.05&df_rep1_res_clearly_as_sensitive_vs_more_sensitive$log2FoldChange > 0.585 , yes = "Up-regulated",no = "Not Significant"))
ggplot(df_rep1_res_clearly_as_sensitive_vs_more_sensitive, aes(x = log2FoldChange, y = -log10(padj), color = Differences)) +
  geom_point(alpha = 0.8, size = 2) +  # Add points
  scale_color_manual(values = c("Down-regulated" = "#E41A1C", "Up-regulated" = "#377EB8","Not Significant" = "grey")) +  # Customize colors
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", size = 0.5) +  # Add threshold line
  geom_vline(xintercept = c(-0.585, 0.585), linetype = "dashed", color = "black", size = 0.5) +  # Add fold change thresholds
  labs(
    title = "Volcano Plot",
    x = "Log2 Fold Change",
    y = "-Log10(p-value)"
  ) +
  theme_minimal() +
  theme(
    legend.position = "top",
    plot.title = element_text(hjust = 0.5, size = 14)
  ) + xlim(-3,5) + ylim(0,8) 
ggsave(filename = "Volcano_plot_Pops_more_sensitive_to_H2O2_than_WT_vs_WT_like_pop.svg", path=output_workspace, width = 20, height = 15, units = "cm",device = svg)

#PYK1 differences in 1st replicates (PYK1 Supplementary Figure SI)
df_gd_rp_fpkm <- fpkm_matrix_without_nonvariable_genes[,!(grepl(pattern = "_RNA_2",x = colnames(new_counts_data),fixed = T))] #|grepl(pattern = "972h-_",x = colnames(new_counts_data),fixed = T)|grepl(pattern = "H11_",x = colnames(new_counts_data),fixed = T)|grepl(pattern = "G11_",x = colnames(new_counts_data),fixed = T)|grepl(pattern = "G5",x = colnames(new_counts_data),fixed = T)
df_pyk1 <- data.frame(Population = gsub(pattern = "_.*", "", colnames(df_gd_rp_fpkm[,-1])),pyk1_expression=unlist(unname(df_gd_rp_fpkm[df_gd_rp_fpkm$Geneid=="SPAC4H3.10c",-1])),stringsAsFactors = F)
rownames(df_pyk1) <- df_pyk1$Population
df_pyk1$sensitivity_cluster <- ifelse(test = df_pyk1$Population%in%c("972h-","G5","G9","H9"),yes = "WT-like sensitivity",no = "More sensitive than WT")
ggplot(data = df_pyk1, mapping = aes(x = sensitivity_cluster , y = pyk1_expression)) +
  geom_boxplot(fill = "grey") + geom_point() + 
  geom_text(aes(label = Population)) +
  labs(y = expression(~italic(pyk1)~"expression level (FPKM)")) + xlab("Population H2O2 sensitivity") + #stat_compare_means(method = "wilcox",paired = F,method.args = list(alternative="less")) +
  theme_bw() + theme(axis.title = element_text(size=16),axis.text = element_text(size=16),legend.title = element_text(size=14),legend.text = element_text(size=12),axis.text.x = element_text(size=12),legend.position = "None")
ggsave(filename = "pyk1_Expression_level_by_H2O2_sensitivity_including_WT_and_all_pops.svg", path=output_workspace, width = 20, height = 15, units = "cm",device = svg)

ggplot(data = subset(df_pyk1,!Population%in%c("G5","G11","H11","H12")), mapping = aes(x = factor(sensitivity_cluster,levels = c("WT-like sensitivity","More sensitive than WT")) , y = pyk1_expression)) +
  geom_boxplot(fill = "grey") + geom_point() + 
  geom_text(aes(label = Population)) +
  scale_y_continuous(limits = c(0,3500),breaks = seq(0,3500,500)) +
  labs(y = expression(~italic(pyk1)~"expression level (FPKM)"),title=paste0("pyk1 DESeq2 Wald test results: log2(Fold change) = ", subset(df_rep1_res_clearly_as_sensitive_vs_more_sensitive,Gene=="SPAC4H3.10c")$log2FoldChange,"; FDR = ",subset(df_rep1_res_clearly_as_sensitive_vs_more_sensitive,Gene=="SPAC4H3.10c")$padj)) + xlab("Population H2O2 sensitivity") + #stat_compare_means(method = "wilcox",paired = F,method.args = list(alternative="greater")) +
  theme_bw() + theme(axis.title = element_text(size=16),axis.text = element_text(size=16),legend.title = element_text(size=14),legend.text = element_text(size=12),axis.text.x = element_text(size=12),title = element_text(size=7),legend.position = "None")
ggsave(filename = "pyk1_Expression_level_by_H2O2_sensitivity_excluding_outlier_pops_in_rep1.svg", path=output_workspace, width = 20, height = 15, units = "cm",device = svg)

#plotCounts(rep1_dds_clearly_as_sensitive_vs_more_sensitive, gene="SPAC4H3.10c", intgroup="H2O2_bin_cluster")

#save session
library("session")
save.session(file = paste0(output_workspace,"Differential_expression_analysis_H2O2_sensitivity_clusters.Rda"))
