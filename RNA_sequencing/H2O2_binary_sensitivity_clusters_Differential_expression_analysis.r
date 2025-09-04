#import packages
library( "DESeq2" )
library("ggplot2")
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


#Create DeSeq object
dds <- DESeqDataSetFromMatrix(countData=new_counts_data,colData=df_metadata_sensitivity_binary_clusters, design=~H2O2_bin_cluster, tidy = TRUE)
dds
dds <- DESeq(dds)
res <- results(dds,alpha = 0.1)
res <- res[order(res$padj),]
head(res)
summary(res)

#Volcano plot
#reset par
par(mfrow=c(1,1))
# Make a basic volcano plot
with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-3,3)))

# Add colored points: blue if padj<0.01, red if log2FC>1 and padj<0.05)
with(subset(res, padj< 0.1 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
# with(subset(rep1_res, padj<0.1 & (log2FoldChange <= -1 || log2FoldChange >= 0.585)), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
with(subset(res, padj< 0.1 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
# with(subset(rep1_res, padj<0.1 & (log2FoldChange <= -1 || log2FoldChange >= 0.585)), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))

#raw count PCA
vsdata <- vst(dds, blind=FALSE)
plotPCA(vsdata, intgroup="H2O2_bin_cluster") 

#log2(FPKM+1) PCA
mtx_log2_fpkm_matrix_without_nonvariable_genes <- as.matrix(log2_fpkm_matrix_without_nonvariable_genes[,-1])
rownames(mtx_log2_fpkm_matrix_without_nonvariable_genes) <- log2_fpkm_matrix_without_nonvariable_genes$Geneid
pca <- prcomp(t(mtx_log2_fpkm_matrix_without_nonvariable_genes), scale. = TRUE)
explained_variance <- pca$sdev^2 / sum(pca$sdev^2)
variance_df <- data.frame(
  PC = 1:length(explained_variance),
  Variance = explained_variance
)
  #Scree Plot
ggplot(variance_df, aes(x = PC, y = Variance)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_line(group = 1, color = "darkblue") +
  geom_point(size = 3, color = "darkblue") +
  labs(title = "Scree Plot", x = "Principal Component", y = "Proportion of Variance") +
  theme_minimal()
  #visualize first 2 PCs
pca_df <- as.data.frame(pca$x)
pca_df$Sample <- rownames(pca_df)
pca_df <- merge(pca_df, df_metadata_sensitivity_binary_clusters, by.x = "Sample", by.y = "id")

# Plot
ggplot(pca_df, aes(x = PC1, y = PC2, color = H2O2_bin_cluster, label = Sample)) +
  geom_point(size = 4, alpha = 0.8) +
  geom_text_repel(size = 3) +
  labs(title = "PCA Scatter Plot", x = paste0("PC1 (", round(explained_variance[1] * 100, 1), "%)"),
       y = paste0("PC2 (", round(explained_variance[2] * 100, 1), "%)")) +
  scale_color_brewer(palette = "Set2") +
  theme_minimal()
#View genes loadings
#View(as.data.frame(pca$rotation))

#PYK1 differences
plotCounts(dds, gene="SPAC4H3.10c", intgroup="H2O2_bin_cluster")

#PYK1 vs ACT1
# v_from_sample_sensitivity <- df_metadata_sensitivity_binary_clusters$H2O2_bin_cluster
# names(v_from_sample_sensitivity) <- df_metadata_sensitivity_binary_clusters$id
# subset_new_count <- new_counts_data[new_counts_data$Geneid%in%c("SPAC4H3.10c","SPBC32H8.12c"),-1]
# ggplot(mapping = aes(x=reorder(colnames(subset_new_count),-unname(as.matrix(subset_new_count[1,]/subset_new_count[2,])[1,])),y=unname(as.matrix(subset_new_count[1,]/subset_new_count[2,])[1,]),fill=v_from_sample_sensitivity[colnames(subset_new_count)] )) + geom_col() +theme(axis.text.x = element_text(size=10,angle = 30,hjust=1)) + ylab("PYK1/ACT1") + xlab("sample")

##################################################ONLY FOCUS ON 1ST REPLICATE##############################
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
#visualize first 2 PCs
rep1_pca_df <- as.data.frame(rep1_pca$x)
rep1_pca_df$Sample <- rownames(rep1_pca_df)
rep1_pca_df <- merge(rep1_pca_df, df_metadata_sensitivity_binary_clusters, by.x = "Sample", by.y = "id")

# Plot
ggplot(rep1_pca_df, aes(x = PC1, y = PC2, color = H2O2_bin_cluster, label = Sample)) +
  geom_point(size = 4, alpha = 0.8) +
  geom_text_repel(size = 3) +
  labs(title = "PCA Scatter Plot", x = paste0("PC1 (", round(rep1_explained_variance[1] * 100, 1), "%)"),
       y = paste0("PC2 (", round(rep1_explained_variance[2] * 100, 1), "%)")) +
  scale_color_brewer(palette = "Set2") +
  theme_minimal()
#View genes loadings
#View(as.data.frame(rep1_pca$rotation))

#new count_data with gene names as first column
rep1_new_counts_data <- new_counts_data[,!grepl(pattern = "_RNA_2",x = colnames(new_counts_data),fixed = T)]
#   #focus on the more sensitive vs as sensitive
# rep1_new_counts_data <- rep1_new_counts_data[,!grepl(pattern = "H11_",x = colnames(rep1_new_counts_data),fixed = T)]
# rep1_new_counts_data <- rep1_new_counts_data[,!grepl(pattern = "G11_",x = colnames(rep1_new_counts_data),fixed = T)]
# rep1_new_counts_data <- rep1_new_counts_data[,!grepl(pattern = "G5_",x = colnames(rep1_new_counts_data),fixed = T)]

        ##########################################

#Create DeSeq object for clearly-differentiable "as sensitive" and "more sensitive"
rep1_new_counts_data_clearly_as_sensitive_vs_more_sensitive <- rep1_new_counts_data[,colnames(rep1_new_counts_data)%in% c("Geneid",subset(df_metadata_sensitivity_binary_clusters,(!(grepl(pattern = "_RNA_2",x = id,fixed = T)|grepl(pattern = "H11_",x = id,fixed = T)|grepl(pattern = "G11_",x = id,fixed = T)|grepl(pattern = "G5_",x = id,fixed = T))) )$id )]
rep1_dds_clearly_as_sensitive_vs_more_sensitive <- DESeqDataSetFromMatrix(countData=rep1_new_counts_data_clearly_as_sensitive_vs_more_sensitive,colData=subset(df_metadata_sensitivity_binary_clusters,!(grepl(pattern = "_RNA_2",x = id,fixed = T)|grepl(pattern = "H11_",x = id,fixed = T)|grepl(pattern = "G11_",x = id,fixed = T)|grepl(pattern = "G5_",x = id,fixed = T)) ), design=~H2O2_bin_cluster, tidy = TRUE)
rep1_dds_clearly_as_sensitive_vs_more_sensitive
rep1_dds_clearly_as_sensitive_vs_more_sensitive <- DESeq(rep1_dds_clearly_as_sensitive_vs_more_sensitive)
rep1_res_clearly_as_sensitive_vs_more_sensitive <- results(rep1_dds_clearly_as_sensitive_vs_more_sensitive,alpha = 0.1)
rep1_res_clearly_as_sensitive_vs_more_sensitive <- rep1_res_clearly_as_sensitive_vs_more_sensitive[order(-abs(rep1_res_clearly_as_sensitive_vs_more_sensitive$log2FoldChange)),]
df_rep1_res_clearly_as_sensitive_vs_more_sensitive <- as.data.frame(rep1_res_clearly_as_sensitive_vs_more_sensitive@listData)
df_rep1_res_clearly_as_sensitive_vs_more_sensitive$Gene <- rep1_res_clearly_as_sensitive_vs_more_sensitive@rownames
df_rep1_res_clearly_as_sensitive_vs_more_sensitive$abs_log2FoldChange <- abs(df_rep1_res_clearly_as_sensitive_vs_more_sensitive$log2FoldChange)
View(subset(df_rep1_res_clearly_as_sensitive_vs_more_sensitive,padj<0.1))
summary(rep1_res_clearly_as_sensitive_vs_more_sensitive)

#save downregulated differential genes log2FC<-0.585 (or absolute FC of 1/1.5) and padj<0.1
write.table(x=subset(df_rep1_res_clearly_as_sensitive_vs_more_sensitive,padj<0.1&log2FoldChange < -1)$Gene,file = paste0(output_workspace,"downregulated_differential_genes_clearly_as_sensitive_vs_more_sensitive.csv"),sep = "\t",na = "NA",row.names = FALSE,col.names = FALSE)
#save upregulared differential genes
write.table(x=subset(df_rep1_res_clearly_as_sensitive_vs_more_sensitive,padj<0.1&log2FoldChange>0.585)$Gene,file = paste0(output_workspace,"upregulated_differential_genes_clearly_as_sensitive_vs_more_sensitive.csv"),sep = "\t",na = "NA",row.names = FALSE,col.names = FALSE)
df_upregulated_gene_in_more_sensitive_compared_to_as_sensitive <- subset(df_rep1_res_clearly_as_sensitive_vs_more_sensitive,padj<0.1&log2FoldChange>0.585)[,c("Gene","log2FoldChange","padj")]
colnames(df_upregulated_gene_in_more_sensitive_compared_to_as_sensitive) <- c("Systematic_ID","log2FoldChange","padj")
Upregulated_genes_in_More_sensitive_pops_Compared_to_as_sensitive_pops <- read.delim(paste0(output_workspace,"Upregulated_genes_in_More_sensitive_pops_Compared_to_as_sensitive_pops.tsv"))
df_upregulated_gene_in_more_sensitive_compared_to_as_sensitive <- merge(x = df_upregulated_gene_in_more_sensitive_compared_to_as_sensitive, y = Upregulated_genes_in_More_sensitive_pops_Compared_to_as_sensitive_pops, by = "Systematic_ID")
write.table(x=df_upregulated_gene_in_more_sensitive_compared_to_as_sensitive,file = paste0(output_workspace,"Table_S1_upregulated_gene_in_more_sensitive_compared_to_as_sensitive.csv"),sep = "\t",na = "NA",row.names = FALSE,col.names = TRUE)

#Volcano plot
df_rep1_res_clearly_as_sensitive_vs_more_sensitive$Significance <- ifelse(df_rep1_res_clearly_as_sensitive_vs_more_sensitive$padj<0.1 , yes = "Significant (FDR<0.1)",no = "Not Significant")
df_rep1_res_clearly_as_sensitive_vs_more_sensitive$Differences <- ifelse(df_rep1_res_clearly_as_sensitive_vs_more_sensitive$padj<0.1&df_rep1_res_clearly_as_sensitive_vs_more_sensitive$log2FoldChange < -1 , yes = "Down-regulated",no =ifelse(df_rep1_res_clearly_as_sensitive_vs_more_sensitive$padj<0.1&df_rep1_res_clearly_as_sensitive_vs_more_sensitive$log2FoldChange > 0.585 , yes = "Up-regulated",no = "Not Significant"))
ggplot(df_rep1_res_clearly_as_sensitive_vs_more_sensitive, aes(x = log2FoldChange, y = -log10(pvalue), color = Differences)) +
  geom_point(alpha = 0.8, size = 2) +  # Add points
  scale_color_manual(values = c("Down-regulated" = "#E41A1C", "Up-regulated" = "#377EB8","Not Significant" = "grey")) +  # Customize colors
  geom_hline(yintercept = -log10(0.1), linetype = "dashed", color = "black", size = 0.5) +  # Add threshold line
  geom_vline(xintercept = c(-1, 0.585), linetype = "dashed", color = "black", size = 0.5) +  # Add fold change thresholds
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

# 
# #reset par
# par(mfrow=c(1,1))
# # Make a basic volcano plot
# with(rep1_res_clearly_as_sensitive_vs_more_sensitive, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-3,5),ylim=c(0,8) ))
# 
# # Add colored points: blue if padj<0.1, red if log2FC>=0.585 (increase of >=1.5) or log2FC <= -1 (or absolute decrease of <= 50%) and padj<0.1)
# with(subset(rep1_res_clearly_as_sensitive_vs_more_sensitive, padj< 0.1 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
# with(subset(rep1_res_clearly_as_sensitive_vs_more_sensitive, padj<0.1 & (log2FoldChange <= -1 || log2FoldChange >= 0.585)), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))

#PYK1 differences in 1st replicates
plotCounts(rep1_dds_clearly_as_sensitive_vs_more_sensitive, gene="SPAC4H3.10c", intgroup="H2O2_bin_cluster")
plotCounts(rep1_dds_clearly_as_sensitive_vs_more_sensitive, gene="SPAC869.06c", intgroup="H2O2_bin_cluster")
plotCounts(rep1_dds_clearly_as_sensitive_vs_more_sensitive, gene="SPBC1D7.03", intgroup="H2O2_bin_cluster")

###############################################
#Create DeSeq object for clearly-differentiable "as sensitive" and "way more sensitive"
rep1_new_counts_data_clearly_as_sensitive_vs_way_more_sensitive <- rep1_new_counts_data[,colnames(rep1_new_counts_data)%in% c("Geneid",subset(df_metadata_sensitivity_binary_clusters,(!(grepl(pattern = "_RNA_2",x = id,fixed = T)|grepl(pattern = "G11_",x = id,fixed = T)|grepl(pattern = "G1_",x = id,fixed = T)|grepl(pattern = "G2_",x = id,fixed = T)|grepl(pattern = "G3_",x = id,fixed = T)|grepl(pattern = "G4_",x = id,fixed = T)|grepl(pattern = "G5_",x = id,fixed = T))) )$id )]
rep1_dds_clearly_as_sensitive_vs_way_more_sensitive <- DESeqDataSetFromMatrix(countData=rep1_new_counts_data_clearly_as_sensitive_vs_way_more_sensitive,colData=subset(df_metadata_sensitivity_binary_clusters,!(grepl(pattern = "_RNA_2",x = id,fixed = T)|grepl(pattern = "G11_",x = id,fixed = T)|grepl(pattern = "G1_",x = id,fixed = T)|grepl(pattern = "G2_",x = id,fixed = T)|grepl(pattern = "G3_",x = id,fixed = T)|grepl(pattern = "G4_",x = id,fixed = T)|grepl(pattern = "G5_",x = id,fixed = T)) ), design=~H2O2_bin_cluster, tidy = TRUE)
rep1_dds_clearly_as_sensitive_vs_way_more_sensitive
rep1_dds_clearly_as_sensitive_vs_way_more_sensitive <- DESeq(rep1_dds_clearly_as_sensitive_vs_way_more_sensitive)
rep1_res_clearly_as_sensitive_vs_way_more_sensitive <- results(rep1_dds_clearly_as_sensitive_vs_way_more_sensitive,alpha = 0.1)
rep1_res_clearly_as_sensitive_vs_way_more_sensitive <- rep1_res_clearly_as_sensitive_vs_way_more_sensitive[order(-abs(rep1_res_clearly_as_sensitive_vs_way_more_sensitive$log2FoldChange)),]
df_rep1_res_clearly_as_sensitive_vs_way_more_sensitive <- as.data.frame(rep1_res_clearly_as_sensitive_vs_way_more_sensitive@listData)
df_rep1_res_clearly_as_sensitive_vs_way_more_sensitive$Gene <- rep1_res_clearly_as_sensitive_vs_way_more_sensitive@rownames
df_rep1_res_clearly_as_sensitive_vs_way_more_sensitive$abs_log2FoldChange <- abs(df_rep1_res_clearly_as_sensitive_vs_way_more_sensitive$log2FoldChange)
View(subset(df_rep1_res_clearly_as_sensitive_vs_way_more_sensitive,padj<0.1))
summary(rep1_res_clearly_as_sensitive_vs_way_more_sensitive)

#save downregulated differential genes log2FC<-0.585 (or absolute FC of 1/1.5) and padj<0.1
write.table(x=subset(df_rep1_res_clearly_as_sensitive_vs_way_more_sensitive,log2FoldChange < -1)$Gene,file = paste0(output_workspace,"downregulated_differential_genes_clearly_as_sensitive_vs_way_more_sensitive.csv"),sep = "\t",na = "NA",row.names = FALSE,col.names = FALSE)
#save upregulared differential genes
write.table(x=subset(df_rep1_res_clearly_as_sensitive_vs_way_more_sensitive,log2FoldChange>0.585)$Gene,file = paste0(output_workspace,"upregulated_differential_genes_clearly_as_sensitive_vs_way_more_sensitive.csv"),sep = "\t",na = "NA",row.names = FALSE,col.names = FALSE)

#Volcano plot
df_rep1_res_clearly_as_sensitive_vs_way_more_sensitive$Significance <- ifelse(df_rep1_res_clearly_as_sensitive_vs_way_more_sensitive$padj<0.1 , yes = "Significant (FDR<0.1)",no = "Not Significant")
df_rep1_res_clearly_as_sensitive_vs_way_more_sensitive$Differences <- ifelse(df_rep1_res_clearly_as_sensitive_vs_way_more_sensitive$padj<0.1&df_rep1_res_clearly_as_sensitive_vs_way_more_sensitive$log2FoldChange < -1 , yes = "Down-regulated",no =ifelse(df_rep1_res_clearly_as_sensitive_vs_way_more_sensitive$padj<0.1&df_rep1_res_clearly_as_sensitive_vs_way_more_sensitive$log2FoldChange > 0.585 , yes = "Up-regulated",no = "Not Significant"))
ggplot(df_rep1_res_clearly_as_sensitive_vs_way_more_sensitive, aes(x = log2FoldChange, y = -log10(pvalue), color = Differences)) +
  geom_point(alpha = 0.8, size = 2) +  # Add points
  scale_color_manual(values = c("Down-regulated" = "#E41A1C", "Up-regulated" = "#377EB8","Not Significant" = "grey")) +  # Customize colors
  geom_hline(yintercept = -log10(0.1), linetype = "dashed", color = "black", size = 0.5) +  # Add threshold line
  geom_vline(xintercept = c(-1, 0.585), linetype = "dashed", color = "black", size = 0.5) +  # Add fold change thresholds
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

# 
# #reset par
# par(mfrow=c(1,1))
# # Make a basic volcano plot
# with(rep1_res_clearly_as_sensitive_vs_way_more_sensitive, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-3,5),ylim=c(0,8) ))
# 
# # Add colored points: blue if padj<0.1, red if log2FC>=0.585 (increase of >=1.5) or log2FC <= -1 (or absolute decrease of <= 50%) and padj<0.1)
# with(subset(rep1_res_clearly_as_sensitive_vs_way_more_sensitive, padj< 0.1 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
# with(subset(rep1_res_clearly_as_sensitive_vs_way_more_sensitive, padj<0.1 & (log2FoldChange <= -1 || log2FoldChange >= 0.585)), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))

#PYK1 differences in 1st replicates
plotCounts(rep1_dds_clearly_as_sensitive_vs_way_more_sensitive, gene="SPAC4H3.10c", intgroup="H2O2_bin_cluster")
plotCounts(rep1_dds_clearly_as_sensitive_vs_way_more_sensitive, gene="SPAC869.06c", intgroup="H2O2_bin_cluster")

    ###############################################
#Create DeSeq object for clearly-differentiable "more sensitive" and "way more sensitive"
rep1_new_counts_data_clearly_more_sensitive_vs_way_more_sensitive <- rep1_new_counts_data[,colnames(rep1_new_counts_data)%in% c("Geneid",subset(df_metadata_sensitivity_binary_clusters,(!(grepl(pattern = "_RNA_2",x = id,fixed = T)|grepl(pattern = "972h-_",x = id,fixed = T)|grepl(pattern = "H9_",x = id,fixed = T)|grepl(pattern = "G9_",x = id,fixed = T)|grepl(pattern = "G11_",x = id,fixed = T)|grepl(pattern = "G5_",x = id,fixed = T))) )$id )]
rep1_dds_clearly_more_sensitive_vs_way_more_sensitive <- DESeqDataSetFromMatrix(countData=rep1_new_counts_data_clearly_more_sensitive_vs_way_more_sensitive,colData=subset(df_metadata_sensitivity_binary_clusters,!(grepl(pattern = "_RNA_2",x = id,fixed = T)|grepl(pattern = "972h-_",x = id,fixed = T)|grepl(pattern = "H9_",x = id,fixed = T)|grepl(pattern = "G9_",x = id,fixed = T)|grepl(pattern = "G11_",x = id,fixed = T)|grepl(pattern = "G5_",x = id,fixed = T)) ), design=~H2O2_bin_cluster, tidy = TRUE)
rep1_dds_clearly_more_sensitive_vs_way_more_sensitive
rep1_dds_clearly_more_sensitive_vs_way_more_sensitive <- DESeq(rep1_dds_clearly_more_sensitive_vs_way_more_sensitive)
rep1_res_clearly_more_sensitive_vs_way_more_sensitive <- results(rep1_dds_clearly_more_sensitive_vs_way_more_sensitive,alpha = 0.1)
rep1_res_clearly_more_sensitive_vs_way_more_sensitive <- rep1_res_clearly_more_sensitive_vs_way_more_sensitive[order(-abs(rep1_res_clearly_more_sensitive_vs_way_more_sensitive$log2FoldChange)),]
df_rep1_res_clearly_more_sensitive_vs_way_more_sensitive <- as.data.frame(rep1_res_clearly_more_sensitive_vs_way_more_sensitive@listData)
df_rep1_res_clearly_more_sensitive_vs_way_more_sensitive$Gene <- rep1_res_clearly_more_sensitive_vs_way_more_sensitive@rownames
df_rep1_res_clearly_more_sensitive_vs_way_more_sensitive$abs_log2FoldChange <- abs(df_rep1_res_clearly_more_sensitive_vs_way_more_sensitive$log2FoldChange)
View(subset(df_rep1_res_clearly_more_sensitive_vs_way_more_sensitive,padj<0.1))
summary(rep1_res_clearly_more_sensitive_vs_way_more_sensitive)

#save downregulated differential genes log2FC<-0.585 (or absolute FC of 1/1.5) and padj<0.1
write.table(x=subset(df_rep1_res_clearly_more_sensitive_vs_way_more_sensitive,padj<0.1&log2FoldChange < -1)$Gene,file = paste0(output_workspace,"downregulated_differential_genes_clearly_more_sensitive_vs_way_more_sensitive.csv"),sep = "\t",na = "NA",row.names = FALSE,col.names = FALSE)
#save upregulared differential genes
write.table(x=subset(df_rep1_res_clearly_more_sensitive_vs_way_more_sensitive,padj<0.1&log2FoldChange>0.585)$Gene,file = paste0(output_workspace,"upregulated_differential_genes_clearly_more_sensitive_vs_way_more_sensitive.csv"),sep = "\t",na = "NA",row.names = FALSE,col.names = FALSE)

#Volcano plot
df_rep1_res_clearly_more_sensitive_vs_way_more_sensitive$Significance <- ifelse(df_rep1_res_clearly_more_sensitive_vs_way_more_sensitive$padj<0.1 , yes = "Significant (FDR<0.1)",no = "Not Significant")
df_rep1_res_clearly_more_sensitive_vs_way_more_sensitive$Differences <- ifelse(df_rep1_res_clearly_more_sensitive_vs_way_more_sensitive$padj<0.1&df_rep1_res_clearly_more_sensitive_vs_way_more_sensitive$log2FoldChange < -1 , yes = "Down-regulated",no =ifelse(df_rep1_res_clearly_more_sensitive_vs_way_more_sensitive$padj<0.1&df_rep1_res_clearly_more_sensitive_vs_way_more_sensitive$log2FoldChange > 0.585 , yes = "Up-regulated",no = "Not Significant"))
ggplot(df_rep1_res_clearly_more_sensitive_vs_way_more_sensitive, aes(x = log2FoldChange, y = -log10(pvalue), color = Differences)) +
  geom_point(alpha = 0.8, size = 2) +  # Add points
  scale_color_manual(values = c("Down-regulated" = "#E41A1C", "Up-regulated" = "#377EB8","Not Significant" = "grey")) +  # Customize colors
  geom_hline(yintercept = -log10(0.1), linetype = "dashed", color = "black", size = 0.5) +  # Add threshold line
  geom_vline(xintercept = c(-1, 0.585), linetype = "dashed", color = "black", size = 0.5) +  # Add fold change thresholds
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

# 
# #reset par
# par(mfrow=c(1,1))
# # Make a basic volcano plot
# with(rep1_res_clearly_more_sensitive_vs_way_more_sensitive, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-3,5),ylim=c(0,8) ))
# 
# # Add colored points: blue if padj<0.1, red if log2FC>=0.585 (increase of >=1.5) or log2FC <= -1 (or absolute decrease of <= 50%) and padj<0.1)
# with(subset(rep1_res_clearly_more_sensitive_vs_way_more_sensitive, padj< 0.1 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
# with(subset(rep1_res_clearly_more_sensitive_vs_way_more_sensitive, padj<0.1 & (log2FoldChange <= -1 || log2FoldChange >= 0.585)), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))

#PYK1 differences in 1st replicates
plotCounts(rep1_dds_clearly_more_sensitive_vs_way_more_sensitive, gene="SPAC4H3.10c", intgroup="H2O2_bin_cluster")
plotCounts(rep1_dds_clearly_more_sensitive_vs_way_more_sensitive, gene="SPAC869.06c", intgroup="H2O2_bin_cluster")

###############################################
#Create DeSeq object for "outlier as sensitive" and other "as sensitive"
df_metadata_sensitivity_binary_clusters$lbl_2 <- df_metadata_sensitivity_binary_clusters$H2O2_bin_cluster
df_metadata_sensitivity_binary_clusters$lbl_2[df_metadata_sensitivity_binary_clusters$id=="G5_RNA_1_S359"] <- "Outlier as sensitive"
rep1_new_counts_data_outlier_as_sensitive_vs_other_as_sensitive <- rep1_new_counts_data[,colnames(rep1_new_counts_data)%in% c("Geneid",subset(df_metadata_sensitivity_binary_clusters,(!(grepl(pattern = "_RNA_2",x = id,fixed = T)|grepl(pattern = "H11_",x = id,fixed = T)|grepl(pattern = "G11_",x = id,fixed = T)|grepl(pattern = "G1_",x = id,fixed = T)|grepl(pattern = "G2_",x = id,fixed = T)|grepl(pattern = "G3_",x = id,fixed = T)|grepl(pattern = "G4_",x = id,fixed = T))) )$id )]
rep1_dds_outlier_as_sensitive_vs_other_as_sensitive <- DESeqDataSetFromMatrix(countData=rep1_new_counts_data_outlier_as_sensitive_vs_other_as_sensitive,colData=subset(df_metadata_sensitivity_binary_clusters,!(grepl(pattern = "_RNA_2",x = id,fixed = T)|grepl(pattern = "H11_",x = id,fixed = T)|grepl(pattern = "G11_",x = id,fixed = T)|grepl(pattern = "G1_",x = id,fixed = T)|grepl(pattern = "G2_",x = id,fixed = T)|grepl(pattern = "G3_",x = id,fixed = T)|grepl(pattern = "G4_",x = id,fixed = T)) ), design=~lbl_2, tidy = TRUE)
rep1_dds_outlier_as_sensitive_vs_other_as_sensitive
rep1_dds_outlier_as_sensitive_vs_other_as_sensitive <- DESeq(rep1_dds_outlier_as_sensitive_vs_other_as_sensitive)
rep1_res_outlier_as_sensitive_vs_other_as_sensitive <- results(rep1_dds_outlier_as_sensitive_vs_other_as_sensitive,alpha = 0.1)
rep1_res_outlier_as_sensitive_vs_other_as_sensitive <- rep1_res_outlier_as_sensitive_vs_other_as_sensitive[order(-abs(rep1_res_outlier_as_sensitive_vs_other_as_sensitive$log2FoldChange)),]
df_rep1_res_outlier_as_sensitive_vs_other_as_sensitive <- as.data.frame(rep1_res_outlier_as_sensitive_vs_other_as_sensitive@listData)
df_rep1_res_outlier_as_sensitive_vs_other_as_sensitive$Gene <- rep1_res_outlier_as_sensitive_vs_other_as_sensitive@rownames
df_rep1_res_outlier_as_sensitive_vs_other_as_sensitive$abs_log2FoldChange <- abs(df_rep1_res_outlier_as_sensitive_vs_other_as_sensitive$log2FoldChange)
View(subset(df_rep1_res_outlier_as_sensitive_vs_other_as_sensitive,padj<0.1))
summary(rep1_res_outlier_as_sensitive_vs_other_as_sensitive)

#save downregulated differential genes log2FC<-0.585 (or absolute FC of 1/1.5) and padj<0.1
write.table(x=subset(df_rep1_res_outlier_as_sensitive_vs_other_as_sensitive,padj<0.1&log2FoldChange < -1)$Gene,file = paste0(output_workspace,"downregulated_differential_genes_outlier_as_sensitive_vs_other_as_sensitive.csv"),sep = "\t",na = "NA",row.names = FALSE,col.names = FALSE)
#save upregulared differential genes
write.table(x=subset(df_rep1_res_outlier_as_sensitive_vs_other_as_sensitive,padj<0.1&log2FoldChange>0.585)$Gene,file = paste0(output_workspace,"upregulated_differential_genes_outlier_as_sensitive_vs_other_as_sensitive.csv"),sep = "\t",na = "NA",row.names = FALSE,col.names = FALSE)

#Volcano plot
df_rep1_res_outlier_as_sensitive_vs_other_as_sensitive$Significance <- ifelse(df_rep1_res_outlier_as_sensitive_vs_other_as_sensitive$padj<0.1 , yes = "Significant (FDR<0.1)",no = "Not Significant")
df_rep1_res_outlier_as_sensitive_vs_other_as_sensitive$Differences <- ifelse(df_rep1_res_outlier_as_sensitive_vs_other_as_sensitive$padj<0.1&df_rep1_res_outlier_as_sensitive_vs_other_as_sensitive$log2FoldChange < -1 , yes = "Down-regulated",no =ifelse(df_rep1_res_outlier_as_sensitive_vs_other_as_sensitive$padj<0.1&df_rep1_res_outlier_as_sensitive_vs_other_as_sensitive$log2FoldChange > 0.585 , yes = "Up-regulated",no = "Not Significant"))
ggplot(df_rep1_res_outlier_as_sensitive_vs_other_as_sensitive, aes(x = log2FoldChange, y = -log10(pvalue), color = Differences)) +
  geom_point(alpha = 0.8, size = 2) +  # Add points
  scale_color_manual(values = c("Down-regulated" = "#E41A1C", "Up-regulated" = "#377EB8","Not Significant" = "grey")) +  # Customize colors
  geom_hline(yintercept = -log10(0.1), linetype = "dashed", color = "black", size = 0.5) +  # Add threshold line
  geom_vline(xintercept = c(-1, 0.585), linetype = "dashed", color = "black", size = 0.5) +  # Add fold change thresholds
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

# 
# #reset par
# par(mfrow=c(1,1))
# # Make a basic volcano plot
# with(rep1_res_outlier_as_sensitive_vs_other_as_sensitive, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-3,5),ylim=c(0,8) ))
# 
# # Add colored points: blue if padj<0.1, red if log2FC>=0.585 (increase of >=1.5) or log2FC <= -1 (or absolute decrease of <= 50%) and padj<0.1)
# with(subset(rep1_res_outlier_as_sensitive_vs_other_as_sensitive, padj< 0.1 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
# with(subset(rep1_res_outlier_as_sensitive_vs_other_as_sensitive, padj<0.1 & (log2FoldChange <= -1 || log2FoldChange >= 0.585)), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))

#PYK1 differences in 1st replicates
plotCounts(rep1_dds_outlier_as_sensitive_vs_other_as_sensitive, gene="SPAC4H3.10c", intgroup="lbl_2")
plotCounts(rep1_dds_outlier_as_sensitive_vs_other_as_sensitive, gene="SPAC869.06c", intgroup="lbl_2")

######################################################################REPLICATE2###############################################
#log2(FPKM+1) PCA on replicates #2
rep2_mtx_log2_fpkm_matrix_without_nonvariable_genes <- mtx_log2_fpkm_matrix_without_nonvariable_genes[,!grepl(pattern = "_RNA_1",x = colnames(mtx_log2_fpkm_matrix_without_nonvariable_genes),fixed = T)]
rep2_mtx_log2_fpkm_matrix_without_nonvariable_genes <- rep2_mtx_log2_fpkm_matrix_without_nonvariable_genes[which(rowSds(as.matrix(rep2_mtx_log2_fpkm_matrix_without_nonvariable_genes[,-1]),na.rm = T)!=0),]
rep2_pca <- prcomp(t(rep2_mtx_log2_fpkm_matrix_without_nonvariable_genes), scale. = T,center = T)
rep2_explained_variance <- rep2_pca$sdev^2 / sum(rep2_pca$sdev^2)
rep2_variance_df <- data.frame(
  PC = 1:length(rep2_explained_variance),
  Variance = rep2_explained_variance
)
#Scree Plot
ggplot(rep2_variance_df, aes(x = PC, y = Variance)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_line(group = 1, color = "darkblue") +
  geom_point(size = 3, color = "darkblue") +
  labs(title = "Scree Plot", x = "Principal Component", y = "Proportion of Variance") +
  theme_minimal()
#visualize first 2 PCs
rep2_pca_df <- as.data.frame(rep2_pca$x)
rep2_pca_df$Sample <- rownames(rep2_pca_df)
rep2_pca_df <- merge(rep2_pca_df, df_metadata_sensitivity_binary_clusters, by.x = "Sample", by.y = "id")

# Plot
ggplot(rep2_pca_df, aes(x = PC1, y = PC2, color = H2O2_bin_cluster, label = Sample)) +
  geom_point(size = 4, alpha = 0.8) +
  geom_text_repel(size = 3) +
  labs(title = "PCA Scatter Plot", x = paste0("PC1 (", round(rep2_explained_variance[1] * 100, 1), "%)"),
       y = paste0("PC2 (", round(rep2_explained_variance[2] * 100, 1), "%)")) +
  scale_color_brewer(palette = "Set2") +
  theme_minimal()
#View genes loadings
View(as.data.frame(rep2_pca$rotation))

#new count_data with gene names as first column
rep2_new_counts_data <- new_counts_data[,!grepl(pattern = "_RNA_1",x = colnames(new_counts_data),fixed = T)]

#Create DeSeq object
rep2_dds <- DESeqDataSetFromMatrix(countData=rep2_new_counts_data,colData=subset(df_metadata_sensitivity_binary_clusters,!grepl(pattern = "_RNA_1",x = id,fixed = T)), design=~H2O2_bin_cluster, tidy = TRUE)
rep2_dds
rep2_dds <- DESeq(rep2_dds)
rep2_res <- results(rep2_dds,alpha = 0.1)
rep2_res <- rep2_res[order(-abs(rep2_res$log2FoldChange)),]
df_rep2_res <- as.data.frame(rep2_res@listData)
df_rep2_res$Gene <- rep2_res@rownames
df_rep2_res$abs_log2FoldChange <- abs(df_rep2_res$log2FoldChange)
View(subset(df_rep2_res,padj<0.40))
summary(rep2_res)

#Volcano plot
df_rep2_res$Significance <- ifelse(df_rep2_res$padj<0.1 , yes = "Significant (FDR<0.1)",no = "Not Significant")
ggplot(df_rep2_res, aes(x = log2FoldChange, y = -log10(pvalue), color = Significance)) +
  geom_point(alpha = 0.8, size = 2) +  # Add points
  scale_color_manual(values = c("Significant (FDR<0.1)" = "red", "Not Significant" = "grey")) +  # Customize colors
  geom_hline(yintercept = -log10(0.1), linetype = "dashed", color = "blue", size = 0.5) +  # Add threshold line
  geom_vline(xintercept = c(-1, 0.585), linetype = "dashed", color = "blue", size = 0.5) +  # Add fold change thresholds
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

# 
# #reset par
# par(mfrow=c(1,1))
# # Make a basic volcano plot
# with(rep2_res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-3,5),ylim=c(0,8) ))
# 
# # Add colored points: blue if padj<0.1, red if log2FC>=0.585 (increase of >=1.5) or log2FC <= -1 (or absolute decrease of <= 50%) and padj<0.1)
# with(subset(rep2_res, padj< 0.1 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
# with(subset(rep2_res, padj<0.1 & (log2FoldChange <= -1 || log2FoldChange >= 0.585)), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))

#PYK1 differences in 2nd replicates
plotCounts(rep2_dds, gene="SPAC4H3.10c", intgroup="H2O2_bin_cluster")

#Differential expression H11 vs the rest
dds_h11_vs_others <- DESeqDataSetFromMatrix(countData=new_counts_data,colData=data.frame(is_H11=substr(df_metadata_sensitivity_binary_clusters$id,1,3)=="H11"), design=~is_H11, tidy = TRUE)
dds_h11_vs_others
dds_h11_vs_others <- DESeq(dds_h11_vs_others)
res_h11_vs_others <- results(dds_h11_vs_others,alpha = 0.1)
res_h11_vs_others <- res_h11_vs_others[order(res$padj),]
head(res_h11_vs_others)
summary(res_h11_vs_others)
View(as.data.frame(res_h11_vs_others))

#save session
library("session")
save.session(file = paste0(output_workspace,"H2O2_binary_sensitivity_clusters_Differential_expression_analysis.Rda"))
