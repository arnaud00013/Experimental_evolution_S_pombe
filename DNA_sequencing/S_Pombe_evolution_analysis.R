#@Author=Arnaud NG
#This script analyses S. Pombe samples (Substitution rate, Mutation bias, Mutation fixation time series, etc)

#import libraries
library("ggplot2")
library("seqinr")
library("RColorBrewer")
library("randomcoloR")
library("FD")
library("vegan")
library("gplots")
library("lmPerm")
library("ggpubr")
library("gridExtra")
library("tidyr")
#import script arguments
#workspace with all the .bam and .tab files
output_workspace <- "C:/Users/arnau/Documents/S_Pombe/"
#name of the depth report file 
depth_report_filename <- "common_depth_report.csv"
#name of the reference genome fasta file
fasta_refseq_filename <- "Schizosaccharomyces_pombe_all_chromosomes.fasta"

#Get list of samples
lst_samples <- unname(sapply(X = (sapply(X = read.csv2(file = paste0(output_workspace,"lst_samples.txt"),sep = "\t",header = FALSE,stringsAsFactors = FALSE)[1], FUN=function(x) gsub(pattern = output_workspace,replacement = "",x,fixed = TRUE))),FUN=function(x) gsub(pattern = "_preprocessed_sorted.bam",replacement = "",x,fixed = TRUE)))
lst_samples_original <- lst_samples
#calculate total number of samples
nb_samples <- length(lst_samples)
nb_samples_original <- length(lst_samples_original)
lst_unique_pops <- sort(unique(unname(vapply(X = lst_samples,FUN = function(x) substr(x,start = 1,stop=gregexpr(pattern = "_",text = x,fixed = T)[[1]][1]-1),FUN.VALUE = ""))))
nb_pops <- length(lst_unique_pops)
v_sample_to_unique_pop <- vapply(X = lst_samples,FUN = function(x) substr(x,start = 1,stop=gregexpr(pattern = "_",text = x,fixed = T)[[1]][1]-1),FUN.VALUE = "")
lst_possible_timepoints <- as.character(sort(unique(unname(vapply(X = lst_samples,FUN = function(x) as.integer(substr(x,start = gregexpr(pattern = "_",text = x,fixed = T)[[1]][1]+1,stop=gregexpr(pattern = "_",text = x,fixed = T)[[1]][2]-1)),FUN.VALUE = 0)))))
v_sample_to_timepoints <- vapply(X = lst_samples,FUN = function(x) as.integer(substr(x,start = gregexpr(pattern = "_",text = x,fixed = T)[[1]][1]+1,stop=gregexpr(pattern = "_",text = x,fixed = T)[[1]][2]-1)),FUN.VALUE = 0)
v_next_timepoint <- c(lst_possible_timepoints[2:length(lst_possible_timepoints)],lst_possible_timepoints[length(lst_possible_timepoints)])
names(v_next_timepoint) <- lst_possible_timepoints
v_previous_timepoint <- c(lst_possible_timepoints[1],lst_possible_timepoints[1:(length(lst_possible_timepoints)-1)])
names(v_previous_timepoint) <- lst_possible_timepoints
col_vector <- RColorBrewer::brewer.pal(length(lst_unique_pops),"Paired")
pie(1:length(col_vector), col=col_vector,labels = col_vector)
palette_pops <- col_vector
names(palette_pops) <- lst_unique_pops

#import reference fasta 
genome_refseq <- read.fasta(paste0(output_workspace,fasta_refseq_filename),seqtype = "DNA",as.string = TRUE,forceDNAtolower = FALSE)#seqinr::getSequence(object = toupper(),as.string = TRUE)[[1]]

v_chrs_length <- vapply(X = names(genome_refseq),FUN = function(x) nchar(genome_refseq[[x]]),FUN.VALUE = 0)
v_chrs_length <- v_chrs_length[names(v_chrs_length)!="mitochondrial"]
palette_chrs <- c("I"="red","II"="blue","chr_II_telomeric_gap"="green3","III"="grey","mating_type_region"="purple")
v_pop_order <- c("G1","G2","G3","G4","G5","G9","G11","H9","H11","H12")
v_pops_max_str_timepoint <- as.character(vapply(X = lst_unique_pops,FUN = function(x) max(unname(vapply(X = lst_samples[v_sample_to_unique_pop[lst_samples]==x],FUN = function(y) as.integer(substr(y,start = gregexpr(pattern = "_",text = y,fixed = T)[[1]][1]+1,stop=gregexpr(pattern = "_",text = y,fixed = T)[[1]][2]-1)),FUN.VALUE = 0))) ,FUN.VALUE = 0))
names(v_pops_max_str_timepoint) <- lst_unique_pops
v_pops_min_str_timepoint <- as.character(vapply(X = lst_unique_pops,FUN = function(x) min(unname(vapply(X = lst_samples[v_sample_to_unique_pop[lst_samples]==x],FUN = function(y) as.integer(substr(y,start = gregexpr(pattern = "_",text = y,fixed = T)[[1]][1]+1,stop=gregexpr(pattern = "_",text = y,fixed = T)[[1]][2]-1)),FUN.VALUE = 0))) ,FUN.VALUE = 0))
names(v_pops_min_str_timepoint) <- lst_unique_pops
#Set language as English for date formatting
Sys.setlocale("LC_ALL","English")

#function for plotting linear model
ggplotRegression <- function (fit,ggsave_path,the_filename,xlabl=NA,ylabl=NA) {
  library(ggplot2)
  bool_gg_save <- TRUE
  if(is.na(xlabl)){
    xlabl <- names(fit$model)[2]
  }
  if(is.na(ylabl)){
    ylabl <- names(fit$model)[1]
  }
  adj_r_sq <- formatC(summary(fit)$adj.r.squared, format = "e", digits = 3)
  slope <-formatC(summary(fit)$coefficients[,1][2], format = "e", digits = 3)
  p_val <- ifelse(test = "Iter"%in%colnames(summary(fit)$coefficients),yes=ifelse(test = unname(summary(fit)$coefficients[,3][2])<(2E-4),yes = "<2E-4",no = formatC(unname(summary(fit)$coefficients[,3][2]), format = "e", digits = 3)),no=ifelse(test = unname(summary(fit)$coefficients[,4][2])<(2e-16),yes = "<2e-16",no = formatC(unname(summary(fit)$coefficients[,4][2]), format = "e", digits = 3)))
  tryCatch(expr = {ggplot(fit$model, aes_string(x = names(fit$model)[2], y = names(fit$model)[1])) +
      geom_point() +
      stat_smooth(method = "lm", col = "red") +
      xlab(xlabl)+
      ylab(ylabl)+
      labs(title = paste("Adj R2 = ",adj_r_sq,
                         " Slope =",slope,
                         " P: ",p_val))+ theme(plot.title=element_text(hjust=0,size=12))},error=function(e) bool_gg_save <- FALSE)
  
  if (bool_gg_save){
    ggsave(filename = the_filename, path=ggsave_path, width = 15, height = 10, units = "cm")
  }else{
    print(paste0(the_filename, "won't be created because of it is irrelevant for gene in path ", ggsave_path))
  }
  #return result as the real float numbers
  adj_r_sq <- unname(summary(fit)$adj.r.squared)
  slope <-unname(summary(fit)$coefficients[,1][2])
  p_val <- ifelse(test = "Iter"%in%colnames(summary(fit)$coefficients),yes=ifelse(test = unname(summary(fit)$coefficients[,3][2])<(2E-4),yes = "<2E-4",no = unname(summary(fit)$coefficients[,3][2])),no=ifelse(test = unname(summary(fit)$coefficients[,4][2])<(2e-16),yes = "<2e-16",no = unname(summary(fit)$coefficients[,4][2])))
  return(list(adj_r_sq_current_lm = adj_r_sq,slope_current_lm = slope,p_val_current_lm=p_val))
}

#Create function to annotate sites
df_loci_annotations <- read.csv2(file = paste0(output_workspace,"Schizosaccharomyces_pombe_all_chromosomes.gff3"),sep = "\t",skip = 1,header = FALSE,stringsAsFactors = FALSE)
names(df_loci_annotations) <- c("chr","db","site_type","start","end","V6","strand","V8","annotation")
df_loci_annotations <- subset(df_loci_annotations,chr!="mitochondrial")
df_loci_annotations$the_product <- unname(vapply(X = df_loci_annotations$annotation, FUN = function(the_ann) ifelse(test = grepl(pattern = "Name=",x = the_ann,fixed = T),yes = substr(the_ann,start = gregexpr(pattern = "Name=",text = the_ann,fixed = T)[[1]]+5,stop = nchar(the_ann)),no = "NA"), FUN.VALUE = ""))
df_loci_annotations$the_product <- ifelse(df_loci_annotations$the_product=="NA",NA,df_loci_annotations$the_product)
df_loci_annotations$accession_number <- unname(vapply(X = df_loci_annotations$annotation, FUN = function(the_ann) ifelse(test = grepl(pattern = ";Name=",x = the_ann,fixed = T),yes = substr(the_ann,start = gregexpr(pattern = "ID=",text = the_ann,fixed = T)[[1]]+3,stop = gregexpr(pattern = ";Name=",text = the_ann,fixed = T)[[1]]-1),no = "NA"), FUN.VALUE = ""))
df_loci_annotations$accession_number <- ifelse(df_loci_annotations$the_product=="NA",NA,df_loci_annotations$accession_number)

find_mutation_variant_type2 <- function(the_row_ind_in_df_variants){
  current_variant_annotations_df <- subset(df_loci_annotations,(chr==df_variants$Chrom[the_row_ind_in_df_variants])&(start<= df_variants$Position[the_row_ind_in_df_variants])&(end >= df_variants$Position[the_row_ind_in_df_variants]))
  if (nrow(current_variant_annotations_df)>0){
    if (nrow(subset(current_variant_annotations_df,(grepl(pattern = "SPAC",x = annotation,fixed = T))&(site_type=="CDS")))>0){
      return("protein-coding")
    }else if (nrow(subset(current_variant_annotations_df,site_type=="intron"))>0){
      return("intron")
    }else if (nrow(subset(current_variant_annotations_df,site_type=="three_prime_UTR"))>0){
      return("3' UTR")
    }else if (nrow(subset(current_variant_annotations_df,site_type=="five_prime_UTR"))>0){
      return("5' UTR")
    }else if (nrow(subset(current_variant_annotations_df,site_type=="tRNA"))>0){
      return("tRNA")
    }else if (nrow(subset(current_variant_annotations_df,site_type=="rRNA"))>0){
      return("rRNA")
    }else if (nrow(subset(current_variant_annotations_df,site_type=="lncRNA"))>0){
      return("lncRNA")
    }else if (nrow(subset(current_variant_annotations_df,site_type=="sncRNA"))>0){
      return("sncRNA")
    }else if (nrow(subset(current_variant_annotations_df,site_type=="snRNA"))>0){
      return("snRNA")
    }else if (nrow(subset(current_variant_annotations_df,site_type=="snoRNA"))>0){
      return("snoRNA")
    }else{
      return("NA")
    }
  }else{
    return("Intergenic")
  }
}

find_mutation_comparison2_category <- function(the_row_ind_in_df_variants){
  current_variant_annotations_df <- subset(df_loci_annotations,(chr==df_variants$Chrom[the_row_ind_in_df_variants])&(start<= df_variants$Position[the_row_ind_in_df_variants])&(end >= df_variants$Position[the_row_ind_in_df_variants]))
  if (nrow(current_variant_annotations_df)>0){
    if (nrow(subset(current_variant_annotations_df,site_type=="pseudogenic_transcript"))>0){
      return("pseudogenic transcript")
    }else if (nrow(subset(current_variant_annotations_df,site_type=="long_terminal_repeat"))>0){
      return("long terminal repeat")
    }else if (nrow(subset(current_variant_annotations_df,site_type=="low_complexity_region"))>0){
      return("low complexity region")
    }else if (nrow(subset(current_variant_annotations_df,site_type=="promoter"))>0){
      return("promoter")
    }else if (nrow(subset(current_variant_annotations_df,site_type=="dg_repeat"))>0){
      return("dg repeat")
    }else if (nrow(subset(current_variant_annotations_df,site_type=="regional_centromere"))>0){
      return("regional centromere")
    }else if (nrow(subset(current_variant_annotations_df,site_type=="nuclear_mt_pseudogene"))>0){
      return("nuclear mitochondrial pseudogene")
    }else if (nrow(subset(current_variant_annotations_df,site_type=="gap"))>0){
      return("telomeric gap")
    }else if (nrow(subset(current_variant_annotations_df,site_type=="repeat_region"))>0){
      return("repeat region")
    }else if (nrow(subset(current_variant_annotations_df,site_type=="dh_repeat"))>0){
      return("dh repeat")
    }else if (nrow(subset(current_variant_annotations_df,site_type=="LTR_retrotransposon"))>0){
      return("LTR retrotransposon")
    }else if (nrow(subset(current_variant_annotations_df,site_type=="origin_of_replication"))>0){
      return("origin of replication")
    }else if (nrow(subset(current_variant_annotations_df,site_type=="regional_centromere_inner_repeat_region"))>0){
      return("RCIRR")
    }else if (nrow(subset(current_variant_annotations_df,site_type=="regional_centromere_central_core"))>0){
      return("RCCC")
    }else if (nrow(subset(current_variant_annotations_df,site_type=="TR_box"))>0){
      return("TR box")
    }else if (nrow(subset(current_variant_annotations_df,site_type=="mating_type_region"))>0){
      return("mating type region")
    }else{
      return("NA")
    }
  }else{
    return("Intergenic")
  }
}

find_gene_product_affected_by_mutation <- function(the_row_ind_in_df_variants){
  current_variant_annotations_df <- subset(df_loci_annotations,(chr==df_variants$Chrom[the_row_ind_in_df_variants])&(start<= df_variants$Position[the_row_ind_in_df_variants])&(end >= df_variants$Position[the_row_ind_in_df_variants])&(!is.na(the_product)))
  if (nrow(current_variant_annotations_df)>0){
    return(paste0(unique(current_variant_annotations_df$the_product),collapse = ";"))
  }else{
    return("NA")
  }
}

find_acession_number_gene_product_affected_by_mutation <- function(the_row_ind_in_df_variants){
  current_variant_annotations_df <- subset(df_loci_annotations,(chr==df_variants$Chrom[the_row_ind_in_df_variants])&(start<= df_variants$Position[the_row_ind_in_df_variants])&(end >= df_variants$Position[the_row_ind_in_df_variants])&(!is.na(accession_number)))
  if (nrow(current_variant_annotations_df)>0){
    return(paste0(unique(current_variant_annotations_df$accession_number),collapse = ";"))
  }else{
    return("NA")
  }
}

find_if_mutation_is_cis <- function(the_row_ind_in_df_ALL_variants){
  current_variant_annotations_df <- subset(df_loci_annotations,(!grepl(pattern = "noncoding",x = tolower(df_ALL_variants$variant_type2[the_row_ind_in_df_ALL_variants]),fixed = T))&(chr==df_ALL_variants$Chrom[the_row_ind_in_df_ALL_variants])&(df_ALL_variants$Position[the_row_ind_in_df_ALL_variants]>=start-1000)&(df_variants$Position[the_row_ind_in_df_ALL_variants]<=end+300)&(!is.na(the_product)))
  if (nrow(current_variant_annotations_df)>0){
    return(TRUE)
  }else{
    return(FALSE)
  }
}

find_potentially_cis_regulated_gene <- function(the_row_ind_in_df_ALL_variants){
  current_variant_annotations_df <- subset(df_loci_annotations,(!grepl(pattern = "noncoding",x = tolower(df_ALL_variants$variant_type2[the_row_ind_in_df_ALL_variants]),fixed = T))&(chr==df_ALL_variants$Chrom[the_row_ind_in_df_ALL_variants])&(df_ALL_variants$Position[the_row_ind_in_df_ALL_variants]>=start-1000)&(df_variants$Position[the_row_ind_in_df_ALL_variants]<=end+300)&(!is.na(the_product)))
  if (nrow(current_variant_annotations_df)>0){
    return(paste0(unique(current_variant_annotations_df$the_product),collapse = ";"))
  }else{
    return("NA")
  }
}

#import and concatenate samples variants calling file
df_variants <- read.csv2(file = paste0(output_workspace,"variant_analysis/variants_",lst_samples[1],".tab"),sep = "\t",header = TRUE,stringsAsFactors = FALSE)
if (nrow(df_variants)!=0){
  df_variants$Sample <- lst_samples[1]
}
for (i in 2:nb_samples){
  df_to_add_to_variants_df <- read.csv2(file = paste0(output_workspace,"variant_analysis/variants_",lst_samples[i],".tab"),sep = "\t",header = TRUE,stringsAsFactors = FALSE)
  if (nrow(df_to_add_to_variants_df)==0){
    next()
  }
  df_to_add_to_variants_df$Sample <- lst_samples[i]
  df_variants <- rbind(df_variants,df_to_add_to_variants_df)
  
}
#Do not confuse nucleotide T with the alias for the boolean value TRUE
df_variants$Ref <- gsub(pattern = "TRUE", replacement = "T",x = df_variants$Ref,fixed=TRUE) 
df_variants$VarAllele <- gsub(pattern = "TRUE", replacement = "T",x = df_variants$VarAllele,fixed=TRUE) 
df_variants$VarAllele <- gsub(pattern = "+", replacement = "",x = df_variants$VarAllele,fixed=TRUE) 
df_variants$VarAllele <- gsub(pattern = "-", replacement = "",x = df_variants$VarAllele,fixed=TRUE)
df_variants <- subset(df_variants,nchar(VarAllele)==1)
df_variants <- subset(df_variants,subset = !duplicated(paste0(Position,VarAllele,Sample,sep="")))
#Make sure to eliminate site with no SNP
df_variants <- subset(df_variants,!is.na(VarAllele))
#make sure to select nuclear chromosomes only
df_variants <- subset(df_variants,Chrom!="mitochondrial")
df_variants$VarFreq <- (as.numeric(gsub(pattern = "%",replacement = "",x = df_variants$VarFreq,fixed = TRUE))/100)
df_variants$is_fixed <- df_variants$VarFreq >= 0.75
#Annotate samples (mutation, population, time point, fixation, present at next time point, site_type and gene)
df_variants$mutation <- paste0(df_variants$Chrom,"_",df_variants$Ref,df_variants$Position,df_variants$VarAllele)
df_variants$population <- v_sample_to_unique_pop[df_variants$Sample]
df_variants$pop_mut <- paste0(df_variants$population,"_",df_variants$mutation)
df_variants$str_time <- as.character(v_sample_to_timepoints[df_variants$Sample])
df_variants$int_time <- as.integer(df_variants$str_time)
min_int_time <- min(df_variants$int_time, na.rm = T)
print(length(unique(subset(df_variants,is_fixed)$mutation)))
#make sure that the variant frequency filter is applied (at least 5 reads supports the variant and VarFreq>=0.05)
df_variants <- subset(df_variants, (Reads1+Reads2>=5)&(VarFreq>=0.05))
print(length(unique(subset(df_variants,is_fixed)$mutation)))
#apply strand bias filter on variants (2 reads supporting SNP on each strand)
df_variants <- df_variants[((df_variants$Reads2Plus>=1)&(df_variants$Reads2Minus>=1)),]
print(length(unique(subset(df_variants,is_fixed)$mutation)))
#remove mutations that are initially present at intermediate or high frequency in >=3 populations (THE INITIAL TIMEPOINT IS DIFFERENT FOR EACH POPULATION)
v_substitutions_at_initial_time <- unique(subset(df_variants,(str_time==unname(v_pops_min_str_timepoint[population])|int_time<3000)&VarFreq>=0.4)$mutation)
v_substitution_initial_prevalence <- vapply(X = v_substitutions_at_initial_time,FUN = function(x) length(unique(subset(df_variants,(str_time==unname(v_pops_min_str_timepoint[population]))&is_fixed&(mutation==x))$population)),FUN.VALUE = 0)
lst_substitutions_in_most_pops_at_min_time <- names(v_substitution_initial_prevalence)[v_substitution_initial_prevalence>=floor(nb_pops/3)]
df_variants <- subset(df_variants,!((mutation%in%lst_substitutions_in_most_pops_at_min_time)))
print(length(unique(subset(df_variants,is_fixed)$mutation)))
#remove variants that are absent in the next timepoint of the same population
#Define next and previous timepoint in a population
df_variants$pop_timepoint <- paste0(df_variants$population,"_",df_variants$str_time)
v_next_pop_timepoint <- paste0(rep(v_pop_order,each=length(v_next_timepoint)),"_",rep(unname(v_next_timepoint),times=length(v_pop_order)))
the_bool_next_timepoint <- v_next_pop_timepoint%in%df_variants$pop_timepoint
v_next_pop_timepoint <- v_next_pop_timepoint[the_bool_next_timepoint]
names(v_next_pop_timepoint) <- paste0(rep(v_pop_order,each=length(lst_possible_timepoints)),"_",rep(unname(lst_possible_timepoints),times=length(v_pop_order)))[the_bool_next_timepoint]
v_previous_pop_timepoint <- paste0(rep(v_pop_order,each=length(v_previous_timepoint)),"_",rep(unname(v_previous_timepoint),times=length(v_pop_order)))
the_bool_previous_timepoint <- v_previous_pop_timepoint%in%df_variants$pop_timepoint
v_previous_pop_timepoint <- v_previous_pop_timepoint[the_bool_previous_timepoint]
names(v_previous_pop_timepoint) <- paste0(rep(v_pop_order,each=length(lst_possible_timepoints)),"_",rep(unname(lst_possible_timepoints),times=length(v_pop_order)))[the_bool_previous_timepoint]
df_variants$is_present_at_adjacent_timepoint <- vapply(X = 1:nrow(df_variants),FUN = function(i) nrow(subset(df_variants,(pop_mut==df_variants$pop_mut[i]) & (((pop_timepoint==v_next_pop_timepoint[df_variants$pop_timepoint[i]])|(str_time==v_pops_max_str_timepoint[df_variants$population[i]]))|((pop_timepoint==v_previous_pop_timepoint[df_variants$pop_timepoint[i]])|(str_time==v_pops_min_str_timepoint[df_variants$population[i]])))))>0,FUN.VALUE = T) #
df_variants <- subset(df_variants,is_present_at_adjacent_timepoint)
print(length(unique(subset(df_variants,is_fixed)$mutation)))
#Get mean pop_mut frequency
v_lst_variants_pop_mut_mean_freq <- vapply(X = sort(unique(df_variants$pop_mut)),FUN = function(the_pop_mut) mean(subset(df_variants,pop_mut==the_pop_mut)$VarFreq,na.rm=T),FUN.VALUE = 0.0)
#Get initial pop_mut frequencies
v_lst_variants_pop_mut_init_freq <- vapply(X = sort(unique(df_variants$pop_mut)),FUN = function(the_pop_mut) subset(df_variants,pop_mut==the_pop_mut)$VarFreq[order(subset(df_variants,pop_mut==the_pop_mut)$int_time)][1],FUN.VALUE = 0.0)
#Get pop_muts frequency range
v_lst_variants_pop_mut_freq_range <- vapply(X = sort(unique(df_variants$pop_mut)),FUN = function(the_pop_mut) abs(diff(range(subset(df_variants,pop_mut==the_pop_mut)$VarFreq,na.rm=T))),FUN.VALUE = 0.0)
#Get pop mut that get fixed at some point
v_lst_variants_pop_mut_fixed_at_some_point <- sort(unique(subset(df_variants,is_fixed)$pop_mut))
#Remove pop mutation with a frequency range <=0.1 that never reach fixation
df_variants <- subset(df_variants,pop_mut%in%names(v_lst_variants_pop_mut_freq_range[v_lst_variants_pop_mut_freq_range>0.1])|(pop_mut%in%v_lst_variants_pop_mut_fixed_at_some_point))
print(length(unique(subset(df_variants,is_fixed)$mutation)))
##Measure autocorrelation
#v_lst_variants_pop_mut_autocorrelation <- vapply(X = sort(unique(df_variants$pop_mut)),FUN = function(the_pop_mut) max(acf(subset(df_variants,pop_mut==the_pop_mut)$VarFreq[order(subset(df_variants,pop_mut==the_pop_mut)$int_time)],plot = F)$acf[acf(subset(df_variants,pop_mut==the_pop_mut)$VarFreq[order(subset(df_variants,pop_mut==the_pop_mut)$int_time)],plot = F)$acf<1],na.rm=T),FUN.VALUE = 0.0)
##Only keep variants that have enough lag1 autocorrelation (acf >=0.2)
#df_variants <- subset(df_variants,pop_mut%in%names(v_lst_variants_pop_mut_autocorrelation[v_lst_variants_pop_mut_autocorrelation>=0.2]))
#print(length(unique(subset(df_variants,is_fixed)$mutation)))
#add gene and loci type annotation to variants dataframe
df_variants$gene_product <- vapply(X = 1:nrow(df_variants),FUN = find_gene_product_affected_by_mutation,FUN.VALUE = "")
df_variants$gene_product <- ifelse(df_variants$gene_product=="NA",yes=NA,no=df_variants$gene_product)
df_variants$accession_number <- vapply(X = 1:nrow(df_variants),FUN = find_acession_number_gene_product_affected_by_mutation,FUN.VALUE = "")
df_variants$accession_number <- ifelse(df_variants$accession_number=="NA",yes=NA,no=df_variants$accession_number)
df_variants$variant_type2 <- vapply(X = 1:nrow(df_variants),FUN = find_mutation_variant_type2,FUN.VALUE = "")
df_variants$variant_type2 <- ifelse(df_variants$variant_type2=="NA",yes=NA,no=df_variants$variant_type2)
df_variants$comparison2_category <- vapply(X = 1:nrow(df_variants),FUN = find_mutation_comparison2_category,FUN.VALUE = "")
df_variants$comparison2_category <- ifelse(df_variants$comparison2_category=="NA",yes=NA,no=df_variants$comparison2_category)
df_variants$variant_type <- "SNV"
#determine the type of SNVs
df_variants$base_change <- paste0(df_variants$Ref,"->",df_variants$VarAllele)
v_base_change_to_SNV_category <- c("A->C"="AT->CG","T->G"="AT->CG","A->T"="AT->TA","T->A"="AT->TA","C->G"="CG->GC","G->C"="CG->GC","G->A"="GC->AT","C->T"="GC->AT","G->T"="GC->TA","C->A"="GC->TA","T->C"="TA->CG","A->G"="TA->CG")
df_variants$SNV_category <- v_base_change_to_SNV_category[df_variants$base_change]
df_variants <- subset(df_variants,!is.na(SNV_category))
# #calculate the length of annotated genomic regions
v_length_annotated_genomic_regions <- readRDS(paste0(output_workspace,"v_length_annotated_genomic_regions.rds"))
# v_unique_annotated_sites <- NULL
# v_length_annotated_genomic_regions <- c("protein-coding"=0,"intron"=0,"3' UTR"=0,"5' UTR"=0,"tRNA"=0,"rRNA"=0,"lncRNA"=0,"sncRNA"=0,"snRNA"=0,"snoRNA"=0,"pseudogenic transcript"=0,"long terminal repeat"=0,"low complexity region"=0,"promoter"=0,"dg repeat"=0,"regional centromere"=0,"nuclear mitochondrial pseudogene"=0,"telomeric gap"=0,"repeat region"=0,"dh repeat"=0,"LTR retrotransposon"=0,"origin of replication"=0,"RCIRR"=0,"RCCC"=0,"TR box"=0,"mating type region"=0,"others"=0)
# for (i in 1:nrow(df_loci_annotations)){
#   if ((grepl(pattern = "SPAC",x = df_loci_annotations$annotation[i],fixed = T))&(df_loci_annotations$site_type[i]=="CDS")){
#     v_length_annotated_genomic_regions["protein-coding"] <- v_length_annotated_genomic_regions["protein-coding"]+ (df_loci_annotations$end[i]-df_loci_annotations$start[i]+1)
#   }else if (df_loci_annotations$site_type[i]=="intron"){
#     v_length_annotated_genomic_regions["intron"] <- v_length_annotated_genomic_regions["intron"]+ (df_loci_annotations$end[i]-df_loci_annotations$start[i]+1)
#   }else if (df_loci_annotations$site_type[i]=="three_prime_UTR"){
#     v_length_annotated_genomic_regions["3' UTR"] <- v_length_annotated_genomic_regions["3' UTR"]+ (df_loci_annotations$end[i]-df_loci_annotations$start[i]+1)
#   }else if (df_loci_annotations$site_type[i]=="five_prime_UTR"){
#     v_length_annotated_genomic_regions["5' UTR"] <- v_length_annotated_genomic_regions["5' UTR"]+ (df_loci_annotations$end[i]-df_loci_annotations$start[i]+1)
#   }else if (df_loci_annotations$site_type[i]=="tRNA"){
#     v_length_annotated_genomic_regions["tRNA"] <- v_length_annotated_genomic_regions["tRNA"]+ (df_loci_annotations$end[i]-df_loci_annotations$start[i]+1)
#   }else if (df_loci_annotations$site_type[i]=="rRNA"){
#     v_length_annotated_genomic_regions["rRNA"] <- v_length_annotated_genomic_regions["rRNA"]+ (df_loci_annotations$end[i]-df_loci_annotations$start[i]+1)
#   }else if (df_loci_annotations$site_type[i]=="lncRNA"){
#     v_length_annotated_genomic_regions["lncRNA"] <- v_length_annotated_genomic_regions["lncRNA"]+ (df_loci_annotations$end[i]-df_loci_annotations$start[i]+1)
#   }else if (df_loci_annotations$site_type[i]=="sncRNA"){
#     v_length_annotated_genomic_regions["sncRNA"] <- v_length_annotated_genomic_regions["sncRNA"]+ (df_loci_annotations$end[i]-df_loci_annotations$start[i]+1)
#   }else if (df_loci_annotations$site_type[i]=="snRNA"){
#     v_length_annotated_genomic_regions["snRNA"] <- v_length_annotated_genomic_regions["snRNA"]+ (df_loci_annotations$end[i]-df_loci_annotations$start[i]+1)
#   }else if (df_loci_annotations$site_type[i]=="snoRNA"){
#     v_length_annotated_genomic_regions["snoRNA"] <- v_length_annotated_genomic_regions["snoRNA"]+ (df_loci_annotations$end[i]-df_loci_annotations$start[i]+1)
#   }else if (df_loci_annotations$site_type[i]=="pseudogenic_transcript"){
#     v_length_annotated_genomic_regions["pseudogenic transcript"] <- v_length_annotated_genomic_regions["pseudogenic transcript"]+ (df_loci_annotations$end[i]-df_loci_annotations$start[i]+1)
#   }else if (df_loci_annotations$site_type[i]=="long_terminal_repeat"){
#     v_length_annotated_genomic_regions["long terminal repeat"] <- v_length_annotated_genomic_regions["long terminal repeat"]+ (df_loci_annotations$end[i]-df_loci_annotations$start[i]+1)
#   }else if (df_loci_annotations$site_type[i]=="low_complexity_region"){
#     v_length_annotated_genomic_regions["low complexity region"] <- v_length_annotated_genomic_regions["low complexity region"]+ (df_loci_annotations$end[i]-df_loci_annotations$start[i]+1)
#   }else if (df_loci_annotations$site_type[i]=="promoter"){
#     v_length_annotated_genomic_regions["promoter"] <- v_length_annotated_genomic_regions["promoter"]+ (df_loci_annotations$end[i]-df_loci_annotations$start[i]+1)
#   }else if (df_loci_annotations$site_type[i]=="dg_repeat"){
#     v_length_annotated_genomic_regions["dg repeat"] <- v_length_annotated_genomic_regions["dg repeat"]+ (df_loci_annotations$end[i]-df_loci_annotations$start[i]+1)
#   }else if (df_loci_annotations$site_type[i]=="regional_centromere"){
#     v_length_annotated_genomic_regions["regional centromere"] <- v_length_annotated_genomic_regions["regional centromere"]+ (df_loci_annotations$end[i]-df_loci_annotations$start[i]+1)
#   }else if (df_loci_annotations$site_type[i]=="nuclear_mt_pseudogene"){
#     v_length_annotated_genomic_regions["nuclear mitochondrial pseudogene"] <- v_length_annotated_genomic_regions["nuclear mitochondrial pseudogene"]+ (df_loci_annotations$end[i]-df_loci_annotations$start[i]+1)
#   }else if (df_loci_annotations$site_type[i]=="gap"){
#     v_length_annotated_genomic_regions["telomeric gap"] <- v_length_annotated_genomic_regions["telomeric gap"]+ (df_loci_annotations$end[i]-df_loci_annotations$start[i]+1)
#   }else if (df_loci_annotations$site_type[i]=="repeat_region"){
#     v_length_annotated_genomic_regions["repeat region"] <- v_length_annotated_genomic_regions["repeat region"]+ (df_loci_annotations$end[i]-df_loci_annotations$start[i]+1)
#   }else if (df_loci_annotations$site_type[i]=="dh_repeat"){
#     v_length_annotated_genomic_regions["dh repeat"] <- v_length_annotated_genomic_regions["dh repeat"]+ (df_loci_annotations$end[i]-df_loci_annotations$start[i]+1)
#   }else if (df_loci_annotations$site_type[i]=="LTR_retrotransposon"){
#     v_length_annotated_genomic_regions["LTR retrotransposon"] <- v_length_annotated_genomic_regions["LTR retrotransposon"]+ (df_loci_annotations$end[i]-df_loci_annotations$start[i]+1)
#   }else if (df_loci_annotations$site_type[i]=="origin_of_replication"){
#     v_length_annotated_genomic_regions["origin of replication"] <- v_length_annotated_genomic_regions["origin of replication"]+ (df_loci_annotations$end[i]-df_loci_annotations$start[i]+1)
#   }else if (df_loci_annotations$site_type[i]=="regional_centromere_inner_repeat_region"){
#     v_length_annotated_genomic_regions["RCIRR"] <- v_length_annotated_genomic_regions["RCIRR"]+ (df_loci_annotations$end[i]-df_loci_annotations$start[i]+1)
#   }else if (df_loci_annotations$site_type[i]=="regional_centromere_central_core"){
#     v_length_annotated_genomic_regions["RCCC"] <- v_length_annotated_genomic_regions["RCCC"]+ (df_loci_annotations$end[i]-df_loci_annotations$start[i]+1)
#   }else if (df_loci_annotations$site_type[i]=="TR_box"){
#     v_length_annotated_genomic_regions["TR box"] <- v_length_annotated_genomic_regions["TR box"]+ (df_loci_annotations$end[i]-df_loci_annotations$start[i]+1)
#   }else if (df_loci_annotations$site_type[i]=="mating_type_region"){
#     v_length_annotated_genomic_regions["mating type region"] <- v_length_annotated_genomic_regions["mating type region"]+ (df_loci_annotations$end[i]-df_loci_annotations$start[i]+1)
#   }
#   v_unique_annotated_sites <- unique(c(v_unique_annotated_sites, paste0(df_loci_annotations$chr[i],"_",df_loci_annotations$start[i]:df_loci_annotations$end[i])))
#   if (i%%5000==0){
#     print(i)
#   }
# }
#   #the length of intergenic regions should be the number of unannotated bases (practical definition)
# v_unique_annotated_sites <- unique(v_unique_annotated_sites)
# v_length_annotated_genomic_regions <- c(v_length_annotated_genomic_regions,"Intergenic"=(sum(unname(v_chrs_length))-length(v_unique_annotated_sites)))
# remove(v_unique_annotated_sites)


#get the dataframe of indels df_indel_variants
find_indel_variant_type2 <- function(the_row_ind_in_df_indel_variants){
  if (df_indel_variants$variant_type[the_row_ind_in_df_indel_variants]=="Deletion"){
    v_out <- rep("NA",nchar(df_indel_variants$VarAllele[the_row_ind_in_df_indel_variants])-1)
    k <- 1
    for (the_deleted_pos_in_chr in (df_indel_variants$Position[the_row_ind_in_df_indel_variants]:(df_indel_variants$Position[the_row_ind_in_df_indel_variants]+nchar(df_indel_variants$VarAllele[the_row_ind_in_df_indel_variants])-1-1))){
      current_variant_annotations_df <- subset(df_loci_annotations,(chr==df_indel_variants$Chrom[the_row_ind_in_df_indel_variants])&(start<= the_deleted_pos_in_chr)&(end >= the_deleted_pos_in_chr))
      if (nrow(current_variant_annotations_df)>0){
        if (nrow(subset(current_variant_annotations_df,(grepl(pattern = "SPAC",x = annotation,fixed = T))&(site_type=="CDS")))>0){
          v_out[k] <- "protein-coding"
        }else if (nrow(subset(current_variant_annotations_df,site_type=="intron"))>0){
          v_out[k] <- "intron"
        }else if (nrow(subset(current_variant_annotations_df,site_type=="three_prime_UTR"))>0){
          v_out[k] <- "3' UTR"
        }else if (nrow(subset(current_variant_annotations_df,site_type=="five_prime_UTR"))>0){
          v_out[k] <- "5' UTR"
        }else if (nrow(subset(current_variant_annotations_df,site_type=="tRNA"))>0){
          v_out[k] <- "tRNA"
        }else if (nrow(subset(current_variant_annotations_df,site_type=="rRNA"))>0){
          v_out[k] <- "rRNA"
        }else if (nrow(subset(current_variant_annotations_df,site_type=="lncRNA"))>0){
          v_out[k] <- "lncRNA"
        }else if (nrow(subset(current_variant_annotations_df,site_type=="sncRNA"))>0){
          v_out[k] <- "sncRNA"
        }else if (nrow(subset(current_variant_annotations_df,site_type=="snRNA"))>0){
          v_out[k] <- "snRNA"
        }else if (nrow(subset(current_variant_annotations_df,site_type=="snoRNA"))>0){
          v_out[k] <- "snoRNA"
        }else{
          v_out[k] <- "NA"
        }
      }else{
        v_out[k] <- "Intergenic"
      }
      k <- k + 1
    }
    v_out <- paste0(unique(v_out),collapse = ";")
  }else{
    the_deleted_pos_in_chr <- df_indel_variants$Position[the_row_ind_in_df_indel_variants]
    current_variant_annotations_df <- subset(df_loci_annotations,(chr==df_indel_variants$Chrom[the_row_ind_in_df_indel_variants])&(start<= the_deleted_pos_in_chr)&(end >= the_deleted_pos_in_chr))
    if (nrow(current_variant_annotations_df)>0){
      if (nrow(subset(current_variant_annotations_df,(grepl(pattern = "SPAC",x = annotation,fixed = T))&(site_type=="CDS")))>0){
        v_out <- "protein-coding"
      }else if (nrow(subset(current_variant_annotations_df,site_type=="intron"))>0){
        v_out <- "intron"
      }else if (nrow(subset(current_variant_annotations_df,site_type=="three_prime_UTR"))>0){
        v_out <- "3' UTR"
      }else if (nrow(subset(current_variant_annotations_df,site_type=="five_prime_UTR"))>0){
        v_out <- "5' UTR"
      }else if (nrow(subset(current_variant_annotations_df,site_type=="tRNA"))>0){
        v_out <- "tRNA"
      }else if (nrow(subset(current_variant_annotations_df,site_type=="rRNA"))>0){
        v_out <- "rRNA"
      }else if (nrow(subset(current_variant_annotations_df,site_type=="lncRNA"))>0){
        v_out <- "lncRNA"
      }else if (nrow(subset(current_variant_annotations_df,site_type=="sncRNA"))>0){
        v_out <- "sncRNA"
      }else if (nrow(subset(current_variant_annotations_df,site_type=="snRNA"))>0){
        v_out <- "snRNA"
      }else if (nrow(subset(current_variant_annotations_df,site_type=="snoRNA"))>0){
        v_out <- "snoRNA"
      }else{
        v_out <- "NA"
      }
    }else{
      v_out <- "Intergenic"
    }
  }
  return(v_out)
}

find_indel_comparison2_category <- function(the_row_ind_in_df_indel_variants){
  if (df_indel_variants$variant_type[the_row_ind_in_df_indel_variants]=="Deletion"){
    v_out <- rep("NA",nchar(df_indel_variants$VarAllele[the_row_ind_in_df_indel_variants])-1)
    k <- 1
    for (the_deleted_pos_in_chr in (df_indel_variants$Position[the_row_ind_in_df_indel_variants]:(df_indel_variants$Position[the_row_ind_in_df_indel_variants]+nchar(df_indel_variants$VarAllele[the_row_ind_in_df_indel_variants])-1-1))){
      current_variant_annotations_df <- subset(df_loci_annotations,(chr==df_indel_variants$Chrom[the_row_ind_in_df_indel_variants])&(start<= df_indel_variants$Position[the_row_ind_in_df_indel_variants])&(end >= df_indel_variants$Position[the_row_ind_in_df_indel_variants]))
      if (nrow(current_variant_annotations_df)>0){
        if (nrow(subset(current_variant_annotations_df,site_type=="pseudogenic_transcript"))>0){
          v_out[k] <- "pseudogenic transcript"
        }else if (nrow(subset(current_variant_annotations_df,site_type=="long_terminal_repeat"))>0){
          v_out[k] <- "long terminal repeat"
        }else if (nrow(subset(current_variant_annotations_df,site_type=="low_complexity_region"))>0){
          v_out[k] <- "low complexity region"
        }else if (nrow(subset(current_variant_annotations_df,site_type=="promoter"))>0){
          v_out[k] <- "promoter"
        }else if (nrow(subset(current_variant_annotations_df,site_type=="dg_repeat"))>0){
          v_out[k] <- "dg repeat"
        }else if (nrow(subset(current_variant_annotations_df,site_type=="regional_centromere"))>0){
          v_out[k] <- "regional centromere"
        }else if (nrow(subset(current_variant_annotations_df,site_type=="nuclear_mt_pseudogene"))>0){
          v_out[k] <- "nuclear mitochondrial pseudogene"
        }else if (nrow(subset(current_variant_annotations_df,site_type=="gap"))>0){
          v_out[k] <- "telomeric gap"
        }else if (nrow(subset(current_variant_annotations_df,site_type=="repeat_region"))>0){
          v_out[k] <- "repeat region"
        }else if (nrow(subset(current_variant_annotations_df,site_type=="dh_repeat"))>0){
          v_out[k] <- "dh repeat"
        }else if (nrow(subset(current_variant_annotations_df,site_type=="LTR_retrotransposon"))>0){
          v_out[k] <- "LTR retrotransposon"
        }else if (nrow(subset(current_variant_annotations_df,site_type=="origin_of_replication"))>0){
          v_out[k] <- "origin of replication"
        }else if (nrow(subset(current_variant_annotations_df,site_type=="regional_centromere_inner_repeat_region"))>0){
          v_out[k] <- "RCIRR"
        }else if (nrow(subset(current_variant_annotations_df,site_type=="regional_centromere_central_core"))>0){
          v_out[k] <- "RCCC"
        }else if (nrow(subset(current_variant_annotations_df,site_type=="TR_box"))>0){
          v_out[k] <- "TR box"
        }else if (nrow(subset(current_variant_annotations_df,site_type=="mating_type_region"))>0){
          v_out[k] <- "mating type region"
        }else{
          v_out[k] <- "NA"
        }
      }else{
        v_out[k] <- "Intergenic"
      }
      k <- k + 1
    }
    v_out <- paste0(unique(v_out),collapse = ";")
  }else{
    current_variant_annotations_df <- subset(df_loci_annotations,(chr==df_indel_variants$Chrom[the_row_ind_in_df_indel_variants])&(start<= df_indel_variants$Position[the_row_ind_in_df_indel_variants])&(end >= df_indel_variants$Position[the_row_ind_in_df_indel_variants]))
    if (nrow(current_variant_annotations_df)>0){
      if (nrow(subset(current_variant_annotations_df,site_type=="pseudogenic_transcript"))>0){
        v_out <- "pseudogenic transcript"
      }else if (nrow(subset(current_variant_annotations_df,site_type=="long_terminal_repeat"))>0){
        v_out <- "long terminal repeat"
      }else if (nrow(subset(current_variant_annotations_df,site_type=="low_complexity_region"))>0){
        v_out <- "low complexity region"
      }else if (nrow(subset(current_variant_annotations_df,site_type=="promoter"))>0){
        v_out <- "promoter"
      }else if (nrow(subset(current_variant_annotations_df,site_type=="dg_repeat"))>0){
        v_out <- "dg repeat"
      }else if (nrow(subset(current_variant_annotations_df,site_type=="regional_centromere"))>0){
        v_out <- "regional centromere"
      }else if (nrow(subset(current_variant_annotations_df,site_type=="nuclear_mt_pseudogene"))>0){
        v_out <- "nuclear mitochondrial pseudogene"
      }else if (nrow(subset(current_variant_annotations_df,site_type=="gap"))>0){
        v_out <- "telomeric gap"
      }else if (nrow(subset(current_variant_annotations_df,site_type=="repeat_region"))>0){
        v_out <- "repeat region"
      }else if (nrow(subset(current_variant_annotations_df,site_type=="dh_repeat"))>0){
        v_out <- "dh repeat"
      }else if (nrow(subset(current_variant_annotations_df,site_type=="LTR_retrotransposon"))>0){
        v_out <- "LTR retrotransposon"
      }else if (nrow(subset(current_variant_annotations_df,site_type=="origin_of_replication"))>0){
        v_out <- "origin of replication"
      }else if (nrow(subset(current_variant_annotations_df,site_type=="regional_centromere_inner_repeat_region"))>0){
        v_out <- "RCIRR"
      }else if (nrow(subset(current_variant_annotations_df,site_type=="regional_centromere_central_core"))>0){
        v_out <- "RCCC"
      }else if (nrow(subset(current_variant_annotations_df,site_type=="TR_box"))>0){
        v_out <- "TR box"
      }else if (nrow(subset(current_variant_annotations_df,site_type=="mating_type_region"))>0){
        v_out <- "mating type region"
      }else{
        v_out <- "NA"
      }
    }else{
      v_out <- "Intergenic"
    }
  }
  return(v_out)
}

find_gene_product_affected_by_indel <- function(the_row_ind_in_df_indel_variants){
  if (df_indel_variants$variant_type[the_row_ind_in_df_indel_variants]=="Deletion"){
    v_out <- rep("NA",nchar(df_indel_variants$VarAllele[the_row_ind_in_df_indel_variants])-1)
    k <- 1
    for (the_deleted_pos_in_chr in (df_indel_variants$Position[the_row_ind_in_df_indel_variants]:(df_indel_variants$Position[the_row_ind_in_df_indel_variants]+nchar(df_indel_variants$VarAllele[the_row_ind_in_df_indel_variants])-1-1))){
      current_variant_annotations_df <- subset(df_loci_annotations,(chr==df_indel_variants$Chrom[the_row_ind_in_df_indel_variants])&(start<= df_indel_variants$Position[the_row_ind_in_df_indel_variants])&(end >= df_indel_variants$Position[the_row_ind_in_df_indel_variants])&(!is.na(the_product)))
      if (nrow(current_variant_annotations_df)>0){
        v_out[k] <- paste0(unique(current_variant_annotations_df$the_product),collapse = ";")
      }else{
        v_out[k] <- "NA"
      }
      k <- k + 1
    }
    v_out <- paste0(unique(v_out),collapse = ";")
  }else{
    current_variant_annotations_df <- subset(df_loci_annotations,(chr==df_indel_variants$Chrom[the_row_ind_in_df_indel_variants])&(start<= df_indel_variants$Position[the_row_ind_in_df_indel_variants])&(end >= df_indel_variants$Position[the_row_ind_in_df_indel_variants])&(!is.na(the_product)))
    if (nrow(current_variant_annotations_df)>0){
      v_out <- paste0(unique(current_variant_annotations_df$the_product),collapse = ";")
    }else{
      v_out <- "NA"
    }
  }
  return(v_out)
}

find_acession_number_gene_product_affected_by_indel<- function(the_row_ind_in_df_indel_variants){
  if (df_indel_variants$variant_type[the_row_ind_in_df_indel_variants]=="Deletion"){
    v_out <- rep("NA",nchar(df_indel_variants$VarAllele[the_row_ind_in_df_indel_variants])-1)
    k <- 1
    for (the_deleted_pos_in_chr in (df_indel_variants$Position[the_row_ind_in_df_indel_variants]:(df_indel_variants$Position[the_row_ind_in_df_indel_variants]+nchar(df_indel_variants$VarAllele[the_row_ind_in_df_indel_variants])-1-1))){
      current_variant_annotations_df <- subset(df_loci_annotations,(chr==df_indel_variants$Chrom[the_row_ind_in_df_indel_variants])&(start<= df_indel_variants$Position[the_row_ind_in_df_indel_variants])&(end >= df_indel_variants$Position[the_row_ind_in_df_indel_variants])&(!is.na(the_product)))
      if (nrow(current_variant_annotations_df)>0){
        v_out[k] <- paste0(unique(current_variant_annotations_df$accession_number),collapse = ";")
      }else{
        v_out[k] <- "NA"
      }
      k <- k + 1
    }
    v_out <- paste0(unique(v_out),collapse = ";")
  }else{
    current_variant_annotations_df <- subset(df_loci_annotations,(chr==df_indel_variants$Chrom[the_row_ind_in_df_indel_variants])&(start<= df_indel_variants$Position[the_row_ind_in_df_indel_variants])&(end >= df_indel_variants$Position[the_row_ind_in_df_indel_variants])&(!is.na(the_product)))
    if (nrow(current_variant_annotations_df)>0){
      v_out <- paste0(unique(current_variant_annotations_df$accession_number),collapse = ";")
    }else{
      v_out <- "NA"
    }
  }
  return(v_out)
}

#import and concatenate samples variants calling file
df_indel_variants <- read.csv2(file = paste0(output_workspace,"variant_analysis/indels_variants_",lst_samples[1],".tab"),sep = "\t",header = TRUE,stringsAsFactors = FALSE)
if (nrow(df_indel_variants)!=0){
  df_indel_variants$Sample <- lst_samples[1]
}
for (i in 2:nb_samples){
  df_to_add_to_variants_df <- read.csv2(file = paste0(output_workspace,"variant_analysis/indels_variants_",lst_samples[i],".tab"),sep = "\t",header = TRUE,stringsAsFactors = FALSE)
  if (nrow(df_to_add_to_variants_df)==0){
    next()
  }
  df_to_add_to_variants_df$Sample <- lst_samples[i]
  df_indel_variants <- rbind(df_indel_variants,df_to_add_to_variants_df)
  
}
#Do not confuse nucleotide T with the alias for the boolean value TRUE
df_indel_variants$Ref <- gsub(pattern = "TRUE", replacement = "T",x = df_indel_variants$Ref,fixed=TRUE) 
df_indel_variants <- subset(df_indel_variants,subset = !duplicated(paste0(Position,VarAllele,Sample,sep="")))
#Make sure to eliminate site with no SNP
df_indel_variants <- subset(df_indel_variants,!is.na(VarAllele))
#make sure to select nuclear chromosomes only
df_indel_variants <- subset(df_indel_variants,Chrom!="mitochondrial")
df_indel_variants$VarFreq <- (as.numeric(gsub(pattern = "%",replacement = "",x = df_indel_variants$VarFreq,fixed = TRUE))/100)
df_indel_variants$is_fixed <- df_indel_variants$VarFreq >= 0.75
#Annotate samples (mutation, population, time point, fixation, present at next time point, site_type and gene)
df_indel_variants$mutation <- paste0(df_indel_variants$Chrom,"_",df_indel_variants$Ref,df_indel_variants$Position,df_indel_variants$VarAllele)
df_indel_variants$population <- v_sample_to_unique_pop[df_indel_variants$Sample]
df_indel_variants$pop_mut <- paste0(df_indel_variants$population,"_",df_indel_variants$mutation)
df_indel_variants$str_time <- as.character(v_sample_to_timepoints[df_indel_variants$Sample])
df_indel_variants$int_time <- as.integer(df_indel_variants$str_time)
print(length(unique(subset(df_indel_variants,is_fixed)$mutation)))
#make sure that the variant frequency filter is applied (at least 5 reads support the indel and VarFreq >=5%)
df_indel_variants <- subset(df_indel_variants, (Reads1+Reads2>=5)&(VarFreq>=0.05))
print(length(unique(subset(df_indel_variants,is_fixed)$mutation)))
#apply strand bias filter on variants (2 reads supporting SNP on each strand)
df_indel_variants <- df_indel_variants[((df_indel_variants$Reads2Plus>=1)&(df_indel_variants$Reads2Minus>=1)),]
print(length(unique(subset(df_indel_variants,is_fixed)$mutation)))
#remove indels that are initially present at intermediate or high frequency in >=3 populations (THE INITIAL TIMEPOINT IS DIFFERENT FOR EACH POPULATION)
v_fixed_indels_at_initial_time <- unique(subset(df_indel_variants,(str_time==unname(v_pops_min_str_timepoint[population])|int_time<3000)&VarFreq>=0.4)$mutation)
v_fixed_indel_initial_prevalence <- vapply(X = v_fixed_indels_at_initial_time,FUN = function(x) length(unique(subset(df_indel_variants,(str_time==unname(v_pops_min_str_timepoint[population]))&is_fixed&(mutation==x))$population)),FUN.VALUE = 0)
lst_fixed_indels_in_most_pops_at_min_time <- names(v_fixed_indel_initial_prevalence)[v_fixed_indel_initial_prevalence>=floor(nb_pops/3)]
df_indel_variants <- subset(df_indel_variants,!((mutation%in%lst_fixed_indels_in_most_pops_at_min_time)))
print(length(unique(subset(df_indel_variants,is_fixed)$mutation)))
#remove indels that are absent in the next timepoint of the same population
#Define next and previous timepoint in a population
df_indel_variants$pop_timepoint <- paste0(df_indel_variants$population,"_",df_indel_variants$str_time)
df_indel_variants$is_present_at_adjacent_timepoint <- vapply(X = 1:nrow(df_indel_variants),FUN = function(i) nrow(subset(df_indel_variants,(pop_mut==df_indel_variants$pop_mut[i]) & (((pop_timepoint==v_next_pop_timepoint[df_indel_variants$pop_timepoint[i]])|(str_time==v_pops_max_str_timepoint[df_indel_variants$population[i]]))|((pop_timepoint==v_previous_pop_timepoint[df_indel_variants$pop_timepoint[i]])|(str_time==v_pops_min_str_timepoint[df_indel_variants$population[i]])))))>0,FUN.VALUE = T) #
df_indel_variants <- subset(df_indel_variants,is_present_at_adjacent_timepoint)
print(length(unique(subset(df_indel_variants,is_fixed)$mutation)))
#Get mean pop_mut frequency
v_lst_indels_pop_mut_mean_freq <- vapply(X = sort(unique(df_indel_variants$pop_mut)),FUN = function(the_pop_mut) mean(subset(df_indel_variants,pop_mut==the_pop_mut)$VarFreq,na.rm=T),FUN.VALUE = 0.0)
#Get initial pop_mut frequencies
v_lst_indels_pop_mut_init_freq <- vapply(X = sort(unique(df_indel_variants$pop_mut)),FUN = function(the_pop_mut) subset(df_indel_variants,pop_mut==the_pop_mut)$VarFreq[order(subset(df_indel_variants,pop_mut==the_pop_mut)$int_time)][1],FUN.VALUE = 0.0)
#Get pop_muts frequency range
v_lst_indels_pop_mut_freq_range <- vapply(X = sort(unique(df_indel_variants$pop_mut)),FUN = function(the_pop_mut) abs(diff(range(subset(df_indel_variants,pop_mut==the_pop_mut)$VarFreq,na.rm=T))),FUN.VALUE = 0.0)
#Get pop mut that get fixed at some point
v_lst_pop_indels_fixed_at_some_point <- sort(unique(subset(df_indel_variants,is_fixed)$pop_mut))
#Remove pop indels with a frequency range <=0.1 that never reach fixation
df_indel_variants <- subset(df_indel_variants,pop_mut%in%names(v_lst_indels_pop_mut_freq_range[v_lst_indels_pop_mut_freq_range>0.1])|(pop_mut%in%v_lst_pop_indels_fixed_at_some_point))
print(length(unique(subset(df_indel_variants,is_fixed)$mutation)))
##Measure autocorrelation
#v_lst_indels_pop_mut_autocorrelation <- vapply(X = sort(unique(df_indel_variants$pop_mut)),FUN = function(the_pop_mut) max(acf(subset(df_indel_variants,pop_mut==the_pop_mut)$VarFreq[order(subset(df_indel_variants,pop_mut==the_pop_mut)$int_time)],plot = F)$acf[acf(subset(df_indel_variants,pop_mut==the_pop_mut)$VarFreq[order(subset(df_indel_variants,pop_mut==the_pop_mut)$int_time)],plot = F)$acf<1],na.rm=T),FUN.VALUE = 0.0)
##Only keep indels that have enough lag1 autocorrelation (acf >=0.1)
#df_indel_variants <- subset(df_indel_variants,pop_mut%in%names(v_lst_indels_pop_mut_autocorrelation[v_lst_indels_pop_mut_autocorrelation>=0.1]))
#print(length(unique(subset(df_indel_variants,is_fixed)$mutation)))
#add gene and loci type annotation to variants dataframe
df_indel_variants$variant_type <- ifelse(test = substr(x = df_indel_variants$VarAllele,1,1)=="+",yes = "Insertion",no="Deletion")
df_indel_variants$gene_product <- vapply(X = 1:nrow(df_indel_variants),FUN = find_gene_product_affected_by_indel,FUN.VALUE = "")
df_indel_variants$gene_product <- ifelse(df_indel_variants$gene_product=="NA",yes=NA,no=df_indel_variants$gene_product)
df_indel_variants$accession_number <- vapply(X = 1:nrow(df_indel_variants),FUN = find_acession_number_gene_product_affected_by_indel,FUN.VALUE = "")
df_indel_variants$accession_number <- ifelse(df_indel_variants$accession_number=="NA",yes=NA,no=df_indel_variants$accession_number)
df_indel_variants$variant_type2 <- vapply(X = 1:nrow(df_indel_variants),FUN = find_indel_variant_type2,FUN.VALUE = "")
df_indel_variants$variant_type2 <- ifelse(df_indel_variants$variant_type2=="NA",yes=NA,no=df_indel_variants$variant_type2)
df_indel_variants$comparison2_category <- vapply(X = 1:nrow(df_indel_variants),FUN = find_indel_comparison2_category,FUN.VALUE = "")
df_indel_variants$comparison2_category <- ifelse(df_indel_variants$comparison2_category=="NA",yes=NA,no=df_indel_variants$comparison2_category)
#remove indels in low-complexity regions
df_indel_variants <- subset(df_indel_variants, mutation != "II_G3983331+TA")

#CALCULATE variants prevalence
df_prevalence <- df_ALL_variants %>%
  group_by(mutation,Chrom,Position) %>%
  dplyr::summarise(prevalence=length(unique(population)))

#Create a combined dataframe with SNVs and indels based on columns c("population","Chrom","Position","Ref","VarFreq","VarAllele","is_fixed","variant_type")
df_ALL_variants <- rbind(df_variants[,c("Sample","population","str_time","Chrom","Position","Ref","VarFreq","VarAllele","mutation","is_fixed","variant_type","gene_product","accession_number")],df_indel_variants[,c("Sample","population","str_time","Chrom","Position","Ref","VarFreq","VarAllele","mutation","is_fixed","variant_type","gene_product","accession_number")])
#Evaluate coevolution to eliminate mutations that are uncorrelated to others positively OR negatively
df_ALL_variants$pop_mut <- paste0(df_ALL_variants$population,"_",df_ALL_variants$mutation)
mtx_allele_freq_time_series <- as.matrix(xtabs(VarFreq~pop_mut+str_time, data=df_ALL_variants))
mtx_autocorrelation_pop_muts <- cor(t(mtx_allele_freq_time_series))
v_lst_variants <- rownames(mtx_autocorrelation_pop_muts)
v_pops_muts <- vapply(X = v_lst_variants,FUN = function(x) substr(x,start = 1,stop=gregexpr(pattern = "_",text = x,fixed = T)[[1]][1]-1),FUN.VALUE = "")
v_nb_positively_coevolving_pop_muts <- rep(F,length(v_lst_variants))
names(v_nb_positively_coevolving_pop_muts) <- v_lst_variants
v_nb_negatively_coevolving_pop_muts <- rep(F,length(v_lst_variants))
names(v_nb_negatively_coevolving_pop_muts) <- v_lst_variants
for (current_variant in v_lst_variants){
  pop_current_variant <- substr(current_variant,start = 1,stop=gregexpr(pattern = "_",text = current_variant,fixed = T)[[1]][1]-1)
  v_nb_positively_coevolving_pop_muts[current_variant] <- sum((v_pops_muts==pop_current_variant)&(mtx_autocorrelation_pop_muts[current_variant,] > 0.9),na.rm=T)
  v_nb_negatively_coevolving_pop_muts[current_variant] <- sum((v_pops_muts==pop_current_variant)&(mtx_autocorrelation_pop_muts[current_variant,] < (-0.9) ),na.rm=T)
}
v_lst_pop_muts_coevolving <- ((v_nb_positively_coevolving_pop_muts>1)|(v_nb_negatively_coevolving_pop_muts>1))&(v_nb_positively_coevolving_pop_muts+v_nb_negatively_coevolving_pop_muts>3)
print(summary(v_nb_positively_coevolving_pop_muts))
print(summary(v_nb_negatively_coevolving_pop_muts))
nb_row_before_coevolution_filter <- nrow(df_ALL_variants)
print(paste0("Number of substitutions before co-evolution filter= ",length(unique(subset(df_variants,is_fixed)$mutation))))
print(paste0("Number of fixed indels before co-evolution filter= ",length(unique(subset(df_indel_variants,is_fixed)$mutation))))
df_ALL_variants <- subset(df_ALL_variants,pop_mut %in%names(v_lst_pop_muts_coevolving)[v_lst_pop_muts_coevolving])
nb_row_after_coevolution_filter <- nrow(df_ALL_variants)
df_variants <- subset(df_variants,pop_mut %in%names(v_lst_pop_muts_coevolving)[v_lst_pop_muts_coevolving])
df_indel_variants <- subset(df_indel_variants,pop_mut %in%names(v_lst_pop_muts_coevolving)[v_lst_pop_muts_coevolving])
print(paste0("Number of substitutions after co-evolution filter= ",length(unique(subset(df_variants,is_fixed)$mutation))))
print(paste0("Number of fixed indels after co-evolution filter= ",length(unique(subset(df_indel_variants,is_fixed)$mutation))))
#Remove mutations that are not unique to each population as parrallelism is very rare
df_ALL_variants <- subset(df_ALL_variants,mutation%in%subset(df_prevalence,prevalence==1)$mutation)

#Variant spectrum across populations (insertion, deletion, SNV) with separate grid for non-fixed and fixed variants
df_ALL_variants$concat_group_by_pop_is_fixed_and_vartype <- paste0(df_ALL_variants$population,"_",df_ALL_variants$is_fixed,"_",df_ALL_variants$variant_type)
df_nb_UNIQUE_variants_per_pop <- df_ALL_variants %>%
  group_by(concat_group_by_pop_is_fixed_and_vartype) %>%
  dplyr::summarise(nb_unique_variants = length(unique(mutation))) %>%
  as.data.frame()
v_nb_non_fixed_variants_per_pop <- table(subset(unique(df_ALL_variants[,c("population","mutation","is_fixed","variant_type")]),!is_fixed)$population)
v_nb_fixed_variants_per_pop <- table(subset(unique(df_ALL_variants[,c("population","mutation","is_fixed","variant_type")]),is_fixed)$population)
df_nb_UNIQUE_variants_per_pop$population <- sapply(X = strsplit(df_nb_UNIQUE_variants_per_pop$concat_group_by_pop_is_fixed_and_vartype,split = "_"),FUN = function(x) x[1])
df_nb_UNIQUE_variants_per_pop$is_fixed <- as.logical(sapply(X = strsplit(df_nb_UNIQUE_variants_per_pop$concat_group_by_pop_is_fixed_and_vartype,split = "_"),FUN = function(x) x[2]))
df_nb_UNIQUE_variants_per_pop$var_type <- sapply(X = strsplit(df_nb_UNIQUE_variants_per_pop$concat_group_by_pop_is_fixed_and_vartype,split = "_"),FUN = function(x) x[3])
df_nb_UNIQUE_variants_per_pop <- df_nb_UNIQUE_variants_per_pop[,c("population","is_fixed","var_type","nb_unique_variants")]
df_nb_UNIQUE_variants_per_pop$proportion <- df_nb_UNIQUE_variants_per_pop$nb_unique_variants/ifelse(test = df_nb_UNIQUE_variants_per_pop$is_fixed,yes = v_nb_fixed_variants_per_pop[df_nb_UNIQUE_variants_per_pop$population],no = v_nb_non_fixed_variants_per_pop[df_nb_UNIQUE_variants_per_pop$population])
v_fdr_variant_fixation_bias_across_pops <- p.adjust(vapply(X = v_pop_order, FUN = function(the_pop) chisq.test(apply(X = as.matrix(reshape2::acast(subset(df_nb_UNIQUE_variants_per_pop,population==the_pop), var_type ~ is_fixed,value.var = 'nb_unique_variants', fill = '0')),MARGIN = 2,FUN = as.numeric))$p.value,FUN.VALUE = 0.0),"fdr")
v_label_top_of_cols <- ifelse(test = v_fdr_variant_fixation_bias_across_pops>=0.05,yes = "N.S.",
                              no=ifelse(test = (v_fdr_variant_fixation_bias_across_pops>=0.01)&(v_fdr_variant_fixation_bias_across_pops<0.05),yes = "*",no =  ifelse(test = (v_fdr_variant_fixation_bias_across_pops>=0.001)&(v_fdr_variant_fixation_bias_across_pops<0.01),yes = "**",no =  "***")))
df_nb_UNIQUE_variants_per_pop$fdr_label <- ifelse(test = df_nb_UNIQUE_variants_per_pop$is_fixed,yes="",no=v_label_top_of_cols[df_nb_UNIQUE_variants_per_pop$population])
v_first_occ_to_add <- v_pop_order
for (i in 1:nrow(df_nb_UNIQUE_variants_per_pop)){
  if (df_nb_UNIQUE_variants_per_pop$population[i]%in%v_first_occ_to_add){
    df_nb_UNIQUE_variants_per_pop$fdr_label[i] <- df_nb_UNIQUE_variants_per_pop$fdr_label[i]
    v_first_occ_to_add <- v_first_occ_to_add[v_first_occ_to_add!=df_nb_UNIQUE_variants_per_pop$population[i]]
  }else{
    df_nb_UNIQUE_variants_per_pop$fdr_label[i] <- ""
  }
}

df_nb_UNIQUE_variants_per_pop$lbl_fixed <- factor(ifelse(df_nb_UNIQUE_variants_per_pop$is_fixed,yes = "Fixed (VAF >= 0.75)",no = "Not fixed (VAF < 0.75)"),levels = c("Not fixed (VAF < 0.75)","Fixed (VAF >= 0.75)"))
ggplot(data = df_nb_UNIQUE_variants_per_pop,mapping = aes(x=factor(population,v_pop_order),y=proportion)) + geom_col(mapping = aes(fill=factor(var_type,levels=c("Deletion","Insertion","SNV")))) + geom_text(aes(label=fdr_label), position=position_dodge(width=0.5), vjust=-15,size=6) + xlab("Population") + ylab("Proportion") +
  theme_bw() + theme(strip.text = element_text(size=16),axis.text.x=element_text(size=14),axis.title.x = element_text(size=16),axis.text.y=element_text(size=14),axis.title.y = element_text(size=16),legend.title = element_text(size=16),legend.text = element_text(size=14)) + guides(fill=guide_legend(title="Variant")) + scale_y_continuous(limits = c(0,1.1),breaks = seq(0,1,0.2)) + scale_fill_brewer(palette = "Dark2") + facet_wrap(~lbl_fixed,ncol=1,scales="free_y")
ggsave(filename = "Fixation_bias_by_variant_type_across_populations.eps", path=output_workspace, width = 25/2.54, height = 25/2.54, device = cairo_ps)
ggsave(filename = "Fixation_bias_by_variant_type_across_populations.svg", path=output_workspace, width = 25/2.54, height = 25/2.54, device = svg)
#SNV spectrum across populations (6 categories) with separate grid for non-fixed and fixed variants
df_variants$concat_group_by_pop_is_fixed_and_SNV_category <- paste0(df_variants$population,"_",df_variants$is_fixed,"_",df_variants$SNV_category)
df_nb_UNIQUE_SNVs_per_pop <- df_variants %>%
  group_by(concat_group_by_pop_is_fixed_and_SNV_category) %>%
  dplyr::summarise(nb_unique_SNVs = length(unique(mutation))) %>%
  as.data.frame()
v_nb_non_fixed_SNVs_per_pop <- table(subset(unique(df_variants[,c("population","mutation","is_fixed","SNV_category")]),!is_fixed)$population)
v_nb_fixed_SNVs_per_pop <- table(subset(unique(df_variants[,c("population","mutation","is_fixed","SNV_category")]),is_fixed)$population)
df_nb_UNIQUE_SNVs_per_pop$population <- sapply(X = strsplit(df_nb_UNIQUE_SNVs_per_pop$concat_group_by_pop_is_fixed_and_SNV_category,split = "_"),FUN = function(x) x[1])
df_nb_UNIQUE_SNVs_per_pop$is_fixed <- as.logical(sapply(X = strsplit(df_nb_UNIQUE_SNVs_per_pop$concat_group_by_pop_is_fixed_and_SNV_category,split = "_"),FUN = function(x) x[2]))
df_nb_UNIQUE_SNVs_per_pop$SNV_category <- sapply(X = strsplit(df_nb_UNIQUE_SNVs_per_pop$concat_group_by_pop_is_fixed_and_SNV_category,split = "_"),FUN = function(x) x[3])
df_nb_UNIQUE_SNVs_per_pop <- df_nb_UNIQUE_SNVs_per_pop[,c("population","is_fixed","SNV_category","nb_unique_SNVs")]
df_nb_UNIQUE_SNVs_per_pop$proportion <- df_nb_UNIQUE_SNVs_per_pop$nb_unique_SNVs/ifelse(test = df_nb_UNIQUE_SNVs_per_pop$is_fixed,yes = v_nb_fixed_SNVs_per_pop[df_nb_UNIQUE_SNVs_per_pop$population],no = v_nb_non_fixed_SNVs_per_pop[df_nb_UNIQUE_SNVs_per_pop$population])
v_fdr_variant_fixation_bias_across_pops <- p.adjust(vapply(X = v_pop_order, FUN = function(the_pop) chisq.test(apply(X = as.matrix(reshape2::acast(subset(df_nb_UNIQUE_SNVs_per_pop,population==the_pop), SNV_category ~ is_fixed,value.var = 'nb_unique_SNVs', fill = '0')),MARGIN = 2,FUN = as.numeric))$p.value,FUN.VALUE = 0.0),"fdr")
v_label_top_of_cols <- ifelse(test = v_fdr_variant_fixation_bias_across_pops>=0.05,yes = "N.S.",
                              no=ifelse(test = (v_fdr_variant_fixation_bias_across_pops>=0.01)&(v_fdr_variant_fixation_bias_across_pops<0.05),yes = "*",no =  ifelse(test = (v_fdr_variant_fixation_bias_across_pops>=0.001)&(v_fdr_variant_fixation_bias_across_pops<0.01),yes = "**",no =  "***")))
df_nb_UNIQUE_SNVs_per_pop$fdr_label <- ifelse(test = df_nb_UNIQUE_SNVs_per_pop$is_fixed,yes="",no=v_label_top_of_cols[df_nb_UNIQUE_SNVs_per_pop$population])
v_first_occ_to_add <- v_pop_order
for (i in 1:nrow(df_nb_UNIQUE_SNVs_per_pop)){
  if (df_nb_UNIQUE_SNVs_per_pop$population[i]%in%v_first_occ_to_add){
    df_nb_UNIQUE_SNVs_per_pop$fdr_label[i] <- df_nb_UNIQUE_SNVs_per_pop$fdr_label[i]
    v_first_occ_to_add <- v_first_occ_to_add[v_first_occ_to_add!=df_nb_UNIQUE_SNVs_per_pop$population[i]]
  }else{
    df_nb_UNIQUE_SNVs_per_pop$fdr_label[i] <- ""
  }
}
df_nb_UNIQUE_SNVs_per_pop$lbl_fixed <- factor(ifelse(df_nb_UNIQUE_SNVs_per_pop$is_fixed,yes = "Fixed (VAF >= 0.75)",no = "Not fixed (VAF < 0.75)"),levels = c("Not fixed (VAF < 0.75)","Fixed (VAF >= 0.75)"))
ggplot(data = df_nb_UNIQUE_SNVs_per_pop,mapping = aes(x=factor(population,v_pop_order),y=proportion)) + geom_col(mapping = aes(fill=SNV_category)) + geom_text(aes(label=fdr_label), position=position_dodge(width=0.75), vjust=-15,size=6) + xlab("Population") + ylab("Proportion") +
  theme_bw() + theme(strip.text = element_text(size=16),axis.text.x=element_text(size=14),axis.title.x = element_text(size=16),axis.text.y=element_text(size=14),axis.title.y = element_text(size=16),legend.title = element_text(size=16),legend.text = element_text(size=14)) + guides(fill=guide_legend(title="Variant")) + scale_y_continuous(limits = c(0,1.1),breaks = seq(0,1,0.2)) + scale_fill_brewer(palette = "Dark2") + facet_wrap(~lbl_fixed,ncol=1,scales="free_y")
ggsave(filename = "Fixation_bias_by_SNV_category_across_populations.eps", path=output_workspace, width = 25/2.54, height = 20/2.54, device = cairo_ps)
ggsave(filename = "Fixation_bias_by_SNV_category_across_populations.svg", path=output_workspace, width = 25/2.54, height = 20/2.54, device = svg)

#Variant spectrum across populations and chromosomes (insertion, deletion, SNV)
df_ALL_variants$concat_group_by_pop_chr_is_fixed_and_vartype <- paste0(df_ALL_variants$population,";",df_ALL_variants$Chrom,";",df_ALL_variants$is_fixed,";",df_ALL_variants$variant_type)
df_nb_UNIQUE_variants_per_pop_and_chr <- df_ALL_variants %>%
  group_by(concat_group_by_pop_chr_is_fixed_and_vartype) %>%
  dplyr::summarise(nb_unique_variants = length(unique(mutation))) %>%
  as.data.frame()
df_nb_UNIQUE_variants_per_pop_and_chr$population <- sapply(X = strsplit(df_nb_UNIQUE_variants_per_pop_and_chr$concat_group_by_pop_chr_is_fixed_and_vartype,split = ";"),FUN = function(x) x[1])
df_nb_UNIQUE_variants_per_pop_and_chr$chromosome <- sapply(X = strsplit(df_nb_UNIQUE_variants_per_pop_and_chr$concat_group_by_pop_chr_is_fixed_and_vartype,split = ";"),FUN = function(x) x[2])
df_nb_UNIQUE_variants_per_pop_and_chr$is_fixed <- as.logical(sapply(X = strsplit(df_nb_UNIQUE_variants_per_pop_and_chr$concat_group_by_pop_chr_is_fixed_and_vartype,split = ";"),FUN = function(x) x[3]))
df_nb_UNIQUE_variants_per_pop_and_chr$var_type <- sapply(X = strsplit(df_nb_UNIQUE_variants_per_pop_and_chr$concat_group_by_pop_chr_is_fixed_and_vartype,split = ";"),FUN = function(x) x[4])
df_nb_UNIQUE_variants_per_pop_and_chr <- df_nb_UNIQUE_variants_per_pop_and_chr[,c("population","chromosome","is_fixed","var_type","nb_unique_variants")]
df_nb_UNIQUE_variants_per_pop_and_chr$density <- df_nb_UNIQUE_variants_per_pop_and_chr$nb_unique_variants/v_chrs_length[df_nb_UNIQUE_variants_per_pop_and_chr$chromosome]
#v_fdr_variant_fixation_bias_across_pops <- p.adjust(vapply(X = v_pop_order, FUN = function(the_pop) chisq.test(apply(X = as.matrix(reshape2::acast(subset(df_nb_UNIQUE_variants_per_pop_and_chr,population==the_pop), var_type ~ is_fixed,value.var = 'nb_unique_variants', fill = '0')),MARGIN = 2,FUN = as.numeric))$p.value,FUN.VALUE = 0.0),"fdr")
#v_label_top_of_cols <- ifelse(test = v_fdr_variant_fixation_bias_across_pops>=0.05,yes = "N.S.",
#                              no=ifelse(test = (v_fdr_variant_fixation_bias_across_pops>=0.01)&(v_fdr_variant_fixation_bias_across_pops<0.05),yes = "*",no =  ifelse(test = (v_fdr_variant_fixation_bias_across_pops>=0.001)&(v_fdr_variant_fixation_bias_across_pops<0.01),yes = "**",no =  "***")))
#df_nb_UNIQUE_variants_per_pop_and_chr$fdr_label <- ifelse(test = df_nb_UNIQUE_variants_per_pop_and_chr$is_fixed,yes="",no=v_label_top_of_cols[df_nb_UNIQUE_variants_per_pop_and_chr$population])
#v_pop_occs_to_add <- rep(v_pop_order,each=3)
#v_chr_occs_to_add <- rep(c("I","II","III"),each=3)
#v_occs_to_add <- paste0(v_pop_occs_to_add,";",v_chr_occs_to_add)
#for (i in 1:nrow(df_nb_UNIQUE_variants_per_pop_and_chr)){
#  if (paste0(df_nb_UNIQUE_variants_per_pop_and_chr$population[i],";",df_nb_UNIQUE_variants_per_pop_and_chr$chromosome[i])%in%v_occs_to_add){
#    df_nb_UNIQUE_variants_per_pop_and_chr$fdr_label[i] <- df_nb_UNIQUE_variants_per_pop_and_chr$fdr_label[i]
#    v_occs_to_add <- v_occs_to_add[v_occs_to_add!=paste0(df_nb_UNIQUE_variants_per_pop_and_chr$population[i],";",df_nb_UNIQUE_variants_per_pop_and_chr$chromosome[i])]
#  }else{
#    df_nb_UNIQUE_variants_per_pop_and_chr$fdr_label[i] <- ""
#  }
#}

df_nb_UNIQUE_variants_per_pop_and_chr$lbl_fixed <- factor(ifelse(df_nb_UNIQUE_variants_per_pop_and_chr$is_fixed,yes = "Fixed (VAF >= 0.75)",no = "Not fixed (VAF < 0.75)"),levels = c("Not fixed (VAF < 0.75)","Fixed (VAF >= 0.75)"))
df_nb_UNIQUE_variants_per_pop_and_chr <- subset(df_nb_UNIQUE_variants_per_pop_and_chr,chromosome!="chr_II_telomeric_gap")
ggplot(data = df_nb_UNIQUE_variants_per_pop_and_chr,mapping = aes(x=factor(population,v_pop_order),y=density)) + geom_col(mapping = aes(fill=factor(var_type,levels=c("Deletion","Insertion","SNV")))) + xlab("Population") + ylab("Density (Number of mutations / bp)") +
  theme_bw() + theme(strip.text = element_text(size=16),axis.text.x=element_text(size=14),axis.title.x = element_text(size=16),axis.text.y=element_text(size=14),axis.title.y = element_text(size=16),legend.title = element_text(size=16),legend.text = element_text(size=14)) + guides(fill=guide_legend(title="Variant")) + scale_fill_brewer(palette = "Dark2") + facet_wrap(chromosome~lbl_fixed,scales="free",ncol=2) #+ geom_text(aes(label=fdr_label), position=position_dodge(width=0.75), vjust=-3,size=6)
ggsave(filename = "SNV_Fixation_bias_by_variant_type_across_populations_and_chromosomes.eps", path=output_workspace, width = 30/2.54, height = 30/2.54, device = cairo_ps)
ggsave(filename = "SNV_Fixation_bias_by_variant_type_across_populations_and_chromosomes.svg", path=output_workspace, width = 30/2.54, height = 30/2.54, device = svg)

#mutation density across populations and LOCI of interest with separate grid for non-fixed and fixed variants
df_variants$concat_group_by_pop_is_fixed_and_variant_type2 <- paste0(df_variants$population,"_",df_variants$is_fixed,"_",df_variants$variant_type2)
df_nb_UNIQUE_variants_per_pop <- df_variants %>%
  group_by(concat_group_by_pop_is_fixed_and_variant_type2) %>%
  dplyr::summarise(nb_unique_variants = length(unique(mutation))) %>%
  as.data.frame()
df_nb_UNIQUE_variants_per_pop$population <- sapply(X = strsplit(df_nb_UNIQUE_variants_per_pop$concat_group_by_pop_is_fixed_and_variant_type2,split = "_"),FUN = function(x) x[1])
df_nb_UNIQUE_variants_per_pop$is_fixed <- as.logical(sapply(X = strsplit(df_nb_UNIQUE_variants_per_pop$concat_group_by_pop_is_fixed_and_variant_type2,split = "_"),FUN = function(x) x[2]))
df_nb_UNIQUE_variants_per_pop$variant_type2 <- sapply(X = strsplit(df_nb_UNIQUE_variants_per_pop$concat_group_by_pop_is_fixed_and_variant_type2,split = "_"),FUN = function(x) x[3])
df_nb_UNIQUE_variants_per_pop <- df_nb_UNIQUE_variants_per_pop[,c("population","is_fixed","variant_type2","nb_unique_variants")]
df_nb_UNIQUE_variants_per_pop$density <- df_nb_UNIQUE_variants_per_pop$nb_unique_variants/v_length_annotated_genomic_regions[df_nb_UNIQUE_variants_per_pop$variant_type2]
df_nb_UNIQUE_variants_per_pop <- subset(df_nb_UNIQUE_variants_per_pop,variant_type2!="NA")
v_fdr_variant_fixation_bias_across_pops <- p.adjust(vapply(X = v_pop_order, FUN = function(the_pop) chisq.test(apply(X = as.matrix(reshape2::acast(subset(df_nb_UNIQUE_variants_per_pop,population==the_pop), variant_type2 ~ is_fixed,value.var = 'nb_unique_variants', fill = '0')),MARGIN = 2,FUN = as.numeric))$p.value,FUN.VALUE = 0.0),"fdr")
v_label_top_of_cols <- ifelse(test = v_fdr_variant_fixation_bias_across_pops>=0.05,yes = "N.S.",
                              no=ifelse(test = (v_fdr_variant_fixation_bias_across_pops>=0.01)&(v_fdr_variant_fixation_bias_across_pops<0.05),yes = "*",no =  ifelse(test = (v_fdr_variant_fixation_bias_across_pops>=0.001)&(v_fdr_variant_fixation_bias_across_pops<0.01),yes = "**",no =  "***")))
df_nb_UNIQUE_variants_per_pop$fdr_label <- ifelse(test = df_nb_UNIQUE_variants_per_pop$is_fixed,yes="",no=v_label_top_of_cols[df_nb_UNIQUE_variants_per_pop$population])
v_occs_to_add <- NULL
for (i in 1:nrow(df_nb_UNIQUE_variants_per_pop)){
  if (!paste0(df_nb_UNIQUE_variants_per_pop$population[i],";",df_nb_UNIQUE_variants_per_pop$variant_type2[i])%in%v_occs_to_add){
    df_nb_UNIQUE_variants_per_pop$fdr_label[i] <- df_nb_UNIQUE_variants_per_pop$fdr_label[i]
    v_occs_to_add <- c(v_occs_to_add,paste0(df_nb_UNIQUE_variants_per_pop$population[i],";",df_nb_UNIQUE_variants_per_pop$variant_type2[i]))
  }else{
    df_nb_UNIQUE_variants_per_pop$fdr_label[i] <- ""
  }
}

df_nb_UNIQUE_variants_per_pop$lbl_fixed <- factor(ifelse(df_nb_UNIQUE_variants_per_pop$is_fixed,yes = "Fixed (VAF >= 0.75)",no = "Not fixed (VAF < 0.75)"),levels = c("Not fixed (VAF < 0.75)","Fixed (VAF >= 0.75)"))
df_nb_UNIQUE_variants_per_pop$log10_density <- log10(df_nb_UNIQUE_variants_per_pop$density)
df_nb_UNIQUE_variants_per_pop <- subset(df_nb_UNIQUE_variants_per_pop,!is.na(variant_type2))
ggplot(data = df_nb_UNIQUE_variants_per_pop,mapping = aes(x=variant_type2,y=-log10_density)) + geom_col(mapping = aes(fill = lbl_fixed),position = "dodge") + geom_text(aes(label=fdr_label), position=position_dodge(width=0.75), vjust=-0.5,size=4) + xlab("Locus") + ylab("Density\n-log10 (number of SNVs / bp)") +
  theme_bw() + theme(strip.text = element_text(size=16),axis.text.x=element_text(size=14,angle=60,hjust=1),axis.title.x = element_text(size=16),axis.text.y=element_text(size=14),axis.title.y = element_text(size=16),legend.title = element_text(size=16),legend.text = element_text(size=14)) + guides(fill=guide_legend(title="Variant")) + scale_fill_brewer(palette = "Set1") + scale_y_continuous(limits = c(0,9),breaks=seq(0,8,2)) + facet_wrap(~factor(population,v_pop_order),ncol=2)
ggsave(filename = "SNV_Fixation_bias_by_locus_across_populations_and_loci_of_interest.eps", path=output_workspace, width = 25/2.54, height = 20/2.54, device = cairo_ps)
ggsave(filename = "SNV_Fixation_bias_by_locus_across_populations_and_loci_of_interest.svg", path=output_workspace, width = 25/2.54, height = 20/2.54, device = svg)

#mutation density across populations and other loci with separate grid for non-fixed and fixed variants
df_variants$concat_group_by_pop_is_fixed_and_comparison2_category <- paste0(df_variants$population,"_",df_variants$is_fixed,"_",df_variants$comparison2_category)
df_nb_UNIQUE_variants_per_pop <- df_variants %>%
  group_by(concat_group_by_pop_is_fixed_and_comparison2_category) %>%
  dplyr::summarise(nb_unique_variants = length(unique(mutation))) %>%
  as.data.frame()
df_nb_UNIQUE_variants_per_pop$population <- sapply(X = strsplit(df_nb_UNIQUE_variants_per_pop$concat_group_by_pop_is_fixed_and_comparison2_category,split = "_"),FUN = function(x) x[1])
df_nb_UNIQUE_variants_per_pop$is_fixed <- as.logical(sapply(X = strsplit(df_nb_UNIQUE_variants_per_pop$concat_group_by_pop_is_fixed_and_comparison2_category,split = "_"),FUN = function(x) x[2]))
df_nb_UNIQUE_variants_per_pop$comparison2_category <- sapply(X = strsplit(df_nb_UNIQUE_variants_per_pop$concat_group_by_pop_is_fixed_and_comparison2_category,split = "_"),FUN = function(x) x[3])
df_nb_UNIQUE_variants_per_pop <- df_nb_UNIQUE_variants_per_pop[,c("population","is_fixed","comparison2_category","nb_unique_variants")]
df_nb_UNIQUE_variants_per_pop$density <- df_nb_UNIQUE_variants_per_pop$nb_unique_variants/v_length_annotated_genomic_regions[df_nb_UNIQUE_variants_per_pop$comparison2_category]
df_nb_UNIQUE_variants_per_pop <- subset(df_nb_UNIQUE_variants_per_pop,comparison2_category!="NA")
v_fdr_variant_fixation_bias_across_pops <- p.adjust(vapply(X = v_pop_order, FUN = function(the_pop) chisq.test(apply(X = as.matrix(reshape2::acast(subset(df_nb_UNIQUE_variants_per_pop,population==the_pop), comparison2_category ~ is_fixed,value.var = 'nb_unique_variants', fill = '0')),MARGIN = 2,FUN = as.numeric))$p.value,FUN.VALUE = 0.0),"fdr")
v_label_top_of_cols <- ifelse(test = v_fdr_variant_fixation_bias_across_pops>=0.05,yes = "N.S.",
                              no=ifelse(test = (v_fdr_variant_fixation_bias_across_pops>=0.01)&(v_fdr_variant_fixation_bias_across_pops<0.05),yes = "*",no =  ifelse(test = (v_fdr_variant_fixation_bias_across_pops>=0.001)&(v_fdr_variant_fixation_bias_across_pops<0.01),yes = "**",no =  "***")))
df_nb_UNIQUE_variants_per_pop$fdr_label <- ifelse(test = df_nb_UNIQUE_variants_per_pop$is_fixed,yes="",no=v_label_top_of_cols[df_nb_UNIQUE_variants_per_pop$population])
v_occs_to_add <- NULL
for (i in 1:nrow(df_nb_UNIQUE_variants_per_pop)){
  if (!paste0(df_nb_UNIQUE_variants_per_pop$population[i],";",df_nb_UNIQUE_variants_per_pop$comparison2_category[i])%in%v_occs_to_add){
    df_nb_UNIQUE_variants_per_pop$fdr_label[i] <- df_nb_UNIQUE_variants_per_pop$fdr_label[i]
    v_occs_to_add <- c(v_occs_to_add,paste0(df_nb_UNIQUE_variants_per_pop$population[i],";",df_nb_UNIQUE_variants_per_pop$comparison2_category[i]))
  }else{
    df_nb_UNIQUE_variants_per_pop$fdr_label[i] <- ""
  }
}

df_nb_UNIQUE_variants_per_pop$lbl_fixed <- factor(ifelse(df_nb_UNIQUE_variants_per_pop$is_fixed,yes = "Fixed (VAF >= 0.75)",no = "Not fixed (VAF < 0.75)"),levels = c("Not fixed (VAF < 0.75)","Fixed (VAF >= 0.75)"))
df_nb_UNIQUE_variants_per_pop$log10_density <- log10(df_nb_UNIQUE_variants_per_pop$density)
df_nb_UNIQUE_variants_per_pop <- subset(df_nb_UNIQUE_variants_per_pop,!is.na(comparison2_category))
ggplot(data = df_nb_UNIQUE_variants_per_pop,mapping = aes(x=comparison2_category,y=-log10_density)) + geom_col(mapping = aes(fill = lbl_fixed),position = "dodge") + geom_text(aes(label=fdr_label), position=position_dodge(width=0.75), vjust=-0.5,size=4) + xlab("Locus") + ylab("Density\n-log10 (number of SNVs / bp)") +
  theme_bw() + theme(strip.text = element_text(size=16),axis.text.x=element_text(size=14,angle=60,hjust=1),axis.title.x = element_text(size=16),axis.text.y=element_text(size=14),axis.title.y = element_text(size=16),legend.title = element_text(size=16),legend.text = element_text(size=14)) + guides(fill=guide_legend(title="Variant")) + scale_fill_brewer(palette = "Set1") + scale_y_continuous(limits = c(0,9),breaks=seq(0,8,2)) + facet_wrap(~factor(population,v_pop_order),ncol=2)
ggsave(filename = "SNV_Fixation_bias_by_locus_across_populations_and_other_loci.eps", path=output_workspace, width = 25/2.54, height = 20/2.54, device = cairo_ps)
ggsave(filename = "SNV_Fixation_bias_by_locus_across_populations_and_other_loci.svg", path=output_workspace, width = 25/2.54, height = 20/2.54, device = svg)

#Indels density across populations and LOCI of interest with separate grid for non-fixed and fixed variants
df_indel_variants$concat_group_by_pop_is_fixed_and_variant_type2 <- paste0(df_indel_variants$population,"_",df_indel_variants$is_fixed,"_",df_indel_variants$variant_type2)
df_nb_UNIQUE_variants_per_pop <- df_indel_variants %>%
  group_by(concat_group_by_pop_is_fixed_and_variant_type2) %>%
  dplyr::summarise(nb_unique_variants = length(unique(mutation))) %>%
  as.data.frame()
df_nb_UNIQUE_variants_per_pop$population <- sapply(X = strsplit(df_nb_UNIQUE_variants_per_pop$concat_group_by_pop_is_fixed_and_variant_type2,split = "_"),FUN = function(x) x[1])
df_nb_UNIQUE_variants_per_pop$is_fixed <- as.logical(sapply(X = strsplit(df_nb_UNIQUE_variants_per_pop$concat_group_by_pop_is_fixed_and_variant_type2,split = "_"),FUN = function(x) x[2]))
df_nb_UNIQUE_variants_per_pop$variant_type2 <- sapply(X = strsplit(df_nb_UNIQUE_variants_per_pop$concat_group_by_pop_is_fixed_and_variant_type2,split = "_"),FUN = function(x) x[3])
df_nb_UNIQUE_variants_per_pop <- df_nb_UNIQUE_variants_per_pop[,c("population","is_fixed","variant_type2","nb_unique_variants")]
df_nb_UNIQUE_variants_per_pop$density <- df_nb_UNIQUE_variants_per_pop$nb_unique_variants/v_length_annotated_genomic_regions[df_nb_UNIQUE_variants_per_pop$variant_type2]
df_nb_UNIQUE_variants_per_pop <- subset(df_nb_UNIQUE_variants_per_pop,variant_type2!="NA")
v_fdr_variant_fixation_bias_across_pops <- p.adjust(vapply(X = v_pop_order, FUN = function(the_pop) chisq.test(apply(X = as.matrix(reshape2::acast(subset(df_nb_UNIQUE_variants_per_pop,population==the_pop), variant_type2 ~ is_fixed,value.var = 'nb_unique_variants', fill = '0')),MARGIN = 2,FUN = as.numeric))$p.value,FUN.VALUE = 0.0),"fdr")
v_label_top_of_cols <- ifelse(test = v_fdr_variant_fixation_bias_across_pops>=0.05,yes = "N.S.",
                              no=ifelse(test = (v_fdr_variant_fixation_bias_across_pops>=0.01)&(v_fdr_variant_fixation_bias_across_pops<0.05),yes = "*",no =  ifelse(test = (v_fdr_variant_fixation_bias_across_pops>=0.001)&(v_fdr_variant_fixation_bias_across_pops<0.01),yes = "**",no =  "***")))
df_nb_UNIQUE_variants_per_pop$fdr_label <- ifelse(test = df_nb_UNIQUE_variants_per_pop$is_fixed,yes="",no=v_label_top_of_cols[df_nb_UNIQUE_variants_per_pop$population])
v_occs_to_add <- NULL
for (i in 1:nrow(df_nb_UNIQUE_variants_per_pop)){
  if (!paste0(df_nb_UNIQUE_variants_per_pop$population[i],";",df_nb_UNIQUE_variants_per_pop$variant_type2[i])%in%v_occs_to_add){
    df_nb_UNIQUE_variants_per_pop$fdr_label[i] <- df_nb_UNIQUE_variants_per_pop$fdr_label[i]
    v_occs_to_add <- c(v_occs_to_add,paste0(df_nb_UNIQUE_variants_per_pop$population[i],";",df_nb_UNIQUE_variants_per_pop$variant_type2[i]))
  }else{
    df_nb_UNIQUE_variants_per_pop$fdr_label[i] <- ""
  }
}

df_nb_UNIQUE_variants_per_pop$lbl_fixed <- factor(ifelse(df_nb_UNIQUE_variants_per_pop$is_fixed,yes = "Fixed (VAF >= 0.75)",no = "Not fixed (VAF < 0.75)"),levels = c("Not fixed (VAF < 0.75)","Fixed (VAF >= 0.75)"))
df_nb_UNIQUE_variants_per_pop$log10_density <- log10(df_nb_UNIQUE_variants_per_pop$density)
ggplot(data = df_nb_UNIQUE_variants_per_pop,mapping = aes(x=variant_type2,y=-log10_density)) + geom_col(mapping = aes(fill = lbl_fixed),position = "dodge") + geom_text(aes(label=fdr_label), position=position_dodge(width=0.75), vjust=-0.5,size=4) + xlab("Locus") + ylab("Density\n-log10 (number of indels / bp)") +
  theme_bw() + theme(strip.text = element_text(size=16),axis.text.x=element_text(size=14,angle=60,hjust=1),axis.title.x = element_text(size=16),axis.text.y=element_text(size=14),axis.title.y = element_text(size=16),legend.title = element_text(size=16),legend.text = element_text(size=14)) + guides(fill=guide_legend(title="Variant")) + scale_fill_brewer(palette = "Set1") + scale_y_continuous(limits = c(0,9),breaks=seq(0,8,2)) + facet_wrap(~factor(population,v_pop_order),ncol=2)
ggsave(filename = "Indels_Fixation_bias_by_locus_across_populations_and_loci_of_interest.eps", path=output_workspace, width = 25/2.54, height = 20/2.54, device = cairo_ps)
ggsave(filename = "Indels_Fixation_bias_by_locus_across_populations_and_loci_of_interest.svg", path=output_workspace, width = 25/2.54, height = 20/2.54, device = svg)


#function that find the original and mutated codons of a variant
get_ref_and_mutated_codon <- function(the_chromosome,the_position,ref_nucl,new_nucl,cds_start){
  pos_in_codon <- ((the_position - cds_start + 1)%%3)+(3*as.integer(((the_position - cds_start + 1)%%3)==0))
  if (pos_in_codon==1){
    the_ref_codon <- paste0(ref_nucl,substr(x = as.character(genome_refseq[[the_chromosome]]),start = the_position+1,stop = the_position+2),sep="")
    the_mut_codon <- paste0(new_nucl,substr(x = as.character(genome_refseq[[the_chromosome]]),start = the_position+1,stop = the_position+2),sep="")
  }else if (pos_in_codon==2){
    the_ref_codon <- paste0(substr(x = as.character(genome_refseq[[the_chromosome]]),start = the_position-1,stop = the_position-1),ref_nucl,substr(x = as.character(genome_refseq[[the_chromosome]]),start = the_position+1,stop = the_position+1),sep="")
    the_mut_codon <- paste0(substr(x = as.character(genome_refseq[[the_chromosome]]),start = the_position-1,stop = the_position-1),new_nucl,substr(x = as.character(genome_refseq[[the_chromosome]]),start = the_position+1,stop = the_position+1),sep="")
  }else if (pos_in_codon==3){
    the_ref_codon <- paste0(substr(x = as.character(genome_refseq[[the_chromosome]]),start = the_position-2,stop = the_position-1),ref_nucl,sep="")
    the_mut_codon <- paste0(substr(x = as.character(genome_refseq[[the_chromosome]]),start = the_position-2,stop = the_position-1),new_nucl,sep="")
  }else{
    stop("Codon position must be between 1 and 3!!!")
  }
  return(list(ref_codon=the_ref_codon,mutated_codon=the_mut_codon))
}

#duplicate entries for variant hitting overlapping genes
for (i in df_ALL_variants){
  if ((!is.na(df_ALL_variants$gene_product[i]))&(grepl(pattern = ";",x = df_ALL_variants$gene_product[i],fixed = T))){
    current_gene_product <- df_ALL_variants$gene_product[i]
    current_gene_accession <- df_ALL_variants$accession_number[i]
    v_lst_overlapping_genes <- strsplit(x = current_gene_product,split = ";",fixed = T)[[1]]
    v_lst_overlapping_genes_accession <- strsplit(x = current_gene_product,split = ";",fixed = T)[[1]]
    current_df_ALL_variants_row <- df_ALL_variants[i,]
    for (j in 1:length(v_lst_overlapping_genes)){
      current_overlapping_gene
      new_row_df_ALL_variants <- current_df_ALL_variants_row
      new_row_df_ALL_variants$gene_product <- v_lst_overlapping_genes[j]
      new_row_df_ALL_variants$accession_number <- v_lst_overlapping_genes_accession[j]
      df_ALL_variants <- rbind(df_ALL_variants,new_row_df_ALL_variants)
    }
  } 
}
df_ALL_variants <- subset(df_ALL_variants,(!grepl(pattern = ";",x = gene_product,fixed = T)))

#Determine variant type 2 (sense. missense, nonsense, non-coding or indels)
df_ALL_variants$old_codon <- NA
df_ALL_variants$new_codon <- NA
df_ALL_variants$old_aa <- NA
df_ALL_variants$new_aa <- NA
df_ALL_variants$variant_type2 <- NA
df_ALL_variants$intron_length <- NA
df_ALL_variants$pos_in_intron <- NA
for (i in 1:nrow(df_ALL_variants)){
  if (df_ALL_variants$variant_type[i]!="SNV"){
    current_gene_df_loci_annotations <- subset(df_loci_annotations,(site_type=="gene")&(chr==df_ALL_variants$Chrom[i])&(start<=df_ALL_variants$Position[i])&(end>=df_ALL_variants$Position[i])& (!( (grepl(pattern = "ncrna",x = tolower(annotation),fixed = T)|grepl(pattern = "trna",x = tolower(annotation),fixed = T)) )))
    current_gene_chromosome <- current_gene_df_loci_annotations$chr[1]
    current_gene_min_pos <- min(current_gene_df_loci_annotations$start,na.rm = T)
    current_gene_max_pos <- max(current_gene_df_loci_annotations$end,na.rm = T)
    current_CDS_df_loci_annotations <- subset(df_loci_annotations,(site_type=="CDS")&(chr==df_ALL_variants$Chrom[i])&(start<=df_ALL_variants$Position[i])&(end>=df_ALL_variants$Position[i])& (!( (grepl(pattern = "ncrna",x = tolower(annotation),fixed = T)|grepl(pattern = "trna",x = tolower(annotation),fixed = T))&(!grepl(pattern = ";",x = tolower(annotation),fixed = T)) )))
    if (nrow(current_CDS_df_loci_annotations)>0){
      v_is_cds_prot_coding <- rep(F,nrow(current_CDS_df_loci_annotations))
      for (indx in 1:nrow(current_CDS_df_loci_annotations)){
        v_is_cds_prot_coding[indx] <- (!any(vapply(X = c("lncRNA", "snRNA", "sncRNA", "snoRNA", "rRNA","NCRNA","nRNA","NRNA","NORNA","RRNA"),FUN = function(the_noncoding_site_type) any(grepl(pattern = the_noncoding_site_type,x = current_CDS_df_loci_annotations$annotation[indx],fixed = T)),FUN.VALUE = F)))
      }
      current_CDS_df_loci_annotations <- subset(current_CDS_df_loci_annotations,v_is_cds_prot_coding)
    }
    df_ALL_variants$old_codon[i] <- NA
    df_ALL_variants$new_codon[i] <- NA
    df_ALL_variants$old_aa[i] <- NA
    df_ALL_variants$new_aa[i] <- NA
    if (nrow(current_CDS_df_loci_annotations)>0){
      df_ALL_variants$variant_type2[i] <- "Indel in protein-coding genes"
    }else if (nrow(subset(df_loci_annotations,(site_type=="intron")&(chr==current_gene_chromosome)&(start>=current_gene_min_pos)&(end<=current_gene_max_pos)))>0){
      current_intron_df_loci_annotations <- subset(df_loci_annotations,(site_type=="intron")&(chr==current_gene_chromosome)&(start>=current_gene_min_pos)&(end<=current_gene_max_pos))
      df_ALL_variants$intron_length[i] <- current_intron_df_loci_annotations$end[1] - current_intron_df_loci_annotations$start[1] + 1
      df_ALL_variants$pos_in_intron[i] <- df_ALL_variants$Position[i] - current_intron_df_loci_annotations$start[1] + 1
      df_ALL_variants$accession_number[i] <- substr(x = current_intron_df_loci_annotations$annotation[1],start = gregexpr(pattern = "ID=",text = current_intron_df_loci_annotations$annotation[1],fixed = T)[[1]]+3,stop = gregexpr(pattern = ":",text = current_intron_df_loci_annotations$annotation[1],fixed = T)[[1]]-1)
      df_ALL_variants$variant_type2[i] <- "Indel in introns"
    }else if (nrow(subset(df_loci_annotations,(site_type%in%c("three_prime_UTR","five_prime_UTR"))&(chr==current_gene_chromosome)&(start>=current_gene_min_pos)&(end<=current_gene_max_pos)))>0){
      current_utr_df_loci_annotations <- subset(df_loci_annotations,(site_type%in%c("three_prime_UTR","five_prime_UTR"))&(chr==current_gene_chromosome)&(start>=current_gene_min_pos)&(end<=current_gene_max_pos))
      df_ALL_variants$accession_number[i] <- substr(x = current_utr_df_loci_annotations$annotation[1],start = gregexpr(pattern = "ID=",text = current_utr_df_loci_annotations$annotation[1],fixed = T)[[1]]+3,stop = gregexpr(pattern = ":",text = current_utr_df_loci_annotations$annotation[1],fixed = T)[[1]]-1)
      df_ALL_variants$variant_type2[i] <- "Indel in UTR"
    }else{
      df_ALL_variants$variant_type2[i] <- "Indel in other noncoding regions"
    }
    
  }else if ((df_ALL_variants$variant_type[i]=="SNV") & (!(is.na(df_ALL_variants$accession_number[i]) | (is.na(df_ALL_variants$gene_product[i])))) & (!( (grepl(pattern = "ncrna",x = tolower(df_ALL_variants$accession_number[i]),fixed = T)|grepl(pattern = "trna",x = tolower(df_ALL_variants$accession_number[i]),fixed = T)) )) ){
    current_gene_df_loci_annotations <- subset(df_loci_annotations,(site_type=="gene")&(chr==df_ALL_variants$Chrom[i])&(start<=df_ALL_variants$Position[i])&(end>=df_ALL_variants$Position[i])& (!( (grepl(pattern = "ncrna",x = tolower(annotation),fixed = T)|grepl(pattern = "trna",x = tolower(annotation),fixed = T)) )))
    current_gene_chromosome <- current_gene_df_loci_annotations$chr[1]
    current_gene_min_pos <- min(current_gene_df_loci_annotations$start,na.rm = T)
    current_gene_max_pos <- max(current_gene_df_loci_annotations$end,na.rm = T)
    current_CDS_df_loci_annotations <- subset(df_loci_annotations,(site_type=="CDS")&(chr==df_ALL_variants$Chrom[i])&(start<=df_ALL_variants$Position[i])&(end>=df_ALL_variants$Position[i])& (!( (grepl(pattern = "ncrna",x = tolower(annotation),fixed = T)|grepl(pattern = "trna",x = tolower(annotation),fixed = T))&(!grepl(pattern = ";",x = tolower(annotation),fixed = T)) )))
    if (nrow(current_CDS_df_loci_annotations)>0){
      v_is_cds_prot_coding <- rep(F,nrow(current_CDS_df_loci_annotations))
      for (indx in 1:nrow(current_CDS_df_loci_annotations)){
        v_is_cds_prot_coding[indx] <- (!any(vapply(X = c("lncRNA", "snRNA", "sncRNA", "snoRNA", "rRNA","NCRNA","nRNA","NRNA","NORNA","RRNA"),FUN = function(the_noncoding_site_type) any(grepl(pattern = the_noncoding_site_type,x = current_CDS_df_loci_annotations$annotation[indx],fixed = T)),FUN.VALUE = F)))
      }
      current_CDS_df_loci_annotations <- subset(current_CDS_df_loci_annotations,v_is_cds_prot_coding)
    }
    
    if ((nrow(current_CDS_df_loci_annotations)>0)){
      the_min_start <- ifelse(test = nrow(current_CDS_df_loci_annotations)==0,yes = min(current_gene_df_loci_annotations$start,na.rm = T),no = min(current_CDS_df_loci_annotations$start,na.rm = T))
      res_get_ref_and_mutated_codon <- get_ref_and_mutated_codon(the_chromosome = current_gene_chromosome,the_position = df_ALL_variants$Position[i],ref_nucl = df_ALL_variants$Ref[i],new_nucl = df_ALL_variants$VarAllele[i],cds_start = the_min_start)
      df_ALL_variants$old_codon[i] <- res_get_ref_and_mutated_codon$ref_codon
      df_ALL_variants$new_codon[i] <- res_get_ref_and_mutated_codon$mutated_codon
      df_ALL_variants$old_aa[i] <- seqinr::translate(seq = unlist(strsplit(df_ALL_variants$old_codon[i],"")))
      df_ALL_variants$new_aa[i] <- seqinr::translate(seq = unlist(strsplit(df_ALL_variants$new_codon[i],"")))
      df_ALL_variants$variant_type2[i] <- ifelse(test = df_ALL_variants$old_aa[i]==df_ALL_variants$new_aa[i],yes = "synonymous",no=ifelse(test = df_ALL_variants$new_aa[i]=="*",yes = "nonsense",no="missense"))
    }else if (nrow(subset(df_loci_annotations,(site_type=="intron")&(chr==current_gene_chromosome)&(start>=current_gene_min_pos)&(end<=current_gene_max_pos)))>0){
      df_ALL_variants$old_codon[i] <- NA
      df_ALL_variants$new_codon[i] <- NA
      df_ALL_variants$old_aa[i] <- NA
      df_ALL_variants$new_aa[i] <- NA
      current_intron_df_loci_annotations <- subset(df_loci_annotations,(site_type=="intron")&(chr==current_gene_chromosome)&(start>=current_gene_min_pos)&(end<=current_gene_max_pos))
      df_ALL_variants$intron_length[i] <- current_intron_df_loci_annotations$end[1] - current_intron_df_loci_annotations$start[1] + 1
      df_ALL_variants$pos_in_intron[i] <- df_ALL_variants$Position[i] - current_intron_df_loci_annotations$start[1] + 1
      df_ALL_variants$variant_type2[i] <- "SNV in introns"
    }else{
      df_ALL_variants$old_codon[i] <- NA
      df_ALL_variants$new_codon[i] <- NA
      df_ALL_variants$old_aa[i] <- NA
      df_ALL_variants$new_aa[i] <- NA
      current_intron_df_loci_annotations <- subset(df_loci_annotations,(site_type%in%c("three_prime_UTR","five_prime_UTR"))&(chr==current_gene_chromosome)&(start>=current_gene_min_pos)&(end<=current_gene_max_pos))
      df_ALL_variants$intron_length[i] <- current_intron_df_loci_annotations$end[1] - current_intron_df_loci_annotations$start[1] + 1
      df_ALL_variants$pos_in_intron[i] <- df_ALL_variants$Position[i] - current_intron_df_loci_annotations$start[1] + 1
      df_ALL_variants$variant_type2[i] <- "SNV in UTR"
    }
  }else{
    df_ALL_variants$old_codon[i] <- NA
    df_ALL_variants$new_codon[i] <- NA
    df_ALL_variants$old_aa[i] <- NA
    df_ALL_variants$new_aa[i] <- NA
    df_ALL_variants$variant_type2[i] <- "SNV in other noncoding regions"
  }
  if (i%%1000==0){
    print(i)
  }
}

#Variant spectrum across populations ("nonsense","synonymous","missense","Indel in introns","SNV in introns","Indel in UTR","SNV in UTR","SNV in other noncoding regions","Indel in protein-coding genes" and "Indel in other noncoding regions") with separate grid for non-fixed and fixed variants
df_ALL_variants$concat_group_by_pop_is_fixed_and_vartype <- paste0(df_ALL_variants$population,"_",df_ALL_variants$is_fixed,"_",df_ALL_variants$variant_type2)
df_nb_UNIQUE_variants_per_pop <- df_ALL_variants %>%
  group_by(concat_group_by_pop_is_fixed_and_vartype) %>%
  dplyr::summarise(nb_unique_variants = length(unique(mutation))) %>%
  as.data.frame()
v_nb_non_fixed_variants_per_pop <- table(subset(unique(df_ALL_variants[,c("population","mutation","is_fixed","variant_type2")]),!is_fixed)$population)
v_nb_fixed_variants_per_pop <- table(subset(unique(df_ALL_variants[,c("population","mutation","is_fixed","variant_type2")]),is_fixed)$population)
df_nb_UNIQUE_variants_per_pop$population <- sapply(X = strsplit(df_nb_UNIQUE_variants_per_pop$concat_group_by_pop_is_fixed_and_vartype,split = "_"),FUN = function(x) x[1])
df_nb_UNIQUE_variants_per_pop$is_fixed <- as.logical(sapply(X = strsplit(df_nb_UNIQUE_variants_per_pop$concat_group_by_pop_is_fixed_and_vartype,split = "_"),FUN = function(x) x[2]))
df_nb_UNIQUE_variants_per_pop$var_type2 <- sapply(X = strsplit(df_nb_UNIQUE_variants_per_pop$concat_group_by_pop_is_fixed_and_vartype,split = "_"),FUN = function(x) x[3])
df_nb_UNIQUE_variants_per_pop <- df_nb_UNIQUE_variants_per_pop[,c("population","is_fixed","var_type2","nb_unique_variants")]
df_nb_UNIQUE_variants_per_pop$proportion <- df_nb_UNIQUE_variants_per_pop$nb_unique_variants/ifelse(test = df_nb_UNIQUE_variants_per_pop$is_fixed,yes = v_nb_fixed_variants_per_pop[df_nb_UNIQUE_variants_per_pop$population],no = v_nb_non_fixed_variants_per_pop[df_nb_UNIQUE_variants_per_pop$population])
v_fdr_variant_fixation_bias_across_pops <- p.adjust(vapply(X = v_pop_order, FUN = function(the_pop) chisq.test(apply(X = as.matrix(reshape2::acast(subset(df_nb_UNIQUE_variants_per_pop,population==the_pop), var_type2 ~ is_fixed,value.var = 'nb_unique_variants', fill = '0')),MARGIN = 2,FUN = as.numeric))$p.value,FUN.VALUE = 0.0),"fdr")
v_label_top_of_cols <- ifelse(test = v_fdr_variant_fixation_bias_across_pops>=0.05,yes = "N.S.",
                              no=ifelse(test = (v_fdr_variant_fixation_bias_across_pops>=0.01)&(v_fdr_variant_fixation_bias_across_pops<0.05),yes = "*",no =  ifelse(test = (v_fdr_variant_fixation_bias_across_pops>=0.001)&(v_fdr_variant_fixation_bias_across_pops<0.01),yes = "**",no =  "***")))
df_nb_UNIQUE_variants_per_pop$fdr_label <- ifelse(test = df_nb_UNIQUE_variants_per_pop$is_fixed,yes="",no=v_label_top_of_cols[df_nb_UNIQUE_variants_per_pop$population])
v_first_occ_to_add <- v_pop_order
for (i in 1:nrow(df_nb_UNIQUE_variants_per_pop)){
  if (df_nb_UNIQUE_variants_per_pop$population[i]%in%v_first_occ_to_add){
    df_nb_UNIQUE_variants_per_pop$fdr_label[i] <- df_nb_UNIQUE_variants_per_pop$fdr_label[i]
    v_first_occ_to_add <- v_first_occ_to_add[v_first_occ_to_add!=df_nb_UNIQUE_variants_per_pop$population[i]]
  }else{
    df_nb_UNIQUE_variants_per_pop$fdr_label[i] <- ""
  }
}

df_nb_UNIQUE_variants_per_pop$lbl_fixed <- factor(ifelse(df_nb_UNIQUE_variants_per_pop$is_fixed,yes = "Fixed (VAF >= 0.75)",no = "Not fixed (VAF < 0.75)"),levels = c("Not fixed (VAF < 0.75)","Fixed (VAF >= 0.75)"))
palette_variant_type2 <- rev(c(brewer.pal(n=8,name = "Dark2"),"#0C7BDC","#E41A1C"))
names(palette_variant_type2) <- c("nonsense","synonymous","missense","Indel in introns","SNV in introns","Indel in UTR","SNV in UTR","SNV in other noncoding regions","Indel in protein-coding genes","Indel in other noncoding regions")
ggplot(data = df_nb_UNIQUE_variants_per_pop,mapping = aes(x=factor(population,v_pop_order),y=proportion)) + geom_col(mapping = aes(fill=factor(var_type2,levels=c("nonsense","synonymous","missense","Indel in introns","SNV in introns","Indel in UTR","SNV in UTR","SNV in other noncoding regions","Indel in protein-coding genes","Indel in other noncoding regions")))) + geom_text(aes(label=fdr_label), position=position_dodge(width=0.75), vjust=-14,size=6) + xlab("Population") + ylab("Proportion") +
  theme_bw() + theme(strip.text = element_text(size=16),axis.text.x=element_text(size=14),axis.title.x = element_text(size=16),axis.text.y=element_text(size=14),axis.title.y = element_text(size=16),legend.title = element_text(size=16),legend.text = element_text(size=14)) + guides(fill=guide_legend(title="Variant")) + scale_y_continuous(limits = c(0,1.1),breaks = seq(0,1,0.2)) + scale_fill_manual(values = palette_variant_type2) + facet_wrap(~lbl_fixed,ncol=1,scales="free_y")
ggsave(filename = "Fixation_bias_by_variant_type2_across_populations.eps", path=output_workspace, width = 25/2.54, height = 20/2.54, device = cairo_ps)
ggsave(filename = "Fixation_bias_by_variant_type2_across_populations.svg", path=output_workspace, width = 25/2.54, height = 20/2.54, device = svg)

#Recreate the  figure in Bergström et al., 2014
palette_variant_type_Bergstrom_et_al_2014 <- c("#CD9462","#CD7BC5","#00A473","#E68B08","#0073B4")
df_nb_UNIQUE_variants_per_pop$var_type_Bergstrom_et_al <- ifelse(test = df_nb_UNIQUE_variants_per_pop$var_type2 %in% c("SNV in introns","SNV in UTR","SNV in other noncoding regions","Indel in other noncoding regions","Indel in UTR","Indel in introns"),yes = "noncoding",no = ifelse(test = df_nb_UNIQUE_variants_per_pop$var_type2 %in% c("Indel in protein-coding genes"),yes = "indel in coding regions",no = df_nb_UNIQUE_variants_per_pop$var_type2))
ggplot(data = df_nb_UNIQUE_variants_per_pop,mapping = aes(x=factor(population,v_pop_order),y=proportion)) + geom_col(mapping = aes(fill=factor(var_type_Bergstrom_et_al,levels=c("noncoding","synonymous","missense","nonsense","indel in coding regions")))) + xlab("Population") + ylab("Proportion") +
  theme_bw() + theme(strip.text = element_text(size=16),axis.text.x=element_text(size=14),axis.title.x = element_text(size=16),axis.text.y=element_text(size=14),axis.title.y = element_text(size=16),legend.title = element_text(size=16),legend.text = element_text(size=14)) + guides(fill=guide_legend(title="Variant")) + scale_y_continuous(limits = c(0,1.1),breaks = seq(0,1,0.2)) + scale_fill_manual(values = palette_variant_type_Bergstrom_et_al_2014) + facet_wrap(~lbl_fixed,ncol=1,scales="free_y")
ggsave(filename = "Fixation_bias_by_variant_type_Bergstrom_et_al_2014.eps", path=output_workspace, width = 25/2.54, height = 20/2.54, device = cairo_ps)
ggsave(filename = "Fixation_bias_by_variant_type_Bergstrom_et_al_2014.svg", path=output_workspace, width = 25/2.54, height = 20/2.54, device = svg)

#add the null model expectations in the dataframe
  #for pombe
pombe_genome_size <- sum(vapply(X = 1:length(genome_refseq),FUN = function(the_chr) nchar(genome_refseq[[the_chr]]),FUN.VALUE = 0))
the_nb_substitutions_expected_in_pombe_null_model <- round(the_mean_substitution_rate*pombe_genome_size*(1E4))
the_nb_fixed_indels_expected_in_pombe_null_model <- round(the_mean_indel_fixation_rate*pombe_genome_size*(1E4))
the_nb_fixed_total_variants_in_pombe_null_model <- the_nb_substitutions_expected_in_pombe_null_model + the_nb_fixed_indels_expected_in_pombe_null_model
the_nb_fixed_noncoding_variants_expected_in_pombe_null_model <- ((the_nb_fixed_indels_expected_in_pombe_null_model/the_nb_fixed_total_variants_in_pombe_null_model)*0.4*the_nb_fixed_total_variants_in_pombe_null_model)+((the_nb_substitutions_expected_in_pombe_null_model/the_nb_fixed_total_variants_in_pombe_null_model)*0.4*the_nb_fixed_total_variants_in_pombe_null_model)
the_nb_fixed_synonymous_substitutions_expected_in_pombe_null_model <- (the_nb_substitutions_expected_in_pombe_null_model/the_nb_fixed_total_variants_in_pombe_null_model)*0.6*(138/576)*the_nb_fixed_total_variants_in_pombe_null_model
the_nb_fixed_missense_substitutions_expected_in_pombe_null_model <- (the_nb_substitutions_expected_in_pombe_null_model/the_nb_fixed_total_variants_in_pombe_null_model)*0.6*(415/576)*the_nb_fixed_total_variants_in_pombe_null_model
the_nb_fixed_nonsense_substitutions_expected_in_pombe_null_model <- (the_nb_substitutions_expected_in_pombe_null_model/the_nb_fixed_total_variants_in_pombe_null_model)*0.6*(23/576)*the_nb_fixed_total_variants_in_pombe_null_model
the_nb_fixed_indels_in_coding_regions_expected_in_pombe_null_model <- (the_nb_fixed_indels_expected_in_pombe_null_model/the_nb_fixed_total_variants_in_pombe_null_model)*0.6*the_nb_fixed_total_variants_in_pombe_null_model
df_nb_UNIQUE_variants_per_pop <- rbind(df_nb_UNIQUE_variants_per_pop,data.frame(population="Null\nmodel",is_fixed=T,var_type2=NA,
                                                                                nb_unique_variants= c(the_nb_fixed_noncoding_variants_expected_in_pombe_null_model,the_nb_fixed_synonymous_substitutions_expected_in_pombe_null_model,the_nb_fixed_missense_substitutions_expected_in_pombe_null_model,the_nb_fixed_nonsense_substitutions_expected_in_pombe_null_model,the_nb_fixed_indels_in_coding_regions_expected_in_pombe_null_model),
                                                                                proportion=c(the_nb_fixed_noncoding_variants_expected_in_pombe_null_model,the_nb_fixed_synonymous_substitutions_expected_in_pombe_null_model,the_nb_fixed_missense_substitutions_expected_in_pombe_null_model,the_nb_fixed_nonsense_substitutions_expected_in_pombe_null_model,the_nb_fixed_indels_in_coding_regions_expected_in_pombe_null_model)/the_nb_fixed_total_variants_in_pombe_null_model,
                                                                                fdr_label=NA,lbl_fixed="Fixed (VAF >= 0.75)",var_type_Bergstrom_et_al=c("noncoding","synonymous","missense","nonsense","indel in coding regions"),stringsAsFactors = F))

ggplot(data = subset(df_nb_UNIQUE_variants_per_pop,lbl_fixed!="Not fixed (VAF < 0.75)"),mapping = aes(x=factor(population,c("Null\nmodel",v_pop_order) ),y=proportion)) + geom_col(mapping = aes(fill=factor(var_type_Bergstrom_et_al,levels=c("noncoding","synonymous","missense","nonsense","indel in coding regions")))) + xlab("Population") + ylab("Proportion") +
  theme_bw() + theme(strip.text = element_text(size=16),axis.text.x=element_text(size=18),axis.title.x = element_text(size=24),axis.text.y=element_text(size=18),axis.title.y = element_text(size=24),legend.title = element_text(size=24),legend.text = element_text(size=24)) + guides(fill=guide_legend(title="Variant")) + scale_y_continuous(limits = c(0,1.01),breaks = seq(0,1,0.2)) + scale_fill_manual(values = palette_variant_type_Bergstrom_et_al_2014)
ggsave(filename = "Substitution_spectrum_Pombe_by_variant_type_Bergstrom_et_al_2014.eps", path=output_workspace, width = 25/2.54, height = 12/2.54, device = cairo_ps)
ggsave(filename = "Substitution_spectrum_Pombe_by_variant_type_Bergstrom_et_al_2014.svg", path=output_workspace, width = 25/2.54, height = 12/2.54, device = svg)
ggsave(filename = "Substitution_spectrum_Pombe_by_variant_type_Bergstrom_et_al_2014.png", path=output_workspace, width = 25/2.54, height = 12/2.54, device = png)

  #for cerevisiae
cerevisiae_genome_size <- 12.07E6
the_nb_substitutions_expected_in_cerevisiae_null_model <- round((2.5E-10)*cerevisiae_genome_size*(1E4)) 
the_nb_fixed_indels_expected_in_cerevisiae_null_model <- round((5E-11)*cerevisiae_genome_size*(1E4))
the_nb_fixed_total_variants_in_cerevisiae_null_model <- the_nb_substitutions_expected_in_cerevisiae_null_model + the_nb_fixed_indels_expected_in_cerevisiae_null_model
the_nb_fixed_noncoding_variants_expected_in_cerevisiae_null_model <- ((the_nb_fixed_indels_expected_in_cerevisiae_null_model/the_nb_fixed_total_variants_in_cerevisiae_null_model)*0.3*the_nb_fixed_total_variants_in_cerevisiae_null_model)+((the_nb_substitutions_expected_in_cerevisiae_null_model/the_nb_fixed_total_variants_in_cerevisiae_null_model)*0.3*the_nb_fixed_total_variants_in_cerevisiae_null_model)
the_nb_fixed_synonymous_substitutions_expected_in_cerevisiae_null_model <- (the_nb_substitutions_expected_in_cerevisiae_null_model/the_nb_fixed_total_variants_in_cerevisiae_null_model)*0.7*(138/576)*the_nb_fixed_total_variants_in_cerevisiae_null_model
the_nb_fixed_missense_substitutions_expected_in_cerevisiae_null_model <- (the_nb_substitutions_expected_in_cerevisiae_null_model/the_nb_fixed_total_variants_in_cerevisiae_null_model)*0.7*(415/576)*the_nb_fixed_total_variants_in_cerevisiae_null_model
the_nb_fixed_nonsense_substitutions_expected_in_cerevisiae_null_model <- (the_nb_substitutions_expected_in_cerevisiae_null_model/the_nb_fixed_total_variants_in_cerevisiae_null_model)*0.7*(23/576)*the_nb_fixed_total_variants_in_cerevisiae_null_model
the_nb_fixed_indels_in_coding_regions_expected_in_cerevisiae_null_model <- (the_nb_fixed_indels_expected_in_cerevisiae_null_model/the_nb_fixed_total_variants_in_cerevisiae_null_model)*0.7*the_nb_fixed_total_variants_in_cerevisiae_null_model
df_nb_UNIQUE_variants_per_cerevisiae_pop <- data.frame(population="Null\nmodel",is_fixed=T,var_type2=NA,
                                                                                nb_unique_variants= c(the_nb_fixed_noncoding_variants_expected_in_cerevisiae_null_model,the_nb_fixed_synonymous_substitutions_expected_in_cerevisiae_null_model,the_nb_fixed_missense_substitutions_expected_in_cerevisiae_null_model,the_nb_fixed_nonsense_substitutions_expected_in_cerevisiae_null_model,the_nb_fixed_indels_in_coding_regions_expected_in_cerevisiae_null_model),
                                                                                proportion=c(the_nb_fixed_noncoding_variants_expected_in_cerevisiae_null_model,the_nb_fixed_synonymous_substitutions_expected_in_cerevisiae_null_model,the_nb_fixed_missense_substitutions_expected_in_cerevisiae_null_model,the_nb_fixed_nonsense_substitutions_expected_in_cerevisiae_null_model,the_nb_fixed_indels_in_coding_regions_expected_in_cerevisiae_null_model)/the_nb_fixed_total_variants_in_cerevisiae_null_model,
                                                                                fdr_label=NA,lbl_fixed="Fixed (VAF >= 0.75)",var_type_Bergstrom_et_al=c("noncoding","synonymous","missense","nonsense","indel in coding regions"),stringsAsFactors = F)

ggplot(data = subset(df_nb_UNIQUE_variants_per_cerevisiae_pop,lbl_fixed!="Not fixed (VAF < 0.75)"),mapping = aes(x=factor(population,c("Null\nmodel",v_pop_order) ),y=proportion)) + geom_col(mapping = aes(fill=factor(var_type_Bergstrom_et_al,levels=c("noncoding","synonymous","missense","nonsense","indel in coding regions")))) + xlab("Population") + ylab("Proportion") +
  theme_bw() + theme(strip.text = element_text(size=16),axis.text.x=element_text(size=18),axis.title.x = element_text(size=24),axis.text.y=element_text(size=18),axis.title.y = element_text(size=24),legend.title = element_text(size=24),legend.text = element_text(size=24)) + guides(fill=guide_legend(title="Variant")) + scale_y_continuous(limits = c(0,1.01),breaks = seq(0,1,0.2)) + scale_fill_manual(values = palette_variant_type_Bergstrom_et_al_2014)
ggsave(filename = "Substitution_spectrum_cerevisiae_by_variant_type_Bergstrom_et_al_2014.eps", path=output_workspace, width = 25/2.54, height = 12/2.54, device = cairo_ps)
ggsave(filename = "Substitution_spectrum_cerevisiae_by_variant_type_Bergstrom_et_al_2014.svg", path=output_workspace, width = 25/2.54, height = 12/2.54, device = svg)
ggsave(filename = "Substitution_spectrum_cerevisiae_by_variant_type_Bergstrom_et_al_2014.png", path=output_workspace, width = 25/2.54, height = 12/2.54, device = png)


#Time series of the number of substitutions per population + estimation of the substitution rate
df_substitutions <- subset(df_variants,is_fixed)
#only keep substitutions that are also fixed at subsequent time point and NOT at max time point (because there is no next time point at max)
df_substitutions$is_fixed_at_next_timepoint <- vapply(X = 1:nrow(df_substitutions),FUN = function(i) nrow(subset(df_substitutions,(pop_mut==df_substitutions$pop_mut[i])&(pop_timepoint==v_next_pop_timepoint[df_substitutions$pop_timepoint[i]])))>0,FUN.VALUE = T)
#|(nrow(subset(df_substitutions,(str_time==v_pops_max_str_timepoint[df_substitutions$population[i]]) & (nrow(subset(df_substitutions,(population==df_substitutions$population[i])&(mutation==df_substitutions$mutation[i])&(VarFreq>0.75)&(str_time==v_previous_timepoint[df_substitutions$str_time[i]])))>0)))>0))
#|((str_time==v_previous_timepoint[df_substitutions$str_time[i]])|(str_time==v_pops_min_str_timepoint[df_substitutions$population[i]]))
df_substitutions <- subset(df_substitutions,is_fixed_at_next_timepoint)
#df_substitutions <- subset(df_substitutions,str_time!=unname(v_pops_max_str_timepoint[df_substitutions$population])) 
df_nb_UNIQUE_variants_per_pop <- df_substitutions %>%
  group_by(population,int_time,.add=TRUE) %>%
  dplyr::summarise(nb_unique_variants = length(unique(mutation))) %>%
  as.data.frame()
#only keep populations with at least 5 timepoints and constant increase
v_nb_timepoints_per_population <- vapply(X = v_pop_order,FUN = function(the_pop) length(unique(subset(df_substitutions,population==the_pop)$str_time )),FUN.VALUE = 0)
v_lst_pops_used_for_substitution_rate_estimation <- names(v_nb_timepoints_per_population)[v_nb_timepoints_per_population>=5]
df_nb_UNIQUE_variants_per_pop <- subset(df_nb_UNIQUE_variants_per_pop,population%in%v_lst_pops_used_for_substitution_rate_estimation)
df_nb_UNIQUE_variants_per_pop <- subset(df_nb_UNIQUE_variants_per_pop,!population%in%c("H9","G11"))
v_lst_pops_used_for_substitution_rate_estimation <- sort(unique(df_nb_UNIQUE_variants_per_pop$population))
df_nb_UNIQUE_variants_per_pop$density <- df_nb_UNIQUE_variants_per_pop$nb_unique_variants/(sum(unname(v_chrs_length)))
df_nb_UNIQUE_variants_per_pop$log10_density <- log10(df_nb_UNIQUE_variants_per_pop$density)
#start populations at 0 fixed indels
#for (current_pop in v_pop_order){
#  df_nb_UNIQUE_variants_per_pop <- rbind(df_nb_UNIQUE_variants_per_pop,data.frame(population=current_pop,int_time=0,nb_unique_variants=0,density=0,log10_density=0,stringsAsFactors = F))
#}
fit_for_subst_rate <- lm(data = df_nb_UNIQUE_variants_per_pop,formula = density~int_time)
#remove outliers
v_cook_dist <- unname(cooks.distance(fit_for_subst_rate))
id_row_outliers <- which(v_cook_dist>=0.5)
if (length(id_row_outliers)>0){
  df_nb_UNIQUE_variants_per_pop <- df_nb_UNIQUE_variants_per_pop[-id_row_outliers,]
}else{
  df_nb_UNIQUE_variants_per_pop <- df_nb_UNIQUE_variants_per_pop
}
fit_for_subst_rate <- lm(data = df_nb_UNIQUE_variants_per_pop,formula = density~int_time)
adj_r_sq <- formatC(summary(fit_for_subst_rate)$adj.r.squared, format = "e", digits = 3)
slope <-formatC(summary(fit_for_subst_rate)$coefficients[,1][2], format = "e", digits = 3)
p_val <- ifelse(test = "Iter"%in%colnames(summary(fit_for_subst_rate)$coefficients),yes=ifelse(test = unname(summary(fit_for_subst_rate)$coefficients[,3][2])<(2E-4),yes = "<2E-4",no = formatC(unname(summary(fit_for_subst_rate)$coefficients[,3][2]), format = "e", digits = 3)),no=ifelse(test = unname(summary(fit_for_subst_rate)$coefficients[,4][2])<(2e-16),yes = "<2e-16",no = formatC(unname(summary(fit_for_subst_rate)$coefficients[,4][2]), format = "e", digits = 3)))

ggplot(data = subset(df_nb_UNIQUE_variants_per_pop,!population%in%c("H9","G11")),mapping = aes(x=int_time,y=density)) + 
  geom_line(aes(col=factor(population,levels = v_pop_order) )) +
  stat_smooth(method = "lm", col = "red") +
  labs(title = paste("Adj R2 = ",adj_r_sq,
                     " Slope =",slope,
                     " P: ",p_val))+ theme(plot.title=element_text(hjust=0,size=12))+
  xlab("Time (Generation)") + ylab("Density\n(number of substitutions / bp)") +
  theme_bw() + theme(strip.text = element_text(size=16),axis.text.x=element_text(size=14,angle=60,hjust=1),axis.title.x = element_text(size=16),axis.text.y=element_text(size=14),axis.title.y = element_text(size=16),legend.title = element_text(size=16),legend.text = element_text(size=14)) +  scale_x_continuous(limits = c(0,10000),breaks = seq(0,10000,1000)) + guides(col=guide_legend(title="Population")) + scale_color_brewer(palette = "Paired") 
ggsave(filename = "Substitution_rate_across_populations.eps", path=output_workspace, width = 25/2.54, height = 20/2.54, device = cairo_ps)
ggsave(filename = "Substitution_rate_across_populations.svg", path=output_workspace, width = 25/2.54, height = 20/2.54, device = svg)
v_per_pop_adj_r_sq <- as.numeric(vapply(X = v_lst_pops_used_for_substitution_rate_estimation,FUN = function(the_current_pop) formatC(summary(lm(data = subset(df_nb_UNIQUE_variants_per_pop,population==the_current_pop),formula = density~int_time))$adj.r.squared, format = "e", digits = 3),FUN.VALUE = ""))
v_per_pop_slope <-as.numeric(vapply(X = v_lst_pops_used_for_substitution_rate_estimation,FUN = function(the_current_pop) formatC(summary(lm(data = subset(df_nb_UNIQUE_variants_per_pop,population==the_current_pop),formula = density~int_time))$coefficients[,1][2], format = "e", digits = 3),FUN.VALUE = ""))
v_per_pop_p_val <- as.numeric(vapply(X = v_lst_pops_used_for_substitution_rate_estimation,FUN = function(the_current_pop) ifelse(test = "Iter"%in%colnames(summary(lm(data = subset(df_nb_UNIQUE_variants_per_pop,population==the_current_pop),formula = density~int_time))$coefficients),yes=ifelse(test = unname(summary(lm(data = subset(df_nb_UNIQUE_variants_per_pop,population==the_current_pop),formula = density~int_time))$coefficients[,3][2])<(2E-4),yes = "<2E-4",no = formatC(unname(summary(lm(data = subset(df_nb_UNIQUE_variants_per_pop,population==the_current_pop),formula = density~int_time))$coefficients[,3][2]), format = "e", digits = 3)),no=ifelse(test = unname(summary(lm(data = subset(df_nb_UNIQUE_variants_per_pop,population==the_current_pop),formula = density~int_time))$coefficients[,4][2])<(2e-16),yes = "<2e-16",no = formatC(unname(summary(lm(data = subset(df_nb_UNIQUE_variants_per_pop,population==the_current_pop),formula = density~int_time))$coefficients[,4][2]), format = "e", digits = 3))),FUN.VALUE = ""))
print(paste0("The minimum substitution rate is:",min(v_per_pop_slope,na.rm=T)))
print(paste0("The median substitution rate is:",median(v_per_pop_slope,na.rm=T)))
print(paste0("The mean substitution rate is:",mean(v_per_pop_slope,na.rm=T)))
the_mean_substitution_rate <- mean(v_per_pop_slope,na.rm=T)
print(paste0("The maximum substitution rate is:",max(v_per_pop_slope,na.rm=T)))
print(paste0("The standard deviation of the substitution rate is:",sd(v_per_pop_slope,na.rm=T)))

#Time series of the number of fixed indels per population + estimation of the indel fixation rate
df_fixed_indels <- subset(df_indel_variants,is_fixed)
#only keep fixed indels that are also fixed at subsequent time point and NOT at max time point (because there is no next time point at max)
df_fixed_indels$is_fixed_at_next_timepoint <- vapply(X = 1:nrow(df_fixed_indels),FUN = function(i) nrow(subset(df_fixed_indels,(pop_mut==df_fixed_indels$pop_mut[i])&(pop_timepoint==v_next_pop_timepoint[df_fixed_indels$pop_timepoint[i]])))>0,FUN.VALUE = T)
#|(nrow(subset(df_fixed_indels,(str_time==v_pops_max_str_timepoint[df_fixed_indels$population[i]]) & (nrow(subset(df_fixed_indels,(population==df_fixed_indels$population[i])&(mutation==df_fixed_indels$mutation[i])&(VarFreq>0.75)&(str_time==v_previous_timepoint[df_fixed_indels$str_time[i]])))>0)))>0))
#|((str_time==v_previous_timepoint[df_fixed_indels$str_time[i]])|(str_time==v_pops_min_str_timepoint[df_fixed_indels$population[i]]))
df_fixed_indels <- subset(df_fixed_indels,is_fixed_at_next_timepoint)
#df_fixed_indels <- subset(df_fixed_indels,str_time!=unname(v_pops_max_str_timepoint[df_fixed_indels$population])) 
df_nb_UNIQUE_variants_per_pop <- df_fixed_indels %>%
  group_by(population,int_time,.add=TRUE) %>%
  dplyr::summarise(nb_unique_variants = length(unique(mutation))) %>%
  as.data.frame()
#only keep populations with at least 5 timepoints
v_nb_timepoints_per_population <- vapply(X = v_pop_order,FUN = function(the_pop) length(unique(subset(df_fixed_indels,population==the_pop)$str_time )),FUN.VALUE = 0)
v_lst_pops_used_for_fixed_indel_rate_estimation <- names(v_nb_timepoints_per_population)[v_nb_timepoints_per_population>=5]
df_nb_UNIQUE_variants_per_pop <- subset(df_nb_UNIQUE_variants_per_pop,population%in%v_lst_pops_used_for_fixed_indel_rate_estimation)
df_nb_UNIQUE_variants_per_pop$density <- df_nb_UNIQUE_variants_per_pop$nb_unique_variants/(sum(unname(v_chrs_length)))
df_nb_UNIQUE_variants_per_pop$log10_density <- log10(df_nb_UNIQUE_variants_per_pop$density)
#start populations at 0 fixed indels
#for (current_pop in v_pop_order){
#  df_nb_UNIQUE_variants_per_pop <- rbind(df_nb_UNIQUE_variants_per_pop,data.frame(population=current_pop,int_time=0,nb_unique_variants=0,density=0,log10_density=0,stringsAsFactors = F))
#}
fit_for_indel_fixation_rate <- lm(data = df_nb_UNIQUE_variants_per_pop,formula = density~int_time)
#remove outliers
v_cook_dist <- unname(cooks.distance(fit_for_indel_fixation_rate))
id_row_outliers <- which(v_cook_dist>=0.5)
if (length(id_row_outliers)>0){
  df_nb_UNIQUE_variants_per_pop <- df_nb_UNIQUE_variants_per_pop[-id_row_outliers,]
}else{
  df_nb_UNIQUE_variants_per_pop <- df_nb_UNIQUE_variants_per_pop
}
fit_for_indel_fixation_rate <- lm(data = df_nb_UNIQUE_variants_per_pop,formula = density~int_time)
adj_r_sq <- formatC(summary(fit_for_indel_fixation_rate)$adj.r.squared, format = "e", digits = 3)
slope <-formatC(summary(fit_for_indel_fixation_rate)$coefficients[,1][2], format = "e", digits = 3)
p_val <- ifelse(test = "Iter"%in%colnames(summary(fit_for_indel_fixation_rate)$coefficients),yes=ifelse(test = unname(summary(fit_for_indel_fixation_rate)$coefficients[,3][2])<(2E-4),yes = "<2E-4",no = formatC(unname(summary(fit_for_indel_fixation_rate)$coefficients[,3][2]), format = "e", digits = 3)),no=ifelse(test = unname(summary(fit_for_indel_fixation_rate)$coefficients[,4][2])<(2e-16),yes = "<2e-16",no = formatC(unname(summary(fit_for_indel_fixation_rate)$coefficients[,4][2]), format = "e", digits = 3)))

ggplot(data = df_nb_UNIQUE_variants_per_pop,mapping = aes(x=int_time,y=density)) + 
  geom_line(aes(col=factor(population,levels = v_pop_order) )) +
  stat_smooth(method = "lm", col = "red") +
  labs(title = paste("Adj R2 = ",adj_r_sq,
                     " Slope =",slope,
                     " P: ",p_val))+ theme(plot.title=element_text(hjust=0,size=12))+
  xlab("Time (Generation)") + ylab("Density\n(number of fixed indels / bp)") +
  theme_bw() + theme(strip.text = element_text(size=16),axis.text.x=element_text(size=14,angle=60,hjust=1),axis.title.x = element_text(size=16),axis.text.y=element_text(size=14),axis.title.y = element_text(size=16),legend.title = element_text(size=16),legend.text = element_text(size=14)) +  scale_x_continuous(limits = c(0,10000),breaks = seq(0,10000,1000)) + guides(col=guide_legend(title="Population")) + scale_color_brewer(palette = "Paired") 

ggsave(filename = "Indel_fixation_rate_across_populations.eps", path=output_workspace, width = 25/2.54, height = 20/2.54, device = cairo_ps)
ggsave(filename = "Indel_fixation_rate_across_populations.svg", path=output_workspace, width = 25/2.54, height = 20/2.54, device = svg)
v_per_pop_adj_r_sq <- as.numeric(vapply(X = v_lst_pops_used_for_fixed_indel_rate_estimation,FUN = function(the_current_pop) formatC(summary(lm(data = subset(df_nb_UNIQUE_variants_per_pop,population==the_current_pop),formula = density~int_time))$adj.r.squared, format = "e", digits = 3),FUN.VALUE = ""))
v_per_pop_slope <-as.numeric(vapply(X = v_lst_pops_used_for_fixed_indel_rate_estimation,FUN = function(the_current_pop) formatC(summary(lm(data = subset(df_nb_UNIQUE_variants_per_pop,population==the_current_pop),formula = density~int_time))$coefficients[,1][2], format = "e", digits = 3),FUN.VALUE = ""))
v_per_pop_p_val <- as.numeric(vapply(X = v_lst_pops_used_for_fixed_indel_rate_estimation,FUN = function(the_current_pop) ifelse(test = "Iter"%in%colnames(summary(lm(data = subset(df_nb_UNIQUE_variants_per_pop,population==the_current_pop),formula = density~int_time))$coefficients),yes=ifelse(test = unname(summary(lm(data = subset(df_nb_UNIQUE_variants_per_pop,population==the_current_pop),formula = density~int_time))$coefficients[,3][2])<(2E-4),yes = "<2E-4",no = formatC(unname(summary(lm(data = subset(df_nb_UNIQUE_variants_per_pop,population==the_current_pop),formula = density~int_time))$coefficients[,3][2]), format = "e", digits = 3)),no=ifelse(test = unname(summary(lm(data = subset(df_nb_UNIQUE_variants_per_pop,population==the_current_pop),formula = density~int_time))$coefficients[,4][2])<(2e-16),yes = "<2e-16",no = formatC(unname(summary(lm(data = subset(df_nb_UNIQUE_variants_per_pop,population==the_current_pop),formula = density~int_time))$coefficients[,4][2]), format = "e", digits = 3))),FUN.VALUE = ""))
print(paste0("The minimum indel fixation rate is:",min(v_per_pop_slope,na.rm=T)))
print(paste0("The median indel fixation rate is:",median(v_per_pop_slope,na.rm=T)))
print(paste0("The mean indel fixation rate is:",mean(v_per_pop_slope,na.rm=T)))
the_mean_indel_fixation_rate <- mean(v_per_pop_slope,na.rm=T)
print(paste0("The maximum indel fixation rate is:",max(v_per_pop_slope,na.rm=T)))
print(paste0("The standard deviation of the indel fixation rate is:",sd(v_per_pop_slope,na.rm=T)))

# ggplot(data = df_nb_UNIQUE_variants_per_pop,mapping = aes(x=int_time,y=density)) + 
#   geom_point(aes(col=factor(population,levels = v_pop_order) )) +
#   stat_smooth(mapping = aes(col=factor(population,levels = v_pop_order)), method = "lm",se=F) +
#   labs(title = paste("Adj R2 = ",adj_r_sq,
#                      " Slope =",slope,
#                      " P: ",p_val))+ theme(plot.title=element_text(hjust=0,size=12))+
#   xlab("Time (Generation)") + ylab("Density\n(number of substitutions / bp)") +
#   theme_bw() + theme(strip.text = element_text(size=16),axis.text.x=element_text(size=14,angle=60,hjust=1),axis.title.x = element_text(size=16),axis.text.y=element_text(size=14),axis.title.y = element_text(size=16),legend.title = element_text(size=16),legend.text = element_text(size=14)) + guides(col=guide_legend(title="Population")) + scale_color_brewer(palette = "Paired") 

#Find genomic regions length by accession_number
v_lst_genomic_regions_accession_number <- sort(unique(subset(df_loci_annotations,!is.na(accession_number))$accession_number))
v_genomic_region_length <- rep(NA,length(v_lst_genomic_regions_accession_number))
names(v_genomic_region_length) <- v_lst_genomic_regions_accession_number
for (current_accession_nb in names(v_lst_genomic_regions_accession_number)){
  current_df_loci_annotation <- subset(df_loci_annotations,(!is.na(accession_number))&accession_number==current_accession_nb)
  if (nrow(current_df_loci_annotation)==1){
    v_genomic_region_length[current_accession_nb] <- current_df_loci_annotation$end - current_df_loci_annotation$start + 1
  }
}

#Find genes' length
v_lst_genes <- sort(unique(subset(df_loci_annotations,!is.na(the_product))$the_product))
v_genes_length <- rep(NA,length(v_lst_genes))
names(v_genes_length) <- v_lst_genes
for (current_gene in names(v_genes_length)){
  current_df_loci_annotation <- subset(df_loci_annotations,(!is.na(the_product))&the_product==current_gene)
  if (nrow(current_df_loci_annotation)==1){
    v_genes_length[current_gene] <- current_df_loci_annotation$end - current_df_loci_annotation$start + 1
  }
}

#Analysis of parallelism (fixed missense SNVs or indels)
df_ALL_variants$str_time <- as.character(v_sample_to_timepoints[df_ALL_variants$Sample])
df_ALL_variants$int_time <- as.integer(df_ALL_variants$str_time)
df_ALL_genes_multiplicity <- subset(df_ALL_variants,(!is.na(gene_product))&(variant_type2%in%c("missense","nonsense")|grepl(pattern = "indel",x = tolower(variant_type2),fixed = TRUE)|grepl(pattern = "intron",x = tolower(variant_type2),fixed = TRUE)|grepl(pattern = "in utr",x = tolower(variant_type2),fixed = TRUE))&is_fixed) %>%
  group_by(mutation,gene_product,accession_number,variant_type2) %>%
  dplyr::summarise(PH=length(unique(population)),multiplicity=length(unique(Position))/v_genes_length[gene_product])
df_ALL_genes_multiplicity <- subset(df_ALL_genes_multiplicity,!is.na(gene_product))
#find genes accession number
for (i in 1:nrow(df_ALL_genes_multiplicity)){
  the_ann <- subset(df_loci_annotations,the_product==df_ALL_genes_multiplicity$gene_product[i])$annotation[1]
  df_ALL_genes_multiplicity$accession_number[i] <- ifelse(test = grepl(pattern = ";Name=",x = the_ann,fixed = T),yes = substr(the_ann,start = gregexpr(pattern = "ID=",text = the_ann,fixed = T)[[1]]+3,stop = gregexpr(pattern = ";Name=",text = the_ann,fixed = T)[[1]]-1))
}
#Gene non-recurrent hit density
df_nonRecuHit_genes_density <- subset(df_ALL_variants,is_fixed&(!is.na(gene_product))&(variant_type2%in%c("missense","nonsense")|grepl(pattern = "indel",x = tolower(variant_type2),fixed = TRUE)|grepl(pattern = "intron",x = tolower(variant_type2),fixed = TRUE)|grepl(pattern = "in utr",x = tolower(variant_type2),fixed = TRUE)))%>%
  group_by(population,gene_product,.add = T)%>%
  dplyr::summarise(density = length(unique(mutation))/v_genes_length[gene_product],nb_non_rec_hit=length(unique(mutation)))
v_lst_nonRecuHit_genes <- sort(unique(df_nonRecuHit_genes_density$gene_product))
df_nonRecuHit_genes_density$minus_log10_density <- -log10(df_nonRecuHit_genes_density$density)
df_nonRecuHit_genes_density <- unique(df_nonRecuHit_genes_density)
mtx_density_nonRecuHitgenes_per_pop <- (as.matrix((reshape2::acast(df_nonRecuHit_genes_density, gene_product~population, value.var="minus_log10_density"))))
mtx_density_nonRecuHitgenes_per_pop[is.na(mtx_density_nonRecuHitgenes_per_pop)] <- 1e-6
v_lst_nonRecuHit_genes_hit_in_multiple_pops <- names(table(df_nonRecuHit_genes_density$gene_product)[table(df_nonRecuHit_genes_density$gene_product)>1])
write.table(x=mtx_density_nonRecuHitgenes_per_pop,file = paste0(output_workspace,"mtx_density_nonRecuHitgenes_per_pop.csv"),sep = ",",na = "NA",row.names = TRUE,col.names = TRUE)
write.table(x=v_lst_nonRecuHit_genes,file = paste0(output_workspace,"lst_nonRecuHit_genes.txt"),sep = ",",na = "NA",row.names = F,col.names = F)

#Multihit genes (Non-recurrent hits in protein-coding regions in different populations)
v_nb_pops_in_which_gene_is_hit <- rowSums(mtx_density_nonRecuHitgenes_per_pop>1e-6)
mtx_nb_nonRecuHits_per_pop_in_multihit_genes_prot_coding_regions <- (as.matrix((reshape2::acast(df_nonRecuHit_genes_density, gene_product~population, value.var="nb_non_rec_hit"))))
mtx_nb_nonRecuHits_per_pop_in_multihit_genes_prot_coding_regions[is.na(mtx_nb_nonRecuHits_per_pop_in_multihit_genes_prot_coding_regions)] <- 0
write.table(x=mtx_nb_nonRecuHits_per_pop_in_multihit_genes_prot_coding_regions,file = paste0(output_workspace,"mtx_nb_nonRecuHits_per_pop_in_multihit_genes.csv"),sep = ",",na = "NA",row.names = TRUE,col.names = TRUE)
v_lst_multihit_genes_in_protein_coding_regions <- names(v_nb_pops_in_which_gene_is_hit[v_nb_pops_in_which_gene_is_hit>1])
write.table(x=v_lst_multihit_genes_in_protein_coding_regions,file = paste0(output_workspace,"lst_multihit_genes.txt"),sep = ",",na = "NA",row.names = F,col.names = F)
#heatmap density
svg(paste0(output_workspace,"Heatmap_per_pop_hit_density_for_Genes_with_muliple_hits.svg"),width=25/2.54,height=20/2.54)
the_hmp <- heatmap.2(x = log10(mtx_density_nonRecuHitgenes_per_pop[rownames(mtx_density_nonRecuHitgenes_per_pop)%in%v_lst_multihit_genes_in_protein_coding_regions,]),  key.xlab = "log10(#hits/gene length)", distfun = function(x) vegan::vegdist(x,method = "eucl"),hclustfun = function(x) hclust(x,method="ward.D2"),cex=6, lwid=c(1.5,4),density.info = "none",na.color = "black",dendrogram = "both",margins = c(6, 5),cexCol = 2,cexRow = 1,xlab = "Population", ylab = "Genes",trace="none",scale="none", Rowv = T, Colv = T,col = colorRampPalette(brewer.pal(name = "Blues",n=9), space = "rgb")(36),symm=F,symkey=F,symbreaks=F)
heatmap.2(x = log10(mtx_density_nonRecuHitgenes_per_pop[rownames(mtx_density_nonRecuHitgenes_per_pop)%in%v_lst_multihit_genes_in_protein_coding_regions,]),  key.xlab = "log10(#hits/gene length)", distfun = function(x) vegan::vegdist(x,method = "eucl"),hclustfun = function(x) hclust(x,method="ward.D2"),cex=8, lwid=c(1.5,4),density.info = "none",na.color = "black",dendrogram = "both",margins = c(6, 5),cexCol = 1.7,cexRow = 1.7,xlab = "Population", ylab = "Genes",trace="none",scale="none", Rowv = T, Colv = T,col = colorRampPalette(brewer.pal(name = "Blues",n=9), space = "rgb")(36),symm=F,symkey=F,symbreaks=F,labRow = rownames(mtx_density_nonRecuHitgenes_per_pop[rownames(mtx_density_nonRecuHitgenes_per_pop)%in%v_lst_multihit_genes_in_protein_coding_regions,]))
dev.off()
remove(the_hmp)
#heatmap presence_absence
svg(paste0(output_workspace,"Heatmap_presence_absence_multihit_Genes.svg"),width=25/2.54,height=20/2.54)
the_hmp <- heatmap.2(x = mtx_density_nonRecuHitgenes_per_pop[rownames(mtx_density_nonRecuHitgenes_per_pop)%in%v_lst_multihit_genes_in_protein_coding_regions,],  key=F, distfun = function(x) vegan::vegdist(x,method = "jaccard"),hclustfun = function(x) hclust(x,method="ward.D2"),cex=6, lwid=c(1.5,4),density.info = "none",na.color = "black",dendrogram = "both",margins = c(6, 5),cexCol = 2,cexRow = 1,xlab = "Population", ylab = "Genes",trace="none",scale="none", Rowv = T, Colv = T,col = colorRampPalette(brewer.pal(name = "Blues",n=9), space = "rgb")(36),symm=F,symkey=F,symbreaks=F)
heatmap.2(x = log10(mtx_density_nonRecuHitgenes_per_pop[rownames(mtx_density_nonRecuHitgenes_per_pop)%in%v_lst_multihit_genes_in_protein_coding_regions,]),  key=F, distfun = function(x) vegan::vegdist(x,method = "jaccard"),hclustfun = function(x) hclust(x,method="ward.D2"), lwid=c(1.5,4),density.info = "none",na.color = "black",dendrogram = "none",margins = c(6, 5),cexCol = 1.7,cexRow = 1.7,xlab = "Population", ylab = "Genes",trace="none",scale="none", Rowv = T, Colv = T,col = colorRampPalette(brewer.pal(name = "Blues",n=9), space = "rgb")(36),symm=F,symkey=F,symbreaks=F,labRow = rownames(mtx_density_nonRecuHitgenes_per_pop[rownames(mtx_density_nonRecuHitgenes_per_pop)%in%v_lst_multihit_genes_in_protein_coding_regions,]))
dev.off()
remove(the_hmp)
#Multihit genes (Non-recurrent hits in noncoding regions in different populations)
df_multiHit_noncoding_regions_nb_hits <- subset(df_ALL_variants,is_fixed&(!is.na(accession_number))&!(variant_type2%in%c("missense","nonsense")|grepl(pattern = "indel",x = tolower(variant_type2),fixed = TRUE)|grepl(pattern = "intron",x = tolower(variant_type2),fixed = TRUE)|grepl(pattern = "in utr",x = tolower(variant_type2),fixed = TRUE)))%>%
  group_by(population,accession_number,.add = T)%>%
  dplyr::summarise(nb_non_rec_hit=length(unique(mutation)))

mtx_nb_nonRecuHits_per_pop_in_multihit_genes_noncoding_regions <- (as.matrix((reshape2::acast(df_multiHit_noncoding_regions_nb_hits, accession_number~population, value.var="nb_non_rec_hit"))))
mtx_nb_nonRecuHits_per_pop_in_multihit_genes_noncoding_regions[is.na(mtx_nb_nonRecuHits_per_pop_in_multihit_genes_noncoding_regions)] <- 0
write.table(x=mtx_nb_nonRecuHits_per_pop_in_multihit_genes_noncoding_regions,file = paste0(output_workspace,"mtx_nb_nonRecuHits_per_pop_in_multihit_genes_noncoding_regions.csv"),sep = ",",na = "NA",row.names = TRUE,col.names = TRUE)
v_nb_pops_in_which_gene_NC_regions_are_hit <- rowSums(mtx_nb_nonRecuHits_per_pop_in_multihit_genes_noncoding_regions>0)
v_lst_multihit_genes_in_noncoding_regions <- names(v_nb_pops_in_which_gene_NC_regions_are_hit[v_nb_pops_in_which_gene_NC_regions_are_hit>1])
write.table(x=v_lst_multihit_genes_in_noncoding_regions,file = paste0(output_workspace,"lst_multihit_genes_in_noncoding_regions.txt"),sep = ",",na = "NA",row.names = F,col.names = F)
write.table(x=rownames(mtx_nb_nonRecuHits_per_pop_in_multihit_genes_noncoding_regions),file = paste0(output_workspace,"lst_all_genes_with_hits_in_noncoding_regions.txt"),sep = ",",na = "NA",row.names = F,col.names = F)


# #Get list of multi-hit genes (with fixed mutations in multiple pop; mutation, gene, accession_number,length(unique(pop)) and subset for count >=2)
# df_multihit_genes <- subset(df_ALL_genes_multiplicity,PH>=2)
# v_lst_multihit_genes <- sort(unique(df_multihit_genes$gene_product))
# v_lst_multihit_genes_accession <- sort(unique(df_multihit_genes$accession_number))
# #***Run enrichment analysis with list of accession numbers on AnGeLi
# #Between-population frequency of multi-hit genes hits (fixed mutations)
# df_ALL_genes_multiplicity$prevalence <- df_ALL_genes_multiplicity$PH/length(v_pop_order)
# #list of hits
# v_lst_hits <- sort(unique(subset(df_ALL_genes_multiplicity,PH>1)$mutation))
# #Heatmap of multi-hit genes number of hits (fixed mutations in multiple population)
# df_multihit_genes_density <- subset(df_ALL_variants,is_fixed&(mutation%in%v_lst_hits)&(!is.na(gene_product))&(variant_type2%in%c("missense","nonsense")|grepl(pattern = "indel",x = tolower(variant_type2),fixed = TRUE)|grepl(pattern = "intron",x = tolower(variant_type2),fixed = TRUE)|grepl(pattern = "in utr",x = tolower(variant_type2),fixed = TRUE))&(gene_product%in%v_lst_multihit_genes))%>%
#   group_by(population,gene_product,.add = T)%>%
#   dplyr::summarise(density = length(unique(mutation))/v_genes_length[gene_product])
# df_multihit_genes_density <- unique(df_multihit_genes_density)
# df_multihit_genes_density$minus_log10_density <- -log10(df_multihit_genes_density$density)
# mtx_density_multihitgenes_per_pop <- (as.matrix((reshape2::acast(df_multihit_genes_density, gene_product~population, value.var="density"))))
# mtx_density_multihitgenes_per_pop[is.na(mtx_density_multihitgenes_per_pop)] <- 1e-6
# write.table(x=mtx_density_multihitgenes_per_pop,file = paste0(output_workspace,"mtx_density_multihitgenes_per_pop.csv"),sep = ",",na = "NA",row.names = TRUE,col.names = TRUE)
# write.table(x=v_lst_multihit_genes,file = paste0(output_workspace,"lst_multihit_genes.txt"),sep = ",",na = "NA",row.names = F,col.names = F)
# # #extract genes most significantly enriched biological process
# # df_GO_enriched_bio_process_multihit_genes <- read.csv2(file = paste0(output_workspace,"GO_enriched_biological_process_for_multihit_genes.csv"),sep = "\t",header = T,stringsAsFactors = FALSE)
# # df_GO_enriched_bio_process_multihit_genes$Corrected_pvalue <- as.numeric(df_GO_enriched_bio_process_multihit_genes$Corrected_pvalue)
# # v_genes_to_accession <- unique(df_ALL_genes_multiplicity$accession_number)
# # names(v_genes_to_accession) <- unique(df_ALL_genes_multiplicity$gene_product)
# # v_genes_most_enriched_bio_process <- rep(NA,nrow(mtx_density_multihitgenes_per_pop))
# # names(v_genes_most_enriched_bio_process) <- rownames(mtx_density_multihitgenes_per_pop)
# # for (current_gene in rownames(mtx_density_multihitgenes_per_pop)){
# #   current_df_GO <- subset(df_GO_enriched_bio_process_multihit_genes,grepl(pattern = v_genes_to_accession[current_gene],x = Genes_In_Common.Interacting_Gene_Pair,fixed = T))
# #   if (nrow(current_df_GO)>0){
# #     v_genes_most_enriched_bio_process[current_gene] <- current_df_GO$GeneSet_Name[which(current_df_GO$Corrected_pvalue==min(current_df_GO$Corrected_pvalue,na.rm = T))[1]]
# #   }
# # }
# # palette_GO <- RColorBrewer::brewer.pal(8,"Dark2")
# # names(palette_GO) <- sort(unique(v_genes_most_enriched_bio_process))
# # v_genes_most_enriched_bio_process[is.na(v_genes_most_enriched_bio_process)] <- "Others"
# # palette_GO["Others"] <- "white"
# #Heatmap_breadth_of_cov_per_sample_for_each_amplicon
# svg(paste0(output_workspace,"Heatmap_multihit_genes_mutation_density_per_pop.svg"),width=25/2.54,height=20/2.54)
# #heatmap.2(x = log10(mtx_density_multihitgenes_per_pop), RowSideColors = unname(palette_GO[v_genes_most_enriched_bio_process[rownames(mtx_density_multihitgenes_per_pop)]]), key.xlab = "log10(#mutations/gene length)", distfun = function(x) vegan::vegdist(x,method = "eucl"),hclustfun = function(x) hclust(x,method="ward.D2"),cex=6, lwid=c(1.5,4),density.info = "none",na.color = "black",dendrogram = "both",margins = c(6, 5),cexCol = 2,cexRow = 0.1,xlab = "Population", ylab = "Genes",trace="none",scale="none", Rowv = T, Colv = T,col = colorRampPalette(brewer.pal(name = "Blues",n=9), space = "rgb")(36),symm=F,symkey=F,symbreaks=F)
# the_hmp <- heatmap.2(x = log10(mtx_density_multihitgenes_per_pop),  key.xlab = "log10(#hits/gene length)", distfun = function(x) vegan::vegdist(x,method = "eucl"),hclustfun = function(x) hclust(x,method="ward.D2"),cex=6, lwid=c(1.5,4),density.info = "none",na.color = "black",dendrogram = "both",margins = c(6, 5),cexCol = 2,cexRow = 0.1,xlab = "Population", ylab = "Genes",trace="none",scale="none", Rowv = T, Colv = T,col = colorRampPalette(brewer.pal(name = "Blues",n=9), space = "rgb")(36),symm=F,symkey=F,symbreaks=F)
# heatmap.2(x = log10(mtx_density_multihitgenes_per_pop),  key.xlab = "log10(#hits/gene length)", distfun = function(x) vegan::vegdist(x,method = "eucl"),hclustfun = function(x) hclust(x,method="ward.D2"),cex=6, lwid=c(1.5,4),density.info = "none",na.color = "black",dendrogram = "both",margins = c(6, 5),cexCol = 1,cexRow = 0.5,xlab = "Population", ylab = "Genes",trace="none",scale="none", Rowv = T, Colv = T,col = colorRampPalette(brewer.pal(name = "Blues",n=9), space = "rgb")(36),symm=F,symkey=F,symbreaks=F,labRow = rownames(mtx_density_multihitgenes_per_pop))
# dev.off()
# remove(the_hmp)
##pie(sort(table(unname(palette_GO[v_genes_most_enriched_bio_process[rownames(mtx_density_multihitgenes_per_pop)]])),dec=T), col=names(sort(table(unname(palette_GO[v_genes_most_enriched_bio_process[rownames(mtx_density_multihitgenes_per_pop)]])),dec=T)),labels = names(palette_GO)[unname(vapply(X = names(sort(table(unname(palette_GO[v_genes_most_enriched_bio_process[rownames(mtx_density_multihitgenes_per_pop)]])),dec=T)), FUN= function(x) which(palette_GO==x), FUN.VALUE = 0))])


#Time series for the allele frequency of the substitutions
ggplot(data = subset(df_variants,pop_mut%in%df_substitutions$pop_mut),mapping = aes(x=int_time,y=VarFreq,group=mutation)) + 
  geom_line(aes(col=mutation))+
  theme(plot.title=element_text(hjust=0,size=12))+
  xlab("Time (Generation)") + ylab("Allele Frequency") +
  theme_bw() + theme(strip.text = element_text(size=16),axis.text.x=element_text(size=14,angle=60,hjust=1),axis.title.x = element_text(size=16),axis.text.y=element_text(size=14),axis.title.y = element_text(size=16),legend.position = "none") +  scale_x_continuous(limits = c(0,10000),breaks = seq(0,10000,1000)) + facet_wrap(~population,ncol = 2) #+ guides(col=guide_legend(title="Substitution"))

ggsave(filename = "Substitutions_Allele_Freq_time_series_across_populations.eps", path=output_workspace, width = 20/2.54, height = 30/2.54, device = cairo_ps)
ggsave(filename = "Substitutions_Allele_Freq_time_series_across_populations.svg", path=output_workspace, width = 20/2.54, height = 30/2.54, device = svg)

#Time series of the mutations allele frequency in protein-coding regions
df_variants_time_series <- df_ALL_variants #subset(df_ALL_variants, (!is.na(gene_product))&(!grepl(pattern = "noncoding",x = tolower(variant_type2),fixed = T)))
df_variants_time_series <- subset(df_variants_time_series,!is.na(Chrom))
df_variants_time_series$mutation <- paste0(df_variants_time_series$Chrom,"_",df_variants_time_series$Ref,df_variants_time_series$Position,df_variants_time_series$VarAllele)
df_variants_time_series$population <- v_sample_to_unique_pop[df_variants_time_series$Sample]
df_variants_time_series$pop_mut <- paste0(df_variants_time_series$population,"_",df_variants_time_series$mutation)
df_variants_time_series$str_time <- as.character(v_sample_to_timepoints[df_variants_time_series$Sample])
df_variants_time_series$int_time <- as.integer(df_variants_time_series$str_time)
df_variants_time_series$is_fixed <- df_variants_time_series$VarFreq >= 0.75
df_variants_time_series$is_multihit_gene <- (df_variants_time_series$gene_product%in%v_lst_nonRecuHit_genes_hit_in_multiple_pops)&(!is.na(df_variants_time_series$gene_product))
df_variants_time_series$Context_5bp_before_and_after <- NA
for (the_row_indx in 1:nrow(df_variants_time_series)){
  df_variants_time_series$Context_5bp_before_and_after[the_row_indx] <- paste0(substr(genome_refseq[[df_variants_time_series$Chrom[the_row_indx]]],start = df_variants_time_series$Position[the_row_indx]-5,stop = df_variants_time_series$Position[the_row_indx]-1),substr(genome_refseq[[df_variants_time_series$Chrom[the_row_indx]]],start = df_variants_time_series$Position[the_row_indx]+1,stop = df_variants_time_series$Position[the_row_indx]+5))
}
df_variants_time_series$is_in_homopolymeric_region <- vapply(X = df_variants_time_series$Context_5bp_before_and_after,FUN = function(the_label_context) grepl(pattern = "AAAA",x = toupper(the_label_context),fixed = T)|grepl(pattern = "TTTT",x = toupper(the_label_context),fixed = T)|grepl(pattern = "GGGG",x = toupper(the_label_context),fixed = T)|grepl(pattern = "CCCC",x = toupper(the_label_context),fixed = T),FUN.VALUE = F)
print(epitools::oddsratio(table(unique(df_variants_time_series[,c("pop_mut","variant_type","is_in_homopolymeric_region")])$variant_type,
                                unique(df_variants_time_series[,c("pop_mut","variant_type","is_in_homopolymeric_region")])$is_in_homopolymeric_region)))

#list of potential multihit genes with evidence for transcription regulation
v_lst_multihit_genes_that_are_evidence_for_transcription_variation <- c("seb1","deb1","atd1","ppk30","ctf18","rsd1","iss9")
print(v_lst_multihit_genes_that_are_evidence_for_transcription_variation)

df_variants_time_series$VAF_time_series_line_color <- ifelse(test = (df_variants_time_series$gene_product%in%v_lst_multihit_genes_that_are_evidence_for_transcription_variation) & (df_variants_time_series$variant_type2!="synonymous"),yes = "orange",no = ifelse(test = df_variants_time_series$variant_type2!="synonymous",yes = ifelse(test = df_variants_time_series$is_multihit_gene,yes = "green4",no = "black"),no = "grey"))
df_variants_time_series$VAF_time_series_line_size <- ifelse(test = (df_variants_time_series$gene_product%in%v_lst_multihit_genes_that_are_evidence_for_transcription_variation) & (df_variants_time_series$variant_type2!="synonymous"),yes = 0.25,no = 0.5)
df_variants_time_series$VAF_time_series_line_type <- ifelse(test = df_variants_time_series$variant_type2=="synonymous",yes = 3,no = 1) # ifelse(test = df_variants_time_series$is_in_homopolymeric_region,yes = 5,no = 1)
df_variants_time_series$VAF_time_series_line_alpha <- ifelse(test = df_variants_time_series$VAF_time_series_line_color=="grey",yes = 0.4,no = 1)

#get the list of sample_generation associations
lst_pop_gentime <- sort(unique(paste0(df_ALL_variants$population,"_",df_ALL_variants$str_time) ))
#put VAF to 0 when it is absent during a time point AND never existed before
for (current_variant in v_lst_variants){
  for (current_str_time in names(v_next_timepoint)){
    the_subset <- subset(df_variants_time_series, (pop_mut==current_variant)&(str_time==current_str_time))
    the_subset_of_prev_gen <- subset(df_variants_time_series, (VarFreq>0)&(pop_mut==current_variant)&(int_time<as.integer(current_str_time)))
    if ( (current_str_time!="10000") & (nrow(the_subset)+nrow(the_subset_of_prev_gen)==0) ){
      the_new_row <- subset(df_variants_time_series, (pop_mut==current_variant))[1,]
      the_new_row$str_time <- current_str_time
      the_new_row$int_time <- as.integer(the_new_row$str_time)
      the_new_row$VarFreq <- 0
      the_new_row$Sample <- ifelse(test = !is.na(the_new_row$population),yes = lst_samples[grepl(pattern = paste0(the_new_row$population,"_",the_new_row$str_time,"_"),x = lst_samples,fixed = T)],no = NA) 
      df_variants_time_series <- rbind(df_variants_time_series,the_new_row)
      #print("Added 0 VAF")
    }
  }
}
#Start the time series at frequency 0
for (current_variant in v_lst_variants){
  current_str_time <- "0"
  the_new_row <- subset(df_variants_time_series, (pop_mut==current_variant))[1,]
  the_new_row$str_time <- current_str_time
  the_new_row$int_time <- as.integer(the_new_row$str_time)
  the_new_row$VarFreq <- 0
  the_new_row$Sample <- ifelse(test = !is.na(the_new_row$population),yes = lst_samples[grepl(pattern = paste0(the_new_row$population,"_",the_new_row$str_time,"_"),x = lst_samples,fixed = T)],no = NA) 
  df_variants_time_series <- rbind(df_variants_time_series,the_new_row)
}
#remove entries without a population (coding artifact)
df_variants_time_series <- subset(df_variants_time_series,!is.na(population))
#make sure that the timepoints that were not there initially have 0 entries
df_variants_time_series$pop_gentime <- paste0(df_variants_time_series$population,"_",df_variants_time_series$str_time)
df_variants_time_series <- subset(df_variants_time_series,pop_gentime%in%lst_pop_gentime)
#remove variants that are not present in more than 4 timepoints and never reached fixation and are absent from generation 10000
nb_timepoints_per_variant <- sort(vapply(unique(df_variants_time_series$pop_mut), function(x) length(unique((subset(df_variants_time_series,pop_mut==x & VarFreq>0)$str_time))), FUN.VALUE = 0),dec=T)
max_freq_of_variants <- sort(vapply(unique(df_variants_time_series$pop_mut), function(x) max(subset(df_variants_time_series,pop_mut==x & VarFreq>0)$VarFreq,na.rm=T), FUN.VALUE = 0.0),dec=T)

lst_pop_muts_present_in_less_than_5_timepoints_never_fixed_and_absent_at_the_end <- intersect(setdiff(names(nb_timepoints_per_variant[nb_timepoints_per_variant<=4]), subset(df_variants_time_series,str_time=="10000")$pop_mut ),names(max_freq_of_variants[max_freq_of_variants<0.75]))
df_variants_time_series <- subset(df_variants_time_series,!pop_mut%in%lst_pop_muts_present_in_less_than_5_timepoints_never_fixed_and_absent_at_the_end)
df_variants_time_series <- subset(df_variants_time_series,Sample%in%lst_samples)
df_variants_time_series$opacity_based_on_is_multihit_genes <- ifelse(df_variants_time_series$is_multihit_gene&!is.na(df_variants_time_series$gene_product),yes = 1,no = 0.8)
print(nrow(df_variants_time_series))
#Number of populations in which multihit genes are hit
v_nb_pops_per_multihit_gene <- sort(rowSums(1*as.matrix(mtx_nb_nonRecuHits_per_pop_in_multihit_genes_prot_coding_regions>0))[rowSums(1*as.matrix(mtx_nb_nonRecuHits_per_pop_in_multihit_genes_prot_coding_regions>0))>1&(names(rowSums(1*as.matrix(mtx_nb_nonRecuHits_per_pop_in_multihit_genes_prot_coding_regions>0)))%in%v_lst_multihit_genes_in_protein_coding_regions)])
#Number of multihit genes in each population
v_nb_multihit_genes_per_pop <- sort(colSums(1*as.matrix(mtx_nb_nonRecuHits_per_pop_in_multihit_genes_prot_coding_regions>0)[rowSums(1*as.matrix(mtx_nb_nonRecuHits_per_pop_in_multihit_genes_prot_coding_regions>0))>1&(names(rowSums(1*as.matrix(mtx_nb_nonRecuHits_per_pop_in_multihit_genes_prot_coding_regions>0)))%in%v_lst_multihit_genes_in_protein_coding_regions),]))
#Number of hits in multihit genes per population
v_nb_hits_in_multihit_genes_per_pop <- vapply(X = names(v_nb_multihit_genes_per_pop),FUN = function(the_pop) length(unique(subset(df_variants_time_series,is_multihit_gene&(!is.na(gene_product))&is_fixed&(variant_type2!="synonymous")&(!is.na(variant_type2))&(population==the_pop))$pop_mut)),FUN.VALUE = 0)
#import ic50 data
v_h2o2_ic50_evolved_pops <- c(1.3368808, 1.2767600, 0.9929333, 1.2583261, 1.9882318, 1.8095702, 1.3376845, 1.8315572,0.5468524, 0.5200863)
names(v_h2o2_ic50_evolved_pops) <- c(paste0("G",c(1:5,9,11)),"H9","H11","H12")
#Evolved populations H2O2 sensistivity (1/IC50) vs Number of multihit genes per population
sqrt(summary(lm(unname(1/v_h2o2_ic50_evolved_pops)~unname(v_nb_multihit_genes_per_pop[names(v_h2o2_ic50_evolved_pops)])))$r.squared) 
#Evolved populations H2O2 sensitivity (1/IC50) vs Number of HITS in multihit genes per population
#sqrt(summary(lm(unname(1/v_h2o2_ic50_evolved_pops)~unname(v_nb_hits_in_multihit_genes_per_pop[names(v_h2o2_ic50_evolved_pops)])))$r.squared) 

#make sure to only retain samples that really exist
#df_variants_time_series <- subset(df_variants_time_series,Sample%in%lst_samples)

#make sure to subset the line categories and draw them in the correct order
df_synonymous <- subset(df_variants_time_series,(VAF_time_series_line_color=="grey")&(VAF_time_series_line_type==3))
df_missense_not_in_multigenes <- subset(df_variants_time_series,(VAF_time_series_line_color=="black")&(VAF_time_series_line_type==1))
df_hits_in_multigenes <- subset(df_variants_time_series,VAF_time_series_line_color=="green4")
df_hits_evidence_of_transcription_modification <- subset(df_variants_time_series,VAF_time_series_line_color=="orange")

for (current_pop in v_pop_order){
  ggplot() + 
    geom_line(mapping = aes(x=subset(df_synonymous,population==current_pop)$int_time,y=subset(df_synonymous,population==current_pop)$VarFreq,group=subset(df_synonymous,population==current_pop)$pop_mut), size= subset(df_synonymous,population==current_pop)$VAF_time_series_line_size,col= subset(df_synonymous,population==current_pop)$VAF_time_series_line_color,alpha= subset(df_synonymous,population==current_pop)$VAF_time_series_line_alpha,lty=subset(df_synonymous,population==current_pop)$VAF_time_series_line_type)+
    geom_line(mapping = aes(x=subset(df_missense_not_in_multigenes,population==current_pop)$int_time,y=subset(df_missense_not_in_multigenes,population==current_pop)$VarFreq,group=subset(df_missense_not_in_multigenes,population==current_pop)$pop_mut), size= subset(df_missense_not_in_multigenes,population==current_pop)$VAF_time_series_line_size,col= subset(df_missense_not_in_multigenes,population==current_pop)$VAF_time_series_line_color,alpha= subset(df_missense_not_in_multigenes,population==current_pop)$VAF_time_series_line_alpha,lty=subset(df_missense_not_in_multigenes,population==current_pop)$VAF_time_series_line_type)+
    geom_line(mapping = aes(x=subset(df_hits_in_multigenes,population==current_pop)$int_time,y=subset(df_hits_in_multigenes,population==current_pop)$VarFreq,group=subset(df_hits_in_multigenes,population==current_pop)$pop_mut), size= subset(df_hits_in_multigenes,population==current_pop)$VAF_time_series_line_size,col= subset(df_hits_in_multigenes,population==current_pop)$VAF_time_series_line_color,alpha= subset(df_hits_in_multigenes,population==current_pop)$VAF_time_series_line_alpha,lty=subset(df_hits_in_multigenes,population==current_pop)$VAF_time_series_line_type)+
    geom_line(mapping = aes(x=subset(df_hits_evidence_of_transcription_modification,population==current_pop)$int_time,y=subset(df_hits_evidence_of_transcription_modification,population==current_pop)$VarFreq,group=subset(df_hits_evidence_of_transcription_modification,population==current_pop)$pop_mut), size= subset(df_hits_evidence_of_transcription_modification,population==current_pop)$VAF_time_series_line_size,col= subset(df_hits_evidence_of_transcription_modification,population==current_pop)$VAF_time_series_line_color,alpha= subset(df_hits_evidence_of_transcription_modification,population==current_pop)$VAF_time_series_line_alpha,lty=subset(df_hits_evidence_of_transcription_modification,population==current_pop)$VAF_time_series_line_type)+
    theme(plot.title=element_text(hjust=0,size=12))+
    xlab("Time (Generation)") + ylab("Allele Frequency") +
    theme_bw() + theme(strip.text = element_text(size=16),axis.text.x=element_text(size=14,angle=60,hjust=1),axis.title.x = element_text(size=16),axis.text.y=element_text(size=14),axis.title.y = element_text(size=16),legend.position = "none") +  scale_x_continuous(limits = c(0,10000),breaks = seq(0,10000,1000)) #+ facet_wrap(~factor(df_variants_time_series$population,levels = v_pop_order),ncol = 2) #+ guides(col=guide_legend(title="Substitution"))
  
  ggsave(filename = paste0(current_pop,"_Mutations_Allele_Freq_time_series_across_populations.svg"), path=output_workspace, width = 15/2.54, height = 10/2.54, device = svg)
  
  ggplot() + 
    geom_line(mapping = aes(x=subset(df_synonymous,population==current_pop&(!is.na(gene_product))&(!grepl(pattern = "noncoding",x = tolower(variant_type2),fixed = T)))$int_time,y=subset(df_synonymous,population==current_pop&(!is.na(gene_product))&(!grepl(pattern = "noncoding",x = tolower(variant_type2),fixed = T)))$VarFreq,group=subset(df_synonymous,population==current_pop&(!is.na(gene_product))&(!grepl(pattern = "noncoding",x = tolower(variant_type2),fixed = T)))$pop_mut), size= subset(df_synonymous,population==current_pop&(!is.na(gene_product))&(!grepl(pattern = "noncoding",x = tolower(variant_type2),fixed = T)))$VAF_time_series_line_size,col= subset(df_synonymous,population==current_pop&(!is.na(gene_product))&(!grepl(pattern = "noncoding",x = tolower(variant_type2),fixed = T)))$VAF_time_series_line_color,alpha= subset(df_synonymous,population==current_pop&(!is.na(gene_product))&(!grepl(pattern = "noncoding",x = tolower(variant_type2),fixed = T)))$VAF_time_series_line_alpha,lty=subset(df_synonymous,population==current_pop&(!is.na(gene_product))&(!grepl(pattern = "noncoding",x = tolower(variant_type2),fixed = T)))$VAF_time_series_line_type)+
    geom_line(mapping = aes(x=subset(df_missense_not_in_multigenes,population==current_pop&(!is.na(gene_product))&(!grepl(pattern = "noncoding",x = tolower(variant_type2),fixed = T)))$int_time,y=subset(df_missense_not_in_multigenes,population==current_pop&(!is.na(gene_product))&(!grepl(pattern = "noncoding",x = tolower(variant_type2),fixed = T)))$VarFreq,group=subset(df_missense_not_in_multigenes,population==current_pop&(!is.na(gene_product))&(!grepl(pattern = "noncoding",x = tolower(variant_type2),fixed = T)))$pop_mut), size= subset(df_missense_not_in_multigenes,population==current_pop&(!is.na(gene_product))&(!grepl(pattern = "noncoding",x = tolower(variant_type2),fixed = T)))$VAF_time_series_line_size,col= subset(df_missense_not_in_multigenes,population==current_pop&(!is.na(gene_product))&(!grepl(pattern = "noncoding",x = tolower(variant_type2),fixed = T)))$VAF_time_series_line_color,alpha= subset(df_missense_not_in_multigenes,population==current_pop&(!is.na(gene_product))&(!grepl(pattern = "noncoding",x = tolower(variant_type2),fixed = T)))$VAF_time_series_line_alpha,lty=subset(df_missense_not_in_multigenes,population==current_pop&(!is.na(gene_product))&(!grepl(pattern = "noncoding",x = tolower(variant_type2),fixed = T)))$VAF_time_series_line_type)+
    geom_line(mapping = aes(x=subset(df_hits_in_multigenes,population==current_pop&(!is.na(gene_product))&(!grepl(pattern = "noncoding",x = tolower(variant_type2),fixed = T)))$int_time,y=subset(df_hits_in_multigenes,population==current_pop&(!is.na(gene_product))&(!grepl(pattern = "noncoding",x = tolower(variant_type2),fixed = T)))$VarFreq,group=subset(df_hits_in_multigenes,population==current_pop&(!is.na(gene_product))&(!grepl(pattern = "noncoding",x = tolower(variant_type2),fixed = T)))$pop_mut), size= subset(df_hits_in_multigenes,population==current_pop&(!is.na(gene_product))&(!grepl(pattern = "noncoding",x = tolower(variant_type2),fixed = T)))$VAF_time_series_line_size,col= subset(df_hits_in_multigenes,population==current_pop&(!is.na(gene_product))&(!grepl(pattern = "noncoding",x = tolower(variant_type2),fixed = T)))$VAF_time_series_line_color,alpha= subset(df_hits_in_multigenes,population==current_pop&(!is.na(gene_product))&(!grepl(pattern = "noncoding",x = tolower(variant_type2),fixed = T)))$VAF_time_series_line_alpha,lty=subset(df_hits_in_multigenes,population==current_pop&(!is.na(gene_product))&(!grepl(pattern = "noncoding",x = tolower(variant_type2),fixed = T)))$VAF_time_series_line_type)+
    geom_line(mapping = aes(x=subset(df_hits_evidence_of_transcription_modification,population==current_pop&(!is.na(gene_product))&(!grepl(pattern = "noncoding",x = tolower(variant_type2),fixed = T)))$int_time,y=subset(df_hits_evidence_of_transcription_modification,population==current_pop&(!is.na(gene_product))&(!grepl(pattern = "noncoding",x = tolower(variant_type2),fixed = T)))$VarFreq,group=subset(df_hits_evidence_of_transcription_modification,population==current_pop&(!is.na(gene_product))&(!grepl(pattern = "noncoding",x = tolower(variant_type2),fixed = T)))$pop_mut), size= subset(df_hits_evidence_of_transcription_modification,population==current_pop&(!is.na(gene_product))&(!grepl(pattern = "noncoding",x = tolower(variant_type2),fixed = T)))$VAF_time_series_line_size,col= subset(df_hits_evidence_of_transcription_modification,population==current_pop&(!is.na(gene_product))&(!grepl(pattern = "noncoding",x = tolower(variant_type2),fixed = T)))$VAF_time_series_line_color,alpha= subset(df_hits_evidence_of_transcription_modification,population==current_pop&(!is.na(gene_product))&(!grepl(pattern = "noncoding",x = tolower(variant_type2),fixed = T)))$VAF_time_series_line_alpha,lty=subset(df_hits_evidence_of_transcription_modification,population==current_pop&(!is.na(gene_product))&(!grepl(pattern = "noncoding",x = tolower(variant_type2),fixed = T)))$VAF_time_series_line_type)+
    theme(plot.title=element_text(hjust=0,size=12))+
    xlab("Time (Generation)") + ylab("Allele Frequency") +
    theme_bw() + theme(strip.text = element_text(size=16),axis.text.x=element_text(size=14,angle=60,hjust=1),axis.title.x = element_text(size=16),axis.text.y=element_text(size=14),axis.title.y = element_text(size=16),legend.position = "none") +  scale_x_continuous(limits = c(0,10000),breaks = seq(0,10000,1000)) #+ facet_wrap(~factor(df_variants_time_series$population,levels = v_pop_order),ncol = 2) #+ guides(col=guide_legend(title="Substitution"))
  
  ggsave(filename = paste0(current_pop,"_coding_regions_mutations_Allele_Freq_time_series_across_populations.svg"), path=output_workspace, width = 15/2.54, height = 10/2.54, device = svg)
  
}

#create G1-G4 allele time series by mutation types and distinguish multi-hit genes
palette_mh_genes <- c("#1F78C8","#ff0000","#33a02c","#6A33C2","#ff7f00","#565656",
           "#FFD700","#a6cee3","#FB6496","#b2df8a","#CAB2D6","#FDBF6F",
           "#999999","#EEE685","#C8308C","#FF83FA","#C814FA","#0000FF",
           "#36648B","#00E2E5","#00FF00","#778B00","#BEBE00","#8B3B00",
           "#A52A3C")
v_lst_hits_pop_mut <- subset(df_variants_time_series, is_fixed&(variant_type2!="synonymous")&(!is.na(variant_type2)))$pop_mut
v_multihit_genes_to_pch <- 0:(length(v_lst_nonRecuHit_genes_hit_in_multiple_pops)-1)
v_multihit_genes_to_pch[v_multihit_genes_to_pch==11] <- length(v_lst_nonRecuHit_genes_hit_in_multiple_pops)+1
names(v_multihit_genes_to_pch) <- v_lst_nonRecuHit_genes_hit_in_multiple_pops
df_variants_time_series$`Multi-hit gene` <- ifelse(test = df_variants_time_series$is_multihit_gene,yes = df_variants_time_series$gene_product,no = NA)
df_variants_time_series$`Multi-hit gene` <- as.factor(df_variants_time_series$`Multi-hit gene`)
df_variants_time_series$pch <- ifelse(test = df_variants_time_series$is_multihit_gene,yes = v_multihit_genes_to_pch[df_variants_time_series$gene_product],no = 20)
  #G1-G4
ggplot(data=subset(df_variants_time_series,population%in%v_pop_order[1:4]&(!is.na(gene_product))&(!is.na(pop_mut))&(pop_mut%in%v_lst_hits_pop_mut) )) + 
  geom_point(mapping = aes(x=int_time,y=VarFreq,shape = variant_type2,col=`Multi-hit gene`),size=ifelse(test = subset(df_variants_time_series,population%in%v_pop_order[1:4]&(!is.na(gene_product))&(!is.na(pop_mut))&(pop_mut%in%v_lst_hits_pop_mut))$is_multihit_gene,yes = 3,no = NA)) + 
  geom_line(mapping = aes(x=int_time,y=VarFreq,group=pop_mut,col=`Multi-hit gene`,alpha=opacity_based_on_is_multihit_genes))+
  theme(plot.title=element_text(hjust=0,size=12))+ guides(alpha=F)+scale_shape_manual(values = unname(v_multihit_genes_to_pch)) + 
  xlab("Time (Generation)") + ylab("Allele Frequency") + guides(shape=guide_legend(ncol=2))+
  theme_bw() + theme(strip.text = element_text(size=16),axis.text.x=element_text(size=14,angle=60,hjust=1),axis.title.x = element_text(size=16),axis.text.y=element_text(size=14),axis.title.y = element_text(size=16),legend.position = "right",legend.title = element_text(size=14),legend.text = element_text(size=10) ) + scale_color_manual(values = palette_mh_genes) + scale_x_continuous(limits = c(0,10000),breaks = seq(0,10000,1000)) + facet_wrap(~factor(population,levels = v_pop_order),ncol = 2) #+ guides(col=guide_legend(title="Substitution"))
ggsave(filename = paste0("G1-G4_Hits_Allele_Freq_by_mutation_type_time_series.svg"), path=output_workspace, width = 33/2.54, height = 20/2.54, device = svg)
  #ALL populations
ggplot(data=subset(df_variants_time_series,(!is.na(gene_product))&(!is.na(pop_mut))&(pop_mut%in%v_lst_hits_pop_mut) )) + 
  geom_point(mapping = aes(x=int_time,y=VarFreq,shape = variant_type2,col=`Multi-hit gene`),size=ifelse(test = subset(df_variants_time_series,(!is.na(gene_product))&(!is.na(pop_mut))&(pop_mut%in%v_lst_hits_pop_mut))$is_multihit_gene,yes = 3,no = NA)) + 
  geom_line(mapping = aes(x=int_time,y=VarFreq,group=pop_mut,col=`Multi-hit gene`,alpha=opacity_based_on_is_multihit_genes))+
  theme(plot.title=element_text(hjust=0,size=12))+ guides(alpha=F)+scale_shape_manual(values = unname(v_multihit_genes_to_pch)) + 
  xlab("Time (Generation)") + ylab("Allele Frequency") + guides(shape=guide_legend(ncol=1)) +
  theme_bw() + theme(strip.text = element_text(size=16),axis.text.x=element_text(size=14,angle=60,hjust=1),axis.title.x = element_text(size=16),axis.text.y=element_text(size=14),axis.title.y = element_text(size=16),legend.position = "right",legend.title = element_text(size=14),legend.text = element_text(size=10) ) + scale_color_manual(values = palette_mh_genes) + scale_x_continuous(limits = c(0,10000),breaks = seq(0,10000,1000)) + facet_wrap(~factor(population,levels = v_pop_order),ncol = 2) #+ guides(col=guide_legend(title="Substitution"))
ggsave(filename = paste0("ALL_POPS_Hits_Allele_Freq_by_mutation_type_time_series.svg"), path=output_workspace, width = 33/2.54, height = 35/2.54, device = svg)
ggsave(filename = paste0("ALL_POPS_Hits_Allele_Freq_by_mutation_type_time_series.png"), path=output_workspace, width = 33/2.54, height = 35/2.54)

# nadf_variants_time_series <- subset(df_variants_time_series, !is.na(population))
# nadf_variants_time_series$VarFreq[nadf_variants_time_series$VarFreq==0] <- NA
# write.table(x=nadf_variants_time_series,file = paste0(output_workspace,"NATable_variants_time_series.csv"),sep = ",",na = "NA",row.names = F,col.names = T)

#write.table(x=df_annotated_variants,file = paste0(output_workspace,"Table_annotated_variants_VAF.csv"),sep = ",",na = "NA",row.names = F,col.names = T)
#write.table(x=df_coding_regions_variants,file = paste0(output_workspace,"Table_VAF_in_coding_regions.csv"),sep = ",",na = "NA",row.names = F,col.names = T)
#write.table(x=df_variants_time_series,file = paste0(output_workspace,"Table_variants_time_series.csv"),sep = ",",na = "NA",row.names = F,col.names = T)

#Save the list of genes hit by SNVs or indels to perform a GO entichment analysis suggesting evidence for transcriptional innovation
write.table(x=unique(df_ALL_variants$gene_product),file = paste0(output_workspace,"lst_mutated_genes.txt"),sep = ",",na = "NA",row.names = F,col.names = F)

# #Find further evidence of transcriptomic innovation by looking at potential cis-eQTL
# df_ALL_variants$is_potential_cis_regulator <- vapply(X = 1:nrow(df_ALL_variants),FUN = find_if_mutation_is_cis,FUN.VALUE = F)
# df_ALL_variants$potentially_cis_regulated_gene <- vapply(X = 1:nrow(df_ALL_variants),FUN = find_potentially_cis_regulated_gene,FUN.VALUE = "")
# df_ALL_variants$potentially_cis_regulated_gene <- ifelse(df_ALL_variants$potentially_cis_regulated_gene=="NA",yes=NA,no=df_ALL_variants$potentially_cis_regulated_gene)
# df_ALL_variants$locus <- paste0(df_ALL_variants$Chrom,"_",df_ALL_variants$Position)
# print(paste0("There are ",length(unique(subset(df_ALL_variants,is_potential_cis_regulator)$gene_product))," genes that are potentially cis-regulated!"))
# print(paste0("There are ",length(unique(subset(df_ALL_variants,is_potential_cis_regulator)$locus))," loci that are potentially cis-eQTL!"))

#Plot contamination results
df_contamination <- read.csv2(file = paste0(output_workspace,"Table_contamination.csv"),sep = "\t",header = TRUE,stringsAsFactors = FALSE)
df_contamination$nb_reads <- as.numeric(df_contamination$nb_reads)
df_contamination$population <- v_sample_to_unique_pop[df_contamination$sample]
df_contamination$str_time <- as.character(v_sample_to_timepoints[df_contamination$sample])
df_contamination$int_time <- as.integer(df_contamination$str_time)
df_contamination$label_species <- c("C_glabrata"="*C. glabrata* CBS138","K_lactis_NRRL-Y-1140"="*K. lactis* NRRL-Y-1140","Pichia_pastoris"="*P. pastoris* EAMORG09","S_cerevisiae_YJM978"="*S. cerevisiae* YJM978","S_paradoxus_IFO1804"="*S. paradoxus* IFO1804","Schizosaccharomyces_pombe"="*S. pombe* 972h-")[df_contamination$species]
df_proportion_contamination <- df_contamination %>%
  group_by(sample) %>%
  dplyr::summarise(population=population,species=species,label_species=label_species,str_time=str_time,int_time=int_time,nb_reads=nb_reads,sum_nb_reads_in_sample=sum(nb_reads),proportion = nb_reads/sum(nb_reads)) %>%
  as.data.frame()

ggplot(data = df_proportion_contamination,mapping = aes(x=factor(str_time,as.character(sort(unique(df_variants$int_time)))),y=proportion)) + geom_col(mapping = aes(fill=factor(label_species,levels=c("*P. pastoris* EAMORG09","*K. lactis* NRRL-Y-1140","*C. glabrata* CBS138","*S. paradoxus* IFO1804","*S. cerevisiae* YJM978","*S. pombe* 972h-")))) + xlab("Time (Generation)") + ylab("Proportion") +
  theme_bw() + theme(strip.text = element_text(size=16),axis.text.x=element_text(size=14,angle=60,hjust=1),axis.title.x = element_text(size=16),axis.text.y=element_text(size=14),axis.title.y = element_text(size=16),legend.title = element_text(size=16),legend.text = ggtext::element_markdown(size=14)) + guides(fill=guide_legend(title="Strains")) + scale_y_continuous(limits = c(0,1.1),breaks = seq(0,1,0.2)) + scale_fill_brewer(palette = "Dark2",direction = -1) + facet_wrap(~factor(population,v_pop_order),ncol=2,scales="fixed")
ggsave(filename = "Contamination_level_across_S_pombe_populations.eps", path=output_workspace, width = 30/2.54, height = 25/2.54, device = cairo_ps)
ggsave(filename = "Contamination_level_across_S_pombe_populations.svg", path=output_workspace, width = 30/2.54, height = 25/2.54, device = svg)
ggsave(filename = "Contamination_level_across_S_pombe_populations.png", path=output_workspace, width = 30/2.54, height = 25/2.54, device = png)

#Find noncoding regions length
v_lst_noncoding_regions <- sort(unique(subset(df_loci_annotations,((is.na(the_product)&is.na(accession_number))|((!is.na(accession_number))& grepl(pattern = "rna",x = tolower(accession_number),fixed = T) )))$accession_number))
v_noncoding_regions_length <- rep(NA,length(v_lst_noncoding_regions))
names(v_noncoding_regions_length) <- v_lst_noncoding_regions
for (current_noncoding_region in names(v_noncoding_regions_length)){
  current_df_loci_annotation <- subset(df_loci_annotations,(accession_number==current_noncoding_region))
  if (nrow(current_df_loci_annotation)>=1){
    v_noncoding_regions_length[current_noncoding_region] <- max(current_df_loci_annotation$end - current_df_loci_annotation$start + 1)
  }
}

#Estimate the coding density in S. pombe genome
subset_df_CDS_annotations <- subset(df_loci_annotations,(site_type=="CDS") & (!( (grepl(pattern = "ncrna",x = tolower(annotation),fixed = T)|grepl(pattern = "trna",x = tolower(annotation),fixed = T)) )))
the_CDS_length <- sum(ifelse(test = subset_df_CDS_annotations$end > subset_df_CDS_annotations$start,yes=(subset_df_CDS_annotations$end - subset_df_CDS_annotations$start + 1),no=(subset_df_CDS_annotations$start - subset_df_CDS_annotations$end + 1) ))
proportion_of_genome_spanned_by_annoted_CDS <- 100*the_CDS_length/(12E6)
proportion_of_genome_spanned_by_annoted_CDS

#Find multi-hit noncoding regions
df_ALL_noncoding_regions_multiplicity <- subset(df_ALL_variants,((is.na(gene_product)&is.na(accession_number))|(grepl(pattern = "noncoding",x = tolower(variant_type2),fixed = T)))&is_fixed) %>%
  group_by(mutation,accession_number) %>%
  dplyr::summarise(gene_product=gene_product, variant_type2=variant_type2, PH=length(unique(population)),multiplicity=length(unique(Position))/v_noncoding_regions_length[gene_product])
df_ALL_noncoding_regions_multiplicity <- subset(df_ALL_noncoding_regions_multiplicity,((is.na(gene_product)&is.na(accession_number))|(grepl(pattern = "noncoding",x = tolower(variant_type2),fixed = T))))

#Get list of multi-hit noncoding regions (with fixed mutations in multiple pop; mutation, gene, accession_number,length(unique(pop)) and subset for count >=2)
df_multihit_noncoding_regions <- subset(df_ALL_noncoding_regions_multiplicity,PH>=2)
v_lst_multihit_noncoding_regions_accession <- sort(unique(df_multihit_noncoding_regions$accession_number))
#***Run enrichment analysis with list of accession numbers on AnGeLi
#Between-population frequency of multi-hit noncoding regions hits (fixed mutations)
df_ALL_noncoding_regions_multiplicity$prevalence <- df_ALL_noncoding_regions_multiplicity$PH/length(v_pop_order)

#multihit UTR
print(v_lst_multihit_noncoding_regions_accession)
#recurrent noncoding regions mutations (not necessarily fixed)
print(sort(table(subset(df_ALL_variants,grepl(pattern = "noncoding",x = tolower(variant_type2),fixed = T)&grepl(pattern = "ncrna",x = tolower(accession_number),fixed = T))$mutation)))

#Odd ratio test for fixation bias of var_type_2
df_ALL_variants$is_fixed_in_pop_at_some_point <- df_ALL_variants$pop_mut %in% unique(subset(df_ALL_variants,is_fixed)$pop_mut)
print(epitools::oddsratio(table(factor(unique(df_ALL_variants[,c("pop_mut","variant_type2","is_fixed_in_pop_at_some_point")])$variant_type2,
                                       levels=c("synonymous","nonsense","missense","Indel in introns","SNV in introns","Indel in UTR","SNV in UTR","SNV in other noncoding regions","Indel in protein-coding genes","Indel in other noncoding regions")),
                                unique(df_ALL_variants[,c("pop_mut","variant_type2","is_fixed_in_pop_at_some_point")])$is_fixed_in_pop_at_some_point)))

# #Compare the density of parallel hits in coding vs noncoding regions
#   # assuming chromosome lengths are already known and stored in a vector "v_chrs_length"
# the_window_size <- 10000
# v_starts <- NULL
# for (the_chr in names(v_chrs_length)){
#   v_starts <- c(v_starts,unique(seq(1, v_chrs_length[the_chr], by = the_window_size)))
# }
# v_ends <- NULL
# for (the_chr in names(v_chrs_length)){
#   v_ends <- c(v_ends,unique(c(seq(the_window_size, v_chrs_length[the_chr], by = the_window_size),v_chrs_length[the_chr])))
# }
# df_per_site_type_windows <- data.frame(chr = rep(names(v_chrs_length), times = ceiling(v_chrs_length / the_window_size)), start = v_starts, end = v_ends)
# df_per_site_type_windows$concatenated_genome_start <- NA
# df_per_site_type_windows$concatenated_genome_end <- NA
# df_per_site_type_windows$concatenated_genome_start[1] <- 1
# df_per_site_type_windows$concatenated_genome_end[1] <- df_per_site_type_windows$concatenated_genome_start[1] + df_per_site_type_windows$end[1] - df_per_site_type_windows$start[1]
# for (i in 2:nrow(df_per_site_type_windows)){
#   df_per_site_type_windows$concatenated_genome_start[i] <- df_per_site_type_windows$concatenated_genome_end[i-1] + 1
#   df_per_site_type_windows$concatenated_genome_end[i] <- df_per_site_type_windows$concatenated_genome_start[i] + df_per_site_type_windows$end[i] - df_per_site_type_windows$start[i]
# }
# #Determine the end position of each chromosome in the concatenated genome
# df_end_chrs_in_concat_genome <- df_per_site_type_windows[,c("chr","concatenated_genome_end")] %>%
#   group_by(chr) %>%
#   dplyr::summarise(concat_end = max(concatenated_genome_end,na.rm=T))
# #Function to determine if a variant is in a coding sequence of a gene
# is_variant_in_a_prot_coding_gene_CDS <- function(the_chromosome,the_pos_in_chr){
#   subset_cds <- subset(df_loci_annotations,site_type=="CDS"&chr==the_chromosome&the_pos_in_chr>=start&the_pos_in_chr<=end)
#   return(nrow(subset(subset_cds, (!any(vapply(X = c("lncRNA", "snRNA", "sncRNA", "snoRNA", "rRNA","NCRNA","nRNA","NRNA","NORNA","RRNA"),FUN = function(the_noncoding_site_type) any(grepl(pattern = the_noncoding_site_type,x = annotation,fixed = T)),FUN.VALUE = F)))))>0)
# }
# #Determine the number of parallel hits (fixed mutations in at least 2 populations) in coding and noncoding regions
#   #variant recurrence
# df_recurrence_ALL_variants <- subset(df_ALL_variants, is_fixed) %>%
#   group_by(mutation,Chrom,Position,gene_product,accession_number,variant_type2) %>%
#   dplyr::summarise(recurrence=length(unique(population)))
# df_recurrence_ALL_variants$is_in_prot_gene_CDS <- FALSE
# for (i in 1:nrow(df_recurrence_ALL_variants)){
#   if (is_variant_in_a_prot_coding_gene_CDS(the_chromosome = df_recurrence_ALL_variants$Chrom[i],the_pos_in_chr = df_recurrence_ALL_variants$Position[i])){
#     df_recurrence_ALL_variants$is_in_prot_gene_CDS[i] <- TRUE
#   }
# }
# df_recurrence_ALL_variants$is_in_gene_UTR <- grepl(pattern = "UTR",x = df_recurrence_ALL_variants$variant_type2,ignore.case = T)
# df_recurrence_ALL_variants$is_in_gene_intron <- grepl(pattern = "intron",x = df_recurrence_ALL_variants$variant_type2,ignore.case = T)
# #number of recurrent variant of each type in each sliding window
# df_per_site_type_windows$non_partitioned_window_id <- 1:nrow(df_per_site_type_windows)
# df_per_site_type_windows$non_partitioned_window_id_2 <- paste0(df_per_site_type_windows$chr,"_",df_per_site_type_windows$start,"-",df_per_site_type_windows$end)
# nb_windows <- nrow(df_per_site_type_windows)
# df_per_site_type_windows$site_type <- "protein-coding"
# df_per_site_type_windows$nb_recurrent_mutations <- 0
# df_per_site_type_windows$length_of_site_type <- NA
# df_per_site_type_windows$average_mutation_recurrence <- 0
# for (i in 1:nb_windows){
#   current_line_df <- df_per_site_type_windows[i,]
#   #Calculating protein-coding regions stats
#   subset_prot_coding_gene_CDS_df_loci_annotations <- subset(df_loci_annotations,(site_type=="CDS")&(chr==df_per_site_type_windows$chr[i])&((start<=df_per_site_type_windows$start[i]&end>=df_per_site_type_windows$start[i])|(start>=df_per_site_type_windows$start[i]&start<=df_per_site_type_windows$end[i])) &
#                                               (!( (grepl(pattern = "ncrna",x = tolower(annotation),fixed = T)|grepl(pattern = "trna",x = tolower(annotation),fixed = T))&(!grepl(pattern = ";",x = tolower(annotation),fixed = T)) )))
#   
#   if (nrow(subset_prot_coding_gene_CDS_df_loci_annotations)>0){
#     v_is_cds_prot_coding <- rep(F,nrow(subset_prot_coding_gene_CDS_df_loci_annotations))
#     for (indx in 1:nrow(subset_prot_coding_gene_CDS_df_loci_annotations)){
#       v_is_cds_prot_coding[indx] <- (!any(vapply(X = c("lncRNA", "snRNA", "sncRNA", "snoRNA", "rRNA","NCRNA","nRNA","NRNA","NORNA","RRNA"),FUN = function(the_noncoding_site_type) any(grepl(pattern = the_noncoding_site_type,x = subset_prot_coding_gene_CDS_df_loci_annotations$annotation[indx],fixed = T)),FUN.VALUE = F)))
#     }
#     subset_prot_coding_gene_CDS_df_loci_annotations <- subset(subset_prot_coding_gene_CDS_df_loci_annotations,v_is_cds_prot_coding)
#   }
#   v_pos_prot_cod_CDS <- NULL
#   if (nrow(subset_prot_coding_gene_CDS_df_loci_annotations)>0){
#     for (j in 1:nrow(subset_prot_coding_gene_CDS_df_loci_annotations)){
#       v_pos_prot_cod_CDS <- c(v_pos_prot_cod_CDS,intersect(df_per_site_type_windows$start[i]:df_per_site_type_windows$end[i],subset_prot_coding_gene_CDS_df_loci_annotations$start[j]:subset_prot_coding_gene_CDS_df_loci_annotations$end[j]))
#     }
#   }
#   df_per_site_type_windows$length_of_site_type[i] <- length(v_pos_prot_cod_CDS)
#   remove(v_pos_prot_cod_CDS)
#   df_per_site_type_windows$nb_recurrent_mutations[i] <- length(unique(subset(df_recurrence_ALL_variants,Chrom==df_per_site_type_windows$chr[i]&Position>=df_per_site_type_windows$start[i]&Position<=df_per_site_type_windows$end[i]&recurrence>1&is_in_prot_gene_CDS)$mutation))
#   df_per_site_type_windows$average_mutation_recurrence[i] <- mean(subset(df_recurrence_ALL_variants,Chrom==df_per_site_type_windows$chr[i]&Position>=df_per_site_type_windows$start[i]&Position<=df_per_site_type_windows$end[i]&recurrence>1&is_in_prot_gene_CDS)$recurrence,na.rm = T)
#   
#   #Calculating UTRs stats
#   utr_new_line_df <- current_line_df
#   utr_new_line_df$site_type <- "UTR"
#   subset_UTR_df_loci_annotations <- subset(df_loci_annotations,site_type%in%c("three_prime_UTR","five_prime_UTR"))
#   v_pos_UTR <- NULL
#   if (nrow(subset_UTR_df_loci_annotations)>0){
#     for (j in 1:nrow(subset_UTR_df_loci_annotations)){
#       v_pos_UTR <- c(v_pos_UTR,intersect(df_per_site_type_windows$start[i]:df_per_site_type_windows$end[i],subset_UTR_df_loci_annotations$start[j]:subset_UTR_df_loci_annotations$end[j]))
#     }
#   }
#   utr_new_line_df$length_of_site_type[1] <- length(v_pos_UTR)
#   remove(v_pos_UTR)
#   utr_new_line_df$nb_recurrent_mutations[1] <- length(unique(subset(df_recurrence_ALL_variants,Chrom==df_per_site_type_windows$chr[i]&Position>=df_per_site_type_windows$start[i]&Position<=df_per_site_type_windows$end[i]&recurrence>1&is_in_gene_UTR)$mutation))
#   utr_new_line_df$average_mutation_recurrence[1] <- mean(subset(df_recurrence_ALL_variants,Chrom==df_per_site_type_windows$chr[i]&Position>=df_per_site_type_windows$start[i]&Position<=df_per_site_type_windows$end[i]&recurrence>1&is_in_gene_UTR)$recurrence,na.rm=T)
#   
#   #Calculating Introns stats
#   intron_new_line_df <- current_line_df
#   intron_new_line_df$site_type <- "Intron"
#   subset_intron_df_loci_annotations  <- subset(df_loci_annotations,site_type=="intron")
#   v_pos_intron <- NULL
#   if (nrow(subset_intron_df_loci_annotations)>0){
#     for (j in 1:nrow(subset_intron_df_loci_annotations)){
#       v_pos_intron <- c(v_pos_intron,intersect(df_per_site_type_windows$start[i]:df_per_site_type_windows$end[i],subset_intron_df_loci_annotations$start[j]:subset_intron_df_loci_annotations$end[j]))
#     }
#   }
#   intron_new_line_df$length_of_site_type[1] <- length(v_pos_intron)
#   remove(v_pos_intron)
#   intron_new_line_df$nb_recurrent_mutations[1] <- length(unique(subset(df_recurrence_ALL_variants,Chrom==df_per_site_type_windows$chr[i]&Position>=df_per_site_type_windows$start[i]&Position<=df_per_site_type_windows$end[i]&recurrence>1&is_in_gene_intron)$mutation))
#   intron_new_line_df$average_mutation_recurrence[1] <- mean(subset(df_recurrence_ALL_variants,Chrom==df_per_site_type_windows$chr[i]&Position>=df_per_site_type_windows$start[i]&Position<=df_per_site_type_windows$end[i]&recurrence>1&is_in_gene_intron)$recurrence,na.rm=T)
#   
#   #Calculating Other noncoding regions stats
#   Other_noncoding_regions_new_line_df <- current_line_df
#   Other_noncoding_regions_new_line_df$site_type[1] <- "Other noncoding regions"
#   Other_noncoding_regions_new_line_df$length_of_site_type[1] <- the_window_size - (df_per_site_type_windows$length_of_site_type[i]+utr_new_line_df$length_of_site_type[1]+intron_new_line_df$length_of_site_type[1])
#   Other_noncoding_regions_new_line_df$nb_recurrent_mutations[1] <- length(unique(subset(df_recurrence_ALL_variants,Chrom==df_per_site_type_windows$chr[i]&Position>=df_per_site_type_windows$start[i]&Position<=df_per_site_type_windows$end[i]&recurrence>1&(!(is_in_prot_gene_CDS|is_in_gene_intron|is_in_gene_UTR)))$mutation))
#   Other_noncoding_regions_new_line_df$average_mutation_recurrence[1] <- mean(subset(df_recurrence_ALL_variants,Chrom==df_per_site_type_windows$chr[i]&Position>=df_per_site_type_windows$start[i]&Position<=df_per_site_type_windows$end[i]&recurrence>1&(!(is_in_prot_gene_CDS|is_in_gene_intron|is_in_gene_UTR)))$recurrence,na.rm=T)
#   
#   #Calculating All noncoding regions stats
#   All_noncoding_regions_new_line_df <- current_line_df
#   All_noncoding_regions_new_line_df$site_type[1] <- "noncoding"
#   All_noncoding_regions_new_line_df$length_of_site_type[1] <- the_window_size - df_per_site_type_windows$length_of_site_type[i]
#   All_noncoding_regions_new_line_df$nb_recurrent_mutations[1] <- length(unique(subset(df_recurrence_ALL_variants,Chrom==df_per_site_type_windows$chr[i]&Position>=df_per_site_type_windows$start[i]&Position<=df_per_site_type_windows$end[i]&recurrence>1&(!is_in_prot_gene_CDS))$mutation))
#   All_noncoding_regions_new_line_df$average_mutation_recurrence[1] <- mean(subset(df_recurrence_ALL_variants,Chrom==df_per_site_type_windows$chr[i]&Position>=df_per_site_type_windows$start[i]&Position<=df_per_site_type_windows$end[i]&recurrence>1&(!is_in_prot_gene_CDS))$recurrence,na.rm=T)
#   
#   #add new lines
#   df_per_site_type_windows <- rbind(df_per_site_type_windows,utr_new_line_df)
#   df_per_site_type_windows <- rbind(df_per_site_type_windows,intron_new_line_df)
#   df_per_site_type_windows <- rbind(df_per_site_type_windows,Other_noncoding_regions_new_line_df)
#   df_per_site_type_windows <- rbind(df_per_site_type_windows,All_noncoding_regions_new_line_df)
# }
# 
# df_per_site_type_windows$density_parallel_hits <- df_per_site_type_windows$nb_recurrent_mutations/df_per_site_type_windows$length_of_site_type
# df_per_site_type_windows$log10_density_parallel_hits <- log10((df_per_site_type_windows$nb_recurrent_mutations/df_per_site_type_windows$length_of_site_type)+1e-7)
# 
# #Sanity check about the number of recurrent hits in the windows map
# print(sum(subset(df_per_site_type_windows,site_type=="protein-coding")$nb_recurrent_mutations))
# print(sum(subset(df_per_site_type_windows,site_type=="noncoding")$nb_recurrent_mutations))
# 
# #Summarise window statistics
# v_lst_windows <- sort(unique(df_per_site_type_windows$non_partitioned_window_id_2))
# #In each window, determine:
#   #whether or not there is more recurrent hits in noncoding regions compared to protein-coding regions (solely based on the number of recurrent hits SO IT IS NOT WEIGHTED FOR THE EFFECT OF THESE MUTATIONS NOR THE LEVEL OF RECURRENCE OF THESE MUTATIONS)
#   #the label of the regions ranking based on the number of recurrent hits at resolution type 2
#   #whether or not the recurrent hits aveage prevalence is higher in noncoding regions compared to protein-coding regions (solely based on the level of recurrence not the number of recurrent hits)
#   #the label of the regions ranking based on the number of recurrent hits
# v_label_noncoding_regions <- c("Intron","UTR","ONC")
# v_windows_to_location_with_highest_nb_of_recurrent_hits_resolution_1 <- NULL
# v_windows_to_regions_ranking_based_on_highest_nb_recurrent_hits_resolution_2 <- NULL
# v_windows_to_location_with_highest_Density_of_recurrent_hits_resolution_1 <- NULL
# v_windows_to_regions_ranking_based_on_highest_Density_of_recurrent_hits_resolution_2 <- NULL
# v_windows_to_location_with_highest_recurrent_hits_avg_prevalence_resolution_1 <- NULL
# v_windows_to_regions_ranking_based_on_highest_avg_prevalence_recurrent_hits_resolution_2 <- NULL
# v_windows_to_is_windows_mutated <- NULL
# v_windows_to_is_window_having_only_noncoding_rec_hits <- NULL
# v_windows_to_is_window_having_only_prot_coding_rec_hits <- NULL
# v_windows_to_is_window_having_both_prot_coding_and_noncoding_rec_hits <- NULL
# for (current_window_id_2 in v_lst_windows){
#   current_window_df <- subset(df_per_site_type_windows,non_partitioned_window_id_2==current_window_id_2)
#   current_window_nb_recurrent_hits_in_protein_coding_genes <- subset(current_window_df,site_type=="protein-coding")$nb_recurrent_mutations
#   current_window_nb_recurrent_hits_in_noncoding_regions <- subset(current_window_df,site_type=="noncoding")$nb_recurrent_mutations
#   current_window_nb_recurrent_hits_in_UTRs <- subset(current_window_df,site_type=="UTR")$nb_recurrent_mutations
#   current_window_nb_recurrent_hits_in_Introns <- subset(current_window_df,site_type=="Intron")$nb_recurrent_mutations
#   current_window_nb_recurrent_hits_in_Other_noncoding_regions <- subset(current_window_df,site_type=="Other noncoding regions")$nb_recurrent_mutations
#   
#   current_window_Density_of_recurrent_hits_in_protein_coding_genes <- subset(current_window_df,site_type=="protein-coding")$density_parallel_hits
#   current_window_Density_of_recurrent_hits_in_noncoding_regions <- subset(current_window_df,site_type=="noncoding")$density_parallel_hits
#   
#   current_window_avg_prevalence_in_protein_coding_genes <- subset(current_window_df,site_type=="protein-coding")$average_mutation_recurrence
#   current_window_avg_prevalence_in_protein_coding_genes <- ifelse(test = is.nan(current_window_avg_prevalence_in_protein_coding_genes),yes = 0,no = current_window_avg_prevalence_in_protein_coding_genes)
#   current_window_avg_prevalence_in_noncoding_regions <- subset(current_window_df,site_type=="noncoding")$average_mutation_recurrence
#   current_window_avg_prevalence_in_noncoding_regions <- ifelse(test = is.nan(current_window_avg_prevalence_in_noncoding_regions),yes = 0,no = current_window_avg_prevalence_in_noncoding_regions)
#   current_window_avg_prevalence_in_UTRs <- subset(current_window_df,site_type=="UTR")$average_mutation_recurrence
#   current_window_avg_prevalence_in_UTRs <- ifelse(test = is.nan(current_window_avg_prevalence_in_UTRs),yes = 0,no = current_window_avg_prevalence_in_UTRs)
#   current_window_avg_prevalence_in_Introns <- subset(current_window_df,site_type=="Intron")$average_mutation_recurrence
#   current_window_avg_prevalence_in_Introns <- ifelse(test = is.nan(current_window_avg_prevalence_in_Introns),yes = 0,no = current_window_avg_prevalence_in_Introns)
#   current_window_avg_prevalence_in_Other_noncoding_regions <- subset(current_window_df,site_type=="Other noncoding regions")$average_mutation_recurrence
#   current_window_avg_prevalence_in_Other_noncoding_regions <- ifelse(test = is.nan(current_window_avg_prevalence_in_Other_noncoding_regions),yes = 0,no = current_window_avg_prevalence_in_Other_noncoding_regions)
#   
#   v_windows_to_location_with_highest_nb_of_recurrent_hits_resolution_1 <- c(v_windows_to_location_with_highest_nb_of_recurrent_hits_resolution_1,ifelse(test = current_window_nb_recurrent_hits_in_noncoding_regions==current_window_nb_recurrent_hits_in_protein_coding_genes,yes = "=",no = ifelse(test = current_window_nb_recurrent_hits_in_noncoding_regions>current_window_nb_recurrent_hits_in_protein_coding_genes,yes = "noncoding",no = "protein-coding")))
#   v_decreasing_order_nb_rec_hits_Intron_UTR_ONC <- rev(order(c(current_window_nb_recurrent_hits_in_Introns,current_window_nb_recurrent_hits_in_UTRs,current_window_nb_recurrent_hits_in_Other_noncoding_regions)))
#   if (v_windows_to_location_with_highest_nb_of_recurrent_hits_resolution_1[length(v_windows_to_location_with_highest_nb_of_recurrent_hits_resolution_1)]=="noncoding"){
#     v_windows_to_regions_ranking_based_on_highest_nb_recurrent_hits_resolution_2 <- c(v_windows_to_regions_ranking_based_on_highest_nb_recurrent_hits_resolution_2,paste0("NC>PC & ",v_label_noncoding_regions[v_decreasing_order_nb_rec_hits_Intron_UTR_ONC[1]],">=",v_label_noncoding_regions[v_decreasing_order_nb_rec_hits_Intron_UTR_ONC[2]],">=",v_label_noncoding_regions[v_decreasing_order_nb_rec_hits_Intron_UTR_ONC[3]]))
#   }else if (v_windows_to_location_with_highest_nb_of_recurrent_hits_resolution_1[length(v_windows_to_location_with_highest_nb_of_recurrent_hits_resolution_1)]=="protein-coding"){
#     v_windows_to_regions_ranking_based_on_highest_nb_recurrent_hits_resolution_2 <- c(v_windows_to_regions_ranking_based_on_highest_nb_recurrent_hits_resolution_2,paste0("NC<PC & ",v_label_noncoding_regions[v_decreasing_order_nb_rec_hits_Intron_UTR_ONC[1]],">=",v_label_noncoding_regions[v_decreasing_order_nb_rec_hits_Intron_UTR_ONC[2]],">=",v_label_noncoding_regions[v_decreasing_order_nb_rec_hits_Intron_UTR_ONC[3]]))
#   }else{
#     v_windows_to_regions_ranking_based_on_highest_nb_recurrent_hits_resolution_2 <- c(v_windows_to_regions_ranking_based_on_highest_nb_recurrent_hits_resolution_2,paste0("NC=PC & ",v_label_noncoding_regions[v_decreasing_order_nb_rec_hits_Intron_UTR_ONC[1]],">=",v_label_noncoding_regions[v_decreasing_order_nb_rec_hits_Intron_UTR_ONC[2]],">=",v_label_noncoding_regions[v_decreasing_order_nb_rec_hits_Intron_UTR_ONC[3]]))
#   }
#   
#   v_windows_to_location_with_highest_Density_of_recurrent_hits_resolution_1 <- c(v_windows_to_location_with_highest_Density_of_recurrent_hits_resolution_1,ifelse(test = current_window_nb_recurrent_hits_in_noncoding_regions==current_window_nb_recurrent_hits_in_protein_coding_genes,yes = "=",no = ifelse(test = current_window_nb_recurrent_hits_in_noncoding_regions>current_window_nb_recurrent_hits_in_protein_coding_genes,yes = "noncoding",no = "protein-coding")))
#   v_decreasing_order_nb_rec_hits_Intron_UTR_ONC <- rev(order(c(current_window_nb_recurrent_hits_in_Introns,current_window_nb_recurrent_hits_in_UTRs,current_window_nb_recurrent_hits_in_Other_noncoding_regions)))
#   if (v_windows_to_location_with_highest_Density_of_recurrent_hits_resolution_1[length(v_windows_to_location_with_highest_Density_of_recurrent_hits_resolution_1)]=="noncoding"){
#     v_windows_to_regions_ranking_based_on_highest_Density_of_recurrent_hits_resolution_2 <- c(v_windows_to_regions_ranking_based_on_highest_Density_of_recurrent_hits_resolution_2,paste0("NC>PC & ",v_label_noncoding_regions[v_decreasing_order_nb_rec_hits_Intron_UTR_ONC[1]],">=",v_label_noncoding_regions[v_decreasing_order_nb_rec_hits_Intron_UTR_ONC[2]],">=",v_label_noncoding_regions[v_decreasing_order_nb_rec_hits_Intron_UTR_ONC[3]]))
#   }else if (v_windows_to_location_with_highest_Density_of_recurrent_hits_resolution_1[length(v_windows_to_location_with_highest_Density_of_recurrent_hits_resolution_1)]=="protein-coding"){
#     v_windows_to_regions_ranking_based_on_highest_Density_of_recurrent_hits_resolution_2 <- c(v_windows_to_regions_ranking_based_on_highest_Density_of_recurrent_hits_resolution_2,paste0("NC<PC & ",v_label_noncoding_regions[v_decreasing_order_nb_rec_hits_Intron_UTR_ONC[1]],">=",v_label_noncoding_regions[v_decreasing_order_nb_rec_hits_Intron_UTR_ONC[2]],">=",v_label_noncoding_regions[v_decreasing_order_nb_rec_hits_Intron_UTR_ONC[3]]))
#   }else{
#     v_windows_to_regions_ranking_based_on_highest_Density_of_recurrent_hits_resolution_2 <- c(v_windows_to_regions_ranking_based_on_highest_Density_of_recurrent_hits_resolution_2,paste0("NC=PC & ",v_label_noncoding_regions[v_decreasing_order_nb_rec_hits_Intron_UTR_ONC[1]],">=",v_label_noncoding_regions[v_decreasing_order_nb_rec_hits_Intron_UTR_ONC[2]],">=",v_label_noncoding_regions[v_decreasing_order_nb_rec_hits_Intron_UTR_ONC[3]]))
#   }
#   
#   v_windows_to_location_with_highest_recurrent_hits_avg_prevalence_resolution_1 <- c(v_windows_to_location_with_highest_recurrent_hits_avg_prevalence_resolution_1,ifelse(test = current_window_avg_prevalence_in_noncoding_regions==current_window_avg_prevalence_in_protein_coding_genes,yes = "=",no = ifelse(test = current_window_avg_prevalence_in_noncoding_regions>current_window_avg_prevalence_in_protein_coding_genes,yes = "noncoding",no = "protein-coding")))
#   v_decreasing_order_avg_prevalence_rec_hits_in_Intron_UTR_ONC <- rev(order(c(current_window_avg_prevalence_in_Introns,current_window_avg_prevalence_in_UTRs,current_window_avg_prevalence_in_Other_noncoding_regions)))
#   if (v_windows_to_location_with_highest_recurrent_hits_avg_prevalence_resolution_1[length(v_windows_to_location_with_highest_recurrent_hits_avg_prevalence_resolution_1)]=="noncoding"){
#     v_windows_to_regions_ranking_based_on_highest_avg_prevalence_recurrent_hits_resolution_2 <- c(v_windows_to_regions_ranking_based_on_highest_avg_prevalence_recurrent_hits_resolution_2,paste0("NC>PC & ",v_label_noncoding_regions[v_decreasing_order_avg_prevalence_rec_hits_in_Intron_UTR_ONC[1]],">=",v_label_noncoding_regions[v_decreasing_order_avg_prevalence_rec_hits_in_Intron_UTR_ONC[2]],">=",v_label_noncoding_regions[v_decreasing_order_avg_prevalence_rec_hits_in_Intron_UTR_ONC[3]]))
#   }else if (v_windows_to_location_with_highest_recurrent_hits_avg_prevalence_resolution_1[length(v_windows_to_location_with_highest_recurrent_hits_avg_prevalence_resolution_1)]=="protein-coding"){
#     v_windows_to_regions_ranking_based_on_highest_avg_prevalence_recurrent_hits_resolution_2 <- c(v_windows_to_regions_ranking_based_on_highest_avg_prevalence_recurrent_hits_resolution_2,paste0("NC<PC & ",v_label_noncoding_regions[v_decreasing_order_avg_prevalence_rec_hits_in_Intron_UTR_ONC[1]],">=",v_label_noncoding_regions[v_decreasing_order_avg_prevalence_rec_hits_in_Intron_UTR_ONC[2]],">=",v_label_noncoding_regions[v_decreasing_order_avg_prevalence_rec_hits_in_Intron_UTR_ONC[3]]))
#   }else{
#     v_windows_to_regions_ranking_based_on_highest_avg_prevalence_recurrent_hits_resolution_2 <- c(v_windows_to_regions_ranking_based_on_highest_avg_prevalence_recurrent_hits_resolution_2,paste0("NC=PC & ",v_label_noncoding_regions[v_decreasing_order_avg_prevalence_rec_hits_in_Intron_UTR_ONC[1]],">=",v_label_noncoding_regions[v_decreasing_order_avg_prevalence_rec_hits_in_Intron_UTR_ONC[2]],">=",v_label_noncoding_regions[v_decreasing_order_avg_prevalence_rec_hits_in_Intron_UTR_ONC[3]]))
#   }
#   v_windows_to_is_windows_mutated <- c(v_windows_to_is_windows_mutated, (current_window_nb_recurrent_hits_in_protein_coding_genes+current_window_nb_recurrent_hits_in_noncoding_regions)>0)
#   v_windows_to_is_window_having_only_noncoding_rec_hits <- c(v_windows_to_is_window_having_only_noncoding_rec_hits, current_window_nb_recurrent_hits_in_protein_coding_genes==0&current_window_nb_recurrent_hits_in_noncoding_regions>0)
#   v_windows_to_is_window_having_only_prot_coding_rec_hits <- c(v_windows_to_is_window_having_only_prot_coding_rec_hits, current_window_nb_recurrent_hits_in_protein_coding_genes>0&current_window_nb_recurrent_hits_in_noncoding_regions==0)
#   v_windows_to_is_window_having_both_prot_coding_and_noncoding_rec_hits <- c(v_windows_to_is_window_having_both_prot_coding_and_noncoding_rec_hits, current_window_nb_recurrent_hits_in_protein_coding_genes>0&current_window_nb_recurrent_hits_in_noncoding_regions>0)
#   
# }
# names(v_windows_to_location_with_highest_nb_of_recurrent_hits_resolution_1) <- v_lst_windows
# names(v_windows_to_regions_ranking_based_on_highest_nb_recurrent_hits_resolution_2) <- v_lst_windows
# names(v_windows_to_location_with_highest_Density_of_recurrent_hits_resolution_1) <- v_lst_windows
# names(v_windows_to_regions_ranking_based_on_highest_Density_of_recurrent_hits_resolution_2) <- v_lst_windows
# names(v_windows_to_location_with_highest_recurrent_hits_avg_prevalence_resolution_1) <- v_lst_windows
# names(v_windows_to_regions_ranking_based_on_highest_avg_prevalence_recurrent_hits_resolution_2) <- v_lst_windows
# names(v_windows_to_is_windows_mutated) <- v_lst_windows
# names(v_windows_to_is_window_having_only_noncoding_rec_hits) <- v_lst_windows
# names(v_windows_to_is_window_having_only_prot_coding_rec_hits) <- v_lst_windows
# names(v_windows_to_is_window_having_both_prot_coding_and_noncoding_rec_hits) <- v_lst_windows
# 
# palette_window_type_resolution_1 <- c("protein-coding"="#E6AB02","noncoding "="#6A3D9A","="="")
# palette_window_type_resolution_2 <- RColorBrewer::brewer.pal(12,"Paired")[RColorBrewer::brewer.pal(12,"Paired")!= "#6A3D9A"]
# names(palette_window_type_resolution_2) <- unique(unname( c(v_windows_to_regions_ranking_based_on_highest_nb_recurrent_hits_resolution_2,v_windows_to_regions_ranking_based_on_highest_avg_prevalence_recurrent_hits_resolution_2) ))
# df_per_site_type_windows$color_edge_nb_rec_hits_resolution_1 <- v_windows_to_location_with_highest_nb_of_recurrent_hits_resolution_1[df_per_site_type_windows$non_partitioned_window_id_2]
# df_per_site_type_windows$color_edge_nb_rec_hits_resolution_2 <- v_windows_to_regions_ranking_based_on_highest_nb_recurrent_hits_resolution_2[df_per_site_type_windows$non_partitioned_window_id_2]
# df_per_site_type_windows$color_edge_Density_of_rec_hits_resolution_1 <- v_windows_to_location_with_highest_Density_of_recurrent_hits_resolution_1[df_per_site_type_windows$non_partitioned_window_id_2]
# df_per_site_type_windows$color_edge_Density_of_rec_hits_resolution_2 <- v_windows_to_regions_ranking_based_on_highest_Density_of_recurrent_hits_resolution_2[df_per_site_type_windows$non_partitioned_window_id_2]
# df_per_site_type_windows$color_edge_avg_prevalence_resolution_1 <- v_windows_to_location_with_highest_recurrent_hits_avg_prevalence_resolution_1[df_per_site_type_windows$non_partitioned_window_id_2]
# df_per_site_type_windows$color_edge_avg_prevalence_resolution_2 <- v_windows_to_regions_ranking_based_on_highest_avg_prevalence_recurrent_hits_resolution_2[df_per_site_type_windows$non_partitioned_window_id_2]
# #nb windows that only have noncoding recurrent hits
# print(sum(v_windows_to_is_window_having_only_noncoding_rec_hits))
# #nb windows that only have protein-coding recurrent hits
# print(sum(v_windows_to_is_window_having_only_prot_coding_rec_hits))
# #nb windows that have both noncoding and protein-coding recurrent hits
# print(sum(v_windows_to_is_window_having_both_prot_coding_and_noncoding_rec_hits))
# #nb windows that WITHOUT recurrent hits
# print(length(v_lst_windows)-sum(c(v_windows_to_is_window_having_only_noncoding_rec_hits,v_windows_to_is_window_having_only_prot_coding_rec_hits,v_windows_to_is_window_having_both_prot_coding_and_noncoding_rec_hits)))
# #FOCUS ON WINDOWS WITH MUTATIONS
# df_per_site_type_windows$length_of_site_type <- abs(df_per_site_type_windows$length_of_site_type)
# df_per_site_type_windows$density_parallel_hits <- abs(df_per_site_type_windows$density_parallel_hits)
# df_per_site_type_windows$log10_density_parallel_hits <- log10(df_per_site_type_windows$density_parallel_hits + 1e-7)
# df_per_site_type_windows <- subset(df_per_site_type_windows,length_of_site_type<=10000)
# df_per_site_type_mutated_windows <- subset(df_per_site_type_windows,v_windows_to_is_windows_mutated[non_partitioned_window_id_2])
# 
# #Compare noncoding regions to protein-coding regions for the density of recurrent hits and color lines
# ggplot(data = subset(df_per_site_type_mutated_windows,site_type%in%rev(c("protein-coding","noncoding"))),mapping=aes(x=site_type,y=log10_density_parallel_hits)) + geom_violin(aes(fill=site_type)) + geom_point() + xlab("Site type") + ylab("Density of recurrent hits\nlog10(Number of recurrent hits/bp + 1e-7)") + scale_fill_viridis_d() + theme_bw() +theme(axis.title = element_text(size=22),axis.text = element_text(size=22),legend.title = element_text(size=12),legend.text = element_text(size=8),axis.text.x = element_text(size=20),legend.position = "none")  + stat_compare_means(method = "wilcox") + geom_line(mapping=aes(group=non_partitioned_window_id_2,col=))
# ggsave(filename = "Comparison_of_recurrent_hits_density_at_resolution_1.svg", path=output_workspace, width = 20, height = 16, units = "cm", device= svg)
# ggsave(filename = "Comparison_of_recurrent_hits_density_at_resolution_1.png", path=output_workspace, width = 20, height = 16, units = "cm")
# 
# #Compare noncoding regions to protein-coding regions for the average recurrence OR prevalence of recurrent hits and color lines
# ggplot(data =  subset(df_per_site_type_mutated_windows,site_type%in%rev(c("protein-coding","noncoding"))),mapping=aes(x=reorder(site_type,average_mutation_recurrence),y=average_mutation_recurrence)) + geom_violin(aes(fill=reorder(site_type,average_mutation_recurrence))) + geom_point() + xlab("Site type") + ylab("Number of recurrent hits") + scale_fill_viridis_d() + theme_bw() +theme(axis.title = element_text(size=22),axis.text = element_text(size=22),legend.title = element_text(size=12),legend.text = element_text(size=8),axis.text.x = element_text(size=20),legend.position = "none")  + stat_compare_means(method = "wilcox") + geom_line(mapping=aes(group=non_partitioned_window_id_2)) 
# ggsave(filename = "Comparison_of_recurrent_hits_avg_prevalence_at_resolution_2.svg", path=output_workspace, width = 20, height = 15, units = "cm", device= svg)
# ggsave(filename = "Comparison_o_recurrent_hits_avg_prevalence_at_resolution_2.png", path=output_workspace, width = 20, height = 15, units = "cm")
# 
# #Piechart of the location with the highest number of recurrent hits with window type resolution 1 (noncoding vs protein-coding)
# library(dplyr)
# data_location_with_highest_nb_rec_hits_at_resolution_1 <- data.frame(location=factor(unique(df_per_site_type_mutated_windows[,c("non_partitioned_window_id_2","color_edge_nb_rec_hits_resolution_1")])$color_edge_nb_rec_hits_resolution_1,levels = c("noncoding","protein-coding","=")),stringsAsFactors = FALSE)
# df_location_with_highest_nb_rec_hits_at_resolution_1 <- data_location_with_highest_nb_rec_hits_at_resolution_1 %>% 
#   group_by(location) %>% # Variable to be transformed
#   count() %>% 
#   ungroup() %>% 
#   mutate(perc = `n` / sum(`n`)) %>% 
#   arrange(perc) %>%
#   mutate(labels = scales::percent(perc))
# 
# ggplot(df_location_with_highest_nb_rec_hits_at_resolution_1, aes(x = "", y = perc, fill = factor(location,levels=c("noncoding","=","protein-coding")) )) +
#   geom_col(color = "black") +
#   geom_label(aes(label = labels), color = c(1, "white", 1),
#              position = position_stack(vjust = 0.5),
#              show.legend = FALSE) +
#   guides(fill = guide_legend(title = "Proportion of genomic windows where\nthe sites with the highest number\nof recurrent hits are:")) +
#   scale_fill_viridis_d() +
#   coord_polar(theta = "y") + 
#   theme_void()
# ggsave(filename = "Piechart_proportion_windows_in_which_site_type1_with_highest_nb_rec_hits_are.svg", path=output_workspace, width = 20, height = 15, units = "cm", device= svg)
# ggsave(filename = "Piechart_proportion_windows_in_which_site_type1_with_highest_nb_rec_hits_are.png", path=output_workspace, width = 20, height = 15, units = "cm")
# 
# #Piechart of the regions ranking based on the highest number of recurrent hits at window type resolution 2
# data_location_with_highest_nb_rec_hits_at_resolution_2 <- data.frame(location=(unique(df_per_site_type_mutated_windows[,c("non_partitioned_window_id_2","color_edge_nb_rec_hits_resolution_2")])$color_edge_nb_rec_hits_resolution_2),stringsAsFactors = FALSE)
# df_location_with_highest_nb_rec_hits_at_resolution_2 <- data_location_with_highest_nb_rec_hits_at_resolution_2 %>% 
#   group_by(location) %>% #Variable to be transformed
#   count() %>% 
#   ungroup() %>% 
#   mutate(perc = `n` / sum(`n`)) %>% 
#   arrange(perc) %>%
#   mutate(labels = scales::percent(perc))
# 
# ggplot(df_location_with_highest_nb_rec_hits_at_resolution_2, aes(x = "", y = perc, fill = location )) +
#   geom_col(color = "black") +
#   geom_label(aes(label = labels), color = 1,
#              position = position_stack(vjust = 0.5),
#              show.legend = FALSE) +
#   guides(fill = guide_legend(title = "Proportion of genomic windows where\nthe site ranking is:")) +
#   scale_fill_manual(values = palette_window_type_resolution_2[intersect(names(palette_window_type_resolution_2),as.character(df_location_with_highest_nb_rec_hits_at_resolution_2$location))]) +
#   coord_polar(theta = "y") + 
#   theme_void()
# ggsave(filename = "Piechart_proportion_windows_where_site_ranking_based_on_nb_rec_hits_are.svg", path=output_workspace, width = 20, height = 15, units = "cm", device= svg)
# 
# #Piechart of the location with the highest Density of recurrent hits with window type resolution 1 (noncoding vs protein-coding)
# library(dplyr)
# data_location_with_highest_nb_rec_hits_at_resolution_1 <- data.frame(location=factor(unique(df_per_site_type_mutated_windows[,c("non_partitioned_window_id_2","color_edge_nb_rec_hits_resolution_1")])$color_edge_nb_rec_hits_resolution_1,levels = c("noncoding","protein-coding","=")),stringsAsFactors = FALSE)
# df_location_with_highest_nb_rec_hits_at_resolution_1 <- data_location_with_highest_nb_rec_hits_at_resolution_1 %>% 
#   group_by(location) %>% # Variable to be transformed
#   count() %>% 
#   ungroup() %>% 
#   mutate(perc = `n` / sum(`n`)) %>% 
#   arrange(perc) %>%
#   mutate(labels = scales::percent(perc))
# 
# ggplot(df_location_with_highest_nb_rec_hits_at_resolution_1, aes(x = "", y = perc, fill = factor(location,levels=c("noncoding","=","protein-coding")) )) +
#   geom_col(color = "black") +
#   geom_label(aes(label = labels), color = c(1, "white", 1),
#              position = position_stack(vjust = 0.5),
#              show.legend = FALSE) +
#   guides(fill = guide_legend(title = "Proportion of genomic windows where\nthe sites with the highest density\nof recurrent hits are:")) +
#   scale_fill_viridis_d() +
#   coord_polar(theta = "y") + 
#   theme_void()
# ggsave(filename = "Piechart_proportion_windows_in_which_site_type1_with_highest_Density_of_rec_hits_are.svg", path=output_workspace, width = 20, height = 15, units = "cm", device= svg)
# ggsave(filename = "Piechart_proportion_windows_in_which_site_type1_with_highest_Density_of_rec_hits_are.png", path=output_workspace, width = 20, height = 15, units = "cm")
# 
# #Piechart of the regions ranking based on the highest Density of recurrent hits at window type resolution 2
# data_location_with_highest_nb_rec_hits_at_resolution_2 <- data.frame(location=(unique(df_per_site_type_mutated_windows[,c("non_partitioned_window_id_2","color_edge_nb_rec_hits_resolution_2")])$color_edge_nb_rec_hits_resolution_2),stringsAsFactors = FALSE)
# df_location_with_highest_nb_rec_hits_at_resolution_2 <- data_location_with_highest_nb_rec_hits_at_resolution_2 %>% 
#   group_by(location) %>% #Variable to be transformed
#   count() %>% 
#   ungroup() %>% 
#   mutate(perc = `n` / sum(`n`)) %>% 
#   arrange(perc) %>%
#   mutate(labels = scales::percent(perc))
# 
# ggplot(df_location_with_highest_nb_rec_hits_at_resolution_2, aes(x = "", y = perc, fill = location )) +
#   geom_col(color = "black") +
#   geom_label(aes(label = labels), color = 1,
#              position = position_stack(vjust = 0.5),
#              show.legend = FALSE) +
#   guides(fill = guide_legend(title = "Proportion of genomic windows where\nthe site ranking based on the\n density of recurrent hit is:")) +
#   scale_fill_manual(values = palette_window_type_resolution_2[intersect(names(palette_window_type_resolution_2),as.character(df_location_with_highest_nb_rec_hits_at_resolution_2$location))]) +
#   coord_polar(theta = "y") + 
#   theme_void()
# ggsave(filename = "Piechart_proportion_windows_where_site_ranking_based_on_Density_of_rec_hits_are.svg", path=output_workspace, width = 20, height = 15, units = "cm", device= svg)
# 
# #Piechart of the location with the highest recurrent hits average prevalence with window type resolution 1 (noncoding vs protein-coding)
# data_location_with_highest_rec_hits_avg_prevalence_at_resolution_1 <- data.frame(location=factor(unique(df_per_site_type_mutated_windows[,c("non_partitioned_window_id_2","color_edge_avg_prevalence_resolution_1")])$color_edge_avg_prevalence_resolution_1,levels = c("noncoding","protein-coding","=")),stringsAsFactors = FALSE)
# df_location_with_highest_rec_hits_avg_prevalence_at_resolution_1 <- data_location_with_highest_rec_hits_avg_prevalence_at_resolution_1 %>% 
#   group_by(location) %>% # Variable to be transformed
#   count() %>% 
#   ungroup() %>% 
#   mutate(perc = `n` / sum(`n`)) %>% 
#   arrange(perc) %>%
#   mutate(labels = scales::percent(perc))
# 
# ggplot(df_location_with_highest_rec_hits_avg_prevalence_at_resolution_1, aes(x = "", y = perc, fill = factor(location,levels=c("noncoding","=","protein-coding")) )) +
#   geom_col(color = "black") +
#   geom_label(aes(label = labels), color = c(1, "white", 1),
#              position = position_stack(vjust = 0.5),
#              show.legend = FALSE) +
#   guides(fill = guide_legend(title = "Proportion of genomic windows where\nthe sites with the highest average\nrecurrent hits prevalence are:")) +
#   scale_fill_viridis_d() +
#   coord_polar(theta = "y") + 
#   theme_void()
# ggsave(filename = "Piechart_proportion_windows_in_which_site_type1_with_highest_rec_hits_avg_prevalence_are.svg", path=output_workspace, width = 20, height = 15, units = "cm", device= svg)
# ggsave(filename = "Piechart_proportion_windows_in_which_site_type1_with_highest_rec_hits_avg_prevalence_are.png", path=output_workspace, width = 20, height = 15, units = "cm")
# 
# #Piechart of the regions ranking based on the superior mean prevalence at window type resolution 2 
# data_location_with_highest_rec_hits_avg_prevalence_at_resolution_2 <- data.frame(location=(unique(df_per_site_type_mutated_windows[,c("non_partitioned_window_id_2","color_edge_avg_prevalence_resolution_2")])$color_edge_avg_prevalence_resolution_2),stringsAsFactors = FALSE)
# df_location_with_highest_rec_hits_avg_prevalence_at_resolution_2 <- data_location_with_highest_rec_hits_avg_prevalence_at_resolution_2 %>% 
#   group_by(location) %>% #Variable to be transformed
#   count() %>% 
#   ungroup() %>% 
#   mutate(perc = `n` / sum(`n`)) %>% 
#   arrange(perc) %>%
#   mutate(labels = scales::percent(perc))
# 
# ggplot(df_location_with_highest_nb_rec_hits_at_resolution_2, aes(x = "", y = perc, fill = location )) +
#   geom_col(color = "black") +
#   geom_label(aes(label = labels), color = 1,
#              position = position_stack(vjust = 0.5),
#              show.legend = FALSE) +
#   guides(fill = guide_legend(title = "Proportion of genomic windows where\nthe site average hit prevalence ranking is:")) +
#   scale_fill_manual(values = palette_window_type_resolution_2[intersect(names(palette_window_type_resolution_2),as.character(df_location_with_highest_rec_hits_avg_prevalence_at_resolution_2$location))]) +
#   coord_polar(theta = "y") + 
#   theme_void()
# ggsave(filename = "Piechart_proportion_windows_where_site_ranking_based_on_rec_hits_avg_prevalence_are.svg", path=output_workspace, width = 20, height = 15, units = "cm", device= svg)
# 
# #compare the density of parallel hits per window across site types
# ggplot(subset(df_per_site_type_mutated_windows,site_type%in%c("noncoding","protein-coding")),aes(xmin=concatenated_genome_start,xmax=concatenated_genome_end,ymin=0, ymax=density_parallel_hits)) + geom_rect(fill="black") + geom_vline(xintercept = df_end_chrs_in_concat_genome$concat_end, lty = 2, col = "red") + xlab("Locus") + ylab("Density of recurrent hits\n(number of recurrent hits per bp)") + facet_wrap(~site_type,ncol=1) +
#   theme_bw() + theme(axis.text.x=element_blank(),axis.title.x = element_text(size=14),axis.ticks.x = element_blank(),axis.text.y=element_text(size=12),axis.title.y = element_text(size=14))
# ggsave(filename = "Density_of_recurrent_hits_in_each_window_per_site_type.eps", path=output_workspace, width = 20/2.54, height = 12/2.54, device = cairo_ps)
# ggsave(filename = "Density_of_recurrent_hits_in_each_window_per_site_type.png", path=output_workspace, width = 20, height = 12, units = "cm",dpi = 1200)
# #compare the number of parallel hits per window across site types
# ggplot(subset(df_per_site_type_mutated_windows,site_type%in%c("noncoding","protein-coding")),aes(xmin=concatenated_genome_start,xmax=concatenated_genome_end,ymin=0, ymax=nb_recurrent_mutations)) + geom_rect(fill="black") + geom_vline(xintercept = df_end_chrs_in_concat_genome$concat_end, lty = 2, col = "red") + xlab("Locus") + ylab("Number of recurrent hits") + facet_wrap(~site_type,ncol=1) +
#   theme_bw() + theme(axis.text.x=element_blank(),axis.title.x = element_text(size=14),axis.ticks.x = element_blank(),axis.text.y=element_text(size=12),axis.title.y = element_text(size=14))
# ggsave(filename = "Number_of_recurrent_hits_in_each_window_per_site_type.eps", path=output_workspace, width = 20/2.54, height = 12/2.54, device = cairo_ps)
# ggsave(filename = "Number_of_recurrent_hits_in_each_window_per_site_type.png", path=output_workspace, width = 20, height = 12, units = "cm",dpi = 1200)
# #compare the average hits prevalence per window across site types
# ggplot(subset(df_per_site_type_mutated_windows,site_type%in%c("noncoding","protein-coding")),aes(xmin=concatenated_genome_start,xmax=concatenated_genome_end,ymin=0, ymax=average_mutation_recurrence)) + geom_rect(fill="black") + geom_vline(xintercept = df_end_chrs_in_concat_genome$concat_end, lty = 2, col = "red") + xlab("Locus") + ylab("Average hits prevalence") + facet_wrap(~site_type,ncol=1) +
#   theme_bw() + theme(axis.text.x=element_blank(),axis.title.x = element_text(size=14),axis.ticks.x = element_blank(),axis.text.y=element_text(size=12),axis.title.y = element_text(size=14))
# ggsave(filename = "Avg_hits_prevalence_in_each_window_per_site_type.eps", path=output_workspace, width = 20/2.54, height = 12/2.54, device = cairo_ps)
# ggsave(filename = "Avg_hits_prevalence_in_each_window_per_site_type.png", path=output_workspace, width = 20, height = 12, units = "cm",dpi = 1200)
# 
# # #Comparison of recurrent hits density across site types
# # ggplot(data = df_per_site_type_windows,mapping=aes(x=reorder(site_type,log10_density_parallel_hits,function(x) mean(x,na.rm=T) ),y=log10_density_parallel_hits)) + geom_violin(aes(fill=reorder(site_type,-log10_density_parallel_hits,function(x) mean(x,na.rm=T)))) + geom_point() + xlab("Site type") + ylab("Density of recurrent hits\nlog10(Number of recurrent hits/bp + 1e-7)") + scale_fill_manual(values = c("Intron"="#E6AB02","protein-coding"="#666666","UTR"="#66A61E")) + theme_bw() +theme(axis.title = element_text(size=22),axis.text = element_text(size=22),legend.title = element_text(size=12),legend.text = element_text(size=8),axis.text.x = element_text(size=20),legend.position = "none")  + stat_compare_means(method = "wilcox",paired=T, comparisons = list( c("Intron","protein-coding"), c("protein-coding", "UTR"))) + geom_line(mapping=aes(group=non_partitioned_window_id_2)) 
# # ggsave(filename = "Comparison_of_recurrent_hits_density_across_site_types.svg", path=output_workspace, width = 20, height = 15, units = "cm", device= svg)
# # ggsave(filename = "Comparison_of_recurrent_hits_density_across_site_types.png", path=output_workspace, width = 20, height = 15, units = "cm")
# # #Comparison of the number of recurrent hits across site types
# # ggplot(data = df_per_site_type_windows,mapping=aes(x=reorder(site_type,nb_recurrent_mutations),y=nb_recurrent_mutations)) + geom_violin(aes(fill=reorder(site_type,nb_recurrent_mutations))) + geom_point() + xlab("Site type") + ylab("Number of recurrent hits") + scale_fill_manual(values = c("Intron"="#E6AB02","protein-coding"="#666666","UTR"="#66A61E")) + theme_bw() +theme(axis.title = element_text(size=22),axis.text = element_text(size=22),legend.title = element_text(size=12),legend.text = element_text(size=8),axis.text.x = element_text(size=20),legend.position = "none")  + stat_compare_means(method = "wilcox",paired=T, comparisons = list( c("Intron","protein-coding"), c("protein-coding", "UTR"))) + geom_line(mapping=aes(group=non_partitioned_window_id_2)) 
# # ggsave(filename = "Comparison_of_nb_recurrent_hits_across_site_types.svg", path=output_workspace, width = 20, height = 15, units = "cm", device= svg)
# # ggsave(filename = "Comparison_of_nb_recurrent_hits_density_site_types.png", path=output_workspace, width = 20, height = 15, units = "cm")
# 
# #Get the list of accession ID for the genes of the introns and UTR parallel hits
# v_lst_gene_accession_for_intron_and_UTR_recurrent_hits <- sort(unique(subset(df_recurrence_ALL_variants,recurrence>1&(is_in_gene_UTR|is_in_gene_intron)&(!is.na(accession_number)))$accession_number))
# write.table(x=v_lst_gene_accession_for_intron_and_UTR_recurrent_hits,file = paste0(output_workspace,"Table_gene_accession_for_intron_and_UTR_recurrent_hits.txt"),sep = "\t",na = "NA",row.names = FALSE,col.names = FALSE)
# 
# #Get the list of accession ID for the Other noncoding regions exhibiting recurrent hits
# v_lst_gene_accession_for_Other_noncodingregions_recurrent_hits <- sort(unique(subset(df_recurrence_ALL_variants,recurrence>1&(!(is_in_prot_gene_CDS|is_in_gene_UTR|is_in_gene_intron))&(!is.na(accession_number)))$accession_number))
# write.table(x=v_lst_gene_accession_for_intron_and_UTR_recurrent_hits,file = paste0(output_workspace,"Table_accession_for_Other_noncodingregions_exhibiting_recurrent_hits.txt"),sep = "\t",na = "NA",row.names = FALSE,col.names = FALSE)
# 
# #Get the list of accession ID for the genes of the introns and UTR parallel hits
# v_lst_gene_accession_for_intron_and_UTR_top_recurrent_hits <- sort(unique(subset(df_recurrence_ALL_variants,recurrence>4&(is_in_gene_UTR|is_in_gene_intron)&(!is.na(accession_number)))$accession_number))
# write.table(x=v_lst_gene_accession_for_intron_and_UTR_recurrent_hits,file = paste0(output_workspace,"Table_gene_accession_for_intron_and_UTR_top_recurrent_hits.txt"),sep = "\t",na = "NA",row.names = FALSE,col.names = FALSE)
# 
# #Get the list of accession ID for the Other noncoding regions exhibiting recurrent hits
# v_lst_gene_accession_for_Other_noncodingregions_top_recurrent_hits <- sort(unique(subset(df_recurrence_ALL_variants,recurrence>4&(!(is_in_prot_gene_CDS|is_in_gene_UTR|is_in_gene_intron))&(!is.na(accession_number)))$accession_number))
# write.table(x=v_lst_gene_accession_for_intron_and_UTR_recurrent_hits,file = paste0(output_workspace,"Table_accession_for_Other_noncodingregions_exhibiting_top_recurrent_hits.txt"),sep = "\t",na = "NA",row.names = FALSE,col.names = FALSE)

# #Per-chromosome coverage analysis to detect potential duplication events
# df_per_chromosome_cov_stats_across_samples <- read.csv2(file = paste0(output_workspace,"per_chromosome_coverage_stats.csv"),sep = "\t",header = FALSE,stringsAsFactors = FALSE)
# names(df_per_chromosome_cov_stats_across_samples) <- c("Chrom","nb_reads","Sample")
# df_per_chromosome_cov_stats_across_samples$nb_reads_per_bp <- df_per_chromosome_cov_stats_across_samples$nb_reads/v_chrs_length[df_per_chromosome_cov_stats_across_samples$Chrom]
# ggplot(data = df_per_chromosome_cov_stats_across_samples,mapping=aes(x=reorder(Chrom,nb_reads_per_bp,function(x) mean(x,na.rm=T) ),y=nb_reads_per_bp)) + geom_violin(aes(fill=reorder(Chrom,nb_reads_per_bp,function(x) mean(x,na.rm=T) ))) + geom_point() + xlab("Chromosome") + ylab("Number of reads per bp") + scale_fill_viridis_d() + theme_bw() +theme(axis.title = element_text(size=22),axis.text = element_text(size=22),legend.title = element_text(size=12),legend.text = element_text(size=8),axis.text.x = element_text(size=20),legend.position = "none")  + stat_compare_means(mapping = aes(group=Sample), method = "wilcox", comparisons = list( c("I","II"), c("I", "III"), c("II","III") )) + geom_line(mapping=aes(group=Sample)) 
# ggsave(filename = "Per-chromosome_coverage_analysis_to_detect_potential_duplication_events.svg", path=output_workspace, width = 20, height = 15, units = "cm", device= svg)

#Compared the observed distributions of the number of populations hit (for all genes AND multi-hit genes) to the expected null model
df_observed_and_expected_nb_pops_hit_per_gene <- NULL
for (current_gene in names(v_genes_length)){
  #Number of expected PH for current gene 
  nb_exp_PH_for_current_gene <- unname((2E-10)*v_genes_length[current_gene]*10000*10)
  #Number of observed PH for current gene 
  nb_obs_PH_for_current_gene <- length(unique(subset(df_ALL_variants,(is_fixed)&(variant_type2!="synonymous")&(!is.na(variant_type2))&(gene_product==current_gene))$population))
  df_observed_and_expected_nb_pops_hit_per_gene <- rbind(df_observed_and_expected_nb_pops_hit_per_gene, data.frame(gene=current_gene,PH=nb_exp_PH_for_current_gene,Dataset="Expected null model\n(All genes)",stringsAsFactors = T))
  df_observed_and_expected_nb_pops_hit_per_gene <- rbind(df_observed_and_expected_nb_pops_hit_per_gene, data.frame(gene=current_gene,PH=nb_obs_PH_for_current_gene,Dataset="Observed (All genes)",stringsAsFactors = T))
  if (current_gene%in%v_lst_nonRecuHit_genes_hit_in_multiple_pops){
    #Number of observed PH for current gene 
    nb_obs_PH_for_current_multihit_gene <- length(unique(subset(df_ALL_variants,(is_fixed)&(variant_type2!="synonymous")&(!is.na(variant_type2))&(gene_product==current_gene))$population))
    df_observed_and_expected_nb_pops_hit_per_gene <- rbind(df_observed_and_expected_nb_pops_hit_per_gene, data.frame(gene=current_gene,PH=nb_obs_PH_for_current_multihit_gene,Dataset="Observed (Multi-hit genes)",stringsAsFactors = T))
  }
}
#Compare the observed distribution of PH for all genes and multihit genes to the expected distribution for genes with PH>=1
v_lst_genes_with_at_least_one_PH <- sort(unique(subset(df_observed_and_expected_nb_pops_hit_per_gene,PH>=1&Dataset=="Observed (All genes)")$gene))
  #perform a permutation test to get the expected distribution of genes with PH>=1 (average PH over 999 samples of 6 genes)
nb_permutations <- 999
nb_sampled_genes_per_permutation <- 6 #same size as the number of multihit genes

for (i_perm in 1:nb_permutations){
  v_current_perm_sampled_genes <- as.character(sample(x = v_lst_genes_with_at_least_one_PH,size = nb_sampled_genes_per_permutation,replace = F))
  df_observed_and_expected_nb_pops_hit_per_gene <- rbind(df_observed_and_expected_nb_pops_hit_per_gene, data.frame(gene="Average over 6 random\ngenes with PH>=1",PH=mean(subset(df_observed_and_expected_nb_pops_hit_per_gene,PH>=1&Dataset=="Observed (All genes)"&gene%in%v_current_perm_sampled_genes)$PH,na.rm=T),Dataset="Expected null model\n(Genes with PH>=1)",stringsAsFactors = T))
}

ggplot(data = subset(df_observed_and_expected_nb_pops_hit_per_gene,Dataset!="Expected null model\n(All genes)"),mapping=aes(y=PH, x=reorder(Dataset,-PH,median))) + geom_violin(fill="grey",width =1.2) + geom_boxplot(fill="white",width =0.02) + xlab("Gene sets") + ylab("Number of population hits (PH)") + theme_bw() +theme(axis.title = element_text(size=16),axis.text = element_text(size=16),legend.title = element_text(size=12),legend.text = element_text(size=10),axis.text.x = element_text(size=14),legend.position = "none") + stat_compare_means(method = "wilcox.test", paired = F,comparisons = list( c("Expected null model\n(Genes with PH>=1)", "Observed (All genes)"),c("Observed (Multi-hit genes)", "Expected null model\n(Genes with PH>=1)")))
ggsave(filename = "Expected_vs_Observed_PHs.svg", path=output_workspace, width = 20, height = 15, units = "cm", device= svg)
ggsave(filename = "Expected_vs_Observed_PHs.jpg", path=output_workspace, width = 20, height = 15, units = "cm", device= jpeg)
#K-S test expected vs observed PH distributions
  #expected null vs observed all genes
ks.test(subset(df_observed_and_expected_nb_pops_hit_per_gene,Dataset=="Expected null model\n(Genes with PH>=1)")$PH,subset(df_observed_and_expected_nb_pops_hit_per_gene,Dataset=="Observed (All genes)")$PH)$p.value
  #expected null vs observed multi-genes
ks.test(subset(df_observed_and_expected_nb_pops_hit_per_gene,Dataset=="Expected null model\n(Genes with PH>=1)")$PH,subset(df_observed_and_expected_nb_pops_hit_per_gene,Dataset=="Observed (Multi-hit genes)")$PH)$p.value
  #observed multi-genes vs observed all genes
ks.test(subset(df_observed_and_expected_nb_pops_hit_per_gene,Dataset=="Observed (Multi-hit genes)")$PH,subset(df_observed_and_expected_nb_pops_hit_per_gene,Dataset=="Observed (All genes)")$PH)$p.value

#PYK1 data for mapping figure
v_pyk1_fixed_non_rec_mutations_annotation <- sort(unique(subset(df_variants_time_series,gene_product=="pyk1"&VarFreq>=0.75)$Position))
names(v_pyk1_fixed_non_rec_mutations_annotation) <- c("3' UTR", "3' UTR", "3' UTR",  "3' UTR","missense","5' UTR")

###LOOK FOR LOF VARIANTS IN MULTIHIT GENES
  #investigate the insertion in bmc1 (Warning: the gene has an intron)
bmc1_indel_chr <- "II"
bmc1_indel_pos_in_gene <- (2965030-2964878+1)+(2964835-2964227+1) #adding first CDS length
bmc1_WT_dna_seq <- paste0(substr(genome_refseq[[bmc1_indel_chr]],2964227,2964835)[1],substr(genome_refseq[[bmc1_indel_chr]],2964878,2965075)[1])
bmc1_WT_prot_seq <- paste0(seqinr::translate(strsplit(bmc1_WT_dna_seq,split = "")[[1]]),collapse = "")
print(bmc1_WT_prot_seq)
bmc1_dna_seq_AFTER_INDEL <- paste0(substr(genome_refseq[[bmc1_indel_chr]],2964227,2964835)[1],substr(genome_refseq[[bmc1_indel_chr]],2964878,2965029)[1],"A",substr(genome_refseq[[bmc1_indel_chr]],2965030,2965075)[1],"--")
bmc1_prot_seq_AFTER_INDEL <- paste0(seqinr::translate(strsplit(bmc1_dna_seq_AFTER_INDEL,split = "")[[1]]),collapse = "")
print(bmc1_prot_seq_AFTER_INDEL)
#print hypothesis about the indel
if (gregexpr(pattern = "*",text = bmc1_prot_seq_AFTER_INDEL,fixed = T)[[1]][1]!=nchar(bmc1_prot_seq_AFTER_INDEL)){
  print(paste0("This indel in bmc1 causes a premature stop codon and truncates ",nchar(bmc1_prot_seq_AFTER_INDEL)-(gregexpr(pattern = "*",text = bmc1_prot_seq_AFTER_INDEL,fixed = T)[[1]][1])," out of ",nchar(bmc1_WT_prot_seq)," amino acids from the protein! This affects all the domains of the protein"))
}
if (length(gregexpr(pattern = "*",text = bmc1_prot_seq_AFTER_INDEL,fixed = T)[[1]])>0){
  print("This indel in bmc1 probably causes a loss of function!")
}

  #investigate the insertion in bop1 (the gene is on the - strand)
bop1_indel_chr <- "I"
bop1_indel_pos_in_gene <- abs(3408209-3409267)+1
bop1_WT_dna_seq <- seqinr::comp(strsplit(paste0(stringi::stri_reverse(substr(genome_refseq[[bop1_indel_chr]],3409427,3409525)[1]),stringi::stri_reverse(substr(genome_refseq[[bop1_indel_chr]],3407825,3409267)[1])),split="" )[[1]])
bop1_WT_prot_seq <- paste0(seqinr::translate(bop1_WT_dna_seq),collapse = "")
print(bop1_WT_prot_seq)
bop1_dna_seq_AFTER_INDEL <- seqinr::comp(strsplit(paste0(stringi::stri_reverse(substr(genome_refseq[[bop1_indel_chr]],3409427,3409525)[1]),stringi::stri_reverse(substr(genome_refseq[[bop1_indel_chr]],3407825,3408209)[1]),"G",stringi::stri_reverse(substr(genome_refseq[[bop1_indel_chr]],3408209,3409267)[1])),split="" )[[1]])
bop1_prot_seq_AFTER_INDEL <- paste0(seqinr::translate(bop1_dna_seq_AFTER_INDEL),collapse = "")
print(bop1_prot_seq_AFTER_INDEL)
#print hypothesis about the indel
if (gregexpr(pattern = "*",text = bop1_prot_seq_AFTER_INDEL,fixed = T)[[1]][1]!=nchar(bop1_prot_seq_AFTER_INDEL)){
  print(paste0("This indel in bop1 causes a premature stop codon and truncates ",nchar(bop1_prot_seq_AFTER_INDEL)-(gregexpr(pattern = "*",text = bop1_prot_seq_AFTER_INDEL,fixed = T)[[1]][1])," out of ",nchar(bop1_WT_prot_seq)," amino acids from the protein! This affects all the domains of the protein"))
}
if (length(gregexpr(pattern = "*",text = bop1_prot_seq_AFTER_INDEL,fixed = T)[[1]])>0){
  print("This indel in bop1 probably causes a loss of function!")
}

  #investigate the nonsense mutation in cmr2
cmr2_nonsense_mutation_chr <- "I"
start_current_gene_CDSs <- 1124167
end_current_gene_CDSs <- 1128720
cmr2_nonsense_mutation_pos_in_gene <- 1128193-start_current_gene_CDSs+1
cmr2_WT_dna_seq <- substr(genome_refseq[[cmr2_nonsense_mutation_chr]],start_current_gene_CDSs,end_current_gene_CDSs)[1]
cmr2_WT_prot_seq <- paste0(seqinr::translate(strsplit(cmr2_WT_dna_seq,split = "")[[1]]),collapse = "")
print(cmr2_WT_prot_seq)
cmr2_dna_seq_AFTER_NONSENSE_MUTATION <- paste0(substr(cmr2_WT_dna_seq,1,cmr2_nonsense_mutation_pos_in_gene-1),"T",substr(cmr2_WT_dna_seq,cmr2_nonsense_mutation_pos_in_gene+1,nchar(cmr2_WT_dna_seq)))
cmr2_prot_seq_AFTER_NONSENSE_MUTATION <- paste0(seqinr::translate(strsplit(cmr2_dna_seq_AFTER_NONSENSE_MUTATION,split = "")[[1]]),collapse = "")
print(cmr2_prot_seq_AFTER_NONSENSE_MUTATION)
#print the truncation size of the nonsense mutation and characterize the trunctation
if (gregexpr(pattern = "*",text = cmr2_prot_seq_AFTER_NONSENSE_MUTATION,fixed = T)[[1]][1]!=nchar(cmr2_prot_seq_AFTER_NONSENSE_MUTATION)){
  print(paste0("This nonsense mutation in cmr2 causes a truncation of ",nchar(cmr2_prot_seq_AFTER_NONSENSE_MUTATION)-(gregexpr(pattern = "*",text = cmr2_prot_seq_AFTER_NONSENSE_MUTATION,fixed = T)[[1]][1])," amino acids that affects the Response regulator receiver domain and disordered regions in C-terminal!"))
}

  #investigate the insertion in pfl5
pfl5_indel_chr <- "II"
pfl5_indel_pos_in_gene <- 4423656-4421413+1
pfl5_WT_dna_seq <- paste0(substr(genome_refseq[[pfl5_indel_chr]],4421413,4423655)[1],substr(genome_refseq[[pfl5_indel_chr]],4423656,4425264)[1])
pfl5_WT_prot_seq <- paste0(seqinr::translate(strsplit(pfl5_WT_dna_seq,split = "")[[1]]),collapse = "")
print(pfl5_WT_prot_seq)
pfl5_dna_seq_AFTER_INDEL <- paste0(substr(genome_refseq[[pfl5_indel_chr]],4421413,4423655)[1],"AAACT",substr(genome_refseq[[pfl5_indel_chr]],4423656,4425264)[1])
pfl5_prot_seq_AFTER_INDEL <- paste0(seqinr::translate(strsplit(pfl5_dna_seq_AFTER_INDEL,split = "")[[1]]),collapse = "")
print(pfl5_prot_seq_AFTER_INDEL)
#print hypothesis about the indel
if (gregexpr(pattern = "*",text = pfl5_prot_seq_AFTER_INDEL,fixed = T)[[1]][1]!=nchar(pfl5_prot_seq_AFTER_INDEL)){
  print(paste0("This indel in pfl5 causes a premature stop codon and truncates ",nchar(pfl5_prot_seq_AFTER_INDEL)-(gregexpr(pattern = "*",text = pfl5_prot_seq_AFTER_INDEL,fixed = T)[[1]][1])," out of ",nchar(pfl5_WT_prot_seq)," amino acids from the protein! This affects all the domains of the protein"))
}
if (length(gregexpr(pattern = "*",text = pfl5_prot_seq_AFTER_INDEL,fixed = T)[[1]])>0){
  print("This indel in pfl5 probably causes a loss of function!")
}

  #investigate the deletion in pfl5
pfl5_indel_chr <- "II"
pfl5_indel_pos_in_gene <- 4423650-4421413+1
pfl5_WT_dna_seq <- paste0(substr(genome_refseq[[pfl5_indel_chr]],4421413,4423649)[1],substr(genome_refseq[[pfl5_indel_chr]],4423650,4425264)[1])
pfl5_WT_prot_seq <- paste0(seqinr::translate(strsplit(pfl5_WT_dna_seq,split = "")[[1]]),collapse = "")
print(pfl5_WT_prot_seq)
pfl5_dna_seq_AFTER_INDEL <- paste0(substr(genome_refseq[[pfl5_indel_chr]],4421413,4423649)[1],substr(genome_refseq[[pfl5_indel_chr]],4423655,4425264)[1])
pfl5_prot_seq_AFTER_INDEL <- paste0(seqinr::translate(strsplit(pfl5_dna_seq_AFTER_INDEL,split = "")[[1]]),collapse = "")
print(pfl5_prot_seq_AFTER_INDEL)
#print hypothesis about the indel
if (gregexpr(pattern = "*",text = pfl5_prot_seq_AFTER_INDEL,fixed = T)[[1]][1]!=nchar(pfl5_prot_seq_AFTER_INDEL)){
  print(paste0("This indel in pfl5 causes a premature stop codon and truncates ",nchar(pfl5_prot_seq_AFTER_INDEL)-(gregexpr(pattern = "*",text = pfl5_prot_seq_AFTER_INDEL,fixed = T)[[1]][1])," out of ",nchar(pfl5_WT_prot_seq)," amino acids from the protein! This affects all the domains of the protein"))
}
if (length(gregexpr(pattern = "*",text = pfl5_prot_seq_AFTER_INDEL,fixed = T)[[1]])>0){
  print("This indel in pfl5 probably causes a loss of function!")
}

  #investigate the deletion in pof15
pof15_indel_chr <- "I"
pof15_indel_pos_in_gene <- 1891186-1890920+1
pof15_WT_dna_seq <- paste0(substr(genome_refseq[[pof15_indel_chr]],1890920,1891185)[1],substr(genome_refseq[[pof15_indel_chr]],1891186,1891651)[1])
pof15_WT_prot_seq <- paste0(seqinr::translate(strsplit(pof15_WT_dna_seq,split = "")[[1]]),collapse = "")
print(pof15_WT_prot_seq)
pof15_dna_seq_AFTER_INDEL <- paste0(substr(genome_refseq[[pof15_indel_chr]],1890920,1891185)[1],substr(genome_refseq[[pof15_indel_chr]],1891194,1891651)[1])
pof15_prot_seq_AFTER_INDEL <- paste0(seqinr::translate(strsplit(pof15_dna_seq_AFTER_INDEL,split = "")[[1]]),collapse = "")
print(pof15_prot_seq_AFTER_INDEL)
#print hypothesis about the indel
if (gregexpr(pattern = "*",text = pof15_prot_seq_AFTER_INDEL,fixed = T)[[1]][1]!=nchar(pof15_prot_seq_AFTER_INDEL)){
  print(paste0("This indel in pof15 causes a premature stop codon and truncates ",nchar(pof15_prot_seq_AFTER_INDEL)-(gregexpr(pattern = "*",text = pof15_prot_seq_AFTER_INDEL,fixed = T)[[1]][1])," out of ",nchar(pof15_WT_prot_seq)," amino acids from the protein! This affects all the domains of the protein"))
}
if (length(gregexpr(pattern = "*",text = pof15_prot_seq_AFTER_INDEL,fixed = T)[[1]])>0){
  print("This indel in pof15 probably causes a loss of function!")
}

  #investigate the deletion in pzh1 (the gene is on the - strand)
pzh1_indel_chr <- "I"
pzh1_WT_dna_seq <- seqinr::comp(strsplit(paste0(stringi::stri_reverse(substr(genome_refseq[[pzh1_indel_chr]],1529295,1529511)[1]),stringi::stri_reverse(substr(genome_refseq[[pzh1_indel_chr]],1527964,1529294)[1])),split="" )[[1]])
pzh1_WT_prot_seq <- paste0(seqinr::translate(pzh1_WT_dna_seq),collapse = "")
print(pzh1_WT_prot_seq)
pzh1_dna_seq_AFTER_INDEL <- seqinr::comp(strsplit(paste0(stringi::stri_reverse(substr(genome_refseq[[pzh1_indel_chr]],1529301,1529511)[1]),stringi::stri_reverse(substr(genome_refseq[[pzh1_indel_chr]],1527964,1529294)[1])),split="" )[[1]])
pzh1_prot_seq_AFTER_INDEL <- paste0(seqinr::translate(pzh1_dna_seq_AFTER_INDEL),collapse = "")
print(pzh1_prot_seq_AFTER_INDEL)
#print hypothesis about the indel
if (gregexpr(pattern = "*",text = pzh1_prot_seq_AFTER_INDEL,fixed = T)[[1]][1]!=nchar(pzh1_prot_seq_AFTER_INDEL)){
  print(paste0("This indel in pzh1 causes a premature stop codon and truncates ",nchar(pzh1_prot_seq_AFTER_INDEL)-(gregexpr(pattern = "*",text = pzh1_prot_seq_AFTER_INDEL,fixed = T)[[1]][1])," out of ",nchar(pzh1_WT_prot_seq)," amino acids from the protein! This affects all the domains of the protein"))
}
if (length(gregexpr(pattern = "*",text = pzh1_prot_seq_AFTER_INDEL,fixed = T)[[1]])>1){
  print("This indel in pzh1 probably causes a loss of function!")
}


  #investigate the insertion in prr1
prr1_indel_chr <- "I"
prr1_indel_pos_in_gene <- 3666021-3666007+1
prr1_WT_dna_seq <- paste0(substr(genome_refseq[[prr1_indel_chr]],3666007,3666036)[1],substr(genome_refseq[[prr1_indel_chr]],3666134,3666223)[1],substr(genome_refseq[[prr1_indel_chr]],3666980,3667059)[1],substr(genome_refseq[[prr1_indel_chr]],3667113,3667170)[1],substr(genome_refseq[[prr1_indel_chr]],3667314,3668675)[1])
prr1_WT_prot_seq <- paste0(seqinr::translate(strsplit(prr1_WT_dna_seq,split = "")[[1]]),collapse = "")
print(prr1_WT_prot_seq)
prr1_dna_seq_AFTER_INDEL <- paste0(substr(prr1_WT_dna_seq,1,prr1_indel_pos_in_gene),"G",substr(prr1_WT_dna_seq,prr1_indel_pos_in_gene+1,nchar(prr1_WT_dna_seq)))
prr1_prot_seq_AFTER_INDEL <- paste0(seqinr::translate(strsplit(prr1_dna_seq_AFTER_INDEL,split = "")[[1]]),collapse = "")
print(prr1_prot_seq_AFTER_INDEL)
#print hypothesis about the indel
if (gregexpr(pattern = "*",text = prr1_prot_seq_AFTER_INDEL,fixed = T)[[1]][1]!=nchar(prr1_prot_seq_AFTER_INDEL)){
  print(paste0("This indel in prr1 causes a premature stop codon and truncates ",nchar(prr1_prot_seq_AFTER_INDEL)-(gregexpr(pattern = "*",text = prr1_prot_seq_AFTER_INDEL,fixed = T)[[1]][1])," out of ",nchar(prr1_WT_prot_seq)," amino acids from the protein! This affects all the domains of the protein"))
}
if (length(gregexpr(pattern = "*",text = prr1_prot_seq_AFTER_INDEL,fixed = T)[[1]])>0){
  print("This indel in prr1 probably causes a loss of function!")
}

  #investigate the insertion in ubp9
ubp9_indel_chr <- "II"
ubp9_indel_pos_in_gene <- 2937277-2936510+1
ubp9_WT_dna_seq <- paste0(substr(genome_refseq[[ubp9_indel_chr]],2936211,2936237)[1],substr(genome_refseq[[ubp9_indel_chr]],2936282, 2936389)[1],substr(genome_refseq[[ubp9_indel_chr]],2936439, 2936452)[1],substr(genome_refseq[[ubp9_indel_chr]],2936510, 2937276)[1],substr(genome_refseq[[ubp9_indel_chr]],2937277, 2938118)[1])
ubp9_WT_prot_seq <- paste0(seqinr::translate(strsplit(ubp9_WT_dna_seq,split = "")[[1]]),collapse = "")
print(ubp9_WT_prot_seq)
ubp9_dna_seq_AFTER_INDEL <- paste0(substr(genome_refseq[[ubp9_indel_chr]],2936211,2936237)[1],substr(genome_refseq[[ubp9_indel_chr]],2936282, 2936389)[1],substr(genome_refseq[[ubp9_indel_chr]],2936439, 2936452)[1],substr(genome_refseq[[ubp9_indel_chr]],2936510, 2937276)[1],"A",substr(genome_refseq[[ubp9_indel_chr]],2937277, 2938118)[1])
ubp9_prot_seq_AFTER_INDEL <- paste0(seqinr::translate(strsplit(ubp9_dna_seq_AFTER_INDEL,split = "")[[1]]),collapse = "")
print(ubp9_prot_seq_AFTER_INDEL)
#print hypothesis about the indel
if (gregexpr(pattern = "*",text = ubp9_prot_seq_AFTER_INDEL,fixed = T)[[1]][1]!=nchar(ubp9_prot_seq_AFTER_INDEL)){
  print(paste0("This indel in ubp9 causes a premature stop codon and truncates ",nchar(ubp9_prot_seq_AFTER_INDEL)-(gregexpr(pattern = "*",text = ubp9_prot_seq_AFTER_INDEL,fixed = T)[[1]][1])," out of ",nchar(ubp9_WT_prot_seq)," amino acids from the protein! This affects all the domains of the protein"))
}
if (length(gregexpr(pattern = "*",text = ubp9_prot_seq_AFTER_INDEL,fixed = T)[[1]])>0){
  print("This indel in ubp9 probably causes a loss of function!")
}

#Assessing the selective pressure on indels (FIXED in-frame vs frameshift in coding and non-coding regions)
backup_df_indel_variants <- df_indel_variants
df_indel_variants <- subset(df_ALL_variants,variant_type%in% c("Insertion","Deletion"))
df_indel_variants$length <- nchar(df_indel_variants$VarAllele) - 1
df_indel_variants$is_in_frame <- df_indel_variants$length%%3==0
df_indel_variants$is_coding <- !(is.na(df_indel_variants$variant_type2)|df_indel_variants$variant_type2!="Indel in protein-coding genes")
df_prevalence_indel_types <- df_indel_variants %>%
  group_by(is_coding,length) %>% #group_by(is_coding,population,length)
  dplyr::summarise(prevalence=length(unique(mutation))) %>%
  as.data.frame()
max_length_indels <- max(df_prevalence_indel_types$length,na.rm=TRUE)
  #add prevalence zero for non-observed indel types within the range of length observed
for (current_is_coding in c(TRUE,FALSE)){
  for (current_pop in v_pop_order){
    for (current_length in 1:max_length_indels){
      is_indel_features_comb_observed <- nrow(subset(df_prevalence_indel_types,(is_coding==current_is_coding)&(length==current_length)))>=1 #&(population==current_pop)
      if (!is_indel_features_comb_observed){
        df_prevalence_indel_types <- rbind(df_prevalence_indel_types,data.frame(is_coding=current_is_coding,length=current_length,prevalence=0)) #,population=current_pop
      }
    }
  }
}
df_prevalence_indel_types$is_in_frame <- df_prevalence_indel_types$length%%3==0
df_prevalence_indel_types$total_sample_size_for_coding_type <- vapply(X=1:nrow(df_prevalence_indel_types),FUN = function(the_i) ifelse(test = df_prevalence_indel_types$is_coding[the_i],yes = length(unique(subset(df_indel_variants, is_coding)$mutation)),no = length(unique(subset(df_indel_variants,(!is_coding))$mutation))), FUN.VALUE = 0)
df_prevalence_indel_types$target_size_for_coding_type <- vapply(X=1:nrow(df_prevalence_indel_types),FUN = function(the_i) ifelse(test = df_prevalence_indel_types$is_coding[the_i],yes = 12.6*1000*0.6,no = 12.6*1000*0.4), FUN.VALUE = 0)
#df_prevalence_indel_types$total_sample_size_for_coding_type_in_pop <- vapply(X=1:nrow(df_prevalence_indel_types),FUN = function(the_i) ifelse(test = df_prevalence_indel_types$is_coding[the_i],yes = length(unique(subset(df_indel_variants,(population == df_prevalence_indel_types$population[the_i]) & is_coding)$mutation)),no = length(unique(subset(df_indel_variants,(population == df_prevalence_indel_types$population[the_i]) & (!is_coding))$mutation))), FUN.VALUE = 0)
df_prevalence_indel_types$proportion <- df_prevalence_indel_types$prevalence/df_prevalence_indel_types$total_sample_size_for_coding_type #df_prevalence_indel_types$total_sample_size_for_coding_type_in_pop
df_prevalence_indel_types$indel_rate <- df_prevalence_indel_types$prevalence/(10*df_prevalence_indel_types$target_size_for_coding_type)
# #get mean and sd by (is_coding,is_in_frame,length) across populations
# df_summary_stas_prevalence_by_indel_types <- df_prevalence_indel_types %>%
#   group_by(is_coding,is_in_frame,length) %>%
#   dplyr::summarise(mean_proportion = mean(proportion,na.rm = T), sd_proportion = sd(proportion,na.rm=T), mean_prevalence = mean(prevalence,na.rm = T)) %>%
#   as.data.frame()

df_prevalence_indel_types$sample_size <- df_prevalence_indel_types$prevalence
#Purifying selection against frameshift indels is specific to coding regions
ggplot(data = df_prevalence_indel_types,mapping=aes(x=length, y = proportion, col=is_coding)) + geom_line() + geom_point(mapping = aes(shape=is_in_frame, size = sample_size ))  + xlab("Indel length") + ylab("Proportion") + scale_color_brewer(palette = "Set1") + scale_x_continuous(breaks = seq(0,max_length_indels+3,3),limits = c(0,max_length_indels+3)) + theme_bw() +theme(axis.title = element_text(size=16),axis.text = element_text(size=16),legend.title = element_text(size=12),legend.text = element_text(size=10),axis.text.x = element_text(size=14),legend.position = "right")
ggsave(filename = "Purifying_selection_against_frameshift_indels_in_coding_regions.svg", path=output_workspace, width = 15, height = 10, units = "cm", device= svg)
ggsave(filename = "Purifying_selection_against_frameshift_indels_in_coding_regions.jpg", path=output_workspace, width = 15, height = 10, units = "cm", device= jpeg)
ggplot(data = df_prevalence_indel_types,mapping=aes(x=length, y = log10(proportion), col=is_coding)) + geom_line() + geom_point(mapping = aes(shape=is_in_frame, size = sample_size ))  + xlab("Indel length") + ylab("Log10(Proportion)") + scale_color_brewer(palette = "Set1") + scale_x_continuous(breaks = seq(0,max_length_indels+3,3),limits = c(0,max_length_indels+3)) + theme_bw() +theme(axis.title = element_text(size=16),axis.text = element_text(size=16),legend.title = element_text(size=12),legend.text = element_text(size=10),axis.text.x = element_text(size=14),legend.position = "right")
ggsave(filename = "Log10_scale_Purifying_selection_against_frameshift_indels_in_coding_regions.svg", path=output_workspace, width = 15, height = 10, units = "cm", device= svg)
ggsave(filename = "Log10_scale_Purifying_selection_against_frameshift_indels_in_coding_regions.jpg", path=output_workspace, width = 15, height = 10, units = "cm", device= jpeg)

#Indel rate by coding type and length
ggplot(data = df_prevalence_indel_types,mapping=aes(x=length, y = log10(indel_rate), col=is_coding)) + geom_line() + geom_point(mapping = aes(shape=is_in_frame, size = sample_size ))  + xlab("Indel length") + ylab("Log10(Indel rate (per kb))") + scale_color_brewer(palette = "Set1") + scale_x_continuous(breaks = seq(0,max_length_indels+3,3),limits = c(0,max_length_indels+3)) + theme_bw() +theme(axis.title = element_text(size=16),axis.text = element_text(size=16),legend.title = element_text(size=12),legend.text = element_text(size=10),axis.text.x = element_text(size=14),legend.position = "right")
ggsave(filename = "Log10_Indel_rate_vs_indel_coding_type_and_length.svg", path=output_workspace, width = 15, height = 10, units = "cm", device= svg)
ggsave(filename = "Log10_Indel_rate_vs_indel_coding_type_and_length.jpg", path=output_workspace, width = 15, height = 10, units = "cm", device= jpeg)

#make sure that the indels data make sense. The figure above should show that longer indels should be
#rarer (except if sample size is small) and the next figure should show that coding regions indels at
#the beginning of the protein sequence are less frequent
df_loci_annotations$length <- abs(df_loci_annotations$end - df_loci_annotations$start)+1
df_indel_variants$pct_pos_in_prot <- vapply(X = 1:nrow(df_indel_variants),FUN = function(the_indx_row) ifelse(test = df_indel_variants$is_coding[the_indx_row],yes = min(floor((abs(df_indel_variants$Position[the_indx_row] - subset(df_loci_annotations,(chr==df_indel_variants$Chrom[the_indx_row])&((df_indel_variants$Position[the_indx_row] >= start) & (df_indel_variants$Position[the_indx_row] <= end) )&(!grepl(pattern = "RNA",x = toupper(df_loci_annotations$annotation),fixed = T))&site_type=="CDS")$start)+1)/3)+1,na.rm=T)/sum(subset(df_loci_annotations,(chr==df_indel_variants$Chrom[the_indx_row])&(df_indel_variants$Position[the_indx_row]>= start & df_indel_variants$Position[the_indx_row] <= end)&(!grepl(pattern = "RNA",x = toupper(df_loci_annotations$annotation),fixed = T))&site_type=="CDS")$length/3,na.rm=T) ,no = NA),FUN.VALUE = 0.0)
#vapply(X = 1:nrow(df_indel_variants),FUN = function(the_indx_row) ifelse(test = df_indel_variants$is_coding[the_indx_row]&(!is.na(strsplit(df_indel_variants$accession_number[the_indx_row],split = ";")[[1]][1])),yes = min(floor((abs(df_indel_variants$Position[the_indx_row] - subset(df_loci_annotations,grepl(pattern = strsplit(df_indel_variants$accession_number[the_indx_row],split = ";")[[1]][1],x = df_loci_annotations$annotation,fixed = T)&site_type=="CDS")$start)+1)/3)+1,na.rm=T)/sum(subset(df_loci_annotations,grepl(pattern = strsplit(df_indel_variants$accession_number[the_indx_row],split = ";")[[1]][1],x = df_loci_annotations$annotation,fixed = T)&site_type=="CDS")$length/3,na.rm=T) ,no = NA),FUN.VALUE = 0.0)
# df_indel_variants$is_coding <- ifelse(test = (df_indel_variants$pct_pos_in_prot>1)&df_indel_variants$is_coding,yes=FALSE,no=df_indel_variants$is_coding)
# df_indel_variants$pct_pos_in_prot <- ifelse(test = df_indel_variants$is_coding,yes = df_indel_variants$pct_pos_in_prot,no=NA)
df_indel_variants$lbl_fixed <- factor(ifelse(test=df_indel_variants$is_fixed,yes="Fixed",no="Not fixed"),levels = c("Not fixed","Fixed"))
df_indel_variants$variant_type <- factor(df_indel_variants$variant_type,levels = c("Insertion","Deletion"))
ggplot(data = subset(df_indel_variants,variant_type=="Insertion"),mapping=aes(x=100*pct_pos_in_prot)) + geom_density(mapping = aes(col=is_in_frame)) + xlab("% protein length") + ylab("Density") + scale_color_brewer(palette = "Set1") + theme_bw() +theme(axis.title = element_text(size=16),axis.text = element_text(size=16),legend.title = element_text(size=12),legend.text = element_text(size=10),axis.text.x = element_text(size=14),legend.position = "right") + facet_grid(rows=vars(lbl_fixed),scales = "free_y")
ggsave(filename = "Distribution_of_Insertion_position_percentage.svg", path=output_workspace, width = 15, height = 10, units = "cm", device= svg)
ggsave(filename = "Distribution_of_Insertion_position_percentage.jpg", path=output_workspace, width = 15, height = 10, units = "cm", device= jpeg)

length(unique(subset(df_indel_variants,variant_type=="Insertion"&is_in_frame&is_fixed)$mutation))
length(unique(subset(df_indel_variants,variant_type=="Insertion"&(!is_in_frame)&is_fixed)$mutation))

ggplot(data = subset(df_indel_variants,variant_type=="Deletion"),mapping=aes(x=100*pct_pos_in_prot)) + geom_density(mapping = aes(col=is_in_frame)) + xlab("% protein length") + ylab("Density") + scale_color_brewer(palette = "Set1") + theme_bw() +theme(axis.title = element_text(size=16),axis.text = element_text(size=16),legend.title = element_text(size=12),legend.text = element_text(size=10),axis.text.x = element_text(size=14),legend.position = "right") + facet_grid(rows=vars(lbl_fixed),scales = "free_y")
ggsave(filename = "Distribution_of_Deletion_position_percentage.svg", path=output_workspace, width = 15, height = 10, units = "cm", device= svg)
ggsave(filename = "Distribution_of_Deletion_position_percentage.jpg", path=output_workspace, width = 15, height = 10, units = "cm", device= jpeg)

length(unique(subset(df_indel_variants,variant_type=="Deletion"&is_in_frame&is_fixed)$mutation))
length(unique(subset(df_indel_variants,variant_type=="Deletion"&(!is_in_frame)&is_fixed)$mutation))

# #Gain of function evidence for recurrently fixed indels
# v_fixed_indels_recurrence_as_fixed<- vapply(X = sort(unique(subset(df_indel_variants,is_fixed)$mutation)), FUN = function(the_indl) length(unique(subset(df_indel_variants,is_fixed & mutation==the_indl)$population)), FUN.VALUE = 0)
# v_lst_recurrently_fixed_indels <- names(v_fixed_indels_recurrence_as_fixed)[v_fixed_indels_recurrence_as_fixed>1]
# #print(sort(v_fixed_indels_recurrence_as_fixed[v_fixed_indels_recurrence_as_fixed>1],dec=T))
# df_indel_variants$is_recurrently_fixed <- df_indel_variants$mutation%in%v_lst_recurrently_fixed_indels
# v_is_recurrently_fixed_indels_increasing_in_freq_in_time <- NULL
# for (current_recurrently_fixed_indels in v_lst_recurrently_fixed_indels){
#   v_current_lst_pop <- unique(subset(df_indel_variants,mutation==current_recurrently_fixed_indels)$population)
#   v_is_increasing_in_freq_in_pops <- NULL
#   for (current_pop in v_current_lst_pop){
#     range_time_in_current_pop <- range(subset(df_indel_variants,mutation==current_recurrently_fixed_indels&population==current_pop)$int_time,na.rm = T)
#     min_time_in_current_pop <- range_time_in_current_pop[1]
#     max_time_in_current_pop <- range_time_in_current_pop[2]
#     v_is_increasing_in_freq_in_pops <- c(v_is_increasing_in_freq_in_pops,subset(df_indel_variants,mutation==current_recurrently_fixed_indels&population==current_pop&int_time==max_time_in_current_pop)$VarFreq > subset(df_indel_variants,mutation==current_recurrently_fixed_indels&population==current_pop&int_time==min_time_in_current_pop)$VarFreq)
#   }
#   v_is_recurrently_fixed_indels_increasing_in_freq_in_time <- c(v_is_recurrently_fixed_indels_increasing_in_freq_in_time,sum(v_is_increasing_in_freq_in_pops)>0) #sum(v_is_increasing_in_freq_in_pops)/length(v_is_increasing_in_freq_in_pops)>0.5
# }
# 
# names(v_is_recurrently_fixed_indels_increasing_in_freq_in_time) <- v_lst_recurrently_fixed_indels
# v_lst_recurrently_fixed_indels_increasing_in_freq_in_time <- names(v_is_recurrently_fixed_indels_increasing_in_freq_in_time)[v_is_recurrently_fixed_indels_increasing_in_freq_in_time]
# #v_fixed_indels_recurrence_as_fixed[v_lst_recurrently_fixed_indels_increasing_in_freq_in_time]
# df_indel_variants$is_recurrently_fixed_and_increasing_in_freq_in_multiple_pops <- df_indel_variants$mutation %in% v_lst_recurrently_fixed_indels_increasing_in_freq_in_time
# v_lst_products_where_recurrently_fixed_indels_with_increasing_freq_map <- sort(unique(subset(df_indel_variants,is_recurrently_fixed_and_increasing_in_freq_in_multiple_pops)$gene_product))
#   #find the genes in cis of intergenic fixed indels with increasing freq in time
# 
# #VAF timeseries for recurrent indels
# 
# #Find indels mapping overlapping genes

#Indel bias 
df_per_sample_indel_bias <- subset(df_indel_variants,is_fixed) %>%
  group_by(population,variant_type)  %>%
  dplyr::summarise(prevalence=length(unique(mutation))) %>%
  as.data.frame()
df_per_sample_indel_bias2 <- df_per_sample_indel_bias %>%
  group_by(population)  %>%
  dplyr::summarise(insertion_bias=prevalence[2]/prevalence[1]) %>%
  as.data.frame()

#Visualization of PYK1 and its mutations 
  #Gviz Load libraries
library(genemodel)

# Convert the pyk1 chromosome and the gene to Biostrings objects
the_genome_name <- paste0(output_workspace,"ChrI_Spom.fasta")
the_genome <- Biostrings::readDNAStringSet(the_genome_name)
the_gene_name <- "pyk1"
pyk1_start <- subset(df_loci_annotations,the_product==the_gene_name)$start
pyk1_end <- subset(df_loci_annotations,the_product==the_gene_name)$end
pyk1_mrna_with_utr_total_nt_length <- abs(pyk1_end- pyk1_start) + 1
pyk1_strand <- subset(df_loci_annotations,the_product==the_gene_name)$strand
chrom_name <- subset(df_loci_annotations,the_product==the_gene_name)$chr  # Get chromosome name
chrom_length <- nchar(genome_refseq[[chrom_name]])  # Get chromosome length
the_chromosome <- GRanges(seqnames = chrom_name, 
                      ranges = IRanges(start = 1, end = chrom_length))
dna_string <- DNAStringSet(paste(genome_refseq[[chrom_name]], collapse=""))
pyk1_acc_nb <- subset(df_loci_annotations,the_product=="pyk1")$accession_number

all_pyk1_gff <- subset(df_loci_annotations,grepl(pattern = pyk1_acc_nb,x = annotation,fixed = T))

# Generate the gene model
# Convert GFF3 data to a gene model dataframe with 'coordinates' column
create_gene_model <- function(gff_data) {
  
  # Initialize an empty data frame to store gene model
  gene_model <- data.frame(type = character(0), 
                           coordinates = character(0), 
                           stringsAsFactors = FALSE)
  
  # Iterate through GFF3 data and convert it to gene model format
  for (i in 1:nrow(gff_data)) {
    current_entry <- gff_data[i,]
    # Assign type based on feature type
    v_conv <- c("gene" = "ORF",             # ORF equivalent to gene in gene model
                           "CDS" = "coding_region",    # CDS to Coding Region
                           "five_prime_UTR" = "5' utr", # 5' UTR equivalent
                           "three_prime_UTR" = "3' utr",# 3' UTR equivalent
                           "intron" = "intron") # If no match, ignore
    current_entry$type <- v_conv[current_entry$site_type]
    current_entry$coordinates <- paste(abs(current_entry$start - pyk1_end)+1, abs(current_entry$end - pyk1_end)+1, sep = "-")
    gene_model <- rbind(gene_model, current_entry[,c("type","coordinates")])
  }
  
  return(gene_model)
}

# Apply the function to the GFF3 data
gene_model_df <- create_gene_model(all_pyk1_gff)
gene_model_df <- subset(gene_model_df,!is.na(type))

genemodel::genemodel.plot(model=gene_model_df, start=1, bpstop=abs(pyk1_start - pyk1_end)+1, orientation="reverse", xaxis=T)

#import mutations 
Table_pyk1_unique_pop_muts <- readRDS(paste0(output_workspace,"Table_pyk1_unique_pop_muts.rds")) #unique(subset(df_variants_time_series, gene_product=="pyk1")[,colnames(Table_pyk1_unique_pop_muts)])
Table_pyk1_unique_pop_muts <- subset(Table_pyk1_unique_pop_muts,!is.na(pop_mut))
Table_pyk1_unique_pop_muts$pos_in_gene <- abs(Table_pyk1_unique_pop_muts$Position - pyk1_end) + 1
Table_pyk1_unique_pop_muts$genomic_end <- NA
#add mutations one at a time
for (current_pop_mut_id in 1:nrow(Table_pyk1_unique_pop_muts)){
  if (grepl("+",x = Table_pyk1_unique_pop_muts$variant_type2,fixed = T)){#Insertion
    #genemodel::mutation.plot(abs(Table_pyk1_unique_pop_muts$Position[current_pop_mut_id]-pyk1_end)+1, abs(Table_pyk1_unique_pop_muts$Position[current_pop_mut_id]-pyk1_end)+1+nchar(Table_pyk1_unique_pop_muts$VarAllele[current_pop_mut_id]), text=Table_pyk1_unique_pop_muts$pop_mut[current_pop_mut_id], col="black", drop=-.5, haplotypes=c("red"))
    Table_pyk1_unique_pop_muts$genomic_end[current_pop_mut_id] <- Table_pyk1_unique_pop_muts$Position[current_pop_mut_id]-1
  }else if (grepl("-",x = Table_pyk1_unique_pop_muts$variant_type2,fixed = T)){#Deletion
    #genemodel::mutation.plot(abs(Table_pyk1_unique_pop_muts$Position[current_pop_mut_id]-pyk1_end)+1, abs(Table_pyk1_unique_pop_muts$Position[current_pop_mut_id]-pyk1_end)+1+nchar(Table_pyk1_unique_pop_muts$VarAllele[current_pop_mut_id]), text=Table_pyk1_unique_pop_muts$pop_mut[current_pop_mut_id], col="black", drop=-.5, haplotypes=c("red"))
    Table_pyk1_unique_pop_muts$genomic_end[current_pop_mut_id] <- Table_pyk1_unique_pop_muts$Position[current_pop_mut_id]-nchar(Table_pyk1_unique_pop_muts$VarAllele[current_pop_mut_id])-1
  }else if (grepl("SNV",x = Table_pyk1_unique_pop_muts$variant_type2,fixed = T)){#SNP in noncoding region
    #genemodel::mutation.plot(abs(Table_pyk1_unique_pop_muts$Position[current_pop_mut_id]-pyk1_end)+1, abs(Table_pyk1_unique_pop_muts$Position[current_pop_mut_id]-pyk1_end)+1, text=Table_pyk1_unique_pop_muts$pop_mut[current_pop_mut_id], col="black", drop=-.5, haplotypes=c("red"))
    Table_pyk1_unique_pop_muts$genomic_end[current_pop_mut_id] <- Table_pyk1_unique_pop_muts$Position[current_pop_mut_id]-1
  }else if (grepl("misssense",x = Table_pyk1_unique_pop_muts$variant_type2,fixed = T)){#SNP in noncoding region
    #genemodel::mutation.plot(abs(Table_pyk1_unique_pop_muts$Position[current_pop_mut_id]-pyk1_end)+1, abs(Table_pyk1_unique_pop_muts$Position[current_pop_mut_id]-pyk1_end)+1, text=paste0(Table_pyk1_unique_pop_muts$pop_mut[current_pop_mut_id],"\n",Table_pyk1_unique_pop_muts$old_aa[current_pop_mut_id],"->",Table_pyk1_unique_pop_muts$new_aa[current_pop_mut_id]), col="black", drop=-.5, haplotypes=c("red"))
    Table_pyk1_unique_pop_muts$genomic_end[current_pop_mut_id] <- Table_pyk1_unique_pop_muts$Position[current_pop_mut_id]-1
  }
}

#With gggenomes
library("gggenomes")
gggenomes(paste0(output_workspace,"Schizosaccharomyces_pombe_all_chromosomes.gff3")) + geom_gene() 
#import gene info
#s0 <- read_seqs(paste0(output_workspace,"Schizosaccharomyces_pombe_all_chromosomes.fasta"))
g0 <- read_feats(paste0(output_workspace,"pyk1_S_pombe_.gff3"))
#mutations
v0 <- tibble::tibble(file_id=rep("pyk1_S_~",nrow(Table_pyk1_unique_pop_muts)),
                     seq_id = rep("I",nrow(Table_pyk1_unique_pop_muts)), start = Table_pyk1_unique_pop_muts$Position,
                     end = Table_pyk1_unique_pop_muts$genomic_end, length = abs(Table_pyk1_unique_pop_muts$Position-Table_pyk1_unique_pop_muts$genomic_end),
                     type = c("SNP", "SNP", "SNP", "SNP", "Deletion", "Deletion"),
                     ALT = c(Table_pyk1_unique_pop_muts$VarAllele[1:4], ".", "."),
                     REF = c("T", "C", "G", "C", "G", "T"),
                     pop_mut = Table_pyk1_unique_pop_muts$pop_mut
) 
o0 <- dplyr::bind_rows(g0,v0)
#figure
gggenomes(feats=o0) +
  geom_seq() +         # draw contig/chromosome lines
  geom_seq_label() +   # label each sequence 
  geom_gene() +
  geom_variant(aes(shape = type),offset = -0.1) +
  scale_shape_variant() +
  geom_variant(aes(label = pop_mut), geom="text",offset = -c(0.25,0.35,0.45,0.25,0.5,0.75),angle=0,position = position_dodge(0.1)) +
  geom_bin_label()

# library(Gviz)
# library(GenomicFeatures)
# library(GenomicRanges)
# #create the loci annotation DB
# sp_db <- GenomicFeatures::makeTxDbFromGFF(paste0(output_workspace,"Schizosaccharomyces_pombe_all_chromosomes.gff3"), format = "gff3")
# granges_genome <- rtracklayer::import(paste0(output_workspace,"Schizosaccharomyces_pombe_all_chromosomes.gff3"),format = "gff3")
# mcols(granges_genome) <- NULL
# chrI_granges_genome <- granges_genome[seqnames(granges_genome)=="I"&start(granges_genome) >= pyk1_end - 10000&end(granges_genome) <= pyk1_start + 10000]
# genome(chrI_granges_genome) <- "S_pombe_972h-"
# chr <- as.character(unique(seqnames(chrI_granges_genome)))
# gen <- genome(chrI_granges_genome)
# atrack <- AnnotationTrack(chrI_granges_genome, name = "chrI")
# plotTracks(atrack)
# gtrack <- GenomeAxisTrack()
# itrack <- IdeogramTrack(genome = gen, chromosome = chr)
# 
# #find pyk1 on the db
# genes_data <- genes(sp_db)
# pyk1_gene <- genes_data[genes_data$gene_id == pyk1_acc_nb, ]
# 
# #Create a GRanges object for the chromosome
# chromosome <- GRanges(seqnames = chrom_name,
#                       ranges = IRanges(start = 1, end = chrom_length))
# 
# #pyk1 "gene" loci annotation
# pyk1_gene <- subset(df_loci_annotations,site_type=="gene"&grepl(pattern = pyk1_acc_nb,x = annotation,fixed = T))
# 
# #pyk1 "exon" loci annotation
# pyk1_exons <- subset(df_loci_annotations,site_type=="CDS"&grepl(pattern = pyk1_acc_nb,x = annotation,fixed = T))
# pyk1_exons$exon_nb <- 1:nrow(pyk1_exons)
# #pyk1 "UTRs" loci annotation
# pyk1_utrs <- subset(df_loci_annotations,(site_type=="five_prime_UTR"|site_type=="three_prime_UTR")&grepl(pattern = pyk1_acc_nb,x = annotation,fixed = T))
# 
# #Convert to GRanges objects (Define the gene coordinates on its chromosome)
# pyk1_gr <- GRanges(seqnames = "pyk1", 
#                    ranges = IRanges(start = subset(df_loci_annotations,the_product=="pyk1")$start, end = subset(df_loci_annotations,the_product=="pyk1")$end),  # Example gene location
#                    strand = "-")
# # 
# # pyk1_gr <- GRanges(seqnames = seqnames(pyk1_gene),
# #                    ranges = ranges(pyk1_gene),
# #                    strand = strand(pyk1_gene),
# #                    gene = "pyk1")
# 
# pyk1_exon_gr <- GRanges(seqnames = paste0("Exon",pyk1_exons$exon_nb),
#                         ranges = IRanges(start = pyk1_exons$start, end = pyk1_exons$end),
#                         strand = pyk1_exons$strand)
# 
# pyk1_utr_gr <- GRanges(seqnames = pyk1_utrs$site_type,
#                        ranges = IRanges(start = pyk1_utrs$start, end = pyk1_utrs$end),
#                        strand = pyk1_utrs$strand)

#pyk1 visualization with trackViewer
library(trackViewer)
#promoter position in the pyk1 loci
promoter_start <- abs(3846665- pyk1_end) + 1
promoter_end <- abs(3846724 - pyk1_end) + 1
#CDS position in the pyk1 gene
CDS_start <- abs(subset(df_loci_annotations,grepl(pattern = pyk1_acc_nb,x = annotation,fixed = T)&site_type=="CDS")$end- pyk1_end) + 1
CDS_end <- abs(subset(df_loci_annotations,grepl(pattern = pyk1_acc_nb,x = annotation,fixed = T)&site_type=="CDS")$start - pyk1_end) + 1
#UTR position in the pyk1 loci
five_prime_utr_start <- abs(subset(df_loci_annotations,grepl(pattern = pyk1_acc_nb,x = annotation,fixed = T)&site_type=="five_prime_UTR")$end- pyk1_end) + 1
five_prime_utr_end <- abs(subset(df_loci_annotations,grepl(pattern = pyk1_acc_nb,x = annotation,fixed = T)&site_type=="five_prime_UTR")$start - pyk1_end) + 1
three_prime_utr_start <- abs(subset(df_loci_annotations,grepl(pattern = pyk1_acc_nb,x = annotation,fixed = T)&site_type=="three_prime_UTR")$end- pyk1_end) + 1
three_prime_utr_end <- abs(subset(df_loci_annotations,grepl(pattern = pyk1_acc_nb,x = annotation,fixed = T)&site_type=="three_prime_UTR")$start - pyk1_end) + 1
#Define the SNP Grange
#SNP_Table_pyk1_unique_pop_muts <- subset(Table_pyk1_unique_pop_muts,!grepl(pattern = "Indel",x = Table_pyk1_unique_pop_muts$variant_type2,fixed = TRUE))
sample.gr <- GRanges("chrI", IRanges(Table_pyk1_unique_pop_muts$pos_in_gene, width=1, names=Table_pyk1_unique_pop_muts$pop_mut),strand = "-")
sample.gr$color <- ifelse(test=grepl(pattern = "Indel",x = Table_pyk1_unique_pop_muts$variant_type2,fixed = TRUE),yes="#E41A1CFF",no="#FFFFFF")
sample.gr$label.parameter.rot <- 60
#Define the indel, CDS, promoter and UTR Granges
#Indel_Table_pyk1_unique_pop_muts <- subset(Table_pyk1_unique_pop_muts,grepl(pattern = "Indel",x = Table_pyk1_unique_pop_muts$variant_type2,fixed = TRUE))
drawing_order <- c(1,3,4,2)
features <- GRanges("chrI", IRanges(c(CDS_start,promoter_start,five_prime_utr_start,three_prime_utr_start)[drawing_order], width=c(abs(CDS_end-CDS_start)+1,abs(promoter_end-promoter_start)+1,abs(five_prime_utr_end-five_prime_utr_start) +1,abs(three_prime_utr_end-three_prime_utr_start)+1)[drawing_order], names=c("CDS","promoter","5' UTR","3' UTR")[drawing_order]),strand = "-")
features$fill <- c("#4DAF4AFF","#51C6E6","#DFA32D","#FF8833")[drawing_order]
lolliplot(sample.gr, features, legend="legend")

#save session
library("session")
save.session(file = paste0(output_workspace,"S_Pombe_evo.Rda"))

