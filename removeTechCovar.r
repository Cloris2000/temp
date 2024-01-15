#import packages
library(GEOquery)
library(limma)
library(ggplot2)
library(dplyr)
library(tidyverse)
library(matrixStats)
library(broom)
library(knitr)
library(ggpubr)
library(biomaRt)
library(ggrepel)
library(patchwork)
library(ggsignif)
library(modelr)
library(cowplot)
library(gridExtra)
library(RColorBrewer)
library(DESeq2)
library(factoextra)
library(PCAtools)
library(markerGeneProfile)


#import RNAseq expression data from ROSMAP
#count_cpm_filtered <- readRDS('/external/rprshnas01/netdata_kcni/stlab/Xiaolin/metabrain/data/ROSMAP_DLPFC_count_cpm_filtered.RData')

#load metadata: col names in exp == row names in meta

### meta data manipulation
#combine different batches into a large matrix
batch1 = read_tsv('/external/rprshnas01/external_data/rosmap/gene_expression/RNAseq_Harmonization/Gene Expression (Raw Gene Counts)/Rosmap_Gene_Quantification/Rosmap_Batch1_Stranded/ROSMAP_batch1_gene_all_counts_matrix.txt')
batch2 = read_tsv('/external/rprshnas01/external_data/rosmap/gene_expression/RNAseq_Harmonization/Gene Expression (Raw Gene Counts)/Rosmap_Gene_Quantification/Rosmap_Batch2_Stranded/ROSMAP_batch2_gene_all_counts_matrix.txt')
batch3 = read_tsv('/external/rprshnas01/external_data/rosmap/gene_expression/RNAseq_Harmonization/Gene Expression (Raw Gene Counts)/Rosmap_Gene_Quantification/Rosmap_Batch3_Stranded/ROSMAP_batch3_gene_all_counts_matrix.txt')
batch4 = read_tsv('/external/rprshnas01/external_data/rosmap/gene_expression/RNAseq_Harmonization/Gene Expression (Raw Gene Counts)/Rosmap_Gene_Quantification/Rosmap_Batch4_Stranded/ROSMAP_batch4_gene_all_counts_matrix.txt')


sample_names = c(colnames(batch1), colnames(batch2), colnames(batch3), colnames(batch4))

rosmap_meta = read_csv('/external/rprshnas01/external_data/rosmap/gene_expression/RNAseq_Harmonization/Gene Expression (Raw Gene Counts)/Metadata/RNAseq_Harmonization_ROSMAP_combined_metadata.csv')

rosmap_meta %>% distinct(projid, .keep_all = T) %>% dim #1169 samples 38

batch1_prov = read_csv('/external/rprshnas01/external_data/rosmap/gene_expression/RNAseq_Harmonization/Gene Expression (Raw Gene Counts)/Rosmap_Gene_Quantification/Rosmap_Batch1_Stranded/ROSMAP_batch1_provenance.csv')
batch2_prov = read_csv('/external/rprshnas01/external_data/rosmap/gene_expression/RNAseq_Harmonization/Gene Expression (Raw Gene Counts)/Rosmap_Gene_Quantification/Rosmap_Batch2_Stranded/ROSMAP_batch2_provenance.csv')
batch3_prov = read_csv('/external/rprshnas01/external_data/rosmap/gene_expression/RNAseq_Harmonization/Gene Expression (Raw Gene Counts)/Rosmap_Gene_Quantification/Rosmap_Batch3_Stranded/ROSMAP_batch3_provenance.csv')
batch4_prov = read_csv('/external/rprshnas01/external_data/rosmap/gene_expression/RNAseq_Harmonization/Gene Expression (Raw Gene Counts)/Rosmap_Gene_Quantification/Rosmap_Batch4_Stranded/ROSMAP_batch4_provenance.csv')
synapse_meta = bind_rows(list(batch1_prov, batch2_prov, batch3_prov, batch4_prov))
id_matcher <- synapse_meta[, c("specimenID", "id")] 

# import tech covariates during alignment
batch1_star = read.table('/external/rprshnas01/external_data/rosmap/gene_expression/RNAseq_Harmonization/Gene Expression (Raw Gene Counts)/Rosmap_Gene_Quantification/Rosmap_Batch1_Stranded/ROSMAP_batch1_Star_Log_Merged_clean.txt', header = T, sep = '\t')
batch2_star = read.table('/external/rprshnas01/external_data/rosmap/gene_expression/RNAseq_Harmonization/Gene Expression (Raw Gene Counts)/Rosmap_Gene_Quantification/Rosmap_Batch2_Stranded/ROSMAP_batch2_Star_Log_Merged_clean.txt', header = T, sep = '\t')
batch3_star = read.table('/external/rprshnas01/external_data/rosmap/gene_expression/RNAseq_Harmonization/Gene Expression (Raw Gene Counts)/Rosmap_Gene_Quantification/Rosmap_Batch3_Stranded/ROSMAP_batch3_Star_Log_Merged_clean.txt', header = T, sep = '\t')
batch4_star = read.table('/external/rprshnas01/external_data/rosmap/gene_expression/RNAseq_Harmonization/Gene Expression (Raw Gene Counts)/Rosmap_Gene_Quantification/Rosmap_Batch4_Stranded/ROSMAP_batch4_Star_Log_Merged_clean.txt', header = T, sep = '\t')
star_meta = bind_rows(list(batch1_star, batch2_star, batch3_star, batch4_star))
colnames(star_meta)[1] = "specimenID"

batch1_align = read.table('/external/rprshnas01/external_data/rosmap/gene_expression/RNAseq_Harmonization/Gene Expression (Raw Gene Counts)/Rosmap_Gene_Quantification/Rosmap_Batch1_Stranded/ROSMAP_batch1_Study_all_metrics_matrix_clean.txt', header = T, sep = '\t')
batch2_align = read.table('/external/rprshnas01/external_data/rosmap/gene_expression/RNAseq_Harmonization/Gene Expression (Raw Gene Counts)/Rosmap_Gene_Quantification/Rosmap_Batch2_Stranded/ROSMAP_batch2_Study_all_metrics_matrix_clean.txt', header = T, sep = '\t')
batch3_align = read.table('/external/rprshnas01/external_data/rosmap/gene_expression/RNAseq_Harmonization/Gene Expression (Raw Gene Counts)/Rosmap_Gene_Quantification/Rosmap_Batch3_Stranded/ROSMAP_batch3_Study_all_metrics_matrix_clean.txt', header = T, sep = '\t')
batch4_align = read.table('/external/rprshnas01/external_data/rosmap/gene_expression/RNAseq_Harmonization/Gene Expression (Raw Gene Counts)/Rosmap_Gene_Quantification/Rosmap_Batch4_Stranded/ROSMAP_batch4_Study_all_metrics_matrix_clean.txt', header = T, sep = '\t')
align_meta = bind_rows(list(batch1_align, batch2_align, batch3_align, batch4_align))
colnames(align_meta)[1] = "specimenID"

#combine with tech covar in batch star and batch align
rosmap_meta = left_join(rosmap_meta %>% filter(assay == 'rnaSeq'), star_meta, by = 'specimenID') 
rosmap_meta = left_join(rosmap_meta %>% filter(assay == 'rnaSeq'), align_meta, by = 'specimenID') 


rosmap_meta = left_join(rosmap_meta %>% filter(assay == 'rnaSeq'), id_matcher, by = 'specimenID') 
rosmap_meta$synapseID <- rosmap_meta$id
rosmap_meta <- rosmap_meta[rosmap_meta$synapseID %in% sample_names, ]

rosmap_meta %>% filter(assay == 'rnaSeq') %>% distinct(projid, tissue) %>% group_by(tissue) %>% tally

rosmap_meta$Age_norm <- gsub("90+", "90", rosmap_meta$age_death)
rosmap_meta$Age_norm = as.numeric(gsub("[+]", "", rosmap_meta$Age_norm))


### include the samples in DLPC region for four batches count matrix
use_tissues = c('dorsolateral prefrontal cortex')

use_tissue_synapse_sample_ids = rosmap_meta %>% filter(tissue %in% use_tissues) %>%
   pull(synapseID)  #[1] 1142
rosmap_meta_DLPFC <- rosmap_meta %>% filter(tissue %in% use_tissues)
rownames(rosmap_meta_DLPFC) <- rosmap_meta_DLPFC$synapseID

saveRDS(rosmap_meta_DLPFC, file="/external/rprshnas01/netdata_kcni/stlab/Xiaolin/metabrain_PCA/data/rosmap_meta_DLPFC.RData")
# rosmap_meta_DLPFC <- readRDS('/external/rprshnas01/netdata_kcni/stlab/Xiaolin/metabrain_PCA/data/rosmap_meta_DLPFC.RData')



#combine different batches into a large matrix
combined_rnaseq_mat = cbind(batch1[, intersect(c(use_tissue_synapse_sample_ids), colnames(batch1))], 
                                  batch2[, intersect(c(use_tissue_synapse_sample_ids), colnames(batch2))],
                                         batch3[, intersect(c(use_tissue_synapse_sample_ids), colnames(batch3))], 
                                                batch4[, intersect(c(use_tissue_synapse_sample_ids), colnames(batch4))]) 
# raw count matrix for ROSMAP in DLPFC
combined_rnaseq_mat = combined_rnaseq_mat[5:60607, ] %>% as.data.frame()
rownames(combined_rnaseq_mat) = batch1[5:60607, ]  %>% pull(feature) %>% substr(., 1, 15) %>% make.names(., unique = T)

##################################################################### DESeq ##########################################################################
####RNA-seq alignment and Quality control PCA-based outlier identification
### PCA for RAW count matrix without normalization and covariate correction
#make sure data is matrix, not a data frame
count.matix <- as.matrix(combined_rnaseq_mat) # raw count
 # prepare certain phenotypes: facotr or numeric values
rosmap_meta_DLPFC_cleaned <- rosmap_meta_DLPFC[!is.na(rosmap_meta_DLPFC$RIN), ] #remove the samples with NA RIN

count.matix_cleaned <- count.matix[, !(colnames(count.matix) %in% rosmap_meta_DLPFC[is.na(rosmap_meta_DLPFC$RIN), ]$synapseID)]
dds <- DESeqDataSetFromMatrix(countData=count.matix_cleaned, colData=rosmap_meta_DLPFC_cleaned, design=~1) #no design?

zscore_data <- readRDS('/external/rprshnas01/netdata_kcni/stlab/Xiaolin/metabrain_PCA/data/ROSMAP_DLPFC_zscore_data.RData') #41796  1142
#Prepare metadata for PCATools::pca()
rosmap_meta_DLPFC_cleaned <- read.csv("/external/rprshnas01/netdata_kcni/stlab/Xiaolin/metabrain_PCA/rosmap_meta_DLPFC_cleaned.csv") #[1] 1141  115
#remove tech cov with NA
rosmap_meta_DLPFC_cleaned_noNA <- rosmap_meta_DLPFC_cleaned[complete.cases(t(rosmap_meta_DLPFC_cleaned))] #[1] 1141   87
#remove tech cov with only one value 
# Step 1: Calculate the number of unique values for each column
num_unique_values <- sapply(rosmap_meta_DLPFC_cleaned_noNA, function(x) length(unique(x)))
# Step 2: Create a logical vector to identify columns with more than one unique value
cols_to_keep <- num_unique_values > 1
# Step 3: Subset the DataFrame to include only those columns
rosmap_meta_DLPFC_cleaned_noNA_rmVar <- rosmap_meta_DLPFC_cleaned_noNA[, cols_to_keep] #[1] 1141   73
rownames(rosmap_meta_DLPFC_cleaned_noNA_rmVar) <- rosmap_meta_DLPFC_cleaned_noNA_rmVar$id

tech_covar_list <- colnames(rosmap_meta_DLPFC_cleaned_noNA_rmVar)
# remove columns that are not tech cov
tech_covar_list <- tech_covar_list[!tech_covar_list %in% c("X", "specimenID", "id","synapseID","tissue","assay","organ",
"Started.job.on","Started.mapping.on","Finished.on","Number.of.splices..Total", "Number.of.splices..Annotated..sjdb.", 
"Number.of.splices..GT.AG", "Number.of.splices..GC.AG", "Number.of.splices..AT.AC", "Number.of.splices..Non.canonical",
"Mismatch.rate.per.base...","Deletion.rate.per.base","Deletion.average.length","Insertion.rate.per.base", "Insertion.average.length")] #include Mapping.speed..Million.of.reads.per.hour here?
print(tech_covar_list) # [55]

# prepare certain phenotypes: facotr or numeric values
metadata <- rosmap_meta_DLPFC_cleaned_noNA_rmVar[,tech_covar_list]
for (tech_cov in tech_covar_list){
  if(is.numeric(metadata[[tech_cov]])){
    metadata[[tech_cov]] <- as.numeric(metadata[[tech_cov]])}
  else{
    metadata[[tech_cov]] <- factor(metadata[[tech_cov]])
    print(tech_cov)
  }
}

# double check to remove samples from the pdata that have any NA value
discard <- apply(metadata, 1, function(x) any(is.na(x)))
metadata <- metadata[!discard,] #checked! nothing changed:  1141   73

# filter the expression data to match the samples in our pdata
zscore_data <- zscore_data[,which(colnames(zscore_data) %in% rownames(metadata))]

# check that sample names match exactly between pdata and expression data 
all(colnames(zscore_data) == rownames(metadata)) #[1] 1141


# try PCAtools 
p <- PCAtools::pca(zscore_data, metadata = metadata)


eigencorvalue <- function(pcaobj, components = getComponents(pcaobj, seq_len(10)), metavars,  corUSE = 'pairwise.complete.obs', corFUN = 'pearson') {
  data <- pcaobj$rotated
  metadata <- pcaobj$metadata
  corvals <- list()

  # Issue warning if any columns to use are not numeric
  for (i in seq_len(length(components))) {
    if (!is.numeric(data[, components[i]])) {
      warning(components[i], ' is not numeric - please check the source data')
    }
  }
  for (i in seq_len(length(metavars))) {
    if (!is.numeric(metadata[, metavars[i]])) {
      warning(metavars[i], ' is not numeric - please check the source data')
    }
  }

  # Convert the data for x and y to data matrix
  xvals <- data.matrix(data[, which(colnames(data) %in% components), drop = FALSE])
  yvals <- metadata[, which(colnames(metadata) %in% metavars), drop = FALSE]

  # Ensure that non-numeric variables are coerced to numeric
  chararcter_columns <- unlist(lapply(yvals, is.numeric))
  chararcter_columns <- !chararcter_columns
  chararcter_columns <- names(which(chararcter_columns))
  for (c in chararcter_columns) {
    yvals[, c] = as.numeric(as.factor(yvals[, c]))
  }

  yvals <- data.matrix(yvals)

  # Create correlation table
  cor_matrix <- cor(xvals, yvals, use = corUSE, method = corFUN)

  # Store the correlation values and their corresponding variable names in a list
  corvals$x_names <- colnames(xvals)
  corvals$y_names <- colnames(yvals)
  corvals$correlation_values <- cor_matrix

  # Return the list containing the correlation values and variable names
  return(corvals)
}

all_corvals <- eigencorvalue(pcaobj=p, metavars=tech_covar_list, corFUN = "pearson")

# get the top 20 tech cov (all positive? or both positive and negative?)
# Extract the correlation values
correlation_matrix <- all_corvals$correlation_values

# Calculate the absolute values of correlation coefficients
abs_correlation_values <- abs(correlation_matrix)

# Flatten the matrix into a vector (since you have multiple x and y variables)
cor_vector <- as.vector(abs_correlation_values)

# Find the indices of the top 20 values
top_20_indices <- order(cor_vector, decreasing = TRUE)[1:40]

# Get the corresponding variable names
x_names <- all_corvals$x_names
y_names <- all_corvals$y_names

# Extract the top 20 technical covariates and their names
top_20_covariates <- data.frame(
  x_var = rep(x_names, each = length(y_names)),
  y_var = rep(y_names, times = length(x_names)),
  correlation = cor_vector,
  stringsAsFactors = FALSE
)

# Sort the top 20 covariates by correlation value
top_20_covariates <- top_20_covariates[top_20_indices, ]

# Print or return the top 20 technical covariates
top_20 <- unique(top_20_covariates$y_var)


zscore_data_removedBatchEff <- zscore_data
# for (tech_cov in top_20){
#   zscore_data_removedBatchEff <- removeBatchEffect(zscore_data_removedBatchEff, as.vector(metadata[, tech_cov]))#what should we include as covariates here? 
# }

#zscore data with batch effect removed
zscore_data_removedBatchEff <- removeBatchEffect(x=zscore_data_removedBatchEff, batch=as.vector(metadata[,"sequencingBatch"]), batch2=as.vector(metadata[,"libraryPrep"]))

#zscore data with batch effect and tech covariates removed
## need to make sure the tech cov matrix contains numerical values for every column
metadata_nobatch <- metadata[, !colnames(metadata) %in% c("sequencingBatch", "libraryPrep")]

# Check if every column contains all numeric values
non_numeric_columns <- names(which(!sapply(metadata_nobatch, function(x) all(is.numeric(x)))))
print(non_numeric_columns)

# modify the non-numeric columns to numeric
# 1. For "notes", change data contribution batch 1 to 1 ...
metadata_nobatch[, "notes"] <- str_split(metadata_nobatch[, "notes"], " ") %>% lapply(function(x) x[4]) %>% unlist()
metadata_nobatch[,"Uniquely.mapped.reads.."] <- as.numeric(sub("%", "", metadata_nobatch[,"Uniquely.mapped.reads.."]))
metadata_nobatch[,"X..of.reads.mapped.to.multiple.loci"] <- as.numeric(sub("%", "", metadata_nobatch[,"X..of.reads.mapped.to.multiple.loci"]))
metadata_nobatch[,"X..of.reads.mapped.to.too.many.loci"] <- as.numeric(sub("%", "", metadata_nobatch[,"X..of.reads.mapped.to.too.many.loci"]))
metadata_nobatch[,"X..of.reads.unmapped..too.short"] <- as.numeric(sub("%", "", metadata_nobatch[,"X..of.reads.unmapped..too.short"]))
metadata_nobatch[,"X..of.reads.unmapped..other"] <- as.numeric(sub("%", "", metadata_nobatch[,"X..of.reads.unmapped..other"]))

#make sure all columns numeric now
for (i in seq_len(length(non_numeric_columns))) {
  metadata_nobatch[, non_numeric_columns[i]] <- as.numeric(as.factor(metadata_nobatch[, non_numeric_columns[i]]))
}






############ newly updated 12.2
metadata_nobatch <- metadata[, !colnames(metadata) %in% c("sequencingBatch","libraryPrep")]
metadata_nobatch <- metadata_nobatch[, top_20[top_20 != "libraryPrep"]]
# Check if every column contains all numeric values
non_numeric_columns <- names(which(!sapply(metadata_nobatch, function(x) all(is.numeric(x)))))
print(non_numeric_columns)

# modify the non-numeric columns to numeric
# 1. For "notes", change data contribution batch 1 to 1 ...
metadata_nobatch[, "notes"] <- str_split(metadata_nobatch[, "notes"], " ") %>% lapply(function(x) x[4]) %>% unlist()
metadata_nobatch[,"Uniquely.mapped.reads.."] <- as.numeric(sub("%", "", metadata_nobatch[,"Uniquely.mapped.reads.."]))
metadata_nobatch[,"X..of.reads.mapped.to.multiple.loci"] <- as.numeric(sub("%", "", metadata_nobatch[,"X..of.reads.mapped.to.multiple.loci"]))
metadata_nobatch[,"X..of.reads.mapped.to.too.many.loci"] <- as.numeric(sub("%", "", metadata_nobatch[,"X..of.reads.mapped.to.too.many.loci"]))
metadata_nobatch[,"X..of.reads.unmapped..too.short"] <- as.numeric(sub("%", "", metadata_nobatch[,"X..of.reads.unmapped..too.short"]))
metadata_nobatch[,"X..of.reads.unmapped..other"] <- as.numeric(sub("%", "", metadata_nobatch[,"X..of.reads.unmapped..other"]))

#make sure all columns numeric now
for (i in seq_len(length(non_numeric_columns))) {
  metadata_nobatch[, non_numeric_columns[i]] <- as.numeric(as.factor(metadata_nobatch[, non_numeric_columns[i]]))
}
zscore_data_removedBatchEff_cov <- removeBatchEffect(x=zscore_data_removedBatchEff, batch=as.vector(metadata[,"sequencingBatch"]), batch2=as.vector(metadata[,"libraryPrep"]), covariates=metadata_nobatch) #covariates=metadata[, (!tech_cov %in% c("sequencingBatch", "libraryPrep"))]
#orig data: wired... 53 tech_cov??saveRDS(zscore_data_removedBatchEff_cov, file="/external/rprshnas01/netdata_kcni/stlab/Xiaolin/metabrain_PCA/data/ROSMAP_DLPFC_zscore_data_removedBatchEff_cov.RData")
###############################################################


#orig data: wired... 55 tech_cov??saveRDS(zscore_data_removedBatchEff_cov, file="/external/rprshnas01/netdata_kcni/stlab/Xiaolin/metabrain_PCA/data/ROSMAP_DLPFC_zscore_data_removedBatchEff_cov.RData")
p_removedBatchEff_cov <- PCAtools::pca(zscore_data_removedBatchEff_cov, metadata = metadata)
biplot_pca_removedBatchEff_cov <- biplot(p_removedBatchEff_cov)
#save into png file
# png('/external/rprshnas01/netdata_kcni/stlab/Xiaolin/metabrain_PCA/data/ROSMAP_DLPFC_biplot_pca_removedBatchEff_cov.png')
# biplot_pca_removedBatchEff_cov
# dev.off()
png('/external/rprshnas01/netdata_kcni/stlab/Xiaolin/metabrain_PCA/data/ROSMAP_DLPFC_biplot_pca_removedBatchEff_cov_12_2.png')
biplot_pca_removedBatchEff_cov
dev.off()




# eigencor_plot_removedBatchEff_cov <- eigencorplot(p_removedBatchEff_cov, metavars = tech_covar_list, cexLabY = 0.5, rotLabY = 0.8, corFUN = "pearson", main = "Correlation of PCs with technical covariates in ROSMAP and significancies", titleX = "PCs", titleY = "technical covariates", corMultipleTestCorrection = "hochberg", signifSymbols = c('***', '**', '*', ''), signifCutpoints = c(0, 0.001, 0.01, 0.05, 1), scale=FALSE)
# png('/external/rprshnas01/netdata_kcni/stlab/Xiaolin/metabrain_PCA/data/ROSMAP_DLPFC_eigencor_plot_removedBatchEff_cov.png', width = 15, height = 10, units = 'in', res = 300)
# eigencor_plot_removedBatchEff_cov
# dev.off()


eigencor_plot_removedBatchEff_cov <- eigencorplot(p_removedBatchEff_cov, metavars = tech_covar_list, cexLabY = 0.5, rotLabY = 0.8, corFUN = "pearson", main = "Correlation of PCs with technical covariates in ROSMAP and significancies", titleX = "PCs", titleY = "technical covariates", corMultipleTestCorrection = "hochberg", signifSymbols = c('***', '**', '*', ''), signifCutpoints = c(0, 0.001, 0.01, 0.05, 1), scale=FALSE)
png('/external/rprshnas01/netdata_kcni/stlab/Xiaolin/metabrain_PCA/data/ROSMAP_DLPFC_eigencor_plot_removedBatchEff_cov_12_2.png', width = 15, height = 10, units = 'in', res = 300)
eigencor_plot_removedBatchEff_cov
dev.off()


##########################################################################################################################
#Converting Ensembl ID to Gene Name
# mart = useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host = "https://dec2021.archive.ensembl.org/") 
mart = useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host = "https://useast.ensembl.org")
ensembl = getBM(attributes = c("ensembl_gene_id", "hgnc_symbol", "entrezgene_id","description","gene_biotype","percentage_gene_gc_content"), mart = mart)
ensembl_to_gene = (data.frame(ensembl$ensembl_gene_id, ensembl$hgnc_symbol))
names(ensembl_to_gene) = c("gene_symbol", "gene_name")
#remove duplicates
ensembl_to_gene = ensembl_to_gene[!duplicated(ensembl_to_gene[,1]),]
count_cpm_filtered <- as.data.frame(zscore_data_removedBatchEff_cov)
count_cpm_filtered$gene_symbol <- rownames(count_cpm_filtered)
count_cpm_filtered = merge(x=count_cpm_filtered, y=ensembl_to_gene, by = "gene_symbol", all.x = T)
# some gene names are duplicated after matching (‘’, ‘DUXAP8’, ‘GOLGA8M’, ‘ITFG2-AS1’, ‘LINC01238’, ‘PINX1’, ‘POLR2J4’, ‘RN7SL274P’, ‘SIGLEC5’, ‘TUBB7P’), in order to set them as unqiue row names, we only keep the first occurence of the duplicated gene names
count_cpm_filtered <- count_cpm_filtered[!duplicated(count_cpm_filtered$gene_name), ]
count_cpm_filtered <- count_cpm_filtered[count_cpm_filtered$gene_name != '', ] #remove the empty gene name
count_cpm_filtered <- count_cpm_filtered[!is.na(count_cpm_filtered$gene_name), ]
rownames(count_cpm_filtered) = count_cpm_filtered$gene_name #[1] 29977  1143
#count_cpm_filtered <- count_cpm_filtered[ , !names(count_cpm_filtered) %in% c("gene_symbol", "gene_name")] # [1] 29977 genes 1141 samples





##############################################################################################################################
### Cell-Type Proportion Estimation: Try updated marker list from Micaela's paper
### GET MARKERS FOR MGP ANALYSIS
# note that this is the list of markers from micaela's paper - you get similar but diff results if you use the markers from the aging paper
sonny_markers = read_csv(url('https://raw.githubusercontent.com/sonnyc247/MarkerSelection/master/Data/Outputs/CSVs_and_Tables/Markers/MTG_and_CgG_lfct2/new_MTGnCgG_lfct2.5_Publication.csv'))
colnames(sonny_markers) = colnames(sonny_markers) %>% make.names() %>% tolower()

# I find it helpful to map some gene symbols to ensembl ids manually using mappings from hgnc, you can get those from here: 
# http://ftp.ebi.ac.uk/pub/databases/genenames/hgnc/tsv/hgnc_complete_set.txt

hgnc_mapping = read_tsv('/external/rprshnas01/kcni/dkiss/cell_prop_psychiatry/data/hgnc_complete_set.txt')

# now, this is the list of sonnys markers with entrez ids and ensembl ids where possible
sonny_hgnc_merged_markers = left_join(sonny_markers %>% dplyr::rename(entrez_id = entrez.gene.id), 
                                      hgnc_mapping %>% distinct(entrez_id, .keep_all = T)%>% 
                                        dplyr::select(entrez_id, ensembl_gene_id) %>% 
                                        dplyr::rename(ensembl_id = ensembl_gene_id)) %>% 
  dplyr::select(gene, entrez_id, ensembl_id, -ensembl.gene.id, everything()) %>% 
  group_by(subclass) %>% 
  arrange(subclass, -average.log.fold.change) %>% 
  ungroup()

# get ensembl list of markers
new_markers = sonny_hgnc_merged_markers %>% filter(used.in.mgp == "TRUE")
new_cell_types = new_markers %>% filter(!is.na(subclass)) %>% pull(subclass) %>% unique
new_marker_list  = lapply(new_cell_types, function(cell_type){
  return(new_markers %>% filter(subclass == cell_type, 
                                ensembl_id %in% unique(c(count_cpm_filtered$gene_symbol)),
  ) %>% pull(ensembl_id))
})
names(new_marker_list) = c('Astrocyte', 'Endothelial', 'Exc_IT', 'Exc_L4_IT', 'Exc_L5_ET', 'Exc_L5/6_IT_Car3', 'Exc_L5/6_NP', 'Exc_L6_CT', 'Exc_L6b', 'Inh_LAMP5', 'Microglia', 'Oligodendrocyte', 'OPC', 'Inh_PAX6', 'Pericyte', 'Inh_PVALB', 'Inh_SST', 'Inh_VIP', 'VLMC')



##############################################################################################################################


#Remove ensembl_ID and move gene names to first column 
count_cpm_filtered_gene_symbol <- count_cpm_filtered
rownames(count_cpm_filtered_gene_symbol) = count_cpm_filtered_gene_symbol$gene_name
count_cpm_filtered_gene_symbol = count_cpm_filtered_gene_symbol[, -c(ncol(count_cpm_filtered_gene_symbol))] #remove the last column: hgnc gene_name
#rename the count_cpm_filtered ensembl id to gene name
colnames(count_cpm_filtered_gene_symbol)[1] = "gene_name" #gene_names is ensembl id.. wired
 #[1] 39006  1143
 #the rownames are hgnc gene names; and the first column is ensembl id



rosmap_estimations =  mgpEstimate(
  exprData = count_cpm_filtered_gene_symbol,
  genes = new_marker_list,
  geneColName = 'gene_name',
  outlierSampleRemove = F, # should outlier samples removed. This is done using boxplot stats
  geneTransform = NULL, # this is the default option for geneTransform
  groups = NULL, # if there are experimental groups provide them here. if not desired set to NULL
  seekConsensus = FALSE, # ensures gene rotations are positive in both of the groups
  removeMinority = F)

#Coerce estimations list into data frame 
rosmap_estimations_scaled = rosmap_estimations$estimates %>% as.data.frame() %>% scale() %>% as.data.frame() %>% tibble::rownames_to_column(var = "specimenID")

names(rosmap_estimations_scaled)[names(rosmap_estimations_scaled) == "specimenID"] = "synapseID"
saveRDS(rosmap_estimations_scaled, file="/external/rprshnas01/netdata_kcni/stlab/Xiaolin/metabrain_PCA/data/rosmap_estimations_scaled_techCovRemoved.RData")

### Merge cell type proportions with sample metadata
#Change column name for join
rosmap_estimations_metadata = right_join(rosmap_meta_DLPFC, rosmap_estimations_scaled, by="synapseID") #[1] 1142   59

#Remove '+' from ageDeath for modelling
rosmap_estimations_metadata$age_death = as.numeric(gsub("[+]", "", rosmap_estimations_metadata$age_death))

# change the identifer for each sample from synapseID into Study+projid since the cov.txt uses this identifier
# Remove duplicates for mapping specimenID to projid
matcher <- rosmap_meta_DLPFC[duplicated(rosmap_meta_DLPFC$synapseID) == FALSE, c("synapseID", "projid", "Study")] 
matcher <- matcher[duplicated(matcher$projid) == FALSE,] #only 1119 samples have unique projid, so 1142 samples with unique synapse ID will reduce to 1119 samples
# Combine columns and create a new column
matcher$combined <- paste(matcher$Study, matcher$projid, sep = "") 
rosmap_estimations_metadata_unique_id <- merge(rosmap_estimations_metadata, matcher, by.x = "synapseID", by.y = "synapseID", all.x = TRUE, all.y = TRUE) #[1] 1119   62
rosmap_estimations_metadata_unique_id <- rosmap_estimations_metadata_unique_id[!is.na(rosmap_estimations_metadata_unique_id$combined),] 

mgp <- rosmap_estimations_metadata_unique_id

mgp <- mgp[, ((ncol(mgp)-21):ncol(mgp))]
mgp <- mgp[, -c((ncol(mgp)-1), (ncol(mgp)-2))]
colnames(mgp)[ncol(mgp)] <- "FID"
mgp <- mgp[, c(ncol(mgp), 1:(ncol(mgp) - 1))]
mgp$IID <- mgp$FID
mgp <- mgp[, c("FID", "IID", names(mgp)[-c(1, ncol(mgp))])]

sn_proportions_raw <- read_csv("/external/rprshnas01/netdata_kcni/stlab/cross_cohort_MGPs/rosmap_single_nuc_proportions.csv")
meta_raw <- read_csv("/external/rprshnas01/external_data/rosmap/gene_expression/RNAseq_Harmonization/Gene Expression (Raw Gene Counts)/Metadata/RNAseq_Harmonization_ROSMAP_combined_metadata.csv") 
# Remove duplicates for mapping specimenID to projid
matcher <- meta_raw[duplicated(meta_raw$specimenID) == FALSE, c("specimenID", "projid", "Study", "msex", "age_death")]
matcher <- matcher[duplicated(matcher$projid) == FALSE,]
# Assign specimenID to sn proportions which only has projid
sn_proportions <- merge(sn_proportions_raw, matcher, by.x = "ID", by.y = "projid", all.x = TRUE, all.y = FALSE)
sn_proportions <- sn_proportions[!is.na(sn_proportions$specimenID),] 

names(sn_proportions)[1] <- "projid"
sn_proportions$Sample <- paste0(sn_proportions$Study, sn_proportions$projid)

# sn_proportions_ordered <- sn_proportions[,c(ncol(sn_proportions),3,4,5,6,8,10,11,16,17,18,19,20,21)]
# #re-order/subset df by common cell types
# mgp_DLPFC_df <- mgp
# mgp_DLPFC_ordered <- mgp_DLPFC_df[,c(1,3,4,18,20,19,16,12,8,9,11,10,7,6)]


# colnames(sn_proportions_ordered) <- colnames(mgp_DLPFC_ordered)

# #match FID order in sn_proportions_ordered and mgp_DLPFC_ordered
# sn_proportions_ordered <- sn_proportions_ordered[match(mgp_DLPFC_ordered$FID, sn_proportions_ordered$FID),]

# re-order sn_proportions cell types
sn_proportions_ordered <- sn_proportions[,c(ncol(sn_proportions),3,4,5,6,8,10,11,16,17,18,19,20,21, 28, 29)]
#re-order/subset df by common cell types
mgp_DLPFC_ordered <- mgp_DLPFC_df[,c(1,4,5,19,21,20,17,13,9,10,12,11,8,7)]
#unify the cell types names
colnames(sn_proportions_ordered)[1:(ncol(sn_proportions_ordered)-2)] <- colnames(mgp_DLPFC_ordered)

sn_proportions_ordered$age_death <- gsub("90+", "91", sn_proportions_ordered$age_death)
sn_proportions_ordered$age_death = as.numeric(gsub("[+]", "", sn_proportions_ordered$age_death))

sn_proportions_ordered$msex <- as.numeric(sn_proportions_ordered$msex)

sn_proportions_ordered$age_death_sex <- sn_proportions_ordered$age_death*sn_proportions_ordered$msex
sn_proportions_ordered$age_death2 <- sn_proportions_ordered$age_death*sn_proportions_ordered$age_death
sn_proportions_ordered$age_death2_sex <- sn_proportions_ordered$age_death*sn_proportions_ordered$age_death*sn_proportions_ordered$msex #412 19



### 10.3 update: remove samples with <500 nuclei to remove outliers (ie: one sample has 100% endo) in any cell type
cell_counts <- read.csv("/external/rprshnas01/netdata_kcni/stlab/rosmap_p400/p400_cell_type_counts.tsv", sep = "\t", header=TRUE)
#remove duplicated cell counts: ie-Astro.2
cell_counts_cleaned <- cell_counts[!grepl("\\.", cell_counts$cell_type),]
# print the outlier sample with one cell type takes up too much proportion
unique_sample_length <- length(unique(cell_counts_cleaned$ID))
unique_sample <- unique(cell_counts_cleaned$ID)
for (i in 1:unique_sample_length) {
  sample <- unique_sample[i]
  unique_cell_type <- unique(cell_counts_cleaned[cell_counts_cleaned$ID == sample,]$cell_type)
  total_cell_counts_per_sample <- sum((cell_counts_cleaned[cell_counts_cleaned$ID == sample,])$num_cells)
  for ( cell_type in unique_cell_type){
    #print(unique_cell_type[j])
    if (((cell_counts_cleaned[cell_counts_cleaned$ID == sample & cell_counts_cleaned$cell_type == cell_type,])$num_cells)/total_cell_counts_per_sample > 0.8) { #one cell type takes up too much proportion
      cat("sample: ", sample)
      cat("cell type: ", cell_type)
  }
  }
}
# outlier: cell_counts_cleaned[cell_counts_cleaned$ID == "21232244" & cell_counts_cleaned$cell_type == "Endo",] = 1

#remove outlier sample with nuclei < 500
outlier_sample <- c()
outlier_cell_type <- c()
for (i in 1:unique_sample_length) {
  sample <- unique_sample[i]
  unique_cell_type <- unique(cell_counts_cleaned[cell_counts_cleaned$ID == sample,]$cell_type)
  total_cell_counts_per_sample <- sum((cell_counts_cleaned[cell_counts_cleaned$ID == sample,])$num_cells)
  if (sum(cell_counts_cleaned[cell_counts_cleaned$ID == sample,]$num_cells) < 500) { #the sum of cell_counts per sample is >500
    cat("sample: ", sample)
    cat("cell type: ", cell_type)
    outlier_sample <- c(outlier_sample, sample)
    outlier_cell_type <- c(outlier_cell_type, unique_cell_type)
  }
}
# sample:  21232244cell type:  Endo
cell_counts_cleaned_500 <- cell_counts_cleaned[!(cell_counts_cleaned$ID %in% outlier_sample),] # remove 1, 10616 samples left
sn_proportions_ordered <- sn_proportions_ordered[!(str_extract(sn_proportions_ordered$FID, "\\d+") %in% outlier_sample),] #remove 1, 411 samples left





### regress sex and age for sn proportions
combined_df <- sn_proportions_ordered
cell_types <- colnames(mgp_DLPFC_ordered)[2:ncol(mgp_DLPFC_ordered)]
pheno_lms = lapply(cell_types, function(cell_type){
  lm = paste0("scale(", cell_type, ")", " ~ msex + scale(age_death) + scale(age_death_sex) + scale(age_death2) + scale(age_death2_sex)")
  results = lm(lm, data = combined_df) %>% 
  residuals()
  results <- results[1:nrow(combined_df)]  %>%
  RNOmni::RankNorm() %>% 
  as.data.frame()
  results$cell_type = cell_type
  results$FID = combined_df$Sample[1:nrow(combined_df)]
  return(results)
}) %>% bind_rows()


#format pheno_lms into phenotype dataframe
names(pheno_lms)[1] <- "transformed_residuals"

# Pivot the dataframe
pivot_df <- pivot_wider(data = pheno_lms, 
                        names_from = "cell_type", 
                        values_from = "transformed_residuals") # 411 samples 14 cell types




# 2.0 draw scatter plots for each cell type 
#### MGP estimated cell type prop from raw count matrix V.S. Actual Proportions after regressing out sex and age and removing outliers
predicted_df <- data.frame(mgp_DLPFC_ordered) #767
actual_df <- data.frame(sn_proportions_ordered) #411

# Merge the predicted and actual data frames
merged_df <- merge(predicted_df, actual_df, by = "FID") #[1] 334  27

# Create an empty list to store the scatter plots
plots_list_sv1 <- list()
#accuracy_list <- c()
for (cell_type in colnames(predicted_df)[2:ncol(predicted_df)]) {
  # Calculate the prediction accuracy (e.g., correlation coefficient)
  x <- paste(cell_type, ".x", sep = "")
  y <- paste(cell_type, ".y", sep = "")
  accuracy <- cor(merged_df[[x]], merged_df[[y]], use = "pairwise.complete.obs")

   #accuracy_list <- c(accuracy_list, accuracy)
  #Create the scatter plot with a trend line
  plot <- ggplot(merged_df, aes(x = .data[[y]], y = .data[[x]])) +
    geom_point() +
    geom_smooth(method = "lm", se = FALSE, color = "red") +  # Add a trend line
    labs(x = "Single-nucleus Proportions", y = "MGP from Raw") +
    ggtitle(paste("Estimation Accuracy for", cell_type, ":", round(accuracy, 2)))

  # Store the scatter plot in the list
  plots_list_sv1[[cell_type]] <- plot
}


#create a pdf file storing the results
pdf_file <- "/nethome/kcni/xzhou/GWAS_tut/ROSMAP/MGP_techCovRemoved_vs_snProp_sexAgeAdj.pdf"

# Open the PDF file with specified page size
pdf(pdf_file, width = 10, height = 10)

# Iterate over each pair of plots and arrange them side by side
for (i in 1:length(plots_list_sv1)) {
  # Arrange the i-th plot from each list side by side
  grid.arrange(plots_list_sv1[[i]])
    # Add a page break after each iteration (except the last)
  if (i < length(plots_list_sv1)) {
    cat("\f")  # Page break character
  }
}

# Close the PDF device
dev.off()






#####################################################################
#read ROSMAP metadata for DLPFC
zscore_data <- readRDS('/external/rprshnas01/netdata_kcni/stlab/Xiaolin/metabrain_PCA/data/ROSMAP_DLPFC_zscore_data.RData')
rosmap_meta_DLPFC <- readRDS('/external/rprshnas01/netdata_kcni/stlab/Xiaolin/metabrain_PCA/data/rosmap_meta_DLPFC.RData')

list <- c("RnaSeqMetrics__MEDIAN_CV_COVERAGE", "RnaSeqMetrics__PCT_RIBOSOMAL_BASES", "RnaSeqMetrics__PCT_CODING_BASES", "RnaSeqMetrics__PCT_UTR_BASES",
                           "RnaSeqMetrics__MEDIAN_5PRIME_TO_3PRIME_BIAS", "AlignmentSummaryMetrics__PCT_PF_READS_ALIGNED", "Study", "RnaSeqMetrics__MEDIAN_3PRIME_BIAS", "RnaSeqMetrics__PCT_INTERGENIC_BASES")
# metadata <- rosmap_meta_DLPFC
# for (tech_cov in list){
#   if(is.numeric(metadata[[tech_cov]])){
#     metadata[[tech_cov]] <- as.numeric(metadata[[tech_cov]])}
#   else{
#     metadata[[tech_cov]] <- factor(metadata[[tech_cov]])
#   }
# }

# try Stuart's design matrix for removing tech cov
#tech cov not exist in ROSMAP meta: LOG_ESTIMATED_LIBRARY_SIZE; PERCENT_DUPLICATION;  apoe4d; smoking; alcohol_g_bl; cesdsum_lv
# we mght not need this in our design since we only care about technical differences: apoe4d; smoking; alcohol_g_bl; cesdsum_lv

design3 <- model.matrix( ~ RnaSeqMetrics__MEDIAN_CV_COVERAGE + RnaSeqMetrics__PCT_RIBOSOMAL_BASES + 
                           RnaSeqMetrics__PCT_CODING_BASES + RnaSeqMetrics__PCT_UTR_BASES + 
                           #we dont have this
                           RnaSeqMetrics__MEDIAN_5PRIME_TO_3PRIME_BIAS + 
                           AlignmentSummaryMetrics__PCT_PF_READS_ALIGNED + Study + 
                           RnaSeqMetrics__MEDIAN_3PRIME_BIAS + RnaSeqMetrics__PCT_INTERGENIC_BASES
                           + msex + pmi + age_death + educ, #apoe4d + smoking + alcohol_g_bl + cesdsum_lv
                           data=rosmap_meta_DLPFC)


#removed_sample <- colnames(zscore_data)[!colnames(zscore_data) %in% rownames(design3)] #"syn4213070" "syn4213093"
zscore_data_updated <- zscore_data[, colnames(zscore_data) %in% rownames(design3)]
zscore_data_design3 <- removeBatchEffect(x=zscore_data_removedBatchEff, batch=as.vector(metadata[,"sequencingBatch"]), covariates=metadata_nobatch[,colnames(metadata_nobatch) %in% list]) #, design=design3
#saveRDS(zscore_data_design3, file="/external/rprshnas01/netdata_kcni/stlab/Xiaolin/metabrain_PCA/data/ROSMAP_DLPFC_zscore_data_design3.RData")
#zscore_data_design3 <- readRDS('/external/rprshnas01/netdata_kcni/stlab/Xiaolin/metabrain_PCA/data/ROSMAP_DLPFC_zscore_data_design3.RData')

p_removedBatchEff_cov <- PCAtools::pca(zscore_data_design3, metadata = metadata)
biplot_pca_removedBatchEff_cov <- biplot(p_removedBatchEff_cov)

#save into png file
png('/external/rprshnas01/netdata_kcni/stlab/Xiaolin/metabrain_PCA/data/ROSMAP_DLPFC_biplot_pca_removedBatchEff_cov_stuart.png')
biplot_pca_removedBatchEff_cov
dev.off()

eigencor_plot_removedBatchEff_cov <- eigencorplot(p_removedBatchEff_cov , metavars = tech_covar_list, cexLabY = 0.5, rotLabY = 0.8, corFUN = "pearson", main = "Correlation of PCs with technical covariates in ROSMAP and significancies", titleX = "PCs", titleY = "technical covariates", corMultipleTestCorrection = "hochberg", signifSymbols = c('***', '**', '*', ''), signifCutpoints = c(0, 0.001, 0.01, 0.05, 1), scale=FALSE)
png('/external/rprshnas01/netdata_kcni/stlab/Xiaolin/metabrain_PCA/data/ROSMAP_DLPFC_eigencor_plot_removedBatchEff_cov_stuart.png', width = 15, height = 10, units = 'in', res = 300)
eigencor_plot_removedBatchEff_cov
dev.off()


##########################################################################################################################
# MGP cell type estimation for Stuart design matrix

##########################################################################################################################
#Converting Ensembl ID to Gene Name
# mart = useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host = "https://dec2021.archive.ensembl.org/") 
mart = useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host = "https://useast.ensembl.org")
ensembl = getBM(attributes = c("ensembl_gene_id", "hgnc_symbol", "entrezgene_id","description","gene_biotype","percentage_gene_gc_content"), mart = mart)
ensembl_to_gene = (data.frame(ensembl$ensembl_gene_id, ensembl$hgnc_symbol))
names(ensembl_to_gene) = c("gene_symbol", "gene_name")
#remove duplicates
ensembl_to_gene = ensembl_to_gene[!duplicated(ensembl_to_gene[,1]),]
count_cpm_filtered <- as.data.frame(zscore_data_design3)
count_cpm_filtered$gene_symbol <- rownames(count_cpm_filtered)
count_cpm_filtered = merge(x=count_cpm_filtered, y=ensembl_to_gene, by = "gene_symbol", all.x = T)
# some gene names are duplicated after matching (‘’, ‘DUXAP8’, ‘GOLGA8M’, ‘ITFG2-AS1’, ‘LINC01238’, ‘PINX1’, ‘POLR2J4’, ‘RN7SL274P’, ‘SIGLEC5’, ‘TUBB7P’), in order to set them as unqiue row names, we only keep the first occurence of the duplicated gene names
count_cpm_filtered <- count_cpm_filtered[!duplicated(count_cpm_filtered$gene_name), ]
count_cpm_filtered <- count_cpm_filtered[count_cpm_filtered$gene_name != '', ] #remove the empty gene name
count_cpm_filtered <- count_cpm_filtered[!is.na(count_cpm_filtered$gene_name), ]
rownames(count_cpm_filtered) = count_cpm_filtered$gene_name #[1] 29977  1143

sonny_markers = read_csv(url('https://raw.githubusercontent.com/sonnyc247/MarkerSelection/master/Data/Outputs/CSVs_and_Tables/Markers/MTG_and_CgG_lfct2/new_MTGnCgG_lfct2.5_Publication.csv'))
colnames(sonny_markers) = colnames(sonny_markers) %>% make.names() %>% tolower()

# I find it helpful to map some gene symbols to ensembl ids manually using mappings from hgnc, you can get those from here: 
# http://ftp.ebi.ac.uk/pub/databases/genenames/hgnc/tsv/hgnc_complete_set.txt

hgnc_mapping = read_tsv('/external/rprshnas01/kcni/dkiss/cell_prop_psychiatry/data/hgnc_complete_set.txt')

# now, this is the list of sonnys markers with entrez ids and ensembl ids where possible
sonny_hgnc_merged_markers = left_join(sonny_markers %>% dplyr::rename(entrez_id = entrez.gene.id), 
                                      hgnc_mapping %>% distinct(entrez_id, .keep_all = T)%>% 
                                        dplyr::select(entrez_id, ensembl_gene_id) %>% 
                                        dplyr::rename(ensembl_id = ensembl_gene_id)) %>% 
  dplyr::select(gene, entrez_id, ensembl_id, -ensembl.gene.id, everything()) %>% 
  group_by(subclass) %>% 
  arrange(subclass, -average.log.fold.change) %>% 
  ungroup()

# get ensembl list of markers
new_markers = sonny_hgnc_merged_markers %>% filter(used.in.mgp == "TRUE")
new_cell_types = new_markers %>% filter(!is.na(subclass)) %>% pull(subclass) %>% unique
new_marker_list  = lapply(new_cell_types, function(cell_type){
  return(new_markers %>% filter(subclass == cell_type, 
                                ensembl_id %in% unique(c(count_cpm_filtered$gene_symbol)),
  ) %>% pull(ensembl_id))
})
names(new_marker_list) = c('Astrocyte', 'Endothelial', 'Exc_IT', 'Exc_L4_IT', 'Exc_L5_ET', 'Exc_L5/6_IT_Car3', 'Exc_L5/6_NP', 'Exc_L6_CT', 'Exc_L6b', 'Inh_LAMP5', 'Microglia', 'Oligodendrocyte', 'OPC', 'Inh_PAX6', 'Pericyte', 'Inh_PVALB', 'Inh_SST', 'Inh_VIP', 'VLMC')

##############################################################################################################################


#Remove ensembl_ID and move gene names to first column 
count_cpm_filtered_gene_symbol <- count_cpm_filtered
rownames(count_cpm_filtered_gene_symbol) = count_cpm_filtered_gene_symbol$gene_name
count_cpm_filtered_gene_symbol = count_cpm_filtered_gene_symbol[, -c(ncol(count_cpm_filtered_gene_symbol))] #remove the last column: hgnc gene_name
#rename the count_cpm_filtered ensembl id to gene name
colnames(count_cpm_filtered_gene_symbol)[1] = "gene_name" #gene_names is ensembl id.. wired
 #[1] 39006  1143
 #the rownames are hgnc gene names; and the first column is ensembl id



rosmap_estimations =  mgpEstimate(
  exprData = count_cpm_filtered_gene_symbol,
  genes = new_marker_list,
  geneColName = 'gene_name',
  outlierSampleRemove = F, # should outlier samples removed. This is done using boxplot stats
  geneTransform = NULL, # this is the default option for geneTransform
  groups = NULL, # if there are experimental groups provide them here. if not desired set to NULL
  seekConsensus = FALSE, # ensures gene rotations are positive in both of the groups
  removeMinority = F)

#Coerce estimations list into data frame 
rosmap_estimations_scaled = rosmap_estimations$estimates %>% as.data.frame() %>% scale() %>% as.data.frame() %>% tibble::rownames_to_column(var = "specimenID")

names(rosmap_estimations_scaled)[names(rosmap_estimations_scaled) == "specimenID"] = "synapseID"
#saveRDS(rosmap_estimations_scaled, file="/external/rprshnas01/netdata_kcni/stlab/Xiaolin/metabrain_PCA/data/rosmap_estimations_scaled_stuart.RData")
mgp_est_stuart <- readRDS('/external/rprshnas01/netdata_kcni/stlab/Xiaolin/metabrain_PCA/data/rosmap_estimations_scaled_stuart.RData')

matcher <- rosmap_meta_DLPFC[duplicated(rosmap_meta_DLPFC$synapseID) == FALSE, c("synapseID", "projid", "Study")] 
matcher <- matcher[duplicated(matcher$projid) == FALSE,] #only 1119 samples have unique projid, so 1142 samples with unique synapse ID will reduce to 1119 samples
# Combine columns and create a new column
matcher$combined <- paste(matcher$Study, matcher$projid, sep = "") 
mgp_est_stuart_unique_id <- merge(mgp_est_stuart, matcher, by.x = "synapseID", by.y = "synapseID", all.x = TRUE, all.y = TRUE) #[1] 1119   62
mgp_est_stuart_unique_id <- mgp_est_stuart_unique_id[!is.na(mgp_est_stuart_unique_id$combined),] 

mgp_est_stuart <- mgp_est_stuart_unique_id
mgp_est_stuart$synapseID <- mgp_est_stuart$combined
colnames(mgp_est_stuart)[1] <- "FID"

mgp_est_stuart_ordered <- mgp_est_stuart[,c(1,2,3,17,19,18,15,11,7,8,10,9,6,5,12,14,13)]
rownames(mgp_est_stuart_ordered) <- mgp_est_stuart_ordered$FID
mgp_est_stuart_ordered <- mgp_est_stuart_ordered[,-c(1)]





# re-order sn_proportions cell types
sn_proportions_ordered <- sn_proportions[,c(ncol(sn_proportions),3,4,5,6,8,10,11,16,17,18,19,20,21,23,24,25)]
#re-order/subset df by common cell types
mgp_DLPFC_df <- mgp 
mgp_DLPFC_ordered <- mgp_DLPFC_df[,c(1,3,4,18,20,19,16,12,8,9,11,10,7,6,13,15,14)]
#unify the cell types names
colnames(sn_proportions_ordered)[1:(ncol(sn_proportions_ordered)-2)] <- colnames(mgp_DLPFC_ordered)

sn_proportions_ordered$age_death <- gsub("90+", "91", sn_proportions_ordered$age_death)
sn_proportions_ordered$age_death = as.numeric(gsub("[+]", "", sn_proportions_ordered$age_death))

sn_proportions_ordered$msex <- as.numeric(sn_proportions_ordered$msex)

sn_proportions_ordered$age_death_sex <- sn_proportions_ordered$age_death*sn_proportions_ordered$msex
sn_proportions_ordered$age_death2 <- sn_proportions_ordered$age_death*sn_proportions_ordered$age_death
sn_proportions_ordered$age_death2_sex <- sn_proportions_ordered$age_death*sn_proportions_ordered$age_death*sn_proportions_ordered$msex #412 19



### 10.3 update: remove samples with <500 nuclei to remove outliers (ie: one sample has 100% endo) in any cell type
cell_counts <- read.csv("/external/rprshnas01/netdata_kcni/stlab/rosmap_p400/p400_cell_type_counts.tsv", sep = "\t", header=TRUE)
#remove duplicated cell counts: ie-Astro.2
cell_counts_cleaned <- cell_counts[!grepl("\\.", cell_counts$cell_type),]
# print the outlier sample with one cell type takes up too much proportion
unique_sample_length <- length(unique(cell_counts_cleaned$ID))
unique_sample <- unique(cell_counts_cleaned$ID)
for (i in 1:unique_sample_length) {
  sample <- unique_sample[i]
  unique_cell_type <- unique(cell_counts_cleaned[cell_counts_cleaned$ID == sample,]$cell_type)
  total_cell_counts_per_sample <- sum((cell_counts_cleaned[cell_counts_cleaned$ID == sample,])$num_cells)
  for ( cell_type in unique_cell_type){
    #print(unique_cell_type[j])
    if (((cell_counts_cleaned[cell_counts_cleaned$ID == sample & cell_counts_cleaned$cell_type == cell_type,])$num_cells)/total_cell_counts_per_sample > 0.8) { #one cell type takes up too much proportion
      cat("sample: ", sample)
      cat("cell type: ", cell_type)
  }
  }
}
# outlier: cell_counts_cleaned[cell_counts_cleaned$ID == "21232244" & cell_counts_cleaned$cell_type == "Endo",] = 1

#remove outlier sample with nuclei < 500
outlier_sample <- c()
outlier_cell_type <- c()
for (i in 1:unique_sample_length) {
  sample <- unique_sample[i]
  unique_cell_type <- unique(cell_counts_cleaned[cell_counts_cleaned$ID == sample,]$cell_type)
  total_cell_counts_per_sample <- sum((cell_counts_cleaned[cell_counts_cleaned$ID == sample,])$num_cells)
  if (sum(cell_counts_cleaned[cell_counts_cleaned$ID == sample,]$num_cells) < 500) { #the sum of cell_counts per sample is >500
    cat("sample: ", sample)
    cat("cell type: ", cell_type)
    outlier_sample <- c(outlier_sample, sample)
    outlier_cell_type <- c(outlier_cell_type, unique_cell_type)
  }
}
# sample:  21232244cell type:  Endo
cell_counts_cleaned_500 <- cell_counts_cleaned[!(cell_counts_cleaned$ID %in% outlier_sample),] # remove 1, 10616 samples left
sn_proportions_ordered <- sn_proportions_ordered[!(str_extract(sn_proportions_ordered$FID, "\\d+") %in% outlier_sample),] #remove 1, 411 samples left


###################################################################################################################################
# MGP for top20 tech cov regressed

#Converting Ensembl ID to Gene Name
# mart = useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host = "https://dec2021.archive.ensembl.org/") 
mart = useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host = "https://useast.ensembl.org")
ensembl = getBM(attributes = c("ensembl_gene_id", "hgnc_symbol", "entrezgene_id","description","gene_biotype","percentage_gene_gc_content"), mart = mart)
ensembl_to_gene = (data.frame(ensembl$ensembl_gene_id, ensembl$hgnc_symbol))
names(ensembl_to_gene) = c("gene_symbol", "gene_name")
#remove duplicates
ensembl_to_gene = ensembl_to_gene[!duplicated(ensembl_to_gene[,1]),]
zscore_data_top20 <- as.data.frame(zscore_data_removedBatchEff_cov)
zscore_data_top20$gene_symbol <- rownames(zscore_data_top20)
zscore_data_top20 = merge(x=zscore_data_top20, y=ensembl_to_gene, by = "gene_symbol", all.x = T)
# some gene names are duplicated after matching (‘’, ‘DUXAP8’, ‘GOLGA8M’, ‘ITFG2-AS1’, ‘LINC01238’, ‘PINX1’, ‘POLR2J4’, ‘RN7SL274P’, ‘SIGLEC5’, ‘TUBB7P’), in order to set them as unqiue row names, we only keep the first occurence of the duplicated gene names
zscore_data_top20 <- zscore_data_top20[!duplicated(zscore_data_top20$gene_name), ]
zscore_data_top20 <- zscore_data_top20[zscore_data_top20$gene_name != '', ] #remove the empty gene name
zscore_data_top20 <- zscore_data_top20[!is.na(zscore_data_top20$gene_name), ]
rownames(zscore_data_top20) = zscore_data_top20$gene_name #[1] 29977  1143

sonny_markers = read_csv(url('https://raw.githubusercontent.com/sonnyc247/MarkerSelection/master/Data/Outputs/CSVs_and_Tables/Markers/MTG_and_CgG_lfct2/new_MTGnCgG_lfct2.5_Publication.csv'))
colnames(sonny_markers) = colnames(sonny_markers) %>% make.names() %>% tolower()

# I find it helpful to map some gene symbols to ensembl ids manually using mappings from hgnc, you can get those from here: 
# http://ftp.ebi.ac.uk/pub/databases/genenames/hgnc/tsv/hgnc_complete_set.txt

hgnc_mapping = read_tsv('/external/rprshnas01/kcni/dkiss/cell_prop_psychiatry/data/hgnc_complete_set.txt')

# now, this is the list of sonnys markers with entrez ids and ensembl ids where possible
sonny_hgnc_merged_markers = left_join(sonny_markers %>% dplyr::rename(entrez_id = entrez.gene.id), 
                                      hgnc_mapping %>% distinct(entrez_id, .keep_all = T)%>% 
                                        dplyr::select(entrez_id, ensembl_gene_id) %>% 
                                        dplyr::rename(ensembl_id = ensembl_gene_id)) %>% 
  dplyr::select(gene, entrez_id, ensembl_id, -ensembl.gene.id, everything()) %>% 
  group_by(subclass) %>% 
  arrange(subclass, -average.log.fold.change) %>% 
  ungroup()

# get ensembl list of markers
new_markers = sonny_hgnc_merged_markers %>% filter(used.in.mgp == "TRUE")
new_cell_types = new_markers %>% filter(!is.na(subclass)) %>% pull(subclass) %>% unique
new_marker_list  = lapply(new_cell_types, function(cell_type){
  return(new_markers %>% filter(subclass == cell_type, 
                                ensembl_id %in% unique(c(zscore_data_top20$gene_symbol)),
  ) %>% pull(ensembl_id))
})
names(new_marker_list) = c('Astrocyte', 'Endothelial', 'Exc_IT', 'Exc_L4_IT', 'Exc_L5_ET', 'Exc_L5/6_IT_Car3', 'Exc_L5/6_NP', 'Exc_L6_CT', 'Exc_L6b', 'Inh_LAMP5', 'Microglia', 'Oligodendrocyte', 'OPC', 'Inh_PAX6', 'Pericyte', 'Inh_PVALB', 'Inh_SST', 'Inh_VIP', 'VLMC')

##############################################################################################################################


#Remove ensembl_ID and move gene names to first column 
zscore_data_top20_gene_symbol <- zscore_data_top20
rownames(zscore_data_top20_gene_symbol) = zscore_data_top20_gene_symbol$gene_name
zscore_data_top20_gene_symbol = zscore_data_top20_gene_symbol[, -c(ncol(zscore_data_top20_gene_symbol))] #remove the last column: hgnc gene_name
#rename the zscore_data_top20 ensembl id to gene name
colnames(zscore_data_top20_gene_symbol)[1] = "gene_name" 

rosmap_estimations_top20 =  mgpEstimate(
  exprData = zscore_data_top20_gene_symbol,
  genes = new_marker_list,
  geneColName = 'gene_name',
  outlierSampleRemove = F, # should outlier samples removed. This is done using boxplot stats
  geneTransform = NULL, # this is the default option for geneTransform
  groups = NULL, # if there are experimental groups provide them here. if not desired set to NULL
  seekConsensus = FALSE, # ensures gene rotations are positive in both of the groups
  removeMinority = F)

#Coerce estimations list into data frame 
rosmap_estimations_scaled_top20 = rosmap_estimations_top20$estimates %>% as.data.frame() %>% scale() %>% as.data.frame() %>% tibble::rownames_to_column(var = "specimenID")

names(rosmap_estimations_scaled_top20)[names(rosmap_estimations_scaled_top20) == "specimenID"] = "synapseID"

rosmap_estimations_metadata = right_join(rosmap_meta_DLPFC, rosmap_estimations_scaled_top20, by="synapseID") #[1] 1142   59

#Remove '+' from ageDeath for modelling
rosmap_estimations_metadata$age_death = as.numeric(gsub("[+]", "", rosmap_estimations_metadata$age_death))

# change the identifer for each sample from synapseID into Study+projid since the cov.txt uses this identifier
# Remove duplicates for mapping specimenID to projid
matcher <- rosmap_meta_DLPFC[duplicated(rosmap_meta_DLPFC$synapseID) == FALSE, c("synapseID", "projid", "Study")] 
matcher <- matcher[duplicated(matcher$projid) == FALSE,] #only 1119 samples have unique projid, so 1142 samples with unique synapse ID will reduce to 1119 samples
# Combine columns and create a new column
matcher$combined <- paste(matcher$Study, matcher$projid, sep = "") 
rosmap_estimations_metadata_unique_id <- merge(rosmap_estimations_metadata, matcher, by.x = "synapseID", by.y = "synapseID", all.x = TRUE, all.y = TRUE) #[1] 1119   62
rosmap_estimations_metadata_unique_id <- rosmap_estimations_metadata_unique_id[!is.na(rosmap_estimations_metadata_unique_id$combined),] 



sn_proportions_ordered <- data.frame(sn_proportions[,c(ncol(sn_proportions),3,4,5,6,8,10,11,16,17,18,19,20,21,23,24,25)])
# remove bad seq sample in sn
sn_proportions_ordered <- sn_proportions_ordered[!(str_extract(sn_proportions_ordered$Sample, "\\d+") %in% outlier_sample),] 
rownames(sn_proportions_ordered) <- sn_proportions_ordered$Sample
sn_proportions_ordered <- sn_proportions_ordered[,-c(1)]



mgp_top20 <- rosmap_estimations_metadata_unique_id
mgp_top20 <- mgp_top20[, ((ncol(mgp_top20)-21):ncol(mgp_top20))]
mgp_top20 <- mgp_top20[, -c((ncol(mgp_top20)-1), (ncol(mgp_top20)-2))]
colnames(mgp_top20)[ncol(mgp_top20)] <- "FID"
mgp_top20 <- mgp_top20[, c(ncol(mgp_top20), 1:(ncol(mgp_top20) - 1))]
mgp_top20$IID <- mgp_top20$FID
mgp_top20 <- mgp_top20[, c("FID", "IID", names(mgp_top20)[-c(1, ncol(mgp_top20))])]
mgp_top20_ordered <- mgp_top20[,c(1,3,4,18,20,19,16,12,8,9,11,10,7,6,13,15,14)]
common_samples <- intersect(sn_proportions_ordered$FID, mgp_top20_ordered$FID) #408
rownames(mgp_top20_ordered) <- mgp_top20_ordered$FID
mgp_top20_ordered <- mgp_top20_ordered[,-c(1)]
mgp_top20_ordered <- mgp_top20_ordered[common_samples,]

colnames(sn_proportions_ordered) <- colnames(mgp_top20_ordered)
sn_proportions_ordered <- sn_proportions_ordered[common_samples,] #[1] 408  16




#load cell type estimation from gene counts with minimum cleaning
mgp_est_min <- readRDS("/external/rprshnas01/netdata_kcni/stlab/Xiaolin/cell_deconv_data/data/rosmap_estimations_scaled.RData")
rosmap_estimations_metadata = right_join(rosmap_meta_DLPFC, mgp_est_min, by="synapseID") #[1] 1142   59

#Remove '+' from ageDeath for modelling
rosmap_estimations_metadata$age_death = as.numeric(gsub("[+]", "", rosmap_estimations_metadata$age_death))

# change the identifer for each sample from synapseID into Study+projid since the cov.txt uses this identifier
# Remove duplicates for mapping specimenID to projid
matcher <- rosmap_meta_DLPFC[duplicated(rosmap_meta_DLPFC$synapseID) == FALSE, c("synapseID", "projid", "Study")] 
matcher <- matcher[duplicated(matcher$projid) == FALSE,] #only 1119 samples have unique projid, so 1142 samples with unique synapse ID will reduce to 1119 samples
# Combine columns and create a new column
matcher$combined <- paste(matcher$Study, matcher$projid, sep = "") 
rosmap_estimations_metadata_unique_id <- merge(rosmap_estimations_metadata, matcher, by.x = "synapseID", by.y = "synapseID", all.x = TRUE, all.y = TRUE) #[1] 1119   62
rosmap_estimations_metadata_unique_id <- rosmap_estimations_metadata_unique_id[!is.na(rosmap_estimations_metadata_unique_id$combined),] 

mgp_est_min <- rosmap_estimations_metadata_unique_id
mgp_est_min <- mgp_est_min[, ((ncol(mgp_est_min)-21):ncol(mgp_est_min))]
mgp_est_min <- mgp_est_min[, -c((ncol(mgp_est_min)-1), (ncol(mgp_est_min)-2))]
colnames(mgp_est_min)[ncol(mgp_est_min)] <- "FID"
mgp_est_min <- mgp_est_min[, c(ncol(mgp_est_min), 1:(ncol(mgp_est_min) - 1))]
mgp_est_min$IID <- mgp_est_min$FID
mgp_est_min <- mgp_est_min[, c("FID", "IID", names(mgp_est_min)[-c(1, ncol(mgp_est_min))])]
mgp_est_min_ordered <- mgp_est_min[,c(1,3,4,18,20,19,16,12,8,9,11,10,7,6,13,15,14)]
rownames(mgp_est_min_ordered) <- mgp_est_min_ordered$FID
mgp_est_min_ordered <- mgp_est_min_ordered[,-c(1)]
mgp_est_min_ordered <- mgp_est_min_ordered[common_samples,]

mgp_est_stuart_ordered <- mgp_est_stuart_ordered[common_samples,] #[1] 408  16












############################################ Plotting Predicted V.S. Actual Proportions Starts Here ##########################################
# Create data frames for the predicted and actual cell type proportions
cor_func <- function(estimations_df, actual_prop){
  predicted_df <- data.frame(Sample = rownames(estimations_df), estimations_df)
  actual_df <- data.frame(Sample = rownames(actual_prop), actual_prop)

  # Merge the predicted and actual data frames
  merged_df <- merge(predicted_df, actual_df, by = "Sample")

  accuracy_list <- c()
  for (cell_type in colnames(estimations_df)) {
  # Calculate the prediction accuracy (e.g., correlation coefficient)
    x <- paste(cell_type, ".x", sep = "")
    y <- paste(cell_type, ".y", sep = "")
    #print(head(merged_df[[x]]))
    #print(head(merged_df[[y]]))
    accuracy <- cor(merged_df[[x]], merged_df[[y]], use = "pairwise.complete.obs")
    accuracy_list <- c(accuracy_list, accuracy)
  }
  return(accuracy_list)

}

mgp_min_sn_cor <- cor_func(mgp_est_min_ordered, sn_proportions_ordered)
mgp_top20_sn_cor <- cor_func(mgp_top20_ordered, sn_proportions_ordered)
mgp_stuart_sn_cor <- cor_func(mgp_est_stuart_ordered, sn_proportions_ordered)

accuracy_list <- c(mgp_min_sn_cor, mgp_top20_sn_cor, mgp_stuart_sn_cor)

#####
categories <- c("no_cleaning", "top20_tech_cov", "Stuart_design_matrix")
names <- colnames(mgp_est_stuart_ordered)

category_vector <- rep(categories, each = 16)
name_vector <- rep(names, 3)
df <- data.frame(Name = name_vector, Category = category_vector, Accuracy = accuracy_list)

color <- brewer.pal(4, "Set3")
# Use ggplot2 to create the barplot.
p <- ggplot(df, aes(x = Name, y = Accuracy, fill = Category)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = color <- brewer.pal(5, "Set1")) +
  coord_flip() +
  labs(title = "Pearson Correlaton Coefficient for Cell Type Proportion Prediction by dtangle",
       x = "Cell Type and parameter settings",
       y = "Pearson Correlation") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

print(p)

#save into png file
png('/external/rprshnas01/netdata_kcni/stlab/Xiaolin/metabrain_PCA/data/barplots_for_three_diff_tech_cov_reg.png')
p
dev.off()

# scp xzhou@dev01.camhres.ca:/external/rprshnas01/netdata_kcni/stlab/Xiaolin/metabrain_PCA/data/barplots_for_three_diff_tech_cov_reg.png /Users/songchenning/Desktop
