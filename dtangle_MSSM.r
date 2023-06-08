#import packages
library(GEOquery)
library(dtangle)
library(hgu133plus2.db)
library(AnnotationDbi)
library(limma)
library(ggplot2)
library(reshape2)
library(edgeR)
library(matrixStats)
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
### IMPORTANT REMINDER:
# 1. This script is used to generate the results for comparing the predicted cell type proportions from dtangle and the actual data provided in the CMC MSSM dataset.
# 2. Count matrix and Reference Matrix preparing processes are time-consuming and might lead to killed terminal. If you want to replicate the results, please try to directly load data from the prepared files [#load('/path/to/my/directory/sth.Rdata')].

# read in count matrix and bulk metadata
meta <- read.csv('/external/rprshnas01/netdata_kcni/stlab/cross_cohort_MGPs/psychencode_snCTPs.csv')
count <- read.csv('/external/rprshnas01/kcni/dkiss/cell_prop_psychiatry/data/CMC_count_matrix.csv')
count_cpm <- read.csv('/external/rprshnas01/kcni/dkiss/cell_prop_psychiatry/data/cmc_matrix_cpm.csv')
sinai <- meta[meta$Institution =='MSSM',]

count_MSSM <- count[,str_detect(colnames(count), "MSSM")] # [1] 57820   332
count_MSSM <- cbind(count[,1], count_MSSM)


###################################################### Count Matrix Processing #############################################3

names(count_MSSM) = gsub("X", '', names(count_MSSM))
names(count_MSSM)[1] = "gene_symbol"
#Rename Ensembl id to be compatible with gene names
count_MSSM$gene_symbol = gsub("\\..*", '', count_MSSM$gene_symbol)

#Make EnsemblIDs row names
row.names(count_MSSM) = count_MSSM$gene_symbol
count_MSSM = count_MSSM[,-1]

#Convert matrices to counts per million
count_cpm = cpm(count_MSSM, log = TRUE, prior.count = 0.1)

#Remove genes with low standard deviations
count_sds = rowSds(count_cpm, na.rm = T)
count_cpm_filtered = count_cpm[count_sds > 0.1, ] %>% as.data.frame() %>% rownames_to_column( var = "gene_symbol")  #[1] 46542 genes   333

#Converting Ensembl ID to Gene Name
# mart = useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host = "https://dec2021.archive.ensembl.org/") 
mart = useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host = "https://useast.ensembl.org")
ensembl = getBM(attributes = c("ensembl_gene_id", "hgnc_symbol", "entrezgene_id","description","gene_biotype","percentage_gene_gc_content"), mart = mart)
ensembl_to_gene = (data.frame(ensembl$ensembl_gene_id, ensembl$hgnc_symbol))
names(ensembl_to_gene) = c("gene_symbol", "gene_name")
#remove duplicates
ensembl_to_gene = ensembl_to_gene[!duplicated(ensembl_to_gene[,1]),]
count_cpm_filtered = merge(x=count_cpm_filtered, y=ensembl_to_gene, by = "gene_symbol", all.x = T)
# some gene names are duplicated after matching (‘’, ‘DUXAP8’, ‘GOLGA8M’, ‘ITFG2-AS1’, ‘LINC01238’, ‘PINX1’, ‘POLR2J4’, ‘RN7SL274P’, ‘SIGLEC5’, ‘TUBB7P’), in order to set them as unqiue row names, we only keep the first occurence of the duplicated gene names
count_cpm_filtered <- count_cpm_filtered[!duplicated(count_cpm_filtered$gene_name), ] #[1] 32541  1002
count_cpm_filtered <- count_cpm_filtered[count_cpm_filtered$gene_name != '', ] #remove the empty gene name
count_cpm_filtered <- count_cpm_filtered[!is.na(count_cpm_filtered$gene_name), ]
rownames(count_cpm_filtered) = count_cpm_filtered$gene_name
count_cpm_filtered <- count_cpm_filtered[ , !names(count_cpm_filtered) %in% c("gene_symbol", "gene_name")] # [1] 32582 genes  332 samples
#save(count_cpm_filtered, file='/nethome/kcni/xzhou/cell_deconv/data/count_cpm_filtered.RData')
#load('/nethome/kcni/xzhou/cell_deconv/data/count_cpm_filtered.RData')
#########################################################################################################################################



##################################################### REFERENCE PREPARATION  ###########################################################
# REFERENCE: load the subset of Darmanis et al single cell brain data as our REFERENCE data set.  https://wm1693.box.com/s/c53107rblevygbfnos8ic89zn5pexbc6.]
ref_exp <- read.csv("/external/rprshnas01/netdata_kcni/stlab/Public/AIBS_scRNAseq_2019/human/matrix.csv") #sce:sc_RNAseq gene expression data
ref_anno <- read.csv("/external/rprshnas01/netdata_kcni/stlab/Public/AIBS_scRNAseq_2019/human/metadata.csv")  #anno: cell ID with subcluster label
#save(ref_exp, file='/nethome/kcni/xzhou/cell_deconv/data/ref_exp.RData')
#load('/nethome/kcni/xzhou/cell_deconv/data/ref_exp.RData')
#save(ref_anno, file='/nethome/kcni/xzhou/cell_deconv/data/ref_anno.RData')
#load('/nethome/kcni/xzhou/cell_deconv/data/ref_anno.RData')

#ref data preprocessing and normalization to gene expression cpm
sample_names <- ref_exp[, 1]
ref_update <- ref_exp[, -1] #remove the sample name column
gene_names <- colnames(ref_update)
rownames(ref_update) <- sample_names
ref_update_t <- t(ref_update)
rownames(ref_update_t) <- gene_names
ref_cpm = cpm(ref_update_t, log = TRUE, prior.count = 0.1) 

#Remove genes with low standard deviations
ref_sds = rowSds(ref_cpm, na.rm = T) 
ref_matrix = ref_cpm[ref_sds > 0.1, ]

#Does that solve some of duplicate issues? - Yes, background noise was reduced.
length(unique(rownames(ref_matrix))) # 45478 - No duplication HGNC symbols
dim(ref_matrix) #45478 49417: 45478 genes and 49417 cells
#save(ref_matrix, file='/nethome/kcni/xzhou/cell_deconv/data/ref_matrix.RData')
#load('/nethome/kcni/xzhou/cell_deconv/data/ref_matrix.RData')
############################################## REFERENCE File preparation End ##############################################


########################################################## Dtangle #########################################################
commongenes <- intersect(rownames(count_cpm_filtered), rownames(ref_matrix))
count_final_common <- count_cpm_filtered[pmatch(commongenes, rownames(count_cpm_filtered)), ]
ref_final_common <- ref_matrix[pmatch(commongenes, rownames(ref_matrix)), ]
#join the datasets
y <- cbind(ref_final_common, count_final_common)
y <- normalizeQuantiles(y)
y <- t(y)
#load("/nethome/kcni/xzhou/cell_deconv/data/y_61.RData")

all_cell_type <- unique(ref_anno$subclass_label)
pure_samples <- lapply(1:length(all_cell_type), function(i) {
    which(ref_anno$subclass_label == all_cell_type[i])
})
#all_cell_type[1] <-"non-neuronal (empty)"
names(pure_samples) = all_cell_type

marker_list = find_markers(y,pure_samples=pure_samples,data_type="rna-seq",marker_method='ratio')
#the top 10% of all marker genes for each type
q = .1
quantiles = lapply(marker_list$V,function(x)quantile(x,1-q))
K = length(pure_samples)
n_markers = sapply(1:K,function(i){max(which(marker_list$V[[i]] > quantiles[[i]]))})

#run the deconvolution
marks = marker_list$L
dc <- dtangle(y, pure_samples=pure_samples, n_markers=n_markers, data_type = 'rna-seq', markers = marks)
#save(dc, file='/nethome/kcni/xzhou/cell_deconv/data/dc_61.RData')

#retrieve the estimated cell type proportions from dtangle results
final_est <- dc$estimates[(dim(ref_anno)[1]+1):dim(y)[1],]
colnames(final_est) <- all_cell_type
head(final_est)
#################################################### Dtangle End #########################################################  



######################################## Combine the estimation proportions with meta ############################################
#rename the names of cell types
colnames(final_est) = c('', 'Inh_VIP', 'Inh_LAMP5', 'Exc_IT', 'Inh_PAX6', 'Oligodendrocyte', 'Astrocyte', 'Exc_L5_6_IT_Car3', 'Exc_L5_6_NP', 'Inh_SST', 'Exc_L6_CT', 'OPC', 'Inh_PVALB', 'Exc_L6b', 'Microglia', 'Exc_L5_ET', 'Pericyte', 'Endothelial', 'Exc_L4_IT', 'VLMC')
cmc_estimations = final_est %>% as.data.frame() %>% tibble::rownames_to_column(var = "specimenID")
names(cmc_estimations)[names(cmc_estimations) == "specimenID"] = "individualID"

#combine the estimated cell type proportions with the metadata
CMC_metadata = read.csv("/external/rprshnas01/netdata_kcni/stlab/PsychENCODE/CMC/Metadata/SYNAPSE_TABLE_QUERY_123020650.csv")
names(CMC_metadata)[names(CMC_metadata) == "Individual_ID"] = "individualID"

psychencode_metadata = read.csv(("/external/rprshnas01/netdata_kcni/stlab/PsychENCODE/Metadata/CapstoneCollection_Metadata_Clinical.csv"))
CMC_metadata = CMC_metadata %>% inner_join(psychencode_metadata)

#get the lists of individual id and specimen id for the samples from MSSM
ind_specimen_id <- CMC_metadata[str_detect(CMC_metadata$individualID, "MSSM"),c("individualID", "specimenID")] #[1] 5228 samples with _MSSM_    2
#match the estimated cell type proportions with the individual id in ind_specimen_id
cmc_estimations$specimenID <- ind_specimen_id$specimenID[match(cmc_estimations$individualID, ind_specimen_id$individualID)]
#find the common specimen id between the estimated cell type proportions and the bulk RNA-seq data
common_specimen_id <- intersect(cmc_estimations$specimenID, sinai$bulk_sample_id)
#subset the estimated cell type proportions with the common specimen id
cmc_estimations_common <- cmc_estimations[cmc_estimations$specimenID %in% common_specimen_id,]
# repeat the similar process for the bulk RNA-seq data
sinai_common <- sinai[sinai$bulk_sample_id %in% common_specimen_id,]
sinai_common <- sinai_common[,colnames(sinai_common) %in% c("bulk_sample_id","Astrocyte","Endothelial","LAMP5","Microglia","Oligodendrocyte","OPC","Pericyte","PVALB","SST","VIP")]
rownames(sinai_common) <- sinai_common$bulk_sample_id
sinai_common <- sinai_common[,-11]

#reorder and subset the cell types in the estimated cell type proportions to match meta
new_order <- c(22, 8, 19, 4, 16, 7,13,18,14,11,3)
cmc_estimations_common <- cmc_estimations_common[,new_order]
rownames(cmc_estimations_common) <- cmc_estimations_common$specimenID
cmc_estimations_common <- cmc_estimations_common[,-1]
colnames(cmc_estimations_common) <- colnames(sinai_common)

#makre sure two data frames have the same order of samples
sinai_common = sinai_common[order(rownames(sinai_common), decreasing = T),]
cmc_estimations_common = cmc_estimations_common[order(rownames(cmc_estimations_common), decreasing = T),]
############################################# Combination of Data Ends Here ####################################################


############################################ Plotting Predicted V.S. Actual Proportions Starts Here ##########################################
# Create data frames for the predicted and actual cell type proportions
predicted_df <- data.frame(Sample = rownames(cmc_estimations_common), cmc_estimations_common)
actual_df <- data.frame(Sample = rownames(sinai_common), sinai_common)

# Merge the predicted and actual data frames
merged_df <- merge(predicted_df, actual_df, by = "Sample")

# Create an empty list to store the scatter plots
plots_list <- list()

for (cell_type in colnames(cmc_estimations_common)) {
  # Calculate the prediction accuracy (e.g., correlation coefficient)
  x <- paste(cell_type, ".x", sep = "")
  y <- paste(cell_type, ".y", sep = "")
  accuracy <- cor(merged_df[[x]], merged_df[[y]], use = "pairwise.complete.obs")

  # Create the scatter plot with a trend line
  plot <- ggplot(merged_df, aes(x = .data[[x]], y = .data[[y]])) +
    geom_point() +
    geom_smooth(method = "lm", se = FALSE, color = "red") +  # Add a trend line
    labs(x = "Predicted Proportions", y = "Actual Proportions") +
    ggtitle(paste("Prediction Accuracy for", cell_type, ":", round(accuracy, 2)))

  # Store the scatter plot in the list
  plots_list[[cell_type]] <- plot
}

# Arrange and display the scatter plots in a grid
grid.arrange(grobs = plots_list, ncol = 3)
#############################################################################################################################################