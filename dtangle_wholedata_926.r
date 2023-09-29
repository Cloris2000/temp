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


#read count matrices from Dan's results
GVEX_matrix = read.csv("/external/rprshnas01/kcni/dkiss/cell_prop_psychiatry/data/GVEX_count_matrix.csv")
names(GVEX_matrix) = gsub("X", '', names(GVEX_matrix))
names(GVEX_matrix)[1] = "gene_symbol"
LIBD_matrix = read.csv("/external/rprshnas01/kcni/dkiss/cell_prop_psychiatry/data/LIBD_count_matrix.csv")
names(LIBD_matrix) = gsub("X", '', names(LIBD_matrix))
names(LIBD_matrix)[1] = "gene_symbol"
CMC_matrix = read.csv("/external/rprshnas01/kcni/dkiss/cell_prop_psychiatry/data/CMC_count_matrix.csv")
names(CMC_matrix) = gsub("X", '', names(CMC_matrix))
names(CMC_matrix)[1] = "gene_symbol"

#Rename Ensembl id to be compatible with gene names
GVEX_matrix$gene_symbol = gsub("\\..*", '', GVEX_matrix$gene_symbol)
LIBD_matrix$gene_symbol = gsub("\\..*", '', LIBD_matrix$gene_symbol)
CMC_matrix$gene_symbol = gsub("\\..*", '', CMC_matrix$gene_symbol)

#Make EnsemblIDs row names
row.names(GVEX_matrix) = GVEX_matrix$gene_symbol
GVEX_matrix = GVEX_matrix[,-1]
row.names(LIBD_matrix) = LIBD_matrix$gene_symbol
LIBD_matrix = LIBD_matrix[,-1]
row.names(CMC_matrix) = CMC_matrix$gene_symbol
CMC_matrix = CMC_matrix[,-1]

#import metadata
GVEX_metadata = read.delim("/external/rprshnas01/netdata_kcni/stlab/PsychENCODE/BrainGVEX/RNAseq/SYNAPSE_METADATA_MANIFEST.tsv") %>% subset(dataType == "geneExpression") #%>% select(specimenID, PMI, hemisphere, pH, BrodmannArea, RIN)
LIBD_metadata = read.delim("/external/rprshnas01/netdata_kcni/stlab/PsychENCODE/LIBD__szControl/RNAseq/SYNAPSE_METADATA_MANIFEST.tsv") %>% subset(dataType == "geneExpression") 
# CMC_metadata = read.csv("/external/rprshnas01/netdata_kcni/stlab/PsychENCODE/CMC/Metadata/SYNAPSE_TABLE_QUERY_123020650.csv")%>% subset(dataType == "geneExpression") 
### ^^^ Different from original version !!!
CMC_metadata = read.delim("/external/rprshnas01/netdata_kcni/stlab/PsychENCODE/CMC/Metadata/SYNAPSE_METADATA_MANIFEST.tsv") %>% subset(dataType == "geneExpression")
#[names(CMC_metadata) == "Individual_ID"] = "individualID"

psychencode_metadata = read.csv(("/external/rprshnas01/netdata_kcni/stlab/PsychENCODE/Metadata/CapstoneCollection_Metadata_Clinical.csv"))

GVEX_metadata = GVEX_metadata %>% inner_join(psychencode_metadata)
LIBD_metadata = LIBD_metadata %>% inner_join(psychencode_metadata)
CMC_metadata = CMC_metadata %>% inner_join(psychencode_metadata)

#Convert matrices to counts per million
GVEX_cpm = cpm(GVEX_matrix, log = TRUE, prior.count = 0.1)
LIBD_cpm = cpm(LIBD_matrix, log = TRUE, prior.count = 0.1)
CMC_cpm = cpm(CMC_matrix, log = TRUE, prior.count = 0.1)

#Remove genes with low standard deviations
GVEX_sds = rowSds(GVEX_cpm, na.rm = T) 
LIBD_sds = rowSds(LIBD_cpm, na.rm = T)
CMC_sds = rowSds(CMC_cpm, na.rm = T)
GVEX_matrix = GVEX_cpm[GVEX_sds > 0.1, ] %>% as.data.frame() %>% rownames_to_column( var = "gene_symbol") 
LIBD_matrix = LIBD_cpm[LIBD_sds > 0.1, ] %>% as.data.frame() %>% rownames_to_column( var = "gene_symbol") 
CMC_matrix = CMC_cpm[CMC_sds > 0.1, ] %>% as.data.frame() %>% rownames_to_column( var = "gene_symbol")

#Converting Ensembl ID to Gene Name
mart = useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host = "https://useast.ensembl.org") #, host = "https://dec2021.archive.ensembl.org/"
ensembl = getBM(attributes = c("ensembl_gene_id", "hgnc_symbol", "entrezgene_id","description","gene_biotype","percentage_gene_gc_content"), mart = mart) 
ensembl_to_gene = (data.frame(ensembl$ensembl_gene_id, ensembl$hgnc_symbol))
names(ensembl_to_gene) = c("gene_symbol", "gene_name")
#remove duplicates
ensembl_to_gene = ensembl_to_gene[!duplicated(ensembl_to_gene[,1]),]

GVEX_matrix = merge(x=GVEX_matrix, y=ensembl_to_gene, by = "gene_symbol", all.x = T) #[1] 48296 genes  432 samples
LIBD_matrix = merge(x=LIBD_matrix, y=ensembl_to_gene, by = "gene_symbol", all.x = T) #[1] 47106   497
CMC_matrix = merge(x=CMC_matrix, y=ensembl_to_gene, by = "gene_symbol", all.x = T) #[1] 47580  1002

# some gene names are duplicated after matching, in order to set them as unqiue row names, we only keep the first occurence of the duplicated gene names
name_list <- c("GVEX", "LIBD", "CMC")
for (name in name_list){
    #print(cpm_filtered)
    matrix <- get(paste(name, "_matrix", sep = ""))
    #print(head(matrix))
    cpm_filtered <- matrix[!duplicated(matrix$gene_name), ] # remove duplicated gene names
    cpm_filtered <- cpm_filtered[cpm_filtered$gene_name != '', ] #remove the empty gene name
    cpm_filtered <- cpm_filtered[!is.na(cpm_filtered$gene_name), ] #remove the NA gene name
    rownames(cpm_filtered) = cpm_filtered$gene_name
    cpm_filtered <- cpm_filtered[ , !names(cpm_filtered) %in% c("gene_symbol", "gene_name")]
    #print(head(cpm_filtered))
    output_file_name <- paste0('/nethome/kcni/xzhou/cell_deconv/data/', name, '_count_cpm_filtered.RData')
    saveRDS(cpm_filtered, file=output_file_name)
}

GVEX_cpm_filtered <- readRDS('/nethome/kcni/xzhou/cell_deconv/data/GVEX_count_cpm_filtered.RData') #[1] 33638   430
LIBD_cpm_filtered <- readRDS('/nethome/kcni/xzhou/cell_deconv/data/LIBD_count_cpm_filtered.RData') # [1] 32782   495
CMC_cpm_filtered <- readRDS('/nethome/kcni/xzhou/cell_deconv/data/CMC_count_cpm_filtered.RData') #[1] 33148 genes  1000 samples





##################################################### REFERENCE PREPARATION  ###########################################################
### !!! Please loading the Rdata directly; otherise the reading and calculation time will be extremely long.
# REFERENCE: load the subset of Darmanis et al single cell brain data as our REFERENCE data set.  https://wm1693.box.com/s/c53107rblevygbfnos8ic89zn5pexbc6.]
# ref_exp <- read.csv("/external/rprshnas01/netdata_kcni/stlab/Public/AIBS_scRNAseq_2019/human/matrix.csv") #sce:sc_RNAseq gene expression data
# ref_anno <- read.csv("/external/rprshnas01/netdata_kcni/stlab/Public/AIBS_scRNAseq_2019/human/metadata.csv")  #anno: cell ID with subcluster label
#save(ref_exp, file='/nethome/kcni/xzhou/cell_deconv/data/ref_exp.RData')
load('/nethome/kcni/xzhou/cell_deconv/data/ref_exp.RData')
#save(ref_anno, file='/nethome/kcni/xzhou/cell_deconv/data/ref_anno.RData')
load('/nethome/kcni/xzhou/cell_deconv/data/ref_anno.RData')
#final_est <- dc$estimates[(dim(ref_anno)[1]+1):dim(y)[1],]
#colnames(final_est) <- all_cell_type

#remove the outlier sample failed QC with empty cell type annotation
ref_anno_filtered <- ref_anno[ref_anno$subclass_label != '',]
ref_exp_filtered <- ref_exp[ref_exp$sample_name %in% ref_anno_filtered$sample_name,] #[1] 47432 50282

#ref data preprocessing and normalization to gene expression cpm
sample_names <- ref_exp_filtered[, 1]
ref_update <- ref_exp_filtered[, -1] #remove the sample name column
gene_names <- colnames(ref_update)
rownames(ref_update) <- sample_names
ref_update_t <- t(ref_update)
rownames(ref_update_t) <- gene_names
ref_cpm = edgeR::cpm(ref_update_t, log = TRUE, prior.count = 0.1) 

#Remove genes with low standard deviations
ref_sds = rowSds(ref_cpm, na.rm = T) 
ref_matrix = ref_cpm[ref_sds > 0.1, ]

#Does that solve some of duplicate issues? - Yes, background noise was reduced.
length(unique(rownames(ref_matrix))) # 45478 - No duplication HGNC symbols
dim(ref_matrix) #45091 47432: 45091 genes and 47432 cells
#save(ref_matrix, file='/nethome/kcni/xzhou/cell_deconv/data/ref_matrix_all_annotated.RData')
load('/nethome/kcni/xzhou/cell_deconv/data/ref_matrix_all_annotated.RData')
############################################## REFERENCE File preparation End ##############################################


########################################################## Dtangle #########################################################

### run dtangle (marker_method = 'ratio') on GVEX to estimate the cell type proportions
GVEX_commongenes <- intersect(rownames(GVEX_cpm_filtered), rownames(ref_matrix))
GVEX_final_common <- GVEX_cpm_filtered[pmatch(GVEX_commongenes, rownames(GVEX_cpm_filtered)), ]
ref_final_common <- ref_matrix[pmatch(GVEX_commongenes, rownames(ref_matrix)), ]
### Load y directly from the path given below will significantly increase the runtime.
#join the datasets
GVEX_y <- cbind(ref_final_common, GVEX_final_common)
GVEX_y <- normalizeQuantiles(GVEX_y)
GVEX_y <- t(GVEX_y)
load("/nethome/kcni/xzhou/cell_deconv/data/GVEX_y_926.RData")

all_cell_type <- unique(ref_anno_filtered$subclass_label)
pure_samples <- lapply(1:length(all_cell_type), function(i) {
    which(ref_anno$subclass_label == all_cell_type[i])
})
#all_cell_type[1] <-"non-neuronal (empty)"
names(pure_samples) = all_cell_type

marker_list = find_markers(GVEX_y,pure_samples=pure_samples,data_type="rna-seq",marker_method='ratio')
#the top 10% of all marker genes for each type
q = .1
quantiles = lapply(marker_list$V,function(x)quantile(x,1-q))
K = length(pure_samples)
n_markers = sapply(1:K,function(i){max(which(marker_list$V[[i]] > quantiles[[i]]))})

#run the deconvolution
marks = marker_list$L
dc <- dtangle(GVEX_y, pure_samples=pure_samples, n_markers=n_markers, data_type = 'rna-seq', markers = marks)
#save(dc, file='/nethome/kcni/xzhou/cell_deconv/data/dc_GVEX_ratio.RData')
load('/nethome/kcni/xzhou/cell_deconv/data/dc_GVEX_ratio_926.RData')

#retrieve the estimated cell type proportions from dtangle results
GVEX_final_est <- dc$estimates[(dim(ref_anno_filtered)[1]+1):dim(GVEX_y)[1],]
colnames(GVEX_final_est) <- all_cell_type
head(GVEX_final_est)
#saveRDS(GVEX_final_est, file='/nethome/kcni/xzhou/cell_deconv/data/GVEX_dtangle_est.RData')
GVEX_final_est <- readRDS('/nethome/kcni/xzhou/cell_deconv/data/GVEX_dtangle_est.RData')


### Combine the estimation proportions with meta 
# #rename the names of cell types
# colnames(GVEX_final_est) = c('', 'Inh_VIP', 'Inh_LAMP5', 'Exc_IT', 'Inh_PAX6', 'Oligodendrocyte', 'Astrocyte', 'Exc_L5_6_IT_Car3', 'Exc_L5_6_NP', 'Inh_SST', 'Exc_L6_CT', 'OPC', 'Inh_PVALB', 'Exc_L6b', 'Microglia', 'Exc_L5_ET', 'Pericyte', 'Endothelial', 'Exc_L4_IT', 'VLMC')
# GVEX_estimations = GVEX_final_est %>% as.data.frame() %>% tibble::rownames_to_column(var = "specimenID")
# names(GVEX_estimations)[names(GVEX_estimations) == "specimenID"] = "individualID"

# #combine the estimated cell type proportions with the metadata
# #rename the colnames: sample names to group names
# # CMC_metadata = read.delim("/external/rprshnas01/netdata_kcni/stlab/PsychENCODE/CMC/Metadata/SYNAPSE_METADATA_MANIFEST.tsv") %>% subset(dataType == "geneExpression")
# # psychencode_metadata = read.csv(("/external/rprshnas01/netdata_kcni/stlab/PsychENCODE/Metadata/CapstoneCollection_Metadata_Clinical.csv"))
# # CMC_metadata = CMC_metadata %>% inner_join(psychencode_metadata) #[1] 613  44

# #get the lists of individual id and specimen id for the samples from MSSM
# ind_specimen_id <- GVEX_metadata[,c("individualID", "specimenID")]
# #match the estimated cell type proportions with the individual id in ind_specimen_id
# GVEX_estimations$specimenID <- ind_specimen_id$specimenID[match(GVEX_estimations$individualID, ind_specimen_id$individualID)]
#Coerce estimations list into data frame 
gvex_estimations = GVEX_final_est %>% as.data.frame() %>% tibble::rownames_to_column(var = "specimenID")


#Merge cell type proportions with sample metadata
gvex_estimations_metadata = inner_join(GVEX_metadata %>% mutate(specimenID = make.names(specimenID)), 
                                       gvex_estimations %>% mutate(specimenID = make.names(specimenID)))

#Remove '+' from ageDeath for modelling
gvex_estimations_metadata$ageDeath = as.numeric(gsub("[+]", "", gvex_estimations_metadata$ageDeath))

write.csv(gvex_estimations_metadata, "/nethome/kcni/xzhou/cell_deconv/data/GVEX_estimations_metadata.csv")
#############################################################################################################################

### run dtangle on LIBD to estimate the cell type proportions
LIBD_commongenes <- intersect(rownames(LIBD_cpm_filtered), rownames(ref_matrix))
LIBD_final_common <- LIBD_cpm_filtered[pmatch(LIBD_commongenes, rownames(LIBD_cpm_filtered)), ]
ref_final_common <- ref_matrix[pmatch(LIBD_commongenes, rownames(ref_matrix)), ]
#join the datasets
LIBD_y <- cbind(ref_final_common, LIBD_final_common)
LIBD_y <- normalizeQuantiles(LIBD_y)
LIBD_y <- t(LIBD_y)
load("/nethome/kcni/xzhou/cell_deconv/data/LIBD_y_926.RData")

all_cell_type <- unique(ref_anno_filtered$subclass_label)
pure_samples <- lapply(1:length(all_cell_type), function(i) {
    which(ref_anno_filtered$subclass_label == all_cell_type[i])
})
#all_cell_type[1] <-"non-neuronal (empty)"
names(pure_samples) = all_cell_type

marker_list = find_markers(LIBD_y,pure_samples=pure_samples,data_type="rna-seq",marker_method='ratio')
#the top 10% of all marker genes for each type
q = .1
quantiles = lapply(marker_list$V,function(x)quantile(x,1-q))
K = length(pure_samples)
n_markers = sapply(1:K,function(i){max(which(marker_list$V[[i]] > quantiles[[i]]))})

#run the deconvolution
marks = marker_list$L
dc <- dtangle(LIBD_y, pure_samples=pure_samples, n_markers=n_markers, data_type = 'rna-seq', markers = marks)
#load the dc for LIBD
load('/external/rprshnas01/netdata_kcni/stlab/Xiaolin/wholedata_dtangle/dc_LIBD_926_ratio.RData')

#retrieve the estimated cell type proportions from dtangle results
LIBD_final_est <- dc$estimates[(dim(ref_anno_filtered)[1]+1):dim(LIBD_y)[1],]
colnames(LIBD_final_est) <- all_cell_type
head(LIBD_final_est)


libd_estimations = LIBD_final_est %>% as.data.frame() %>% tibble::rownames_to_column(var = "specimenID")


#Merge cell type proportions with sample metadata
libd_estimations_metadata = inner_join(LIBD_metadata %>% mutate(specimenID = make.names(specimenID)), 
                                       libd_estimations %>% mutate(specimenID = make.names(specimenID)))

#Remove '+' from ageDeath for modelling
libd_estimations_metadata$ageDeath = as.numeric(gsub("[+]", "", libd_estimations_metadata$ageDeath))

write.csv(libd_estimations_metadata, "/nethome/kcni/xzhou/cell_deconv/data/LIBD_estimations_metadata.csv")
#############################################################################################################################

### run dtangle on CMC to estimate the cell type proportions
CMC_commongenes <- intersect(rownames(CMC_cpm_filtered), rownames(ref_matrix))
CMC_final_common <- CMC_cpm_filtered[pmatch(CMC_commongenes, rownames(CMC_cpm_filtered)), ]
ref_final_common <- ref_matrix[pmatch(CMC_commongenes, rownames(ref_matrix)), ]
#join the datasets
CMC_y <- cbind(ref_final_common, CMC_final_common)
CMC_y <- normalizeQuantiles(CMC_y)
CMC_y <- t(CMC_y)
load("/nethome/kcni/xzhou/cell_deconv/data/CMC_y_926.RData")

all_cell_type <- unique(ref_anno_filtered$subclass_label)
pure_samples <- lapply(1:length(all_cell_type), function(i) {
    which(ref_anno_filtered$subclass_label == all_cell_type[i])
})
#all_cell_type[1] <-"non-neuronal (empty)"
names(pure_samples) = all_cell_type

marker_list = find_markers(CMC_y,pure_samples=pure_samples,data_type="rna-seq",marker_method='ratio')
#the top 10% of all marker genes for each type
q = .1
quantiles = lapply(marker_list$V,function(x)quantile(x,1-q))
K = length(pure_samples)
n_markers = sapply(1:K,function(i){max(which(marker_list$V[[i]] > quantiles[[i]]))})

#run the deconvolution
marks = marker_list$L
dc <- dtangle(CMC_y, pure_samples=pure_samples, n_markers=n_markers, data_type = 'rna-seq', markers = marks)
#sbatch -J dc_CMC --mem=4G -N 2 -t 0-1:0 -o /nethome/kcni/xzhou/cell_deconv/dc_CMC_ratio.slurm.log --wrap="Rscript /nethome/kcni/xzhou/cell_deconv/dc_CMC.r" 
#save(dc, file='/nethome/kcni/xzhou/cell_deconv/data/dc_GVEX_ratio.RData')
load('/external/rprshnas01/netdata_kcni/stlab/Xiaolin/wholedata_dtangle/dc_CMC_ratio_926.RData')

#retrieve the estimated cell type proportions from dtangle results
CMC_final_est <- dc$estimates[(dim(ref_anno_filtered)[1]+1):dim(CMC_y)[1],]
colnames(CMC_final_est) <- all_cell_type
head(CMC_final_est)

cmc_estimations = CMC_final_est %>% as.data.frame() %>% tibble::rownames_to_column(var = "specimenID")


#Merge cell type proportions with sample metadata by individualID (changed from previous version due to meta change!!!)
cmc_estimations_metadata = inner_join(CMC_metadata, cmc_estimations %>% mutate(individualID = make.names(specimenID)), by="individualID")

#Remove '+' from ageDeath for modelling
cmc_estimations_metadata$ageDeath = as.numeric(gsub("[+]", "", cmc_estimations_metadata$ageDeath))

write.csv(cmc_estimations_metadata, "/nethome/kcni/xzhou/cell_deconv/data/CMC_estimations_metadata.csv")
#############################################################################################################################################
#load the estmated+meta data
gvex_estimations_metadata <- read.csv("/nethome/kcni/xzhou/cell_deconv/data/GVEX_estimations_metadata.csv")
libd_estimations_metadata <- read.csv("/nethome/kcni/xzhou/cell_deconv/data/LIBD_estimations_metadata.csv")
cmc_estimations_metadata <- read.csv("/nethome/kcni/xzhou/cell_deconv/data/CMC_estimations_metadata.csv")

##############################################################################################################################################
### Modelling cell type proportion based on primary diagnosis
#GVEX - Model cell type proportion based on primaryDiagnosis + covariates
cell_types = c('Inh_VIP', 'Inh_LAMP5', 'Exc_IT', 'Inh_PAX6', 'Oligodendrocyte', 'Astrocyte', 'Exc_L5_6_IT_Car3', 'Exc_L5_6_NP', 'Inh_SST', 'Exc_L6_CT', 'OPC', 'Inh_PVALB', 'Exc_L6b', 'Microglia', 'Exc_L5_ET', 'Pericyte', 'Endothelial', 'Exc_L4_IT', 'VLMC')
colnames(gvex_estimations_metadata)[44:ncol(gvex_estimations_metadata)] <- cell_types
gvex_estimations_metadata$primaryDiagnosis = factor(gvex_estimations_metadata$primaryDiagnosis, levels =c('control', "Bipolar Disorder", "Schizophrenia"))
gvex_lms = lapply(cell_types, function(cell_type){
  lm = paste0(cell_type, " ~ scale(RIN) + scale(ageDeath) + primaryDiagnosis + reportedGender")
  results = lm(lm, data = gvex_estimations_metadata) %>% tidy() %>% as.data.frame()
  results$term = c("Intercept", "RIN", "ageDeath", "bipolar_disorder", "schizophrenia", "genderMale")
  results$cell_type = cell_type
  return(results)
}) %>% bind_rows() %>%
  # adjust for multiple comparisons using the Benjamini-Hochberg method
  mutate(padj = p.adjust(`p.value`, method = 'BH')) %>%
  mutate(class = case_when(
    str_detect(cell_type, "Inh") ~ "Inhibitory",
    str_detect(cell_type, "Exc") ~ "Excitatory",
    TRUE ~ "Non-Neuronal"
  ))
#Save model dataframe
write_csv(gvex_lms, "/nethome/kcni/xzhou/cell_deconv/data/gvex_models_926.csv", col_names = T)
gvex_lms <- read_csv("/nethome/kcni/xzhou/cell_deconv/data/gvex_models_926.csv")

#LIBD - Model cell type proportion based on primaryDiagnosis + covariates
names(libd_estimations_metadata)[41:ncol(libd_estimations_metadata)] = cell_types
libd_estimations_metadata$primaryDiagnosis = factor(libd_estimations_metadata$primaryDiagnosis, levels =c('control', "Schizophrenia"))
libd_lms = lapply(cell_types, function(cell_type){
  lm = paste0("scale(", cell_type, ")", " ~ scale(PMI) + scale(RIN) + scale(ageDeath) + scale(pH) + primaryDiagnosis + reportedGender")
  results = lm(lm, data = libd_estimations_metadata) %>% tidy() %>% as.data.frame()
  results$term = c("Intercept", "PMI", "RIN", "ageDeath", "pH", "schizophrenia", "genderMale")
  results$cell_type = cell_type
  return(results)
}) %>% bind_rows() %>%
  # adjust for multiple comparisons using the Benjamini-Hochberg method
  mutate(padj = p.adjust(`p.value`, method = 'BH')) %>%
  mutate(class = case_when(
    str_detect(cell_type, "Inh") ~ "Inhibitory",
    str_detect(cell_type, "Exc") ~ "Excitatory",
    TRUE ~ "Non-Neuronal"
  ))
#Save model dataframe
write_csv(libd_lms, "/nethome/kcni/xzhou/cell_deconv/data/libd_models_926.csv", col_names = T)
libd_lms <- read_csv("/nethome/kcni/xzhou/cell_deconv/data/libd_models_926.csv")

#CMC - Model cell type proportion based on primaryDiagnosis + covariates
names(cmc_estimations_metadata)[46:ncol(cmc_estimations_metadata)] = cell_types
cmc_estimations_metadata$primaryDiagnosis = factor(cmc_estimations_metadata$primaryDiagnosis, levels =c('control', "Schizophrenia"))
cmc_lms = lapply(cell_types, function(cell_type){
  lm = paste0(cell_type, " ~scale(RIN) + scale(ageDeath) + primaryDiagnosis + reportedGender")
  results = lm(lm, data = cmc_estimations_metadata) %>% tidy() %>% as.data.frame()
  results$term = c("Intercept", "RIN", "ageDeath", "schizophrenia", "genderMale")
  results$cell_type = cell_type
  return(results)
}) %>% bind_rows() %>%
  # adjust for multiple comparisons using the Benjamini-Hochberg method
  mutate(padj = p.adjust(`p.value`, method = 'BH')) %>%
  mutate(class = case_when(
    str_detect(cell_type, "Inh") ~ "Inhibitory",
    str_detect(cell_type, "Exc") ~ "Excitatory",
    TRUE ~ "Non-Neuronal"
  ))
#Save model dataframe
write_csv(cmc_lms, "/nethome/kcni/xzhou/cell_deconv/data/cmc_models_926.csv", col_names = T)
cmc_lms <- read_csv("/nethome/kcni/xzhou/cell_deconv/data/cmc_models_926.csv")

#Add study column for colour coding - combine data
gvex_lms$study = rep("GVEX", nrow(gvex_lms))
libd_lms$study = rep("LIBD", nrow(libd_lms))
cmc_lms$study = rep("CMC", nrow(cmc_lms))
combined_lms = rbind(gvex_lms, libd_lms, cmc_lms) 

#Plotting beta coefficients per cell type for each disorder in each cohort
beta_plot_gvex_scz = gvex_lms %>% 
  filter(term %in% 'schizophrenia') %>% 
  mutate(cell_type = fct_reorder(cell_type, estimate)) %>% 
  ggplot(aes(x = cell_type, y = estimate)) + 
  geom_hline(yintercept = 0) + 
  geom_bar(stat = "identity") + 
  geom_errorbar(aes(ymin = estimate - std.error, ymax = estimate + std.error)) + 
  ggtitle("GVEX: Schizophrenia vs. Controls") +
  ylab('Beta Coefficient') + 
  xlab('Cell Type Proportions') + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  facet_wrap(~class, drop = T, scale = "free")
beta_plot_gvex_scz

#libd
beta_plot_libd_scz = libd_lms %>% 
  filter(term %in% 'schizophrenia') %>% 
  mutate(cell_type = fct_reorder(cell_type, estimate)) %>% 
  ggplot(aes(x = cell_type, y = estimate)) + 
  geom_hline(yintercept = 0) + 
  geom_bar(stat = "identity") + 
  geom_errorbar(aes(ymin = estimate - std.error, ymax = estimate + std.error)) + 
  ggtitle("LIBD: Schizophrenia vs. Controls") +
  ylab('Beta Coefficient') + 
  xlab('Cell Type Proportions') + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  facet_wrap(~class, drop = T, scale = "free")
beta_plot_libd_scz

#cmc
beta_plot_cmc_scz = cmc_lms %>% 
  filter(term %in% 'schizophrenia') %>% 
  mutate(cell_type = fct_reorder(cell_type, estimate)) %>% 
  ggplot(aes(x = cell_type, y = estimate)) + 
  geom_hline(yintercept = 0) + 
  geom_bar(stat = "identity") + 
  geom_errorbar(aes(ymin = estimate - std.error, ymax = estimate + std.error)) + 
  ggtitle("CMC: Schizophrenia vs. Controls") +
  ylab('Beta Coefficient') + 
  xlab('Cell Type Proportions') + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  facet_wrap(~class, drop = T, scale = "free")
beta_plot_cmc_scz