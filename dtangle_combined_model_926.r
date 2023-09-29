library(dplyr)
library(tidyverse)
library(edgeR)
library(markerGeneProfile) 
library(matrixStats)
library(cowplot)
library(broom)
library(knitr)
library(ggpubr)
library(biomaRt)
library(ggrepel)
library(patchwork)
library(ggsignif)
library(modelr)
theme_set(theme_classic2())
#Colour palette
cbPalette <- c("#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7","#000000","#E69F00")

#Read in three count matrices after normalization and filtering - preprocessed from /nethome/kcni/xzhou/cell_deconv/dtanlge_wholedata.r
GVEX_matrix <- readRDS('/nethome/kcni/xzhou/cell_deconv/data/GVEX_count_cpm_filtered.RData') #[1] 33638   430
LIBD_matrix <- readRDS('/nethome/kcni/xzhou/cell_deconv/data/LIBD_count_cpm_filtered.RData') # [1] 32782   495
CMC_matrix <- readRDS('/nethome/kcni/xzhou/cell_deconv/data/CMC_count_cpm_filtered.RData') #[1] 33633  1000

#Metadata for individual cohorts + Psychencode -> Merge for ease of access
GVEX_metadata = read.delim("/external/rprshnas01/netdata_kcni/stlab/PsychENCODE/BrainGVEX/RNAseq/SYNAPSE_METADATA_MANIFEST.tsv") %>% subset(dataType == "geneExpression") #%>% select(specimenID, PMI, hemisphere, pH, BrodmannArea, RIN)
LIBD_metadata = read.delim("/external/rprshnas01/netdata_kcni/stlab/PsychENCODE/LIBD__szControl/RNAseq/SYNAPSE_METADATA_MANIFEST.tsv") %>% subset(dataType == "geneExpression") 
CMC_metadata = read.delim("/external/rprshnas01/netdata_kcni/stlab/PsychENCODE/CMC/Metadata/SYNAPSE_METADATA_MANIFEST.tsv") %>% subset(dataType == "geneExpression")
psychencode_metadata = read.csv(("/external/rprshnas01/netdata_kcni/stlab/PsychENCODE/Metadata/CapstoneCollection_Metadata_Clinical.csv"))
GVEX_metadata = GVEX_metadata %>% inner_join(psychencode_metadata)
LIBD_metadata = LIBD_metadata %>% inner_join(psychencode_metadata)
CMC_metadata = CMC_metadata %>% inner_join(psychencode_metadata)

#LIBD has prenatal samples that can affect MGP estimation - remove before analysis
LIBD_matrix = LIBD_matrix[,-which(names(LIBD_matrix) %in% (LIBD_metadata$specimenID %>% subset(LIBD_metadata$contributingStudy == "LIBD_szControl" & LIBD_metadata$ageDeath <=0)))]

#Cell type proportions for each data set merged with metadata - ready for modelling - preprocessed from /nethome/kcni/xzhou/cell_deconv/dtanlge_wholedata.r
gvex_estimations <- read.csv("/nethome/kcni/xzhou/cell_deconv/data/GVEX_estimations_metadata.csv")
libd_estimations  <- read.csv("/nethome/kcni/xzhou/cell_deconv/data/LIBD_estimations_metadata.csv")
cmc_estimations  <- read.csv("/nethome/kcni/xzhou/cell_deconv/data/CMC_estimations_metadata.csv")
gvex_estimations$dataset = rep("GVEX", nrow(gvex_estimations))
libd_estimations$dataset = rep("LIBD", nrow(libd_estimations))
libd_estimations[which(is.na(libd_estimations$individualIdSource)), 'individualIdSource'] = 'LIBD_szControl'
cmc_estimations$dataset = rep("CMC", nrow(cmc_estimations))


cell_type_names <- c('Inh_VIP', 'Inh_LAMP5', 'Exc_IT', 'Inh_PAX6', 'Oligodendrocyte', 'Astrocyte', 'Exc_L5_6_IT_Car3', 'Exc_L5_6_NP', 'Inh_SST', 'Exc_L6_CT', 'OPC', 'Inh_PVALB', 'Exc_L6b', 'Microglia', 'Exc_L5_ET', 'Pericyte', 'Endothelial', 'Exc_L4_IT', 'VLMC')
names(gvex_estimations)[45:(ncol(gvex_estimations)-1)] = cell_type_names
names(cmc_estimations)[47:(ncol(cmc_estimations)-1)] = cell_type_names
names(libd_estimations)[42:(ncol(libd_estimations)-1)] = cell_type_names



#Load single cell samples 
psychencode_snCTPs = read_csv('/external/rprshnas01/netdata_kcni/stlab/cross_cohort_MGPs/psychencode_snCTPs.csv') 

all_estimations = bind_rows(bind_rows(gvex_estimations, libd_estimations), cmc_estimations)
all_estimations$primaryDiagnosis = factor(all_estimations$primaryDiagnosis, levels =c('control', "Bipolar Disorder", "Schizophrenia"))
all_estimations$reportedGender = factor(all_estimations$reportedGender, levels =c('female','male'))
all_estimations = all_estimations %>% filter(!is.na(individualIdSource), !is.na(primaryDiagnosis),individualIdSource != 'BSHRI', ageDeath >= 20)



#all_estimations = all_estimations %>% mutate(individualIdSource = recode(individualID, NA ="LIBD"))
all_estimations_long = all_estimations %>% pivot_longer(cols = Astrocyte:VLMC, names_to = 'cell_type', values_to = 'rel_prop') 
# It appears that SMRI cohorts have very few subjects after filtering. Let's create a new column (new_study) that collapses SHMI cohorts back into GVEX, whilst keeping all other cohorts on their own. 
new_study = all_estimations_long$individualIdSource 
new_study[grepl("SMRI", new_study)] = "GVEX"
all_estimations_long$newStudy = new_study %>% factor(levels = c("Penn", "MSSM", "LIBD_szControl","Pitt", "NIMH_HBCC", "GVEX"))

#Model cell type proportion based on psychiatric diagnosis
combined_lms = all_estimations_long %>% 
  # group stacked data by cell_type
  group_by(newStudy, cell_type) %>%
  
  # fit all the cell_type_prop data accorting to the model 
  # using the broom package to tidy the results 
  do(tidy(lm(rel_prop ~ scale(RIN) + scale(ageDeath)  +  scale(PMI) +
               reportedGender + primaryDiagnosis ,  data = .))) %>%
  
  # unstack the data and adjust for multiple comparisons using the Benjamini-Hochberg method
  ungroup() %>% 
  mutate(padj = p.adjust(`p.value`, method = 'BH')) %>%
  #add cell class labels
  mutate(class = case_when(
    str_detect(cell_type, "Inh") ~ "Inhibitory",
    str_detect(cell_type, "Exc") ~ "Excitatory",
    TRUE ~ "Non-Neuronal")) %>%
  
  # clean up the names the the term column
  mutate(term = recode(term,
                       `(Intercept)` = "Intercept",
                       `reportedGendermale` = "gender:Male",
                       `primaryDiagnosisSchizophrenia` = "SCZ",
                       `primaryDiagnosisBipolar Disorder` = "BP",
                       `DxMDD` = "MDD",
                       `scale(ageDeath)` = "Age",
                       #`scale(PMI..in.hours.)` = "PMI", 
                       `scale(RIN)` = "RIN"))
#Add study name for colour-coding
combined_lms = merge(combined_lms, unique(all_estimations_long[c('dataset', "newStudy")]), by.x = "newStudy")
save(combined_lms, file = "/nethome/kcni/xzhou/cell_deconv/data/combined_lms_926.RData")



# data visualization
scz_plot_combined = combined_lms %>% 
  filter(term %in% 'SCZ') %>% 
  mutate(cell_type = fct_reorder(cell_type, estimate)) %>% 
  ggplot(aes(x = newStudy, y = estimate, fill = dataset)) + 
  geom_hline(yintercept = 0) + 
  geom_bar(stat = "identity", position = "dodge") + 
  geom_errorbar(aes(ymin = estimate - std.error, ymax = estimate + std.error)) + 
  ggtitle("Schizophrenia vs. Controls") +
  ylab('Beta Coefficient') + 
  xlab('Study') + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  facet_wrap(~cell_type, drop = T, scale = "free")
ggsave(scz_plot_combined, file = "/nethome/kcni/xzhou/cell_deconv/pic/scz_plot_combined_926.png", width = 10, height = 10, units = "in", dpi = 300)