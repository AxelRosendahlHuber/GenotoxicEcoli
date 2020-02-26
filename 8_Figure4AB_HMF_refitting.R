# Code accompanying the publication # "Mutational signature in colorectal cancer induced by genotoxic pks+ E.Coli"
# C. Pleguezuelos-Manzano, J. Puschhof and A. Rosendahl Huber et al.

# Please direct any questions or comments to A.K.M.RosendahlHuber@prinsesmaximacentrum.nl

# NOTE: This script uses patient-level somatic variant and clinical data obtained from the Hartwig Medical Foundation. 

# ---- HMF -----
# NOTE: This script uses patient-level somatic variant and clinical data which have been obtained from the Hartwig Medical Foundation under the data request number DR-084. Somatic variant and clinical data are freely available for academic use from the Hartwig Medical Foundation through standardized procedures. Privacy and publication policies, including co-authorship policies, can be retrieved from: https://www.hartwigmedicalfoundation.nl/en/data-policy/. 
# Data request forms can be downloaded from https://www.hartwigmedicalfoundation.nl/en/applying-for-data/.
# To gain access to the data, this data request form should be emailed to info@hartwigmedicalfoundation.nl., upon which it will be evaluated within 6 weeks by the HMF Scientific Council and an independent Data Access Board.
# When access is granted, the requested data become available through a download link provided by HMF.

# Mutation matrix can be generated using script 6: 
setwd(dir)
source("Utils.R")

mut_mat = read.delim("Mutation Matrix generated from HMF data") %>% as.matrix()
metadata = read.delim("Metadata file obtained from HMF")
metadata$mutmat_cols = colnames(mut_mat)

mut_mat_colon = mut_mat[,metadata$sampleId[metadata$primaryTumorLocation == "Colon/Rectum"]]

working_sigs = # Read in sigProfiler signatures from Alexandrov et al., 2020. https://www.synapse.org/#!Synapse:syn12009743. 
  # All signatures indicated as 'Possible sequencing artefacts' on https://cancer.sanger.ac.uk/cosmic/signatures/SBS/ are removed from this list. 
  
SBS_pks_subtracted = read.delim("Output/SBS_pks_subtracted.txt")
all_sigs = cbind(working_sigs, SBS_pks_subtracted) %>% as.matrix()
all_sigs = prop.table(all_sigs,2)
fit_res = fit_to_signatures(mut_mat, all_sigs)

rel_contribution = prop.table(fit_res$contribution,2)
rel_contribution_snv = rel_contribution
rel_contribution_order = rel_contribution[,order(rel_contribution[48,], decreasing = T)]
rel_contribution_order_snv = rel_contribution[,order(rel_contribution[18,], decreasing = T)]

abs_contribution_order = fit_res$contribution[,order(rel_contribution[48,], decreasing = T)]

# Plot figure 3 SNV contribution of samples 
data_SBS_pks = data.frame(name = colnames(rel_contribution_order[,1:20]), SBS_pks = rel_contribution_order[48,1:20])
data_SBS_pks$name = factor(data_SBS_pks$name, levels = data_SBS_pks$name)
metadata_ordered = metadata[order(rel_contribution[48,], decreasing = T),]
data_SBS_pks$tissue = metadata_ordered$primaryTumorLocation[1:20]
contribution_plot = ggplot(data_SBS_pks, aes(x = name, y = SBS_pks, fill = tissue)) + geom_bar(stat = "identity", alpha = 0.85) + 
  theme_BM() +scale_y_continuous(expand = c(0, 0), limits = c(0, 0.8), breaks = c(0,0.2, 0.4,0.6, 0.8)) + scale_fill_manual(values = c("#d95f02", "#1b9e77", "#7570b3"))
ggsave("Figures/Main/Figure4/FigureA.pdf", contribution_plot, width = 5, height = 4.5)

# ----- Indels 
##########Documentation for the indel and dbs profiles############
library(pracma)
#  Read processed indel file from HMF
name = c(NA, NA, paste0(metadata$sampleId, "_", metadata$primaryTumorLocation))

indel_counts = read.delim("Indel Mutation Matrix generated from HMF data. See folder HMF for script to obtain this matrix")
indel_counts$muttype = factor(x = indel_counts$muttype, levels = c("C_deletion","T_deletion","C_insertion","T_insertion","2bp_deletion","3bp_deletion","4bp_deletion","5+bp_deletion","2bp_insertion","3bp_insertion","4bp_insertion","5+bp_insertion","2bp_deletion_with_microhomology","3bp_deletion_with_microhomology","4bp_deletion_with_microhomology","5+bp_deletion_with_microhomology"))
indel_counts_named = rbind(indel_counts, name)
colnames(indel_counts_named)[3:3670] = indel_counts_named[84,3:3670]
indel_counts_named = indel_counts_named[-84,]
indel_counts_named[,3:3668] = lapply(indel_counts_named[,3:3668], as.integer)

id_sigs = # Read in sigProfiler indel signatures from Alexandrov et al., 2020. https://www.synapse.org/#!Synapse:syn12009743. 
ID_pks = read.delim("Output/ID_pks.txt")
id_sigs$ID_pks = ID_pks$x
row.names(id_sigs) = id_sigs[,1]
id_sigs = as.matrix(id_sigs[,-1])

#Turn the indel counts into a matrix
indel_counts_m = indel_counts_named %>% dplyr::select(-muttype, -muttype_sub) %>% data.matrix() 

####Perform the regular signature re-fitting
refit_out = fit_to_signatures_indels(indel_counts_m, id_sigs)

rel_contribution = prop.table(refit_out$contribution,2)
rel_contribution_indel = rel_contribution
rel_contribution_order = rel_contribution[,order(rel_contribution[18,], decreasing = T)]
rel_contribution_order_indel = rel_contribution[,order(rel_contribution[18,], decreasing = T)]

counts = indel_counts_named[,c("muttype", "muttype_sub",names_top_IDpks)]

# Plot figure 4
data_Indel_pks = data.frame(name = colnames(rel_contribution_order[,1:20]), Indel_pks = rel_contribution_order[18,1:20])
data_Indel_pks$name = factor(data_Indel_pks$name, levels = data_Indel_pks$name)
metadata_ordered = metadata[order(rel_contribution[18,], decreasing = T),]
data_Indel_pks$tissue = metadata_ordered$primaryTumorLocation[1:20]
contribution_plot = ggplot(data_Indel_pks, aes(x = name, y = Indel_pks, fill = tissue)) + geom_bar(stat = "identity", alpha = 0.85) + 
  theme_BM() +scale_y_continuous(expand = c(0, 0), limits = c(0, 0.8), breaks = c(0, 0.2, 0.4, 0.6, 0.8)) + scale_fill_manual(values = c('#d95f02', '#1b9e77','#e7298a','#66a61e','#7570b3'))
contribution_plot 
ggsave("Figures/Main/Figure4/FigureB.pdf", indel_context_per_sample, width = 12, height = 15)



# ----- Analyze correlation of signatures 
tissue = gsub(".*_","", colnames(rel_contribution_indel)) %>% as.factor()
correlation = data.frame(indel = rel_contribution_indel[18,],snv = rel_contribution_snv[48,], tissue = tissue)

rsq <- function (x, y) cor(x, y) ^ 2

correlation$tissue = as.character(correlation$tissue)
correlations = correlation
correlation$sample = gsub("_.*", "", rownames(correlation))
correlation$tissue[!grepl("Colon/Rectum|Urinary tract|Head and neck", correlation$tissue)] = "Other"
unique(correlation$tissue)
rsq(correlation$indel, correlation$snv)
rsq(correlations$indel, correlations$snv)


# Correlation dot-plot for all samples
all_plot = ggplot(correlation, aes(x= indel, y = snv, color = tissue)) + geom_point(alpha = 0.7,size =2, ) + 
  geom_hline(yintercept =  0.05, linetype = "dashed", color = "#828282") + geom_vline(xintercept = 0.05, linetype = "dashed", color = "#828282") + theme_BM() +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 0.8)) +scale_x_continuous(expand = c(0, 0), limits = c(0,0.8))  +
  scale_color_manual(values = c("#d95f02", "#1b9e77",  "#bdbdbd", "#7570b3"))
all_plot

ggsave("Figures/Main/Figure4/FigureC.pdf", all_plot, width = 7, height = 5, useDingbats = F)

rsq(correlation_colon$indel,correlation_colon$snv)

correlation_colon_SBS_ID_positive =correlation_colon[correlation_colon$indel > 0.05 & correlation_colon$snv > 0.05, ]
correlation_colon_n = rownames(correlation_colon_SBS_ID_positive)
correlation_colon_names = gsub("T_.*", "", correlation_colon_n)
correlation_colon_names_grep = paste(correlation_colon_names, collapse = "|")

correlation_36 = correlation[correlation$indel > 0.05 & correlation$snv > 0.05,]
correlation_SBSpkspos = correlation[correlation$snv > 0.05,]
correlation_IDpkspos = correlation[correlation$indel > 0.05,]

correlations_SBSpkspos = correlations[correlation$snv > 0.05,]

table(correlation_SBSpkspos$tissue)
table(correlation_IDpkspos$tissue)

rsq(correlation_colon$indel, correlation_colon$snv)
rsq(correlation$indel, correlation$snv)

lm = lm(formula = indel~snv, data = correlation_colon)
summary(lm)
lm = lm(formula = indel~snv, data = correlation_colon_SBS_ID_positive)
summary(lm)