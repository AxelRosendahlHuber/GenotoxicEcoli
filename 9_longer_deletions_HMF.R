# Code accompanying the publication # "Mutational signature in colorectal cancer induced by genotoxic pks+ E.Coli"
# C. Pleguezuelos-Manzano, J. Puschhof and A. Rosendahl Huber et al.

# Please direct any questions or comments to A.K.M.RosendahlHuber@prinsesmaximacentrum.nl

# NOTE: This script uses patient-level somatic variant and clinical data obtained from the Hartwig Medical Foundation. 

# ---- HMF -----
# NOTE: This script uses patient-level somatic variant and clinical data which have been obtained from the Hartwig Medical Foundation under the data request number DR-084. Somatic variant and clinical data are freely available for academic use from the Hartwig Medical Foundation through standardized procedures. Privacy and publication policies, including co-authorship policies, can be retrieved from: https://www.hartwigmedicalfoundation.nl/en/data-policy/. 
# Data request forms can be downloaded from https://www.hartwigmedicalfoundation.nl/en/applying-for-data/.
# To gain access to the data, this data request form should be emailed to info@hartwigmedicalfoundation.nl., upon which it will be evaluated within 6 weeks by the HMF Scientific Council and an independent Data Access Board.
# When access is granted, the requested data become available through a download link provided by HMF.

# Script to analyze the percentage of indels + 5bp flanking ends containing pks motifs.
# NOTE: This script requires input generated from the HMF source data using scripts in the folder 'HMF_analysis'

# Axel Rosendahl Huber
setwd(dir)
source("Utils.R")
source("8_Figure4AB_HMF_refitting.R")

# ---- Histograms STE72 indels
HMF_indel= read.delim("indel_status_mc.txt")
HMF_indels= HMF_indel
HMF_indels$pks_status = factor(HMF_indels$pks_status, levels = c("pks_at_flank", "no_pks_pattern"))
intervals = c("2bp_deletions", "3-5bp deletions", "6-10bp deletions", "11-15bp deletions", "15bp+ deletions")
HMF_indels$interval = intervals[findInterval(HMF_indels$indel_length, c(0,3, 6, 11, 16))]
HMF_indels$interval = factor(HMF_indels$interval, levels = c("2bp_deletions", "3-5bp deletions", "6-10bp deletions", "11-15bp deletions", "15bp+ deletions"))

# optional = select only the smaller indels 
#HMF_indels = HMF_indels[HMF_indels$interval == "2bp_deletions" | HMF_indels$interval == "3-5bp deletions"  ,]

table_totals = table(HMF_indels[,c(2,3)]) %>% t()
df_totals = data.frame(pks_at_flank = table_totals[,1], no_pks_pattern = table_totals[,2], 
                       total_indels = rowSums(table_totals))
df_totals$percentage_pks = df_totals$pks_at_flank/df_totals$total_indels
df_totals = df_totals[order(df_totals$percentage_pks, decreasing = T),]

correlation$sample = gsub("_.*$", "", rownames(correlation))
correlation$tissue = gsub("^.*_", "", rownames(correlation))
match(rownames(df_totals), correlation$sample)
correlation$multi_indel = df_totals$percentage_pks[match(correlation$sample, rownames(df_totals))]
ggplot(correlation, aes(multi_indel)) + geom_density() + geom_vline(data = correlation_colon, aes(x = multi_indel))

correlation$tissue = as.character(correlation$tissue)
correlations = correlation
correlation$tissue[!grepl("Colon/Rectum|Urinary tract|Head and neck", correlation$tissue)] = "Other"
unique(correlation$tissue)
rsq(correlation$indel, correlation$snv)
rsq(correlations$indel, correlations$snv)


# Correlation dot-plot for all samples
multi_SBS_correlation = ggplot(correlation, aes(x= multi_indel, y = snv, color = tissue)) + geom_point(alpha = 0.7,size =2, ) + 
  #geom_hline(yintercept =  0.05, linetype = "dashed", color = "#828282") + geom_vline(xintercept = 0.05, linetype = "dashed", color = "#828282") 
  theme_BM() +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 0.8), breaks = seq(0, 0.8, 0.2)) +scale_x_continuous(expand = c(0, 0), limits = c(0,0.6), breaks = seq(0, 0.6, 0.2))  +
  scale_color_manual(values = c("#d95f02", "#1b9e77",  "#bdbdbd", "#7570b3")) +
  xlab(expression(paste("Fraction of Indels containing ",  italic("pks"), "-motif"))) + ylab(expression(paste("Relative SBS-",  italic("pks"), " contribution")))
multi_SBS_correlation
ggsave("Plots/HMF_multi_indel_SBSpks_correlation.pdf", multi_SBS_correlation, width = 7, height = 5)

multi_ID_correlation = ggplot(correlation, aes(x= multi_indel, y = indel, color = tissue)) + geom_point(alpha = 0.7,size =2, ) + 
  #geom_hline(yintercept =  0.05, linetype = "dashed", color = "#828282") + geom_vline(xintercept = 0.05, linetype = "dashed", color = "#828282") 
  theme_BM() +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 0.8), breaks = seq(0, 0.8, 0.2)) +scale_x_continuous(expand = c(0, 0), limits = c(0,0.6), breaks = seq(0, 0.6, 0.2))  +
  scale_color_manual(values = c("#d95f02", "#1b9e77",  "#bdbdbd", "#7570b3")) + 
  xlab(expression(paste("Fraction of Indels containing ",  italic("pks"), "-motif"))) + ylab(expression(paste("Relative ID-",  italic("pks"), " contribution")))
multi_ID_correlation
ggsave("Plots/HMF_multi_indel_IDpks_correlation.pdf", multi_ID_correlation, width = 7, height = 5)

correlation_colon_SBS_ID_positive = correlation[correlation$indel >0.05 & correlation$snv > 0.05,]
correlation$double_positive  = "negative"
correlation$double_positive[correlation$sample %in% correlation_colon_SBS_ID_positive$sample] = "positive"

all_plot = ggplot(correlation, aes(x= multi_indel, y = snv, color = double_positive)) + geom_point(alpha = 0.7,size =2, ) + 
  #geom_hline(yintercept =  0.05, linetype = "dashed", color = "#828282") + geom_vline(xintercept = 0.05, linetype = "dashed", color = "#828282") 
  theme_BM() +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 0.7)) +scale_x_continuous(expand = c(0, 0), limits = c(0,0.7))  +
  scale_color_manual(values = c("#d95f02", "#1b9e77"))
all_plot



