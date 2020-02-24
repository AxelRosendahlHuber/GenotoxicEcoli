# Transcriptional strand bias I5 data (no dye) 
# Axel Rosendahl Huber 
# 20-8-2019

# ---- Load STE0072 I5 vcf data ----
library(MutationalPatterns)
library(BSgenome)
library(BSgenome.Hsapiens.UCSC.hg19)
library(reshape2)
library(ggplot2)
source("~/surfdrive/Shared/Boxtel_General/Scripts/pmc_vanboxtel/Axel_BMseq/Utils.R")
ref_genome = "BSgenome.Hsapiens.UCSC.hg19"
setwd("~/surfdrive/Shared/Projects/Axel/Axel_GenoEcoli/I5/STE0072/")

folder = "~/surfdrive/Shared/Boxtel_General/Data/Mutation_data/SNVs/GenoEcoli/I5/"
vcf_list = list.files(folder, pattern = ".vcf", full.names = F)
names = gsub("_72C12TM0_Q50_CGQ10_SGQ99_PASS_20X_VAF0.3_NoY_nonBlacklist_final.vcf", "", vcf_list)
vcf_list = paste0(folder, vcf_list)
vcfs = read_vcfs_as_granges(vcf_files = vcf_list,names, genome = ref_genome)

# Generate mut_mat_s 
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(pracma)
genes_hg19 <- genes(TxDb.Hsapiens.UCSC.hg19.knownGene)

# Load only the adult HSPCs as VCFs 
mut_mat_s <- mut_matrix_stranded2(vcfs, ref_genome, genes_hg19)
strand_counts <- strand_occurrences(mut_mat_s)

# Plot 192-mutational profiles with strand bias for each sample
for (sample in colnames(mut_mat_s)) {
  ts_profile = plot_192_profile2(as.data.frame(mut_mat_s)[sample], ymax = 0.15) + ggtitle(sample)
  filename = paste0("Plots/192_Mutation_profiles/", sample, ".pdf")
  ggsave(filename, ts_profile,  width = 7, height = 3)   
}

# Plot all samples in total: 
ts_mut_mat_totals = data.frame(EcoliWT = rowSums(mut_mat_s[,4:6]), EKO = rowSums(mut_mat_s[,1:3]),  Subtracted = rowSums(mut_mat_s[,4:6])-rowSums(mut_mat_s[,1:3]))
ts_mut_mat_totals$Subtracted[ts_mut_mat_totals$Subtracted < 0] = 0
SBS_pks_substracted_ts = data.frame("SBS_pks_substracted_ts" = ts_mut_mat_totals$Subtracted)
rownames(SBS_pks_substracted_ts) = rownames(ts_mut_mat_totals)
write.table(SBS_pks_substracted_ts, "SBS_pks_subtracted_ts.txt", quote = F, sep = "\t")



plot_ts_totals = plot_192_profile3(ts_mut_mat_totals, ymax = 0.18)
ggsave("Plots/192_Mutation_profiles/I5_ts_totals_consensus.pdf", plot_ts_totals, width = 12, height = 8)

# Make simple barplot depicting T>N and C>N transcriptional strand biases
ts_mut_mat_totals

ts_TN_CN = data.frame(type = colnames(ts_mut_mat_totals[,2:1]),
  CA_transcribed = colSums(ts_mut_mat_totals[seq(2, 96, by = 2),2:1]),
  CA_untranscribed = colSums(ts_mut_mat_totals[seq(1, 96, by = 2),2:1]),
  TN_transcribed = colSums(ts_mut_mat_totals[seq(98, 192, by = 2),2:1]),
  TN_untranscribed = colSums(ts_mut_mat_totals[seq(97, 192, by = 2),2:1]))

ts_TN_CN_m = melt(ts_TN_CN, )
ts_TN_CN_m$ts = rep(c("transcribed","transcribed", "untranscribed", "untranscribed" ), 2)
ts_TN_CN_m$subs = rep(c("C>N", "T>N"), each = 4)
ts_TN_CN_m$type = factor(ts_TN_CN_m$type, levels = c("EKO", "EcoliWT"))
ts_bias_plot_simple = ggplot(ts_TN_CN_m, aes(x= type, y = value, fill = subs, alpha  = ts)) + 
  geom_bar(stat = "identity", position = "dodge") +
  scale_color_manual() + scale_alpha_manual(values = c(1,0.35)) + theme_BM() + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, 750), breaks = c(0, 250, 500, 750))  

ggsave("Plots/Simple_ts_bias_plot.pdf", ts_bias_plot_simple, width = 5, height = 5) 



# # Plot per sample category: 
# transcriptional_strand_EKO = plot_192_profile2(rowSums(mut_mat_s[,1:3]))
# ggsave("Plots/192_Mutation_profiles/I5_Transcriptional_strand_EKO.pdf", transcriptional_strand_EKO, width = 8, height = 5)
# 
# plot_ts_EWT = plot_192_profile2(rowSums(mut_mat_s[,4:6]))
# ggsave("Plots/192_Mutation_profiles/I5_Transcriptional_strand_EWT.pdf", plot_ts_EWT, width = 8, height = 5)
# 
# plot_ts_subtracted = plot_192_profile2(rowSums(mut_mat_s[,4:6])-rowSums(mut_mat_s[,1:3]), ymax = 0.2)
# ggsave("Plots/192_Mutation_profiles/I5_Transcriptional_strand_subtracted_consensus.pdf", plot_ts_subtracted, width = 8, height = 5)
