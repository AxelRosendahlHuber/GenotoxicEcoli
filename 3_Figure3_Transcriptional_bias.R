# Code accompanying the publication # "Mutational signature in colorectal cancer induced by genotoxic pks+ E.Coli"
# C. Pleguezuelos-Manzano, J. Puschhof and A. Rosendahl Huber et al.

# Description: Analysis of the transcriptional strand bias of single base substitutions (Figure 3D-E)
# Please direct any questions or comments to A.K.M.RosendahlHuber@prinsesmaximacentrum.nl
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(pracma)
setwd(dir)
genes_hg19 <- genes(TxDb.Hsapiens.UCSC.hg19.knownGene)

SNV_vcf_folder = "Data/STE0072/I5/SNVs/"
vcf_list = list.files(SNV_vcf_folder, pattern = ".vcf", full.names = F)
names = gsub("_72C12TM0_Q50_CGQ10_SGQ99_PASS_20X_VAF0.3_NoY_nonBlacklist_final.vcf", "", vcf_list)
vcf_list = paste0(SNV_vcf_folder, vcf_list)
vcfs = read_vcfs_as_granges(vcf_files = vcf_list,names, genome = ref_genome)
vcfs = vcfs[4:9]

# Generate mut_mat_s 
mut_mat_s <- mut_matrix_stranded2(vcfs, ref_genome, genes_hg19)
strand_counts <- strand_occurrences(mut_mat_s)

# Plot all samples in total: 
ts_mut_mat_totals = data.frame(EcoliWT = rowSums(mut_mat_s[,4:6]), EKO = rowSums(mut_mat_s[,1:3]),  Subtracted = rowSums(mut_mat_s[,4:6])-rowSums(mut_mat_s[,1:3]))
ts_mut_mat_totals$Subtracted[ts_mut_mat_totals$Subtracted < 0] = 0
SBS_pks_substracted_ts = data.frame("SBS_pks_substracted_ts" = ts_mut_mat_totals$Subtracted)
rownames(SBS_pks_substracted_ts) = rownames(ts_mut_mat_totals)
write.table(SBS_pks_substracted_ts, "Output/SBS_pks_subtracted_ts.txt", quote = F, sep = "\t")
plot_ts_totals = plot_192_profile3(ts_mut_mat_totals, ymax = 0.18)
ggsave("Figures/Main/Figure3/Figure_3E_192_Mutation_profiles.pdf", plot_ts_totals, width = 12, height = 8)

# -----Make barplot using means and indicate standard deviation using error bars ------
subtypes = data.frame(type = colnames(mut_mat_s),
                      CA_transcribed = colSums(mut_mat_s[seq(2, 96, by = 2),]),
                      CA_untranscribed = colSums(mut_mat_s[seq(1, 96, by = 2),]),
                      TN_transcribed = colSums(mut_mat_s[seq(98, 192, by = 2),]),
                      TN_untranscribed = colSums(mut_mat_s[seq(97, 192, by = 2),]))

subtypes_m = melt(subtypes)
subtypes_m$name = ifelse(grepl("EKO", subtypes_m$type), "EKO", "EcoliWT")
subtypes_m$subs = ifelse(grepl("CA", subtypes_m$variable), "CA", "TN")
subtypes_m$ts = ifelse(grepl("untranscribed", subtypes_m$variable), "untranscribed", "transcribed")

sds_EKO = sapply(subtypes[1:3,-1], sd)
sds_EcoliWT = sapply(subtypes[4:6, -1], sd)
sds_total = rbind(sds_EKO, sds_EcoliWT)

means_EKO = sapply(subtypes[1:3,-1], mean)
means_EcoliWT = sapply(subtypes[4:6, -1], mean)
means_total = rbind(means_EKO, means_EcoliWT)

sds_total_m = melt(sds_total)
means_total_m = melt(means_total)
means_total_m$sd  = sds_total_m$value
means_total_m$ts = gsub("^.*_", "", means_total_m$Var2) 
means_total_m$subs = gsub("_.*$", "", means_total_m$Var2)
means_total_m$type = factor(gsub("means_", "", means_total_m$Var1), levels = c("EKO", "EcoliWT"))

ts_bias_plot_simple = ggplot(means_total_m, aes(x= type, y = value, fill = subs, alpha  = ts)) + 
  geom_bar(stat = "identity", position = "dodge") +
  geom_jitter(data = subtypes_m, aes(x = name, y = value, fill = subs, alpha = ts), color = "black", shape = 21,size = 2.3, position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.9)) +
  scale_color_manual() + scale_alpha_manual(values = c(1,0.35)) + theme_BM() + 
  geom_errorbar(aes(ymin = value - sd, ymax = value + sd),width = 0.4, position = position_dodge(0.9), stat = "identity") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 400), breaks = c(0, 100, 200, 300, 400))  + ylab("Number of single base substitutions")
ts_bias_plot_simple

ggsave("Figures/Main/Figure3/Figure_3D_ts_bias_plot_errorbars.pdf", ts_bias_plot_simple, width = 5, height = 5)