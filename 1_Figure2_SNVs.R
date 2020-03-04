# Code accompanying the publication # "Mutational signature in colorectal cancer induced by genotoxic pks+ E.Coli"
# C. Pleguezuelos-Manzano, J. Puschhof and A. Rosendahl Huber et al.

# Description: Set working directory, analyze single base substitutions (Figure 2B-C)
# Please direct any questions or comments to A.K.M.RosendahlHuber@prinsesmaximacentrum.nl

# Enter input directory where .zip file is extracted: 
dir = paste0("", "Genotoxic_Ecoli")

setwd(dir)
source("Utils.R")

vcf_list = list.files("Data/STE0072/I5/SNVs", pattern = ".vcf", full.names = T)
names = gsub("_72C12TM0_Q50_CGQ10_SGQ99_PASS_20X_VAF0.3_NoY_nonBlacklist_final.vcf", "", basename(vcf_list))
vcfs = read_vcfs_as_granges(vcf_files = vcf_list, names, genome = ref_genome)
mut_mat = mut_matrix(vcfs, ref_genome = ref_genome)
plot = plot_96_profile3(mut_mat,ymax = 200, condensed = F) + scale_y_continuous(expand = c(0,0)) + theme(panel.spacing = unit(1.5, "lines")) 
ggsave("Figures/ExtendedData/Extended_data_figure_2/96_profile_all.pdf", plot, width = 7, height = 15)

write.table(mut_mat, "Output/I5_mut_mat.txt", sep = "\t", quote = F)

# Subtract EKO from EWT patterns (PKS - deltaPKS)
total_DYE = rowSums(mut_mat[,1:3])
total_EKO = rowSums(mut_mat[,4:6])
total_EWT = rowSums(mut_mat[,7:9])
Subtracted = total_EWT - total_EKO
Subtracted = as.data.frame(Subtracted)
Subtracted[Subtracted<0] <- 0
Subtracted_total = Subtracted
Subtracted = Subtracted/sum(Subtracted)
write.table(Subtracted, "Output/SBS_pks_subtracted.txt", quote = F, sep = "\t")

# Main figure 2C 
totals = data.frame(total_EcoliWT = total_EWT, total_EKO = total_EKO,  Subtracted = Subtracted_total)
totals_plot = plot_96_profile3(totals, condensed = F, ymax = 400) + scale_y_continuous(expand = c(0, 0)) +theme(panel.spacing = unit(1.5, "lines"))
ggsave("Figures/Main/Figure2/Figure_2C_total_mutational_profiles.pdf", totals_plot, width = 15, height = 7)

subtracted_EKO_DYE = total_EKO-total_DYE
subtracted_EKO_DYE[subtracted_EKO_DYE<0] <- 0 
subtracted_EWT_DYE = total_EWT-total_DYE
subtracted_EWT_DYE[subtracted_EWT_DYE<0] <- 0
totals_subtracted = data.frame(subtracted_EKO_DYE = subtracted_EKO_DYE,
                               subtracted_EWT_DYE = subtracted_EWT_DYE)
Subtracted_toDYE =plot_96_profile3(totals_subtracted, ymax = 400, condensed = F) + scale_y_continuous(expand = c(0, 0)) +theme(panel.spacing = unit(1.5, "lines"))
ggsave("Figures/ExtendedData/Extended_data_figure_2/Total_mutational_profiles_subtracted_DYE.pdf", Subtracted_toDYE, width = 10, height = 5)


totals_df = data.frame(total_DYE = total_DYE,total_EKO = total_EKO, total_EWT= total_EWT)
cosine_matrix = cos_sim_matrix(totals_df, totals_df)
cosine_matrix[upper.tri(cosine_matrix)] <- NA
cosine_heatmap = plot_cosine_heatmap(cosine_matrix, cluster_rows = F, plot_values = T)
ggsave("Figures/ExtendedData/Extended_data_figure_2/Cosine_similarity_SNVs.pdf", cosine_heatmap, width = 5, height = 5)

# Generate barplots C>N T>N (Figure 2B)
mut_mat = mut_mat[,-1:3]
mut_mat_N = data.frame("T-N" = colSums(mut_mat[49:96,]), "C-N" = colSums(mut_mat[1:48,]))
mut_mat_N$name = rownames(mut_mat_N)
mmut_mat_N = melt(mut_mat_N)
mmut_mat_N = mmut_mat_N[!grepl("DYE", mmut_mat_N$name),]
mmut_mat_N$Coli = c(rep("EKO", 3), rep("EWT", 3), rep("EKO", 3), rep("EWT", 3))
mean = data.frame(Coli = c("EKO", "EWT", "EKO", "EWT"), variable = rep(c("T.N", "C.N"), each = 2),
                  means = c(mean(mmut_mat_N$value[1:3]),
                            mean(mmut_mat_N$value[4:6]),
                            mean(mmut_mat_N$value[7:9]),
                            mean(mmut_mat_N$value[10:12])), 
                  sd = c(sd(mmut_mat_N$value[1:3]),
                         sd(mmut_mat_N$value[4:6]),
                         sd(mmut_mat_N$value[7:9]),
                         sd(mmut_mat_N$value[10:12]))
)

mean$min = mean$means - mean$sd
mean$max = mean$means + mean$sd

mean$variable = factor(mean$variable, levels = c("T.N", "C.N"))
mean$max[1:2] = mean$max[1:2] + mean$means[3:4]
mean$min[1:2] = mean$min[1:2] + mean$means[3:4]

# ---- Generate individual datapoints for barplots figure 2B ----- 
individuals_mut_mat_N = melt(mut_mat_N)
individuals_mut_mat_N = individuals_mut_mat_N [!grepl("DYE", individuals_mut_mat_N $name),]
individuals_mut_mat_N $Coli = c(rep("EKO", 3), rep("EWT", 3), rep("EKO", 3), rep("EWT", 3))
individuals_mut_mat_N$value = individuals_mut_mat_N$value + c(rep(mean$means[3],3), rep(mean$means[4], 3), rep(0,6))

condensed_barplot_SNVs = ggplot(mean, aes(x = Coli, y = means, fill = variable )) + 
  geom_bar(stat = "identity") +  
  geom_jitter(data = individuals_mut_mat_N , aes(x = Coli, y = value, fill = variable), color = "black", shape = 21, width = 0.25, size = 2.3) +
  geom_errorbar(aes(min = min, max = max),width = 0.3) +
  theme_BM() +  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  scale_y_continuous(expand = c(0, 0)) + ylab("Single base substitutions") + xlab("") + coord_cartesian(ylim = c(0,1500)) + scale_fill_manual(values = c("#3E3E3E","#D5D5D5"))
condensed_barplot_SNVs

ggsave("Figures/Main/Figure2/Figure_2B_Condensed_barplot.pdf", condensed_barplot_SNVs, width = 2.7, height = 4)