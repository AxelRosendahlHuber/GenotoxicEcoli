# Code accompanying the publication # "Mutational signature in colorectal cancer induced by genotoxic pks+ E.Coli"
# C. Pleguezuelos-Manzano, J. Puschhof and A. Rosendahl Huber et al.

# Description: Set working directory, analyze indels (Figure 2D-E)
# Please direct any questions or comments to A.K.M.RosendahlHuber@prinsesmaximacentrum.nl

##########Documentation for the indel and dbs profiles############
# Define function to plot samples individually
plot_indels <- function(vcfs_grl, path, vcf_fnames) {
  indel_list = list()
  for (i in 1:length(vcfs_grl)) {
    sample_name = strsplit(vcf_fnames[i], split = "_")[[1]][1]
    print(sample_name)
    grl_exp = get_indel_context(vcfs_grl[[i]], ref_genome)
    indel_counts = count_indel_contexts(grl = grl_exp)
    
    main_indel_context_fig = plot_main_indel_contexts(indel_counts) + ggtitle(sample_name)
    ggsave(paste0(path, sample_name, "_main_indels.pdf"), main_indel_context_fig, height = 5, width = 8)
    
    indel_context_fig = plot_indel_contexts(indel_counts) + ggtitle(sample_name)
    ggsave(paste0(path, sample_name, "_indels.pdf"), height = 5, width = 13)
    ggsave(paste0(path, sample_name, "_indels.png"), height = 5, width = 13)
    
    indel_list[[sample_name]] = indel_counts
  }
  return(indel_list)
}

####_____________________Indel___________________####
vcf_fnames = list.files("Data/STE0072/I5/Indels/", recursive = T, full.names = T)
sample_names = gsub("_.*$", "", basename(vcf_fnames))

#Read in indels. (You can't use read_vcfs_as_granges as this function removes indels)
chroms = paste0("chr", c(1:22, "X"))
vcfs_list = lapply(vcf_fnames, function(x){
  vcf = readVcf(x, genome = "hg19")
  gr = granges(vcf)
  gr = remove_snvs(gr)
  seqlevelsStyle(gr) = "UCSC"
  seqlevels(gr, pruning.mode = "coarse") = chroms
  return(gr)
})
vcfs_grl = GRangesList(vcfs_list)
names(vcfs_grl) = sample_names
grl_exp = get_indel_context(vcfs_grl, ref_genome)
indel_counts = count_indel_contexts(grl = grl_exp)
total_indel_context = plot_indel_contexts(indel_counts, same_y = T) + theme(legend.position = "none") + expand_limits( y = 0)
indel_context_fig = plot_main_indel_contexts(indel_counts,same_y = T)

ggsave("Figures/ExtendedData/Extended_data_figure_2/Total_indels_detail_STE0072I5.pdf", total_indel_context, width = 10, height = 15)

# get consensus   
indel_counts$DYE = rowSums(indel_counts[,3:5])
indel_counts$EKO = rowSums(indel_counts[,6:8]) 
indel_counts$EWT = rowSums(indel_counts[,9:11])
indel_counts$subtracted = indel_counts$EWT - indel_counts$EKO
indel_counts$subtracted[indel_counts$subtracted < 0] = 0

indel_subtracted = indel_counts[,9:11]
indel_subs = cbind(indel_counts[,1:2], indel_subtracted)
names(indel_subs)[3:5] = c("deltaPKS", "EcoliWT", "Subtracted")
indel_consensus_reordered = plot_indel_contexts(indel_subs, same_y = T)
ggsave("Figures/Main/Figure2/Figure_2E_Total_consensus_indels_reordered.pdf", indel_consensus_reordered, width = 15, height = 5)

# Get indel subtraction levels: 
indel_counts$EKOminDYE = indel_counts$EKO- indel_counts$DYE
indel_counts$EKOminDYE[indel_counts$EKOminDYE < 0] = 0

indel_counts$EWTminDYE = indel_counts$EWT- indel_counts$DYE
indel_counts$EWTminDYE[indel_counts$EWTminDYE < 0] = 0

indel_subtraction_plot = plot_indel_contexts(indel_counts[,c(1:2, 16, 17)],same_y = T) + ylim(c(0,120))
ggsave("Figures/ExtendedData/Extended_data_figure_2/ED_Fig2E_Total_consensus_indels_subtraction.pdf", indel_subtraction_plot, width = 15, height = 5)

totals_df = data.frame(DYE = indel_counts$DYE, EKO = indel_counts$EKO, EWT = indel_counts$EWT)
cossim_matrix = cos_sim_matrix(totals_df, totals_df)
cossim_matrix[upper.tri(cossim_matrix)] = NA
cosine_plot = plot_cosine_heatmap(cossim_matrix, cluster_rows = F, plot_values = T)
ggsave("Figures/ExtendedData/Extended_data_figure_2/ED_Fig2F.pdf", cosine_plot, width = 5, height = 5)

# Define indel signature 
vcfs_total_EWT = unlist(vcfs_grl[7:9])  %>% GRangesList()
grl_total_EWT = get_indel_context(vcfs_total_EWT, ref_genome)
indel_list_total_EWT = count_indel_contexts(grl = grl_total_EWT)

vcfs_total_EKO = unlist(vcfs_grl[4:6])  %>% GRangesList()
grl_total_EKO = get_indel_context(vcfs_total_EKO, ref_genome)
indel_list_total_EKO = count_indel_contexts(grl = grl_total_EKO)

subtracted_indels = indel_list_total_EWT$counts - indel_list_total_EKO$counts
indel_list_subtracted = indel_list_total_EWT
subtracted_indels[subtracted_indels < 0 ] = 0
indel_list_subtracted$total = subtracted_indels

# Save ID_pks signature 
subtracted_indels_rel = subtracted_indels/sum(subtracted_indels)
write.table(subtracted_indels_rel, "Output/ID_pks.txt", sep = "\t", quote = F)

# Indel subtype barplot
del_counts = indel_counts[,-1:-2]
simple_indel_counts = data.frame(singleTdels_Thomopolymer = colSums(del_counts[7:12,]), other_indels = colSums(del_counts[-7:-12,]) )
simple_indel_counts$name = rownames(simple_indel_counts)
indel_counts_m = melt(simple_indel_counts)
indel_means = indel_counts_m[c(11:12, 26:27),]
indel_samples = indel_counts_m[c(4:9, 19:24),]

indel_means$indel_means = indel_means$value/3
indel_means$sd = c(sd(indel_samples$value[1:3]),
                   sd(indel_samples$value[4:6]),
                   sd(indel_samples$value[7:9]),
                   sd(indel_samples$value[10:12]))

indel_means$min = indel_means$indel_means - indel_means$sd
indel_means$max = indel_means$indel_means + indel_means$sd

indel_means$variable = factor(indel_means$variable, levels = c("singleTdels_Thomopolymer", "other_indels"))
indel_means$max[1:2] = indel_means$max[1:2] + indel_means$indel_means[3:4]
indel_means$min[1:2] = indel_means$min[1:2] + indel_means$indel_means[3:4]

indel_counts_m_added = indel_counts_m[c(4:9, 19:24),]
indel_counts_m_added$names = indel_counts_m_added[,1]
indel_counts_m_added$name = ifelse(grepl("EKO", indel_counts_m_added$names), "EKO", "EWT")
indel_counts_m_added$value = indel_counts_m_added$value + c(rep(indel_means$indel_means[3],3), rep(indel_means$indel_means[4], 3), rep(0,6)) 

condensed_barplot_Indels = ggplot(indel_means, aes(x = name, y = indel_means, fill = variable )) + 
  geom_bar(stat = "identity") +
  geom_jitter(data = indel_counts_m_added, aes(x = name, y = value, fill = variable), color = "black", shape = 21, width = 0.25, size = 2.3) +
  geom_errorbar(aes(min = min, max = max),width = 0.3) +
  theme_BM() +  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_y_continuous(expand = c(0, 0)) + ylab("Number of indels") + xlab("") + coord_cartesian(ylim = c(0,350)) + scale_fill_manual(values = c("#3E3E3E","#D5D5D5"))
condensed_barplot_Indels

ggsave("Figures/Main/Figure2/Figure_2D_indel_barplot.pdf", condensed_barplot_Indels, width = 4, height = 4)

