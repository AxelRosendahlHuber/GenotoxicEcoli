##########Documentation for the indel and dbs profiles############
library(tidyverse)
library(magrittr)
library(VariantAnnotation)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg19)
library(MutationalPatterns)
library(reshape2)
source("~/surfdrive/Shared/Boxtel_General/Scripts/pmc_vanboxtel/Axel_BMseq/Utils.R")
source("~/surfdrive/Shared/Boxtel_General/Scripts/pmc_vanboxtel/Freek_misc/id_context.R")
ref_genome = "BSgenome.Hsapiens.UCSC.hg19"
setwd("~/surfdrive/Shared/Projects/Axel/Axel_GenoEcoli/I5/STE0072/Indels/")

####_____________________Indel___________________####
#List files and patients.
vcf_fnames = list.files("~/surfdrive/Shared/Boxtel_General/Data/Mutation_data/INDELs/GenoEcoli/I5/", recursive = T, full.names = F)
vcf_fnames_full = paste0("~/surfdrive/Shared/Boxtel_General/Data/Mutation_data/INDELs/GenoEcoli/I5/", vcf_fnames)
sample_names = gsub("_.*$", "", vcf_fnames)
head(sample_names)

#Vector used for pooling mutations.
patient = substr(sample_names, start = 1, stop = 9)

#Read in indels. (You can't use read_vcfs_as_granges as this function removes indels)
chroms = paste0("chr", c(1:22, "X"))
vcfs_list = lapply(vcf_fnames_full, function(x){
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

ggsave("Plots/Total_indels_STE0072I5.pdf", indel_context_fig, width = 6, height = 8)
ggsave("Plots/Total_indels_detail_STE0072I5.pdf", total_indel_context, width = 10, height = 15)
ggsave("Plots/Total_indels_detail_STE0072I5.png", total_indel_context, width = 10, height = 15)

# get consensus   
indel_counts$DYE = rowSums(indel_counts[,3:5])
indel_counts$EKO = rowSums(indel_counts[,6:8]) 
indel_counts$EWT = rowSums(indel_counts[,9:11])
indel_counts$subtracted = indel_counts$EWT - indel_counts$EKO
indel_counts$subtracted[indel_counts$subtracted < 0] = 0

indel_subtracted = indel_counts[,9:11]
indel_subs = cbind(indel_counts[,1:2], indel_subtracted)
names(indel_subs)[3:5] = c("deltaPKS", "coliWT", "Subtracted")
indel_consensus_reordered = plot_indel_contexts(indel_subs)
indel_consensus_reordered_same_y = plot_indel_contexts(indel_subs, same_y = T)
ggsave("Plots/Total_consensus_indels_reordered.pdf", indel_consensus_reordered, width = 15, height = 5)
ggsave("Plots/Total_consensus_indels_reordered.pdf", indel_consensus_reordered_same_y, width = 15, height = 5)

# Get indel subtraction levels: 
indel_counts$EKOminDYE = indel_counts$EKO- indel_counts$DYE
indel_counts$EKOminDYE[indel_counts$EKOminDYE < 0] = 0

indel_counts$EWTminDYE = indel_counts$EWT- indel_counts$DYE
indel_counts$EWTminDYE[indel_counts$EWTminDYE < 0] = 0

indel_subtraction_plot = plot_indel_contexts(indel_counts[,c(1:2, 16, 17)],same_y = T) + ylim(c(0,120))
ggsave("Plots/Total_consensus_indels_subtraction_plot.pdf", indel_subtraction_plot, width = 15, height = 5)

totals_df = data.frame(DYE = indel_counts$DYE, EKO = indel_counts$EKO, EWT = indel_counts$EWT)
cossim_matrix = cos_sim_matrix(totals_df, totals_df)
cossim_matrix[upper.tri(cossim_matrix)] = NA
cosine_plot = plot_cosine_heatmap(cossim_matrix, cluster_rows = F, plot_values = T)
ggsave("Plots/Cosine_matrix_indels.pdf", cosine_plot, width = 5, height = 5)

# plot samples in indel rainfall-plot 
plot_rainfall_indels(grl_exp[[4]], chromosomes)
chromosomes = seqnames(get(ref_genome))[1:22] 
for (i in 1:length(grl_exp)){
  sample_name = strsplit(vcf_fnames[i], split = "_")[[1]][1]
  plot = plot_rainfall_indels(grl_exp[[i]], 
                        title = sample_name, 
                        chromosomes = chromosomes, 
                        cex = 1.5, ylim = 1e+09)
  plot
  ggsave(paste0("Indels/Plots/Rainfall_plots/", sample_name, ".pdf"), plot, width = 10, height = 5, useDingbats = F)
}

# Define indel signature 
vcfs_total_EWT = unlist(vcfs_grl[4:6])  %>% GRangesList()
indel_list_total_EWT = plot_indels(vcfs_total_EWT,"/Users/ahuber/surfdrive/Shared/Projects/Axel/Axel_GenoEcoli/I5/STE0072/Indels/Plots/", "total_EWT" )
# ToDO count indels function instead of plot_indels
vcfs_total_EKO = unlist(vcfs_grl[1:3])  %>% GRangesList()
indel_list_total_EKO = plot_indels(vcfs_total_EKO,"/Users/ahuber/surfdrive/Shared/Projects/Axel/Axel_GenoEcoli/I5/STE0072/Indels/Plots/", "total_EKO" )

# vcfs_total_EWT_non_unique = unlist(vcfs_grl[1:2]) %>% GRangesList()
# indel_list_total_EWT = plot_indels(vcfs_total_EWT_non_unique,"/Users/ahuber/surfdrive/Shared/Projects/Axel/Axel_GenoEcoli/I5/STE0072/Indels/Plots/", "totalnon_unique" )

subtracted_indels = indel_list_total_EWT[[1]]$My_sample - indel_list_total_EKO[[1]]$My_sample
indel_list_subtracted = indel_list_total_EWT
subtracted_indels[subtracted_indels < 0 ] = 0
indel_list_subtracted$total$My_sample = subtracted_indels
indel_context_fig = plot_indel_contexts(indel_list_subtracted[[1]]) + ggtitle("indels_subtracted")
ggsave("Plots/indels_subtracted_indels.pdf", height = 5, width = 13)
ggsave("Plots/indels_subtracted_indels.png", height = 5, width = 13)

# Save ID_pks signature 
subtracted_indels_rel = subtracted_indels/sum(subtracted_indels)
write.table(subtracted_indels_rel, "ID_pks.txt", sep = "\t", quote = F)


# Plot condensed barplot
#condensed barplot
del_counts = indel_counts[,-1:-2]
simple_indel_counts = data.frame(singleTdels_Thomopolymer = colSums(del_counts[7:12,]), other_indels = colSums(del_counts[-7:-12,]) )
simple_indel_counts$name = rownames(simple_indel_counts)
indel_counts_m = melt(simple_indel_counts)
indel_means = indel_counts_m[c(7:8, 16:17),]
indel_samples = indel_counts_m[c(1:6, 10:15),]

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

condensed_barplot_Indels = ggplot(indel_means, aes(x = name, y = indel_means, fill = variable )) + 
  geom_bar(stat = "identity") + geom_errorbar(aes(min = min, max = max),width = 0.3) +
  theme_BM() +  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_y_continuous(expand = c(0, 0)) + ylab("Number of indels") + xlab("") + coord_cartesian(ylim = c(0,300)) + scale_fill_manual(values = c("#3E3E3E","#D5D5D5"))
condensed_barplot_SNVs

ggsave("Plots/Condensed_barplot.png", condensed_barplot_Indels, width = 4, height = 4)
ggsave("Plots/Condensed_barplot.pdf", condensed_barplot_Indels, width = 4, height = 4)
 








