# Code accompanying the publication # "Mutational signature in colorectal cancer induced by genotoxic pks+ E.Coli"
# C. Pleguezuelos-Manzano, J. Puschhof and A. Rosendahl Huber et al.

# Description: Analysis of deltaClbQ-E.Coli recomplemented with ClbQ (Extended Data Figure 3E-H)
# Please direct any questions or comments to A.K.M.RosendahlHuber@prinsesmaximacentrum.nl

####_____________________Indel___________________####
#List files
vcf_fnames = list.files("Data/STE0072/I3/Indels", recursive = T, full.names = T)
ID_pks = read.delim("Output/ID_pks.txt")
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
total_indel_context = plot_indel_contexts(indel_counts, same_y = T) + theme(legend.position = "none") + scale_y_continuous(limits = c(0,50))
ggsave("Figures/ExtendedData/Extended_data_figure_3/ED_Fig3g.pdf", total_indel_context, width = 13, height = 7)

# Plot total indel levels and subtract I5-indel background 
indel_counts$total = rowSums(indel_counts[3:5]) %>% prop.table()
indel_counts_background = read.delim("/Users/ahuber/surfdrive/Shared/Projects/Axel/Axel_GenoEcoli/I5/STE0072/Indels/EKO_background_indels.txt")
indel_counts$zsubtracted = indel_counts$total - indel_counts_background$My_sample
indel_counts$zsubtracted[indel_counts$zsubtracted < 0] = 0 

indel_context_fig_totals = plot_indel_contexts(indel_counts[,-3:-5], same_y = T)
indel_context_fig_totals

cos_sim(indel_counts$zsubtracted,ID_pks$x)
ggsave("Figures/ExtendedData/Extended_data_figure_3/ED_Fig3h.pdf", indel_context_fig_totals, width = 15, height = 4)

# ---- SNVs 
# Axel Rosendahl Huber 
SBS_pks = read.delim("Output/SBS_pks_subtracted.txt")
vcf_list = list.files("Data/STE0072/I3/SNVs/", pattern = ".vcf", full.names = T)
names = gsub("_72C12TM0_Q50_CGQ10_SGQ99_PASS_20X_VAF0.3_NoY_nonBlacklist_final.vcf", "", basename(vcf_list))
vcfs = read_vcfs_as_granges(vcf_files = vcf_list,names, genome = ref_genome)
mut_mat = mut_matrix(vcfs, ref_genome = ref_genome)
plot = plot_96_profile3(mut_mat,ymax = 100, condensed = F) + scale_y_continuous(expand = c(0,0)) + theme(panel.spacing = unit(1.5, "lines"))  + ylab("Number of Single base substitutions")
ggsave("Figures/ExtendedData/Extended_data_figure_3/ED_Fig3e.pdf", plot, width = 7, height = 4.5)

mut_mat_rel = prop.table(mut_mat, 2) %>% as.data.frame()

STE0072I5 = read.delim("Output/I5_mut_mat.txt")
STE0072I5$EKO = prop.table(rowSums(STE0072I5[,4:6],2))

mut_mat_rel$total_EWQ = prop.table(rowSums(mut_mat_rel[,1:3]))
mut_mat_rel$subtracted = mut_mat_rel$total_EWQ- STE0072I5$EKO
mut_mat_rel[mut_mat_rel<0] = 0
mut_mat_rel$subtracted = prop.table(mut_mat_rel$subtracted)

plot = plot_96_profile3(mut_mat_rel[,4:5], condensed = F)
cos_sim_matrix(mut_mat_rel, SBS_pks)

ggsave("Figures/ExtendedData/Extended_data_figure_3/ED_Fig3f.pdf", plot, width = 7, height = 4)