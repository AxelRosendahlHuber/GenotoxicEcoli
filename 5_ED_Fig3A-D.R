# Code accompanying the publication # "Mutational signature in colorectal cancer induced by genotoxic pks+ E.Coli"
# C. Pleguezuelos-Manzano, J. Puschhof and A. Rosendahl Huber et al.

# Description: Analyze single base substitutions and indels of second clone (Extended Data Figure 3A-D)
# Please direct any questions or comments to A.K.M.RosendahlHuber@prinsesmaximacentrum.nl

# STE0076 Analysis
setwd(dir)
vcf_list = list.files("Data/STE0076/SNVs/", pattern = ".vcf", full.names = T)
names = gsub("_STE0076TM0_Q50_CGQ10_SGQ99_PASS_20X_VAF0.3_NoY_nonBlacklist_final.vcf", "", basename(vcf_list))
vcfs = read_vcfs_as_granges(vcf_files = vcf_list,names, genome = ref_genome)
mut_mat = mut_matrix(vcfs, ref_genome = ref_genome)
plot = plot_96_profile3(mut_mat,ymax = 100, condensed = F) + scale_y_continuous(expand = c(0,0)) + theme(panel.spacing = unit(1.5, "lines"))  + ylab("Number of Single base substitutions")
ggsave("Figures/ExtendedData/Extended_data_figure_3/ED_Fig3a.pdf", plot, width = 7, height = 7)

# import mutational profiles 
SBS_pks = read.delim("Output/SBS_pks_subtracted.txt")
names(SBS_pks) = "SBS-pks_subtracted"
mut_mat = as.data.frame(mut_mat)
mut_mat$total_76C = rowSums(mut_mat[,4:6])
mut_mat$total_76_EKO = rowSums(mut_mat[,1:3])
mut_mat$subtracted = mut_mat$total_76C - mut_mat$total_76_EKO
mut_mat[mut_mat < 0] = 0 

cossim = cos_sim_matrix(mut_mat,SBS_pks)
plot_cosine_heatmap(cossim, plot_values = T, cluster_rows = F)

# Subtract relative levels of in-vitro signature from the 96-trinucleotide contexts. 
profile_totals_76 = plot_96_profile3(mut_mat[ c(7:9)], condensed = F, ymax = 100)
ggsave("Figures/ExtendedData/Extended_data_figure_3/ED_Fig3b.pdf", profile_totals_76, width = 7, height = 6)

# ---- Indels 
#List files and patients.
vcf_fnames = list.files("Data/STE0076/Indels", recursive = T, full.names = T)
sample_names = gsub("_.*$", "", basename(vcf_fnames))
ID_pks = read.delim("Output/ID_pks.txt")

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

ggsave("Figures/ExtendedData/Extended_data_figure_3/ED_Fig3c.pdf", total_indel_context, width = 10, height = 8)

# Plot total indel levels and subtract I5-indel background 
indel_counts$EKO = rowSums(indel_counts[3:5])/3
indel_counts$EWT = rowSums(indel_counts[6:8])/3
indel_counts$subtracted = indel_counts$EWT - indel_counts$EKO
indel_counts$subtracted[indel_counts$subtracted < 0] = 0 

indel_context_fig_totals = plot_indel_contexts(indel_counts[,c("muttype", "muttype_sub", "EWT", "EKO", "subtracted")], same_y = T)
indel_context_fig_totals

cos_sim(indel_counts$subtracted,ID_pks$x)
ggsave("Figures/ExtendedData/Extended_data_figure_3/ED_Fig3d.pdf", indel_context_fig_totals, width = 15, height = 8)