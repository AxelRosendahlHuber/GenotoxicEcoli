# Code accompanying the publication # "Mutational signature in colorectal cancer induced by genotoxic pks+ E.Coli"
# C. Pleguezuelos-Manzano, J. Puschhof and A. Rosendahl Huber et al.

# Description: Set working directory, analyze > 1bp deletions (Figure 3C)
# Please direct any questions or comments to A.K.M.RosendahlHuber@prinsesmaximacentrum.nl

# Load data
setwd(dir)
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

STE72_I5_dels = get_length_pks_indels_all_grls_mc(grl_exp, 1)
# --- multi_indel plot
# ---- EWT vs EKO 
STE72_I5_dels$pks_status = factor(STE72_I5_dels$pks_status, levels = c("pks_pattern_present", "no_pks_pattern"))
intervals = c("2bp_deletions", "3-5bp deletions", "6-10bp deletions", "11-15bp deletions", "15bp+ deletions")
STE72_I5_dels$interval = intervals[findInterval(STE72_I5_dels$indel_length, c(0,3, 6, 11, 16))]
STE72_I5_dels$interval = factor(STE72_I5_dels$interval, levels = c("2bp_deletions", "3-5bp deletions", "6-10bp deletions", "11-15bp deletions", "15bp+ deletions"))
STE72_I5_dels$PKS_exposed = ifelse(grepl("EWT",STE72_I5_dels$sample), "pks", "delta/dye")
STE72_I5_dels$PKS_exposed = factor(STE72_I5_dels$PKS_exposed, levels = c("pks", "delta/dye"))

labels = c("non_exposed" = paste0("Non-pks_exposed clones n=", length(unique(STE72_I5_dels[STE72_I5_dels$PKS_exposed == "non_exposed", "sample"]))), 
           "exposed" = paste0("pks_exposed clones n=",length(unique(STE72_I5_dels[STE72_I5_dels$PKS_exposed == "exposed", "sample"]))))

STE72_I5_dels_EWT = STE72_I5_dels[STE72_I5_dels$PKS_exposed == "pks",]
STE72_I5_dels_EKO = STE72_I5_dels[grepl("EKO", STE72_I5_dels$sample),]

counts_EWT = table(STE72_I5_dels_EWT[,c(2,4)])
counts_EWT = melt(counts_EWT)
counts_EKO = table(STE72_I5_dels_EKO[,c(2,4)])
counts_EKO = melt(counts_EKO)

counts_EWT$PKS_exposed = "pks"
counts_EKO$PKS_exposed = "delta"
counts_EKO_EWT = rbind(counts_EWT, counts_EKO)

# ------ EWT /EKO /  DYE -----------
STE72_I5_dels_DYE = STE72_I5_dels[grepl("DYE", STE72_I5_dels$sample),]
counts_DYE = table(STE72_I5_dels_DYE[,c(2,4)])
counts_DYE = melt(counts_DYE)
counts_DYE$PKS_exposed = "dye"
counts_STE72_total = rbind(counts_EWT, counts_EKO, counts_DYE)
counts_STE72_total$PKS_exposed = factor(counts_STE72_total$PKS_exposed, levels = c("dye", "delta", "pks"))

# Indel plot with error bars 
row <- apply(STE72_I5_dels[,2:4], 1, paste0, collapse = "_");
counts = STE72_I5_dels[,2:4] %>% group_by(sample, interval, pks_status) %>% summarise(n()) %>% ungroup()
counts_per_sample = dcast(counts, formula = counts$sample ~ counts$interval + pks_status)
counts_per_sample[is.na(counts_per_sample)] = 0

counts_per_sample$type = c("DYE", "DYE", "DYE", "EKO", "EKO", "EKO", "EWT", "EWT", "EWT")
means = aggregate(counts_per_sample[,2:11], by = list(counts_per_sample$type), mean)
rownames(means) = means$Group.1
means = t(means)[-1,]
mode(means) = "numeric"

sds = aggregate(counts_per_sample[,2:11], by = list(counts_per_sample$type), sd)
rownames(sds) = sds$Group.1
sds = t(sds)[-1,]
mode(sds) = "numeric"

sds_max  = means + sds
sds_min = means - sds
sds_max[seq(1, nrow(sds), 2 ),] = sds_max[seq(1, nrow(sds), 2 ),] + means[seq(2, nrow(means), 2 ),]
sds_min[seq(1, nrow(sds), 2 ),] = sds_min[seq(1, nrow(sds), 2 ),] + means[seq(2, nrow(means), 2 ),]

means_m = melt(means)
means_m$max = melt(sds_max)$value
means_m$min = melt(sds_min)$value
means_m$min[means_m$min < 0] = 0
means_m$pks_status = gsub("^.*deletions_", "", means_m$Var1)
means_m$interval = gsub("_pks_pattern_present|_no_pks_pattern", "", means_m$Var1)
means_m$pks_status = factor(means_m$pks_status, levels = c("pks_pattern_present", "no_pks_pattern") )
means_m$interval = factor(means_m$interval, levels = c("2bp_deletions", "3-5bp deletions", "6-10bp deletions", "11-15bp deletions", "15bp+ deletions"))

counts_per_sample_added = counts_per_sample[4:9,]

add_to_EKO = counts_per_sample_added[1:3,grep("present", colnames(counts_per_sample_added))] + matrix(rep(means_m$value[c(12,14,16,18,20)], 3), nrow = 3, byrow = T)
add_to_EWT = counts_per_sample_added[4:6,grep("present", colnames(counts_per_sample_added))] + matrix(rep(means_m$value[c(22,24,26,28,30)], 3), nrow = 3, byrow = T)

counts_per_sample_added[1:3,grep("present", colnames(counts_per_sample_added))] = add_to_EKO
counts_per_sample_added[4:6,grep("present", colnames(counts_per_sample_added))] = add_to_EWT
counts_added_m = melt(counts_per_sample_added)
counts_added_m$pks_status = gsub("^.*deletions_", "", counts_added_m$variable)
counts_added_m$interval = gsub("_pks_pattern_present|_no_pks_pattern", "", counts_added_m$variable)
counts_added_m$pks_status = factor(counts_added_m$pks_status, levels = c("pks_pattern_present", "no_pks_pattern") )
counts_added_m$interval = factor(counts_added_m$interval, levels = c("2bp_deletions", "3-5bp deletions", "6-10bp deletions", "11-15bp deletions", "15bp+ deletions"))
colnames(counts_added_m) = c("Var1", "Var2", "variable", "value", "pks_status", "interval")

# Plot containing error-bars and dots indicating individual samples
means_m_noDYE = means_m[!means_m$Var2 == "DYE",]
errorbar_plot_noDYE = ggplot(means_m_noDYE, aes( y = value, x= Var2, fill  = pks_status)) + 
  geom_bar(stat = "identity", color = "black") + 
  geom_point(data = counts_added_m, mapping = aes(x = Var2, y = value, fill = pks_status), color = "black", shape = 21,size = 2.3,position = position_jitterdodge()) +
  facet_grid(.~interval) + geom_errorbar(aes(ymin = min, ymax = max), size = 0.75, position = position_dodge(0.9)) + 
  scale_fill_manual(values = c("#262626", "#e6e6e6")) + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, 18)) + 
  theme_classic() + ylab("Deletion count") + xlab(element_blank()) + 
  theme(panel.spacing.x=unit(0.05, "lines"),panel.spacing.y=unit(1, "lines"))
# NOTE: Plot will result in an error - this is due to errorbars going into negative values on y-axis. Nothing to worry about
errorbar_plot_noDYE
ggsave("Figures/Main/Figure3/Figure_3C_Larger_deletions.pdf" ,errorbar_plot_noDYE, width = 8, height = 5, useDingbats = F)

