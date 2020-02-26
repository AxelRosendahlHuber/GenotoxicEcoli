# Code accompanying the publication # "Mutational signature in colorectal cancer induced by genotoxic pks+ E.Coli"
# C. Pleguezuelos-Manzano, J. Puschhof and A. Rosendahl Huber et al.

# Description: Analysis of the wider context of SBS and T-deletion mutations (Figure 3A-B)
# Please direct any questions or comments to A.K.M.RosendahlHuber@prinsesmaximacentrum.nl

# ===== LOAD LOCAL FUNCTIONS AND SET VARIABLES =====
# Source script which load all preset functions and variables
library(ggseqlogo)
setwd(dir)
options(stringsAsFactors = F)

#  ------ Context SNVs
# ---- Functions ----- 
get_context = function(vcf, size_context){
  gr = rowRanges(vcf)
  strand = ifelse(ref(vcf) == "G" | ref(vcf) == "A", '-', "+")
  start = start(vcf)
  ref = as.character(gr$REF)
  alt = unlist(CharacterList(gr$ALT))
  type = paste0(ref, "/", alt)
  chromosome = paste0("chr", as.character(seqnames(vcf)))
  context = getSeq(Hsapiens, chromosome,  start = start - size_context,
                   end =  start + size_context, 
                   strand = strand)
  
  context_table = data.frame(chr = chromosome, position = start, type = type, strand = strand, context = as.character(context))
  rownames(context_table) = names(vcf)
  return(context_table)
}

# Load data 
vcfs_path <- get_vcfs("Data/STE0072/I5/SNVs/", pattern = '.vcf')

vcf_list = list()
for (i in 1:nrow(vcfs_path)){
  vcf = readVcf(vcfs_path$vcf_files[i], ref_genome)
  vcf_list[[vcfs_path$sample_names[i]]] = vcf
}

grl_list = lapply(vcf_list, granges)
vcfs_grl = GRangesList(grl_list)


DYE_context = lapply(vcf_list[1:3], get_context, 10)
DYE_context = do.call("rbind", EKO_context)
EKO_context = lapply(vcf_list[4:6], get_context, 10)
EKO_context = do.call("rbind", EKO_context)
EWT_context = lapply(vcf_list[7:9], get_context, 10)
EWT_context = do.call("rbind", EWT_context)

EWT_context_table = EWT_context$context

EWT_TN_context = EWT_context[grep(pattern = "A/|T/", EWT_context$type),]
EKO_TN_context = EKO_context[grep(pattern = "A/|T/", EKO_context$type),]
EKO_TN = ggseqlogo(EKO_TN_context$context) + scale_y_continuous(expand = c(0,0), limits = c(0, 2)) + theme_classic()
EWT_TN = ggseqlogo(EWT_TN_context$context) + scale_y_continuous(expand = c(0,0), limits = c(0, 2)) + theme_classic()
ggsave("Figures/Main/Figure3/Figure_3A_SeqLogo_bits_EWT.pdf", EWT_TN, width = 7, height = 3)
EWT_TN

# ----- Extended Context Indels ---- 
vcf_fnames = list.files("Data/STE0072/I5/Indels/", recursive = T, full.names = F)
vcf_fnames_full = paste0("Data/STE0072/I5/Indels/", vcf_fnames)
sample_names = gsub("_.*$", "", vcf_fnames)
head(sample_names)

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

# update the get_1bp_dels function to make sure no T-homopolymer lengths are missed. 
get_1bp_dels = function(gr, mut_size, ref_genome){
  gr = gr[mut_size == -1]
  if (length(gr) == 0){
    return(gr)
  }
  
  del_bases = gr$REF %>% as.character() %>% substring(2)
  withCallingHandlers({#Flank the granges object, to get a sequence, that can be searched for repeats. This can result in a warning message, when the flanked range extends beyond the chrom lenght. This message is suppressed.
    gr_extended = flank(gr, 300, start = F)
  }, warning = function(w) {
    if (grepl("out-of-bound range located on sequence", conditionMessage(w)))
      invokeRestart("muffleWarning")
  })
  gr_extended %<>% trim() #Trim the ranges that are extended beyond the actual length of the chromosome.
  seq = getSeq(eval(as.symbol(ref_genome)), gr_extended)
  
  seq_z = str_replace_all(as.character(seq), del_bases, rep("Z", length(seq))) #For each mut replace the deleted basetype in the flanking sequence with Zs.
  homopolymer_length = gsub("[^Z].*", "", as.character(seq_z)) %>% nchar() + 1 #Remove all bases after the Zs and count how many bases are left. Add one for the deleted base itself.
  del_bases[del_bases == "A"] = "T"
  del_bases[del_bases == "G"] = "C"
  
  gr$muttype = paste0(del_bases, "_deletion")
  gr$muttype_sub = homopolymer_length
  return(gr)
}

approach2_id_context_table = function(gr_exp) {
  size_context = 10
  # Remove T-homopolymer
  gr_exp = gr_exp[gr_exp$muttype == "T_deletion",]
  homopolymer_length = gr_exp$muttype_sub
  muts = length(gr_exp)
  start = start(gr_exp)
  end = end(gr_exp)
  ref = as.character(gr_exp$REF)
  chr = as.character(seqnames(gr_exp))
  strand = ifelse(substr(ref, 2, 2) == "A","-","+")
  total_context = c()
  for (i in 1:length(strand)){
    print(i)
    if (strand[i] == "+") {
      upstr_context = as.character(getSeq(Hsapiens, chr[i],
                                          start = start[i] - size_context + 1, # add +1 value since start base is the base before the homopolymer
                                          end =  start[i], 
                                          strand = strand[i]))
      downstr_context = as.character(getSeq(Hsapiens, chr[i],
                                            start = end[i],
                                            end =  end[i] + homopolymer_length[i] + size_context, 
                                            strand = strand[i]))
      downstr_context = substr(downstr_context, 1, 11)
      
      context = paste0(upstr_context, downstr_context)
    }
    if (strand[i] == "-") {
      downstr_context = as.character(getSeq(Hsapiens, chr[i],
                                            start = start[i] + homopolymer_length[i] - size_context, 
                                            end =  start[i] + homopolymer_length[i], 
                                            strand = strand[i]))
      downstr_context = substr(downstr_context, 1, 11)
      upstr_context = as.character(getSeq(Hsapiens, chr[i],
                                          start = start[i] + homopolymer_length[i] + 1 ,# add +1 value since start base is the base before the homopolymer
                                          end =  start[i] + homopolymer_length[i]  + size_context, 
                                          strand = strand[i]))
      context = paste0(upstr_context, downstr_context)
    }
    
    total_context = c(total_context, context)
    
  }
  return(total_context)
}

total_context_EWT = lapply(grl_exp[7:9], approach2_id_context_table) 
total_context_EWT = unlist(total_context_EWT)
Logo_plot_EWT = ggseqlogo(total_context_EWT) + theme_BM() +scale_y_continuous(expand = c(0,0), limits = c(0, 2))
Logo_plot_EWT

ggsave("Figures/Main/Figure3/Figure_3B_Seqlogo_indels_non_condensed_EWT.pdf", Logo_plot_EWT,  width = 7, height = 3)

# ---- Plot Extended data figure 4 ------- 
# - approach 2 with different indel sizes: 
grl_exp_EWT = GRangesList(unlist(grl_exp[7:9]))
names(grl_exp_EWT) = "Total_Tdels_EWT"

approach_2_id_context_length = function(grl_exp, length_hp) {
  for (name in names(grl_exp)){
    gr_exp = grl_exp[[name]]
    size_context = 10
    gr_exp = gr_exp[gr_exp$muttype == "T_deletion",]
    gr_exp = gr_exp[gr_exp$muttype_sub == length_hp,]
    homopolymer_length = gr_exp$muttype_sub
    start = as.numeric(start(gr_exp))
    end = as.numeric(end(gr_exp))
    chr = as.character(seqnames(gr_exp))
    ref = gr_exp$REF
    strand = ifelse(substr(ref, 2, 2) == "A","-","+")
    total_context = c()
    for (i in 1:length(strand)){
      if (strand[i] == "+") {
        upstr_context = as.character(getSeq(Hsapiens, chr[i],
                                            start = start[i] - size_context + 1, # add +1 value since start base is the base before the homopolymer
                                            end =  start[i], 
                                            strand = strand[i]))
        downstr_context = as.character(getSeq(Hsapiens, chr[i],
                                              start = end[i],
                                              end =  end[i] + homopolymer_length[i] + size_context, 
                                              strand = strand[i]))
        downstr_context = substr(downstr_context, 1, 11)
        
        context = paste0(upstr_context, downstr_context)
      }
      if (strand[i] == "-") {
        downstr_context = as.character(getSeq(Hsapiens, chr[i],
                                              start = start[i] + homopolymer_length[i] - size_context, # add +1 value since start base is the base before the homopolymer
                                              end =  start[i] + homopolymer_length[i], 
                                              strand = strand[i]))
        downstr_context = substr(downstr_context, 1, 11)
        upstr_context = as.character(getSeq(Hsapiens, chr[i],
                                            start = start[i] + homopolymer_length[i] + 1 ,
                                            end =  start[i] + homopolymer_length[i]  + size_context, 
                                            strand = strand[i]))
        context = paste0(upstr_context, downstr_context)
      }
      
      total_context = c(total_context, context)
    }
    plot_mutation_profile(total_context, paste0("Figures/ExtendedData/Extended_data_figure_4/",name, "T_homopolymer", as.character(length_hp),"approach2"))
  }
}

# Plot the 10bp mutational context for deletions in T-homopolymers of length 1 till 8 bp
approach_2_id_context_length(grl_exp_EWT, 1)
approach_2_id_context_length(grl_exp_EWT, 2)
approach_2_id_context_length(grl_exp_EWT, 3)
approach_2_id_context_length(grl_exp_EWT, 4)
approach_2_id_context_length(grl_exp_EWT, 5)
approach_2_id_context_length(grl_exp_EWT, 6)
approach_2_id_context_length(grl_exp_EWT, 7)
approach_2_id_context_length(grl_exp_EWT, 8)