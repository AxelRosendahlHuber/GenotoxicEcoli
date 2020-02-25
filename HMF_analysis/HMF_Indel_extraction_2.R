##########Documentation for the indel and dbs profiles############
library(tidyverse)
library(magrittr)
library(VariantAnnotation)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg19)
source("/hpc/pmc_vanboxtel/tools/scripts/Axel_misc/functions/id_context.R")
ref_genome = "BSgenome.Hsapiens.UCSC.hg19"
setwd("/hpc/pmc_vanboxtel/projects/Axel_GenoEcoli/HMF_analysis/Indels")

####_____________________Indel___________________####
#List files and patients.
vcf_fnames = list.files("/hpc/pmc_vanboxtel/projects/Axel_GenoEcoli/HMF_analysis/somatics/", recursive = T, pattern = "INDEL.vcf$", full.names = T)
print(length(vcf_fnames))
sample_names = basename(vcf_fnames) %>% gsub(pattern = "\\..*$",replacement =  "")


#Read in indels. (You can't use read_vcfs_as_granges as this function removes indels)
chroms = paste0("chr", c(1:22, "X"))
vcfs_list = lapply(vcf_fnames, function(x){
  print(x)
  vcf = readVcf(x, genome = "hg19")
  gr = granges(vcf)
  gr = gr[gr$FILTER == "PASS",]
 # gr = remove_snvs(gr)
  seqlevelsStyle(gr) = "UCSC"
  seqlevels(gr, pruning.mode = "coarse") = chroms
  return(gr)
})

vcfs_grl = GRangesList(vcfs_list)
names(vcfs_grl) = sample_names
grl_exp = get_indel_context(vcfs_grl, ref_genome)
indel_counts = count_indel_contexts(grl = grl_exp)
write.table(file = "Indel_counts_PASS.txt", x = indel_counts, sep = "\t", quote = F, row.names = T, col.names = T)

