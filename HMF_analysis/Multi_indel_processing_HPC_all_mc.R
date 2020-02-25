# Get Multi-indels for the entire HMF-dataset. 
source("/hpc/pmc_vanboxtel/tools/scripts/Axel_misc/pks_context_selection_functions_v2.R")
source("/hpc/pmc_vanboxtel/tools/scripts/Axel_misc/Multi_bp_indel_pks_classification.R")
source("/hpc/pmc_vanboxtel/tools/scripts/Axel_misc/functions/id_context.R")
library(VariantAnnotation)
library(GenomicRanges)
library(tidyverse)
library(BSgenome)
library(BSgenome.Hsapiens.UCSC.hg19)
library(reshape2)
library(pracma)
library(magrittr)
ref_genome = "BSgenome.Hsapiens.UCSC.hg19"
n_cores = 8

setwd("/hpc/pmc_vanboxtel/projects/Axel_GenoEcoli/HMF_analysis/Multi_indel/")

# 31 CRC double positive samples 
vcf_fnames = list.files("/hpc/pmc_vanboxtel/projects/Axel_GenoEcoli/HMF_analysis/somatics/", pattern = "*somatic.vcf$", recursive = T, full.names = T)
sample_names = basename(vcf_fnames) %>% gsub(pattern = ".purple.somatic.vcf", replacement = "")

#Read in indels. (You can't use read_vcfs_as_granges as this function removes indels)
chroms = paste0("chr", c(1:22, "X"))
indel_status_list = list()

index = 1:length(vcf_fnames)
index = split(index, ceiling(seq_along(index)/48))
names(index) = as.character(1:length(index))
for (name in names(index)) {
  i = index[[name]]
  vcfs_list = mclapply(vcf_fnames[i],  function(x){
    vcf = readVcf(x, genome = "hg19")
    gr = granges(vcf)
    gr = gr[elementNROWS(gr$ALT) ==1]
    gr = gr[gr$FILTER == "PASS"]
    gr = remove_snvs(gr)
    seqlevelsStyle(gr) = "UCSC"
    seqlevels(gr, pruning.mode = "coarse") = chroms
    return(gr)
  }, mc.cores = n_cores)
  vcfs_grl = GRangesList(vcfs_list)
  names(vcfs_grl) = sample_names[i]
  grl_exp = get_indel_context_mc(vcfs_grl, ref_genome, n_cores = n_cores)
  indel_status = get_length_pks_indels_all_grls_mc(grl_exp, n_cores = n_cores)
  indel_status_list[[name]] = indel_status
}

indel_status = do.call("rbind", indel_status_list)
write.table(indel_status, "indel_status_mc.txt", sep = "\t", quote = F)



  
