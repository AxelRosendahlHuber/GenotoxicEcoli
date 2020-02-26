# Code accompanying the publication # "Mutational signature in colorectal cancer induced by genotoxic pks+ E.Coli"
# C. Pleguezuelos-Manzano, J. Puschhof and A. Rosendahl Huber et al.

# Please direct any questions or comments to A.K.M.RosendahlHuber@prinsesmaximacentrum.nl

# NOTE: This script uses patient-level somatic variant and clinical data obtained from the Hartwig Medical Foundation. 


# ---- HMF -----
# NOTE: This script uses patient-level somatic variant and clinical data which have been obtained from the Hartwig Medical Foundation under the data request number DR-084. Somatic variant and clinical data are freely available for academic use from the Hartwig Medical Foundation through standardized procedures. Privacy and publication policies, including co-authorship policies, can be retrieved from: https://www.hartwigmedicalfoundation.nl/en/data-policy/. 
# Data request forms can be downloaded from https://www.hartwigmedicalfoundation.nl/en/applying-for-data/.
# To gain access to the data, this data request form should be emailed to info@hartwigmedicalfoundation.nl., upon which it will be evaluated within 6 weeks by the HMF Scientific Council and an independent Data Access Board.
# When access is granted, the requested data become available through a download link provided by HMF.

# Get Multi-indels for the entire HMF-dataset. 
setwd(dir)
source("Utils.R")
source("/HMF_analysis/pks_context_selection_functions_v2.R")
n_cores = 8


# 31 CRC double positive samples 
vcf_fnames = list.files("HMF vcf folder", pattern = "*somatic.vcf.gz$", recursive = T, full.names = T)
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



  
