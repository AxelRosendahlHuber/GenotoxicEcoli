# Code accompanying the publication # "Mutational signature in colorectal cancer induced by genotoxic pks+ E.Coli"
# C. Pleguezuelos-Manzano, J. Puschhof and A. Rosendahl Huber et al.

# Please direct any questions or comments to A.K.M.RosendahlHuber@prinsesmaximacentrum.nl

# NOTE: This script uses patient-level somatic variant and clinical data obtained from the Hartwig Medical Foundation. 


# ---- HMF -----
# NOTE: This script uses patient-level somatic variant and clinical data which have been obtained from the Hartwig Medical Foundation under the data request number DR-084. Somatic variant and clinical data are freely available for academic use from the Hartwig Medical Foundation through standardized procedures. Privacy and publication policies, including co-authorship policies, can be retrieved from: https://www.hartwigmedicalfoundation.nl/en/data-policy/. 
# Data request forms can be downloaded from https://www.hartwigmedicalfoundation.nl/en/applying-for-data/.
# To gain access to the data, this data request form should be emailed to info@hartwigmedicalfoundation.nl., upon which it will be evaluated within 6 weeks by the HMF Scientific Council and an independent Data Access Board.
# When access is granted, the requested data become available through a download link provided by HMF.


##########Documentation for the indel and dbs profiles############
setwd(dir)
source("Utils.R")
ref_genome = "BSgenome.Hsapiens.UCSC.hg19"

####_____________________Indel___________________####
#List files and patients.
vcf_fnames = list.files("HMF vcf folder", recursive = T, pattern = "somatic.vcf.gz$", full.names = T)
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

