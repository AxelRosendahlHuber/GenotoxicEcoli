# Generate Mut_mat from the HMF data
library(ggplot2)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(pracma)
library(matrixStats)
setwd(dir)
source("Utils.R")
genes_hg19 <- genes(TxDb.Hsapiens.UCSC.hg19.knownGene)

# ---- LOAD FUNCTIONS NEEDED TO OPERATE MULTI-CORE ON THE HPC  ----

mut_96_occurrences2 = function(type_context)
{
    vector = rep(0,96)
    names(vector) = TRIPLETS_96
    
    # if type_context is empty, return vector with zeroes
    if (isEmpty(type_context$types) || isEmpty(type_context$context))
        return(vector)
    
    # for all mutations in this sample
    for (i in 1:length(type_context[[1]]))
    {
        # Find mutation type
        type = which(SUBSTITUTIONS == type_context[[1]][i])

        # Find triplet
        if(type < 4)
            context = which(C_TRIPLETS == type_context[[2]][i])
        else
            context = which(T_TRIPLETS == type_context[[2]][i])

        pos = (type - 1)*16 + context
        vector[pos] = vector[pos] + 1
    }

    return(vector)
}


read_vcfs_as_granges2 = function (vcf_files, sample_names, genome, group = "auto+sex",
                                  check_alleles = TRUE)
{
  if (length(vcf_files) != length(sample_names))
    stop("Please provide the same number of sample names as VCF files")
  ref_genome <- base::get(genome)
  ref_organism <- GenomeInfoDb::organism(ref_genome)
  ref_style <- seqlevelsStyle(ref_genome)
  genome_name <- genome(ref_genome)[[1]]
  if (!(class(ref_genome) == "BSgenome"))
    stop("Please provide the name of a BSgenome object.")
  num_cores = 8
  original_warn_state = getOption("warn")
  options(warn = -1)
  warnings <- NULL
  if (!check_alleles) {
    warnings <- c(warnings, paste("check_alleles is set to FALSE.  Make sure your",
                                  "input VCF does not contain any positions with",
                                  "insertions, deletions or multiple alternative",
                                  "alleles, as these positions cannot be analysed",
                                  "with MutationalPatterns and cause obscure", "errors."))
  }
  vcf_list <- mclapply(seq_along(vcf_files), function(index) {
    file <- vcf_files[index]
    vcf <- rowRanges(readVcf(file, genome_name))
    seqlevelsStyle(vcf) <- ref_style[1]
    groups <- c()
    if (group != "none") {
      if (group == "auto+sex") {
        groups <- c(extractSeqlevelsByGroup(species = ref_organism,
                                            style = ref_style, group = "auto"), extractSeqlevelsByGroup(species = ref_organism,
                                                                                                        style = ref_style, group = "sex"))
        groups_names <- names(groups)
        if (!is.null(groups_names)) {
          unique_names <- unique(groups_names)
          groups <- llply(unique_names, function(x) groups[groups_names ==
                                                             x])
          groups <- llply(groups, unlist, recursive = FALSE)
          groups <- unique(as.vector(groups[[1]]))
        }
      }
      else {
        groups <- extractSeqlevelsByGroup(species = ref_organism,
                                          style = ref_style, group = group)
        groups <- unique(as.vector(t(groups)))
      }
      groups <- intersect(groups, seqlevels(vcf))
      vcf <- keepSeqlevels(vcf, groups, pruning.mode = "tidy")
    }
    if (check_alleles) {
      rem <- which(all(!(!is.na(match(vcf$ALT, DNA_BASES)) &
                           !is.na(match(vcf$REF, DNA_BASES)) & (lengths(vcf$ALT) ==
                                                                  1))))
      if (length(rem) > 0) {
        vcf = vcf[-rem]
        warnings <- c(warnings, paste(length(rem), "position(s) with indels and/or multiple",
                                      "alternative alleles are excluded from", paste(sample_names[[index]],
                                                                                     ".", sep = "")))
      }
    }
    return(list(vcf, warnings))
  }, mc.cores = num_cores)
  options(warn = original_warn_state)
  vcf_list <- lapply(vcf_list, function(item) {
    if (class(item) == "try-error")
      stop(item)
    if (!is.null(item[[2]]))
      for (i in item[[2]]) warning(i)
    return(item[[1]])
  })
  vcf_list <- GRangesList(vcf_list)
  names(vcf_list) <- sample_names
  return(vcf_list)
}

mut_matrix2 = function (vcf_list, ref_genome)
{
  df = data.frame()
  num_cores = 8
  rows <- mclapply(as.list(vcf_list), function(vcf) {
    type_context = type_context(vcf, ref_genome)
    row = mut_96_occurrences2(type_context)
    return(row)
  }, mc.cores = num_cores)
  for (row in rows) {
    if (class(row) == "try-error")
      stop(row)
    df = rbind(df, row)
  }
  names(df) = names(row)
  row.names(df) = names(vcf_list)
  return(t(df))
}


# --- READ SAMPLE NAMES
vcfs_names = list.files("HMF vcf folder", pattern = "somatic.vcf.gz$", recursive = T, full.names = T) # HPC
vcfs_sample_names = gsub(".purple.somatic_SNV.vcf", "", basename(vcfs_names))

mut_mat_list = list()
mut_mat_s_list = list()

index = 1:length(vcfs_names)
index = split(index, ceiling(seq_along(index)/96))
names(index) = as.character(1:length(index))
for (name in names(index)) {
  i = index[[name]]
  print(i)
  grls = read_vcfs_as_granges2(vcfs_names[i], vcfs_sample_names[i], ref_genome)
  grls = lapply(grls, FUN = function(x) x[x$FILTER =="PASS",])
  mut_mat = mut_matrix2(grls, ref_genome)
  # mut_mat_s = mut_matrix_stranded2(grls, ref_genome, genes_hg19)
  mut_mat_list[[name]] = mut_mat
  #  mut_mat_s_list[[i]] = mut_mat_s
}


mut_mat = do.call("cbind", mut_mat_list)
#mut_mat_s = do.call("cbind", mut_mat_s_list)

#write.table(mut_mat_s,"Mut_mats/HMF_Mut_mat_s_noPASS.txt", quote = F, sep = "\t")
write.table(mut_mat, "Output/HMF_mut_mat_noPASS.txt", quote = F, sep = "\t")



