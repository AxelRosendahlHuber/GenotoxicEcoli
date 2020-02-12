get_length_pks_indels = function(gr, context_size = 5, name) {
  gr = gr[grepl("bp_deletion", gr$muttype)]
  
  if (length(gr) == 0) {
    warning("sample does not contain indels")
    output = c(0,0,0)
  } else {
    
    gr = remove_mult_alts(gr)
    chr = seqnames(gr)
    start = start(gr) - context_size %>% as.integer()
    end = start(gr) + context_size %>% as.integer()
    strand = strand(gr) %>% as.character()
    contexts_start = getSeq(Hsapiens, names = chr, start =start , end = end ,strand = strand, as.character = T )
    pks_contexts = ("AAAAT|AAATT|AATTT|ATTTT")
    start_pattern = grepl(pks_contexts, contexts_start)
    
    
    start = end(gr) - context_size %>% as.integer()
    end = end(gr) + context_size %>% as.integer()
    strand = strand(gr) %>% as.character()
    contexts_end = getSeq(Hsapiens, names = chr, start =start , end = end ,strand = strand, as.character = T )
    end_pattern = grepl(pks_contexts, contexts_end)
    
    
    gr$pks_status = ifelse(start_pattern | end_pattern, "pks_at_flank", "no_pks_pattern")
    gr$length = gsub("bp_deletion|bp_deletion_with_microhomology", "", gr$muttype) %>% as.numeric()
    
    indel_table = data.frame(indel_length = gr$length, pks_status = gr$pks_status, sample = name)
    
  }
  return(indel_table)
}

get_length_pks_indels_all = function(gr, context_size = 5, name) {
  gr = gr[grepl("bp_deletion", gr$muttype)]
  
  if (length(gr) == 0) {
    warning("sample does not contain indels")
    output = c(0,0,0)
  } else {
    
    gr = remove_mult_alts(gr)
    chr = seqnames(gr)
    start = start(gr) - context_size %>% as.integer()
    end = end(gr) + context_size %>% as.integer()
    strand = strand(gr) %>% as.character()
    contexts_start = getSeq(Hsapiens, names = chr, start =start , end = end ,strand = strand, as.character = T )
    pks_contexts = ("AAAAT|AAATT|AATTT|ATTTT")
    indels_with_pattern = grepl(pks_contexts, contexts_start)
    
    gr$pks_status = ifelse(indels_with_pattern, "pks_pattern_present", "no_pks_pattern")
    gr$length = gsub("bp_deletion|bp_deletion_with_microhomology", "", gr$muttype) %>% as.numeric()
    
    indel_table = data.frame(indel_length = gr$length, pks_status = gr$pks_status, sample = name)
    
  }
  return(indel_table)
}

get_length_pks_indels_all_grls = function(grl){
  if (class(grl)[[1]] == "CompressedGRangesList"){
    names = names(grl)
    indel_counts = lapply(1:length(grl), function(x) get_length_pks_indels_all(grl[[x]], name = names[x]))
    indel_counts = do.call("rbind", indel_counts)
  } else if (class(grl)[[1]] == "GRanges"){
    indel_counts = get_length_pks_indels(grl, name = "sample")
  }
  return(indel_counts)
}

get_length_pks_indels_all_grls_mc = function(grl, n_cores){
  if (class(grl)[[1]] == "CompressedGRangesList"){
    names = names(grl)
    indel_counts = mclapply(1:length(grl), function(x) get_length_pks_indels_all(grl[[x]], name = names[x]), mc.cores = n_cores)
    indel_counts = do.call("rbind", indel_counts)
  } else if (class(grl)[[1]] == "GRanges"){
    indel_counts = get_length_pks_indels(grl, name = "sample")
  }
  return(indel_counts)
}


get_length_pks_indels_grls = function(grl){
  if (class(grl)[[1]] == "CompressedGRangesList"){
    names = names(grl)
    indel_counts = lapply(1:length(grl), function(x) get_length_pks_indels(grl[[x]], name = names[x]))
    indel_counts = do.call("rbind", indel_counts)
  } else if (class(grl)[[1]] == "GRanges"){
    indel_counts = get_length_pks_indels(grl, name = "sample")
  }
  return(indel_counts)
}

get_length_pks_indels_grls_mc = function(grl, n_cores){
  if (class(grl)[[1]] == "CompressedGRangesList"){
    names = names(grl)
    indel_counts = mclapply(1:length(grl), function(x) get_length_pks_indels(grl[[x]], name = names[x]), mc.cores = n_cores)
    indel_counts = do.call("rbind", indel_counts)
  } else if (class(grl)[[1]] == "GRanges"){
    indel_counts = get_length_pks_indels(grl, name = "sample")
  }
  return(indel_counts)
}

remove_mult_alts = function(gr) {
  
  mult_alts = elementNROWS(gr$ALT) > 1
  nr_mult_alts = sum(mult_alts)
  if (nr_mult_alts > 0){
    gr = gr[!mult_alts]
    warning(paste0("Removed ", nr_mult_alts, " variant(s), because they had multiple alternative alleles."), call. = F)
  }
  return(gr)
}



