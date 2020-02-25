fit_to_signatures_indels = function (mut_matrix, signatures)
{
  n_feats = nrow(mut_matrix)
  n_samples = dim(mut_matrix)[2]
  n_signatures = dim(signatures)[2]
  lsq_contribution = matrix(NA, nrow = n_signatures, ncol = n_samples)
  lsq_reconstructed = matrix(NA, nrow = n_feats, ncol = n_samples)
  for (i in 1:ncol(mut_matrix)) {
    y = mut_matrix[, i]
    lsq = lsqnonneg(signatures, y)
    
    #ratio = sum(y) / sum(lsq$x)
    #lsq$x = lsq$x * ratio
    
    lsq_contribution[, i] = lsq$x
    lsq_reconstructed[, i] = signatures %*% as.matrix(lsq$x)
  }
  sample_names = colnames(mut_matrix)
  signature_names = colnames(signatures)
  mut_type_names = rownames(signatures)
  colnames(lsq_contribution) = sample_names
  rownames(lsq_contribution) = signature_names
  colnames(lsq_reconstructed) = sample_names
  rownames(lsq_reconstructed) = mut_type_names
  res = list(lsq_contribution, lsq_reconstructed)
  names(res) = c("contribution", "reconstructed")
  return(res)
}