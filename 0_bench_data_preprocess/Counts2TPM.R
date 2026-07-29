Counts2TPM <- function(counts, gene_length_bp) {
  # counts: gene x sample matrix
  # gene_length_bp: named numeric vector, bp
  common <- intersect(rownames(counts), names(gene_length_bp))
  counts <- counts[common, , drop = FALSE]
  len <- gene_length_bp[common]

  keep <- !is.na(len) & len > 0 & rowSums(counts) > 0
  counts <- counts[keep, , drop = FALSE]
  len <- len[keep]

  rpk <- sweep(counts, 1, len / 1000, "/")
  tpm <- sweep(rpk, 2, colSums(rpk), "/") * 1e6

  return(tpm)
}
