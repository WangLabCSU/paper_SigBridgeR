fpkm_to_tpm <- function(
  fpkm,
  na_as_zero = TRUE,
  check_input = TRUE,
  verbose = TRUE
) {
  # Input validation
  if (!is.matrix(fpkm)) {
    stop("Input must be a matrix.")
  }

  if (!is.numeric(fpkm)) {
    stop("The input matrix must contain numeric FPKM values.")
  }

  if (is.null(rownames(fpkm))) {
    warning("The input matrix has no row names; gene IDs cannot be preserved.")
  }

  if (check_input) {
    if (any(is.infinite(fpkm), na.rm = TRUE)) {
      stop("The input matrix contains infinite values.")
    }

    if (any(fpkm < 0, na.rm = TRUE)) {
      stop("FPKM values must be non-negative.")
    }
  }

  # Handle missing values
  na_count <- sum(is.na(fpkm))

  if (na_count > 0) {
    if (na_as_zero) {
      if (verbose) {
        message(
          sprintf(
            "Replacing %d missing value(s) with zero.",
            na_count
          )
        )
      }

      # This modifies the matrix copy once and avoids creating another
      # full-size matrix.
      fpkm[is.na(fpkm)] <- 0
    } else {
      stop(
        "The input matrix contains missing values. ",
        "Set na_as_zero = TRUE to replace them with zero."
      )
    }
  }

  # Calculate the FPKM sum for each sample
  sample_sums <- colSums(fpkm)

  zero_sum_samples <- which(sample_sums == 0)

  if (length(zero_sum_samples) > 0) {
    warning(
      sprintf(
        "%d sample(s) have a total FPKM of zero; ",
        length(zero_sum_samples)
      ),
      "their TPM values will be set to zero."
    )
  }

  valid_samples <- which(sample_sums > 0)

  # Normalize column by column to reduce peak memory usage.
  # The input matrix is reused as the output matrix.
  if (length(valid_samples) > 0) {
    for (j in valid_samples) {
      fpkm[, j] <- fpkm[, j] / sample_sums[j] * 1e6
    }
  }

  # Preserve the matrix structure and report completion
  if (verbose) {
    message(
      sprintf(
        "TPM conversion completed: %d genes × %d samples.",
        nrow(fpkm),
        ncol(fpkm)
      )
    )
  }

  return(fpkm)
}
