setwd(file.path(usethis::proj_path(), "5_GSEA/all_viz"))

library(dplyr)

degs <- list.files(
  "..",
  pattern = "degs.*\\.qs$",
  recursive = TRUE,
  full.names = TRUE
)
names(degs) <- tools::file_path_sans_ext(degs) %>% basename()

loaded_degs <- lapply(degs, function(x) {
  qs::qread(x, nthreads = 8L)
})

flatten_loaded_degs <- purrr::imap(loaded_degs, \(tissue, tissue_name) {
  purrr::imap(tissue, \(dataset, dataset_name) {
    tbl <- dplyr::bind_rows(dataset)
    tbl_name <- names(dataset)[!purrr::map_lgl(dataset, base::is.null)]
    dplyr::mutate(tbl, screen_method = tbl_name, dataset_name = dataset_name)
  }) %>%
    dplyr::bind_rows() %>%
    dplyr::mutate(tissue_name = tissue_name)
}) %>%
  dplyr::bind_rows()

data.table::fwrite(flatten_loaded_degs, "degs_all_stats.csv")

library(data.table)

# * Load combined GSEA results (leadingEdge as comma-separated string)
combined <- data.table::fread("gsea_res_all_stats.csv")

# * Load DEGs (wide format: gene columns + meta columns)
flatten_loaded_degs <- data.table::fread("degs_all_stats.csv")

# * Build matching key: combined$tumor -> paste0("degs_", tumor) = tissue_name
combined[, match_tissue := paste0("degs_", tumor)]

# * Identify gene columns in flatten_loaded_degs (exclude meta columns)
meta_cols <- c("screen_method", "dataset_name", "tissue_name")
gene_cols_all <- setdiff(names(flatten_loaded_degs), meta_cols)

# * For each row in combined, match with flatten_loaded_degs and compute stats
result_list <- lapply(seq_len(nrow(combined)), function(i) {
  # Parse leadingEdge genes
  genes_str <- combined$leadingEdge[i]
  genes_vec <- strsplit(genes_str, ",\\s*")[[1]]
  genes_vec <- genes_vec[genes_vec != "" & !is.na(genes_vec)]

  # Find matching row in flatten_loaded_degs
  match_idx <- which(
    flatten_loaded_degs$screen_method == combined$screen_method[i] &
      flatten_loaded_degs$dataset_name == combined$dataset[i] &
      flatten_loaded_degs$tissue_name == combined$match_tissue[i]
  )

  # Default: no match or no genes
  if (length(match_idx) == 0 || length(genes_vec) == 0) {
    return(data.table::data.table(
      screen_method = combined$screen_method[i],
      dataset = combined$dataset[i],
      tumor = combined$tumor[i],
      pathway = combined$pathway[i],
      leadingEdge = genes_str,
      n_genes_leadingEdge = length(genes_vec),
      n_genes_matched = 0L,
      mean_logFC = NA_real_,
      median_logFC = NA_real_,
      sum_logFC = NA_real_,
      logFC_data = list(numeric(0))
    ))
  }

  match_row <- flatten_loaded_degs[match_idx[1]]

  # Find which leadingEdge genes exist in the DEG data
  available_genes <- intersect(genes_vec, gene_cols_all)
  if (length(available_genes) == 0) {
    return(data.table::data.table(
      screen_method = combined$screen_method[i],
      dataset = combined$dataset[i],
      tumor = combined$tumor[i],
      pathway = combined$pathway[i],
      leadingEdge = genes_str,
      n_genes_leadingEdge = length(genes_vec),
      n_genes_matched = 0L,
      mean_logFC = NA_real_,
      median_logFC = NA_real_,
      sum_logFC = NA_real_,
      logFC_data = list(numeric(0))
    ))
  }

  # Extract logFC values for matched genes
  logFC_values <- as.numeric(match_row[, ..available_genes])

  data.table::data.table(
    screen_method = combined$screen_method[i],
    dataset = combined$dataset[i],
    tumor = combined$tumor[i],
    pathway = combined$pathway[i],
    leadingEdge = genes_str,
    n_genes_leadingEdge = length(genes_vec),
    n_genes_matched = length(available_genes),
    mean_logFC = mean(logFC_values, na.rm = TRUE),
    median_logFC = median(logFC_values, na.rm = TRUE),
    sum_logFC = sum(logFC_values, na.rm = TRUE),
    logFC_data = list(logFC_values)
  )
})

# * Bind all results into one data.table
leading_edge_stats <- data.table::rbindlist(result_list)

# * Save results
data.table::fwrite(leading_edge_stats, "leading_edge_stats.csv")
