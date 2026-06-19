# ! include both survival and binary

setwd(file.path(usethis::proj_path(), "5_GSEA/test_ov"))

library(dplyr)

# * Hyperparameter
sur_seurat_path <- "/home/data/sigbridger/benchmark_data/ov/ov"
bi_seurat_path <- "/home/data/sigbridger/benchmark_binary/ov"
tissue <- "ov"
mirai::daemons(3L)
# future::plan(future.mirai::mirai_multisession(workers = 3L))

# ------------------------------------------------------------------------------------------------------------------------------------
# * Load Seurat data
phenotype_class <- c("survival", "binary")

labeled_seurat <- c(
  list.files(
    sur_seurat_path,
    pattern = glue::glue("{phenotype_class[[1L]]}.*\\.qs$"),
    full.names = TRUE
  ),
  list.files(
    bi_seurat_path,
    pattern = glue::glue("{phenotype_class[[2L]]}.*\\.qs$"),
    full.names = TRUE
  )
)
names(labeled_seurat) <- basename(labeled_seurat) %>%
  tools::file_path_sans_ext()

labeled_seurat_loaded <-
  #   furrr::future_map(
  mirai::mirai_map(
    labeled_seurat,
    function(x) {
      qs::qread(x, nthreads = 8L)
    },
    .progress = TRUE
  )
labeled_seurat_loaded <- labeled_seurat_loaded[mirai::.progress]

# * Hallmark gene list
geneset_hallmark <- msigdbr::msigdbr(species = "Homo sapiens", category = "H")
gene_list <- split(
  geneset_hallmark$gene_symbol,
  geneset_hallmark$gs_description
)

# * Find DEG
methods <- c(
  "scissor",
  "scAB",
  "scPAS",
  "scPP",
  "DEGAS",
  "PIPET",
  "SCIPAC",
  "LP_SGL"
)

degs <- lapply(seq_along(labeled_seurat_loaded), function(i) {
  sc_data <- labeled_seurat_loaded[[i]]
  existing_cols <- colnames(sc_data[[]])
  filtered_cols <- existing_cols[existing_cols %in% methods] # available methods

  cli::cli_alert_info(
    "Available methods for {.val {names(labeled_seurat_loaded)[i]}}: {
    filtered_cols}"
  )

  #   res <- furrr::future_map(
  res <- mirai::mirai_map(
    filtered_cols,
    function(method) {
      # method is a character
      sc_data[[method]] <- ifelse(
        sc_data[[method]] == "Positive",
        "Positive",
        "non_Positive"
      )
      cell_counts <- cheapr::table_(sc_data[[method]])

      if (!"Positive" %in% names(cell_counts)) {
        cli::cli_warn("{method}: Positive not in cell_counts")
        return(NULL)
      }

      deg_results <- Seurat::FindMarkers(
        sc_data,
        ident.1 = "Positive",
        ident.2 = "non_Positive",
        group.by = method
      )

      gene_list <- deg_results$avg_log2FC

      # 3. 给向量命名 (基因名)
      names(gene_list) <- rownames(deg_results)

      # 4. 按 log2FC 从高到低 降序排序（这是必须执行的一步）
      gene_list <- sort(gene_list, decreasing = TRUE)
      gene_list
    },
    sc_data = sc_data
    # .progress = TRUE
  )
  res <- res[mirai::.progress]
  names(res) <- filtered_cols
  res
})
names(degs) <- names(labeled_seurat_loaded)

qs::qsave(degs, glue::glue("degs_{tissue}.qs"), nthreads = 4L) # < 5MB


# * GSEA
fgsea_res <- lapply(degs, function(one_seurat_deg) {
  mirai::mirai_map(
    one_seurat_deg,
    function(deg_vec) {
      if (is.null(deg_vec)) {
        return(NULL)
      }
      fgsea::fgsea(
        pathways = gene_list,
        stats = deg_vec,
        eps = 0.0,
      )
      # data.table
    },
    gene_list = gene_list
  )[mirai::.progress]
})

mirai::daemons(0L)
# future::plan(future::sequential())

# * Save results
qs::qsave(fgsea_res, glue::glue("gsea_res_{tissue}.qs"), nthreads = 4L) # < 1MB

cli::cli_h1("All done!")
gc()
