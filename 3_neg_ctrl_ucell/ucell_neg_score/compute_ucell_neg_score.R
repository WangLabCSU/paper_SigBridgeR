# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd(file.path(usethis::proj_path(), "3_neg_ctrl_ucell/ucell_neg_score"))

library(dplyr)
set.seed(123L)

random_markers <- list.files(
  path = "../..",
  pattern = "random20_markers_100rep\\.csv",
  full.names = TRUE,
  recursive = TRUE,
  ignore.case = TRUE
)
names(random_markers) <- basename(random_markers) %>%
  tools::file_path_sans_ext() %>%
  gsub("_.*", "", .)

loaded_random_markers <- purrr::map(
  random_markers,
  ~ data.table::fread(.x) %>% as.data.frame()
)

mirai::daemons(12)

neg_ucell_score <- purrr::imap(
  loaded_random_markers,
  function(marker_df, name) {
    cli::cli_alert_info("Processing {.val {name}}")

    n_col <- ncol(marker_df)
    mirai::mirai_map(
      seq_len(n_col),
      function(i) {
        neg_gene_list <- list(marker_df[, i][[1L]])
        names(neg_gene_list) <- names(marker_df)[i]

        seurat_i <- switch(
          name,
          "her2" = qs::qread(
            "/home/data/sigbridger/benchmark_data/brca/seurat_her2.qs",
            nthreads = 8L
          ),
          "tnbc" = qs::qread(
            "/home/data/sigbridger/benchmark_data/brca/seurat_tnbc.qs",
            nthreads = 8L
          ),
          "luad" = qs::qread(
            "/home/data/sigbridger/benchmark_data/lung/luad_GSE123902_seurat.qs",
            nthreads = 8L
          ),
          "ov" = qs::qread(
            "/home/data/sigbridger/benchmark_data/ov/hgsoc_GSE165897_seurat.qs",
            nthreads = 8L
          ),
          cli::cli_abort("Invalid name: {name}")
        )

        score <- UCell::ScoreSignatures_UCell(
          matrix = SeuratObject::LayerData(seurat_i, layer = "data"),
          features = neg_gene_list,
          precalc.ranks = NULL,
          maxRank = 1500,
          w_neg = 1,
          name = "_UCell",
          chunk.size = 100,
          BPPARAM = NULL,
          ties.method = "average",
          force.gc = FALSE
        )

        gc(verbose = FALSE)

        score
      },
      marker_df = marker_df,
      name = name
    )[
      mirai::.progress
    ]
  }
)

qs::qsave(neg_ucell_score, "neg_ucell_score.qs", nthreads = 8L)

cli::cli_alert_success("Successfully saved neg_ucell_score to qs file")

gc()
