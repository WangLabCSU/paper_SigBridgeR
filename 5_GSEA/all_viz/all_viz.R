setwd(file.path(usethis::proj_path(), "5_GSEA/all_viz"))

library(dplyr)
source("plot_gsea_table.R")

# future::plan(future.mirai::mirai_multisession(workers = 4L))

stats_file <- list.files(
  "..",
  pattern = "^gsea_res.*\\.qs",
  recursive = TRUE,
  full.names = TRUE
)
names(stats_file) <- basename(stats_file) %>%
  tools::file_path_sans_ext() %>%
  gsub("gsea_res_", "", .)

degs_file <- list.files(
  "..",
  pattern = "^degs.*\\.qs",
  recursive = TRUE,
  full.names = TRUE
)
names(degs_file) <- basename(degs_file) %>%
  tools::file_path_sans_ext() %>%
  gsub("degs_", "", .)

# * Load files
loaded_stats_file <- lapply(stats_file, \(x) {
  qs::qread(x, nthreads = 8L)
})
loaded_degs_file <- lapply(degs_file, \(x) {
  qs::qread(x, nthreads = 8L)
})

# * viz
# * Hallmark gene list
gene_list <- if (!file.exists("hallmarks_gene_list.qs")) {
  geneset_hallmark <- msigdbr::msigdbr(species = "Homo sapiens", category = "H")
  gene_list <- split(
    geneset_hallmark$gene_symbol,
    geneset_hallmark$gs_description
  )
  qs::qsave(gene_list, "hallmarks_gene_list.qs", nthreads = 4L) # < 1 MB
  gene_list
} else {
  cli::cli_alert_info("Found existing hallmarks")
  qs::qread("hallmarks_gene_list.qs", nthreads = 4L)
}

ts_cli <- SigBridgeRUtils::CreateTimeStampCliEnv()


dir.create("plot", showWarnings = FALSE)
# ~ 20 mins
gsea_plots <- purrr::map2(
  loaded_stats_file,
  loaded_degs_file,
  function(tissue_type, deg_tissue) {
    # 2 lists
    l_tissue <- purrr::pmap(
      list(tissue_type, deg_tissue, names(tissue_type)),
      function(data_group, deg_group, group_name) {
        cli::cli_h2("handling {group_name}")
        # a list
        l_group <- purrr::pmap(
          list(data_group, deg_group, names(data_group)),
          function(dt, deg_method, method_name) {
            if (is.null(dt)) {
              ts_cli$cli_alert_info(
                "Skipping {method_name} because no data"
              )
              return(NULL)
            }
            ts_cli$cli_alert_info(
              "Drawing {method_name}"
            )
            # a data.table
            plot_gsea_table(
              hallmarks = gene_list,
              fgseaRes = dt,
              deg = deg_method
            )
          }
        )
        names(l_group) <- names(data_group)
        l_group
      }
    )
    names(l_tissue) <- names(tissue_type)
    l_tissue
  }
)
names(gsea_plots) <- names(loaded_stats_file)

# * save
purrr::walk(
  gsea_plots,
  function(x) {
    purrr::iwalk(
      x,
      function(x, data_group) {
        purrr::iwalk(
          x,
          function(x, method) {
            file_name <- glue::glue("{data_group}_{method}.png")
            ts_cli$cli_alert_info("Saving {file_name}")
            ggplot2::ggsave(
              plot = x,
              filename = paste0("plot/", file_name),
              width = 20,
              height = 10,
              dpi = 400
            )
          }
        )
      }
    )
  }
)

# future::plan(future::sequential())

cli::cli_h1("All done!")

gc()
