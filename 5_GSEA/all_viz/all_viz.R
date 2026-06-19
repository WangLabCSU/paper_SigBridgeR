setwd(file.path(usethis::proj_path(), "5_GSEA/all_viz"))

library(dplyr)
source("plot_gsea_table.R")

future::plan(future.mirai::mirai_multisession(workers = 4L))

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
geneset_hallmark <- msigdbr::msigdbr(species = "Homo sapiens", category = "H")
gene_list <- split(
  geneset_hallmark$gene_symbol,
  geneset_hallmark$gs_description
)

# ? furrr is not compatible with data.table syntax
dir.create("plot", showWarnings = FALSE)
# ~ 20 mins
gsea_plots <- purrr::map2(
  loaded_stats_file,
  loaded_degs_file,
  function(tissue_type, deg_tissue) {
    # 2 lists
    purrr::map2(
      tissue_type,
      deg_tissue,
      function(data_group, deg_group) {
        # a list
        purrr::map2(
          data_group,
          deg_group,
          function(dt, deg_method) {
            if (is.null(dt)) {
              return(NULL)
            }
            # a data.table
            plot_gsea_table(
              hallmarks = gene_list,
              fgseaRes = dt,
              deg = deg_method
            )
          }
        )
      }
    )
  },
  .progress = "Drawing"
)

# * save
purrr::iwalk(
  gsea_plots,
  function(x, tissue_type) {
    furrr::future_iwalk(
      x,
      function(x, data_group) {
        purrr::iwalk(
          x,
          function(x, method) {
            file_name <- glue::glue("{tissue_type}_{data_group}_{method}.png")
            ggplot2::ggsave(
              plot = x,
              filename = paste0("plot/", file_name),
              width = 10,
              height = 10,
              dpi = 400
            )
          }
        )
      }
    )
  },
  .progress = "Saving plots"
)

future::plan(future::sequential())
