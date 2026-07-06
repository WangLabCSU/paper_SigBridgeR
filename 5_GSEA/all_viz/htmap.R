setwd(file.path(usethis::proj_path(), "5_GSEA/all_viz"))

library(dplyr)
source("wrap_every_n_chars.R")
source("plot_heatmap.R")
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

# * Load files
loaded_stats_file <- lapply(stats_file, \(x) {
  qs::qread(x, nthreads = 8L)
})

combined_stats_file <- purrr::imap(
  loaded_stats_file,
  \(tumor, tumor_name) {
    dt_tumor <- purrr::imap(tumor, \(dataset, dataset_name) {
      dt_dataset <- purrr::imap(dataset, \(method, method_name) {
        method$screen_method <- method_name
        method
      }) %>%
        dplyr::bind_rows()

      dt_dataset$dataset <- dataset_name
      dt_dataset
    }) %>%
      dplyr::bind_rows()

    dt_tumor$tumor <- tumor_name
    dt_tumor
  }
) %>%
  dplyr::bind_rows()

combined_stats_file$leadingEdge <- lapply(
  combined_stats_file$leadingEdge,
  \(x) unlist(x) %>% toString()
)
data.table::fwrite(combined_stats_file, "gsea_res_all_stats.csv")

# ------------------------------------------------------------------------------

combined_stats_file <- data.table::fread("gsea_res_all_stats.csv")


# * binary and survival - padj
plot_heatmap2(
  combined_stats_file = combined_stats_file,
  metrics = "padj",
  col_fun = circlize::colorRamp2(c(0, 80), c("white", "red")),
  filename = "gsea_res_all_stats_padj.png"
)
# * only survival - padj
plot_heatmap3(
  combined_stats_file = combined_stats_file[
    stringr::str_detect(combined_stats_file$dataset, "survival"),
  ],
  metrics = "padj",
  col_fun = circlize::colorRamp2(c(0, 80), c("white", "red")),
  filename = "gsea_res_all_stats_padj_survival.png",
  chr_width = 45L
)
# * only binary - padj
plot_heatmap3(
  combined_stats_file = combined_stats_file[
    stringr::str_detect(combined_stats_file$dataset, "binary"),
  ],
  metrics = "padj",
  col_fun = circlize::colorRamp2(c(0, 80), c("white", "red")),
  filename = "gsea_res_all_stats_padj_binary.png",
  chr_width = 45L
)

# * binary and survival - NES
plot_heatmap2(
  combined_stats_file = combined_stats_file,
  metrics = "NES",
  col_fun = circlize::colorRamp2(c(-5, 0, 5), c("blue", "white", "red")),
  filename = "gsea_res_all_stats_nes.png"
)
# * only survival - NES
plot_heatmap3(
  combined_stats_file = combined_stats_file[
    stringr::str_detect(combined_stats_file$dataset, "survival"),
  ],
  metrics = "NES",
  col_fun = circlize::colorRamp2(c(-5, 0, 5), c("blue", "white", "red")),
  filename = "gsea_res_all_stats_nes_survival.png",
  chr_width = 45L
)
# * only binary - NES
plot_heatmap3(
  combined_stats_file = combined_stats_file[
    stringr::str_detect(combined_stats_file$dataset, "binary"),
  ],
  metrics = "NES",
  col_fun = circlize::colorRamp2(c(-5, 0, 5), c("blue", "white", "red")),
  filename = "gsea_res_all_stats_nes_binary.png",
  chr_width = 45L
)

# * only survival - padj - significance
plot_heatmap4(
  combined_stats_file = combined_stats_file[
    stringr::str_detect(combined_stats_file$dataset, "survival"),
  ],
  metrics = "padj",
  filename = "gsea_res_all_stats_padj_signif_survival.png",
  chr_width = 45L
)
# * only binary - padj - significance
plot_heatmap4(
  combined_stats_file = combined_stats_file[
    stringr::str_detect(combined_stats_file$dataset, "binary"),
  ],
  metrics = "padj",
  filename = "gsea_res_all_stats_padj_signif_binary.png",
  chr_width = 45L
)
