plot_heatmap <- function(
  combined_stats_file,
  metrics = "pval",
  col_fun = circlize::colorRamp2(c(0, 100), c("white", "red"))
) {
  # * Metrics Conversion
  if (metrics == "pval") {
    combined_stats_file <- dplyr::mutate(
      combined_stats_file,
      neg_log10_pval = -log10(pval)
    )
    metrics <- "neg_log10_pval"
  } else if (metrics == "padj") {
    combined_stats_file <- dplyr::mutate(
      combined_stats_file,
      neg_log10_padj = -log10(padj)
    )
    metrics <- "neg_log10_padj"
  }

  mat_data <- dplyr::mutate(
    combined_stats_file,
    pheno = gsub("(binary|survival)_.*", "\\1", combined_stats_file$dataset),
    bulk = gsub(
      "[^_]*_[^_]*_(.*)_merged_seurat",
      "\\1",
      combined_stats_file$dataset
    ),
    tumor = toupper(tumor)
  ) %>%
    dplyr::select(
      pathway, # x - axis (列名)
      !!dplyr::sym(metrics), # 热图填充值
      tumor,
      bulk,
      pheno,
      screen_method
    )

  grouped_mat_data <- mat_data %>%
    dplyr::group_by(tumor, bulk, pheno) %>%
    dplyr::group_split(.keep = TRUE) %>%
    purrr::map(
      ~ tidyr::pivot_wider(
        .x,
        names_from = pathway,
        values_from = !!dplyr::sym(metrics)
      )
    )

  names(grouped_mat_data) <- purrr::map_chr(grouped_mat_data, \(tbl) {
    paste(
      unique(tbl$tumor),
      unique(tbl$pheno),
      unique(tbl$bulk),
      sep = "_"
    )
  })

  grouped_mat_metrics <- purrr::map(grouped_mat_data, \(tbl) {
    tbl$tumor <- NULL
    tbl$pheno <- NULL
    tbl$bulk <- NULL
    tbl <- tibble::column_to_rownames(tbl, "screen_method")
    if ("NA" %in% colnames(tbl)) {
      tbl <- dplyr::select(tbl, -"NA")
    }
    tbl
  })

  # * 构建分面热图列表（每个 tumor_pheno_bulk 一个分面）
  n <- length(grouped_mat_metrics)
  heatmap_list <- vector("list", n)

  for (i in seq_len(n)) {
    mat <- grouped_mat_metrics[[i]]
    grp_name <- names(grouped_mat_metrics)[i]

    heatmap_list[[i]] <- ComplexHeatmap::Heatmap(
      na_col = "#FFFFFF",
      matrix = as.matrix(mat),
      name = metrics,
      col = col_fun,

      cluster_rows = FALSE,
      cluster_columns = FALSE,

      row_title = grp_name,
      row_title_gp = grid::gpar(fontsize = 10, fill = "#FFFFFF"),
      row_names_rot = 0,
      row_title_rot = 0,
      row_names_gp = grid::gpar(fontsize = 8),
      row_gap = grid::unit(1, "mm"),

      show_row_names = TRUE,
      show_column_names = (i == n),

      column_names_rot = 60,
      column_names_gp = grid::gpar(fontsize = 8),
      column_title_side = "bottom",

      heatmap_width = grid::unit(0.7, "npc"),
      heatmap_height = grid::unit(0.8, "npc")
    )
  }

  # * 纵向拼接并绘制
  ht_list <- Reduce(ComplexHeatmap::`%v%`, heatmap_list)

  CairoPNG(
    "gsea_res_all_stats.png",
    width = 3500,
    height = 10000,
    dpi = 400
  )
  ComplexHeatmap::draw(
    ht_list,
    merge_legends = TRUE,
    ht_gap = grid::unit(2, "mm"),
    column_title = "Hallmarks Gene Sets",
    column_title_side = "bottom"
  )
  dev.off()
}
