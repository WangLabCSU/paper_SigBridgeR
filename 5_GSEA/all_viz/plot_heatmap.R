# # * Group by tumor/Bulk label, containing NA rows
# plot_heatmap2 <- function(
#   combined_stats_file,
#   metrics = "pval",
#   col_fun = circlize::colorRamp2(c(0, 100), c("white", "red")),
#   filename = "gsea_res_all_stats.png"
# ) {
#   # * Metrics Conversion
#   if (metrics == "pval") {
#     combined_stats_file <- dplyr::mutate(
#       combined_stats_file,
#       neg_log10_pval = -log10(pval)
#     )
#     metrics <- "neg_log10_pval"
#   } else if (metrics == "padj") {
#     combined_stats_file <- dplyr::mutate(
#       combined_stats_file,
#       neg_log10_padj = -log10(padj)
#     )
#     metrics <- "neg_log10_padj"
#   }

#   mat_data <- dplyr::mutate(
#     combined_stats_file,
#     pheno = gsub("(binary|survival)_.*", "\\1", combined_stats_file$dataset),
#     bulk = gsub(
#       "[^_]*_[^_]*_(.*)_merged_seurat",
#       "\\1",
#       combined_stats_file$dataset
#     ),
#     tumor = toupper(tumor)
#   ) %>%
#     dplyr::select(
#       pathway, # x - axis (列名)
#       !!dplyr::sym(metrics), # 热图填充值
#       tumor,
#       bulk,
#       pheno,
#       screen_method
#     ) %>%
#     dplyr::mutate(
#       y_axis = paste(tumor, pheno, bulk, screen_method, sep = "_")
#     )

#   wide_mat_data <- mat_data %>%
#     dplyr::select(-tumor, -pheno, -bulk, -screen_method) %>%
#     tidyr::pivot_wider(
#       names_from = y_axis,
#       values_from = !!dplyr::sym(metrics)
#     ) %>%
#     t()

#   colnames(wide_mat_data) <- wide_mat_data[1, ]
#   wide_mat_data <- wide_mat_data[-1, -1]

#   wide_mat_data <- as.data.frame(wide_mat_data) %>%
#     dplyr::mutate(dplyr::across(dplyr::everything(), base::as.numeric)) %>%
#     as.matrix()

#   #
#   # 唯一区分标签
#   dataset_label <- rownames(wide_mat_data)

#   palette_tumor <- c(
#     "BRCA_HER2" = "#773232",
#     "TNBC" = "#12298f",
#     "LUNG" = "#23964f",
#     "OV" = "#5a229b"
#   )
#   palette_pheno <- c(
#     "survival" = "#ddcbaf",
#     "binary" = "#b3cad4"
#   )
#   palette_bulk <- setNames(
#     c(
#       "#B8D9A0",
#       "#E8E0DB",
#       "#A8D9D8",
#       "#F5C0A8",
#       "#D5C9E8",
#       "#9BBBD9",
#       "#B8A8E5",
#       "#D5E5C0",
#       "#E8D8E8",
#       "#88D5B0"
#     ),
#     unique(mat_data$bulk)
#   )

#   row_anno_tissue <- ComplexHeatmap::rowAnnotation(
#     "Phenotype" = stringr::str_extract(dataset_label, "binary|survival"),
#     "Bulk Data" = stringr::str_extract(dataset_label, "TCGA_[^_]*|GSE[0-9]*"),
#     annotation_name_gp = grid::gpar(fontface = "bold"),
#     annotation_name_rot = 60,
#     col = list(
#       "Phenotype" = palette_pheno,
#       "Bulk Data" = palette_bulk
#     ),
#     gap = grid::unit(0.6, "mm")
#   )

#   htmap <- ComplexHeatmap::Heatmap(
#     matrix = wide_mat_data,
#     name = if (metrics == "neg_log10_pval") {
#       "-log10(p-value)"
#     } else if (metrics == "neg_log10_padj") {
#       "-log10(p-adjusted)"
#     } else {
#       metrics
#     },
#     col = col_fun,
#     na_col = "#DDDDDD",

#     cluster_rows = FALSE,
#     cluster_columns = TRUE,

#     top_annotation = ComplexHeatmap::columnAnnotation(
#       Distribution = ComplexHeatmap::anno_boxplot(
#         wide_mat_data,
#         gp = grid::gpar(
#           fill = c(
#             "#B8D9A0",
#             "#E8E0DB",
#             "#A8D9D8",
#             "#F5C0A8",
#             "#D5C9E8",
#             "#9BBBD9",
#             "#B8A8E5",
#             "#D5E5C0",
#             "#E8D8E8",
#             "#88D5B0",
#             "#E5B8E5",
#             "#88E0D5",
#             "#A8D5A8",
#             "#E5A8D0",
#             "#F5C8E5",
#             "#F0E890",
#             "#F5E8B8",
#             "#B8BCA0",
#             "#C8E5D0",
#             "#B888D5",
#             "#C5B8E8",
#             "#E5C090",
#             "#D5A0A8",
#             "#88A8D5",
#             "#F5D890",
#             "#D8E5B0",
#             "#B8D5D0",
#             "#E0E8D8",
#             "#E0A8E5",
#             "#A888D5",
#             "#C5B0B8",
#             "#E5A8E8",
#             "#E8A8D5",
#             "#D0C5B8",
#             "#F5A098",
#             "#B0E5D0",
#             "#C0D8E5",
#             "#90C8D5",
#             "#E8A0B0",
#             "#F5E0B8"
#           )
#         ),
#       ),
#       height = grid::unit(5, "cm")
#     ),
#     left_annotation = c(row_anno_tissue),

#     # column_title = "Hallmarks Gene Set",
#     # column_title_side = "bottom",
#     column_names_gp = grid::gpar(fontsize = 10),
#     row_names_gp = grid::gpar(fontsize = 10),

#     row_labels = gsub(".*_", "", dataset_label) %>%
#       stringr::str_replace("SGL", "LP_SGL"),
#     row_split = stringr::str_replace(rownames(wide_mat_data), "_SGL", "SGL") %>%
#       gsub("_.*", "", .) %>%
#       stringr::str_replace("BRCA", "BRCA HER2"),
#     column_labels = wrap_every_n_chars(colnames(wide_mat_data), width = 60),
#     column_names_rot = 60,

#     heatmap_height = grid::unit(0.1, "npc"),
#     height = grid::unit(0.8, "npc"),
#     width = grid::unit(0.9, "npc"),
#   )

#   Cairo::CairoPNG(
#     filename = filename,
#     dpi = 400,
#     width = 6000,
#     height = 10000
#   )
#   ComplexHeatmap::draw(htmap)
#   dev.off()
# }

# * NA rows not included, cluster rows & cols
plot_heatmap3 <- function(
  combined_stats_file,
  metrics = "pval",
  col_fun = circlize::colorRamp2(c(0, 100), c("white", "red")),
  filename = "gsea_res_all_stats.png",
  chr_width = 30L
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
    ) %>%
    dplyr::mutate(
      y_axis = paste(tumor, pheno, bulk, screen_method, sep = "_")
    )

  wide_mat_data <- mat_data %>%
    dplyr::select(-tumor, -pheno, -bulk, -screen_method) %>%
    tidyr::pivot_wider(
      names_from = y_axis,
      values_from = !!dplyr::sym(metrics)
    ) %>%
    t()

  colnames(wide_mat_data) <- wide_mat_data[1, ]
  wide_mat_data <- wide_mat_data[-1, -1]

  wide_mat_data <- as.data.frame(wide_mat_data) %>%
    dplyr::mutate(dplyr::across(dplyr::everything(), base::as.numeric)) %>%
    dplyr::filter(!dplyr::if_all(dplyr::everything(), base::is.na)) %>%
    as.matrix()
  # 唯一区分标签
  dataset_label <- rownames(wide_mat_data)

  #   palette_tumor <- c(
  #     "BRCA_HER2" = "#773232",
  #     "TNBC" = "#12298f",
  #     "LUNG" = "#23964f",
  #     "OV" = "#5a229b"
  #   )

  bulk_names <- unique(mat_data$bulk)
  n_bulk <- length(bulk_names)
  palette_bulk <- setNames(
    c(
      "#B8D9A0",
      "#E8E0DB",
      "#A8D9D8",
      "#F5C0A8",
      "#D5C9E8",
      "#88D5B0",
      "#B8A8E5",
      "#D5E5C0",
      "#E8D8E8",
      "#9BBBD9"
    )[seq_len(n_bulk)],
    unique(mat_data$bulk)
  )

  row_anno <- ComplexHeatmap::rowAnnotation(
    # "Tumor" = stringr::str_extract(dataset_label, "BRCA_HER2|TNBC|LUNG|OV"),
    "Bulk Data" = stringr::str_extract(dataset_label, "TCGA_[^_]*|GSE[0-9]*"),
    annotation_name_gp = grid::gpar(fontface = "bold"),
    annotation_name_rot = 60,
    col = list(
      #   "Tumor" = palette_tumor,
      "Bulk Data" = palette_bulk
    ),
    gap = grid::unit(0.6, "mm")
  )

  htmap <- ComplexHeatmap::Heatmap(
    matrix = wide_mat_data,
    name = if (metrics == "neg_log10_pval") {
      "-log10(p-value)"
    } else if (metrics == "neg_log10_padj") {
      "-log10(p-adjusted)"
    } else {
      metrics
    },
    col = col_fun,
    na_col = "#DDDDDD",

    cluster_rows = TRUE,
    cluster_columns = TRUE,
    column_dend_height = grid::unit(2.25, "cm"),
    row_dend_width = grid::unit(2.25, "cm"),

    top_annotation = ComplexHeatmap::columnAnnotation(
      Distribution = ComplexHeatmap::anno_boxplot(
        wide_mat_data,
        gp = grid::gpar(
          fill = c(
            "#B8D9A0",
            "#E8E0DB",
            "#A8D9D8",
            "#F5C0A8",
            "#D5C9E8",
            "#9BBBD9",
            "#B8A8E5",
            "#D5E5C0",
            "#E8D8E8",
            "#88D5B0",
            "#E5B8E5",
            "#88E0D5",
            "#A8D5A8",
            "#E5A8D0",
            "#F5C8E5",
            "#F0E890",
            "#F5E8B8",
            "#B8BCA0",
            "#C8E5D0",
            "#B888D5",
            "#C5B8E8",
            "#E5C090",
            "#D5A0A8",
            "#88A8D5",
            "#F5D890",
            "#D8E5B0",
            "#B8D5D0",
            "#E0E8D8",
            "#E0A8E5",
            "#A888D5",
            "#C5B0B8",
            "#E5A8E8",
            "#E8A8D5",
            "#D0C5B8",
            "#F5A098",
            "#B0E5D0",
            "#C0D8E5",
            "#90C8D5",
            "#E8A0B0",
            "#F5E0B8"
          )
        ),
      ),
      height = grid::unit(5, "cm")
    ),
    left_annotation = row_anno,

    # column_title = "Hallmarks Gene Set",
    # column_title_side = "bottom",
    column_names_gp = grid::gpar(fontsize = 14),
    row_names_gp = grid::gpar(fontsize = 14),

    row_labels = gsub(".*_", "", dataset_label) %>%
      stringr::str_replace("SGL", "LP_SGL"),
    row_split = stringr::str_replace(rownames(wide_mat_data), "_SGL", "SGL") %>%
      gsub("_.*", "", .) %>%
      stringr::str_replace("BRCA", "BRCA HER2"),
    show_parent_dend_line = FALSE,
    column_labels = wrap_every_n_chars(
      colnames(wide_mat_data),
      width = chr_width
    ),
    column_names_rot = 60,
    column_order = sort(colnames(wide_mat_data)),

    heatmap_height = grid::unit(0.1, "npc"),
    height = grid::unit(0.7, "npc"),
    width = grid::unit(0.85, "npc"),
  )

  Cairo::CairoPNG(
    filename = filename,
    dpi = 400,
    width = 10000,
    height = 10000
  )
  htmap <- ComplexHeatmap::draw(htmap)
  dev.off()

  htmap
}

# * padj only, significance
plot_heatmap4 <- function(
  combined_stats_file,
  metrics = "pval",
  # Color mapping matching discretized significance levels
  # -log10 thresholds: 0.05→1.301, 0.01→2, 0.001→3
  col_fun = circlize::colorRamp2(
    c(0, 1.30103, 2, 3),
    c("#E8E8E8", "#FDB462", "#FB8072", "#B2182B")
  ),
  filename = "gsea_res_all_stats.png",
  chr_width = 30L
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
    ) %>%
    dplyr::mutate(
      y_axis = paste(tumor, pheno, bulk, screen_method, sep = "_")
    )

  wide_mat_data <- mat_data %>%
    dplyr::select(-tumor, -pheno, -bulk, -screen_method) %>%
    tidyr::pivot_wider(
      names_from = y_axis,
      values_from = !!dplyr::sym(metrics)
    ) %>%
    t()

  colnames(wide_mat_data) <- wide_mat_data[1, ]
  wide_mat_data <- wide_mat_data[-1, -1]

  wide_mat_data <- as.data.frame(wide_mat_data) %>%
    dplyr::mutate(dplyr::across(dplyr::everything(), base::as.numeric)) %>%
    dplyr::filter(!dplyr::if_all(dplyr::everything(), base::is.na)) %>%
    as.matrix()

  # 唯一区分标签
  dataset_label <- rownames(wide_mat_data)

  #   palette_tumor <- c(
  #     "BRCA_HER2" = "#773232",
  #     "TNBC" = "#12298f",
  #     "LUNG" = "#23964f",
  #     "OV" = "#5a229b"
  #   )

  bulk_names <- unique(mat_data$bulk)
  n_bulk <- length(bulk_names)
  palette_bulk <- setNames(
    c(
      "#B8D9A0",
      "#E8E0DB",
      "#A8D9D8",
      "#F5C0A8",
      "#D5C9E8",
      "#88D5B0",
      "#B8A8E5",
      "#D5E5C0",
      "#E8D8E8",
      "#9BBBD9"
    )[seq_len(n_bulk)],
    unique(mat_data$bulk)
  )

  row_anno <- ComplexHeatmap::rowAnnotation(
    # "Tumor" = stringr::str_extract(dataset_label, "BRCA_HER2|TNBC|LUNG|OV"),
    "Bulk Data" = stringr::str_extract(dataset_label, "TCGA_[^_]*|GSE[0-9]*"),
    annotation_name_gp = grid::gpar(fontface = "bold"),
    annotation_name_rot = 60,
    col = list(
      #   "Tumor" = palette_tumor,
      "Bulk Data" = palette_bulk
    ),
    gap = grid::unit(0.6, "mm")
  )

  # Cap values > 3 (corresponding to padj < 0.001, -log10(0.001) = 3)
  wide_mat_data[wide_mat_data > 3] <- 3

  # significance legend
  lgd <- ComplexHeatmap::Legend(
    labels = c(
      ">= 0.05",
      "< 0.05",
      "< 0.01",
      "< 0.001"
    ),
    title = "padj",
    legend_gp = grid::gpar(
      fill = c("#DDDDDD", "#FDB462", "#FB8072", "#962835")
    )
  )

  htmap <- ComplexHeatmap::Heatmap(
    matrix = wide_mat_data,
    name = "-log10(padj)",
    col = col_fun,
    na_col = "#BBBBBB",

    cluster_rows = TRUE,
    cluster_columns = FALSE,
    column_dend_height = grid::unit(3, "cm"),
    row_dend_width = grid::unit(3, "cm"),

    top_annotation = ComplexHeatmap::columnAnnotation(
      Distribution = ComplexHeatmap::anno_boxplot(
        wide_mat_data,
        gp = grid::gpar(
          fill = c(
            "#B8D9A0",
            "#E8E0DB",
            "#A8D9D8",
            "#F5C0A8",
            "#D5C9E8",
            "#9BBBD9",
            "#B8A8E5",
            "#D5E5C0",
            "#E8D8E8",
            "#88D5B0",
            "#E5B8E5",
            "#88E0D5",
            "#A8D5A8",
            "#E5A8D0",
            "#F5C8E5",
            "#F0E890",
            "#F5E8B8",
            "#B8BCA0",
            "#C8E5D0",
            "#B888D5",
            "#C5B8E8",
            "#E5C090",
            "#D5A0A8",
            "#88A8D5",
            "#F5D890",
            "#D8E5B0",
            "#B8D5D0",
            "#E0E8D8",
            "#E0A8E5",
            "#A888D5",
            "#C5B0B8",
            "#E5A8E8",
            "#E8A8D5",
            "#D0C5B8",
            "#F5A098",
            "#B0E5D0",
            "#C0D8E5",
            "#90C8D5",
            "#E8A0B0",
            "#F5E0B8"
          )
        ),
      ),
      height = grid::unit(5, "cm")
    ),
    left_annotation = c(row_anno),

    # column_title = "Hallmarks Gene Set",
    # column_title_side = "bottom",
    column_names_gp = grid::gpar(fontsize = 14),
    row_names_gp = grid::gpar(fontsize = 14),

    row_labels = gsub(".*_", "", dataset_label) %>%
      stringr::str_replace("SGL", "LP_SGL"),
    row_split = stringr::str_replace(rownames(wide_mat_data), "_SGL", "SGL") %>%
      gsub("_.*", "", .) %>%
      stringr::str_replace("BRCA", "BRCA HER2"),
    show_parent_dend_line = FALSE,

    column_labels = wrap_every_n_chars(
      colnames(wide_mat_data),
      width = chr_width
    ),
    column_names_rot = 60,

    heatmap_height = grid::unit(0.1, "npc"),
    height = grid::unit(0.7, "npc"),
    width = grid::unit(0.85, "npc"),

    show_heatmap_legend = FALSE
  )

  Cairo::CairoPNG(
    filename = filename,
    dpi = 400,
    width = 10000,
    height = 10000
  )
  ComplexHeatmap::draw(htmap, annotation_legend_list = list(lgd))
  dev.off()
}

# * direction - padj
plot_heatmap_red_blue_padj <- function(
  combined_stats_file,
  metrics = "neg_log10_padj",
  # Color mapping matching discretized significance levels
  # -log10 thresholds: 0.05→1.301, 0.01→2, 0.001→3
  col_fun = circlize::colorRamp2(c(-80, 0, 80), c("blue", "white", "red")),
  filename = "gsea_res_all_stats.png",
  chr_width = 30L,
  col_order = NULL,
  ...
) {
  combined_stats_file <- dplyr::mutate(
    combined_stats_file,
    neg_log10_padj = -log10(padj),
    neg_log10_padj = dplyr::case_when(
      NES < 0 ~ -neg_log10_padj,
      NES > 0 ~ neg_log10_padj,
      TRUE ~ NA
    )
  )

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
    ) %>%
    dplyr::mutate(
      y_axis = paste(tumor, pheno, bulk, screen_method, sep = "_")
    )

  wide_mat_data <- mat_data %>%
    dplyr::select(-tumor, -pheno, -bulk, -screen_method) %>%
    tidyr::pivot_wider(
      names_from = y_axis,
      values_from = !!dplyr::sym(metrics)
    ) %>%
    t()

  colnames(wide_mat_data) <- wide_mat_data[1, ]
  wide_mat_data <- wide_mat_data[-1, -1]

  wide_mat_data <- as.data.frame(wide_mat_data) %>%
    dplyr::mutate(dplyr::across(dplyr::everything(), base::as.numeric)) %>%
    dplyr::filter(!dplyr::if_all(dplyr::everything(), base::is.na)) %>%
    as.matrix()

  # 唯一区分标签
  dataset_label <- rownames(wide_mat_data)

  #   palette_tumor <- c(
  #     "BRCA_HER2" = "#773232",
  #     "TNBC" = "#12298f",
  #     "LUNG" = "#23964f",
  #     "OV" = "#5a229b"
  #   )

  bulk_names <- unique(mat_data$bulk)
  n_bulk <- length(bulk_names)
  palette_bulk <- setNames(
    c(
      "#B8D9A0",
      "#E8E0DB",
      "#A8D9D8",
      "#F5C0A8",
      "#D5C9E8",
      "#88D5B0",
      "#B8A8E5",
      "#D5E5C0",
      "#E8D8E8",
      "#9BBBD9"
    )[seq_len(n_bulk)],
    unique(mat_data$bulk)
  )

  row_anno <- ComplexHeatmap::rowAnnotation(
    # "Tumor" = stringr::str_extract(dataset_label, "BRCA_HER2|TNBC|LUNG|OV"),
    "Bulk Data" = stringr::str_extract(dataset_label, "TCGA_[^_]*|GSE[0-9]*"),
    annotation_name_gp = grid::gpar(fontface = "bold"),
    annotation_name_rot = 60,
    col = list(
      #   "Tumor" = palette_tumor,
      "Bulk Data" = palette_bulk
    ),
    gap = grid::unit(0.6, "mm")
  )

  lgd <- ComplexHeatmap::Legend(
    labels = c(
      "80 & dn",
      "0",
      "80 & up"
    ),
    title = ComplexHeatmap::gt_render(
      "<span> log10(p-adjusted) <br> with direction</span>"
    ),
    at = c(-80, 0, 80),
    col_fun = col_fun
  )

  htmap <- ComplexHeatmap::Heatmap(
    matrix = wide_mat_data,
    name = "-log10(padj)",
    col = col_fun,
    na_col = "#BBBBBB",

    cluster_rows = TRUE,
    cluster_columns = FALSE,
    column_dend_height = grid::unit(3, "cm"),
    row_dend_width = grid::unit(3, "cm"),

    # top_annotation = ComplexHeatmap::columnAnnotation(
    #   Distribution = ComplexHeatmap::anno_boxplot(
    #     wide_mat_data,
    #     gp = grid::gpar(
    #       fill = c(
    #         "#B8D9A0",
    #         "#E8E0DB",
    #         "#A8D9D8",
    #         "#F5C0A8",
    #         "#D5C9E8",
    #         "#9BBBD9",
    #         "#B8A8E5",
    #         "#D5E5C0",
    #         "#E8D8E8",
    #         "#88D5B0",
    #         "#E5B8E5",
    #         "#88E0D5",
    #         "#A8D5A8",
    #         "#E5A8D0",
    #         "#F5C8E5",
    #         "#F0E890",
    #         "#F5E8B8",
    #         "#B8BCA0",
    #         "#C8E5D0",
    #         "#B888D5",
    #         "#C5B8E8",
    #         "#E5C090",
    #         "#D5A0A8",
    #         "#88A8D5",
    #         "#F5D890",
    #         "#D8E5B0",
    #         "#B8D5D0",
    #         "#E0E8D8",
    #         "#E0A8E5",
    #         "#A888D5",
    #         "#C5B0B8",
    #         "#E5A8E8",
    #         "#E8A8D5",
    #         "#D0C5B8",
    #         "#F5A098",
    #         "#B0E5D0",
    #         "#C0D8E5",
    #         "#90C8D5",
    #         "#E8A0B0",
    #         "#F5E0B8"
    #       )
    #     ),
    #   ),
    #   height = grid::unit(5, "cm")
    # ),
    left_annotation = c(row_anno),

    # column_title = "Hallmarks Gene Set",
    # column_title_side = "bottom",
    column_names_gp = grid::gpar(fontsize = 14),
    row_names_gp = grid::gpar(fontsize = 14),

    row_labels = gsub(".*_", "", dataset_label) %>%
      stringr::str_replace("SGL", "LP_SGL"),
    row_split = stringr::str_replace(rownames(wide_mat_data), "_SGL", "SGL") %>%
      gsub("_.*", "", .) %>%
      stringr::str_replace("BRCA", "BRCA HER2"),
    show_parent_dend_line = FALSE,

    column_labels = wrap_every_n_chars(
      colnames(wide_mat_data),
      width = chr_width
    ),
    column_names_rot = 60,
    column_order = col_order,

    heatmap_height = grid::unit(0.1, "npc"),
    height = grid::unit(0.7, "npc"),
    width = grid::unit(0.85, "npc"),

    show_heatmap_legend = FALSE
  )

  Cairo::CairoPNG(
    filename = filename,
    dpi = 400,
    width = 10000,
    height = 10000
  )
  ComplexHeatmap::draw(htmap, annotation_legend_list = list(lgd))
  dev.off()
}

# * direction - padj significance
plot_heatmap_red_blue_padj_signif <- function(
  combined_stats_file,
  metrics = "neg_log10_padj",
  # Color mapping matching discretized significance levels
  # -log10 thresholds: 0.05→1.301, 0.01→2, 0.001→3
  # Negative direction (NES < 0): down-regulated, blue palette
  # Positive direction (NES > 0): up-regulated, red palette
  col_fun = circlize::colorRamp2(
    c(-3, -2, -1.30103, 0, 1.30103, 2, 3),
    c(
      "#2166AC",
      "#67A9CF",
      "#D1E5F0",
      "#E8E8E8",
      "#FDDBC7",
      "#EF8A62",
      "#B2182B"
    )
  ),
  filename = "gsea_res_all_stats.png",
  chr_width = 30L,
  col_order = NULL,
  ...
) {
  combined_stats_file <- dplyr::mutate(
    combined_stats_file,
    neg_log10_padj = -log10(padj),
    neg_log10_padj = dplyr::case_when(
      NES < 0 ~ -neg_log10_padj,
      NES > 0 ~ neg_log10_padj,
      TRUE ~ NA
    )
  )

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
    ) %>%
    dplyr::mutate(
      y_axis = paste(tumor, pheno, bulk, screen_method, sep = "_")
    )

  wide_mat_data <- mat_data %>%
    dplyr::select(-tumor, -pheno, -bulk, -screen_method) %>%
    tidyr::pivot_wider(
      names_from = y_axis,
      values_from = !!dplyr::sym(metrics)
    ) %>%
    t()

  colnames(wide_mat_data) <- wide_mat_data[1, ]
  wide_mat_data <- wide_mat_data[-1, -1]

  wide_mat_data <- as.data.frame(wide_mat_data) %>%
    dplyr::mutate(dplyr::across(dplyr::everything(), base::as.numeric)) %>%
    dplyr::filter(!dplyr::if_all(dplyr::everything(), base::is.na)) %>%
    as.matrix()

  # 唯一区分标签
  dataset_label <- rownames(wide_mat_data)

  #   palette_tumor <- c(
  #     "BRCA_HER2" = "#773232",
  #     "TNBC" = "#12298f",
  #     "LUNG" = "#23964f",
  #     "OV" = "#5a229b"
  #   )

  bulk_names <- unique(mat_data$bulk)
  n_bulk <- length(bulk_names)
  palette_bulk <- setNames(
    c(
      "#B8D9A0",
      "#E8E0DB",
      "#A8D9D8",
      "#F5C0A8",
      "#D5C9E8",
      "#88D5B0",
      "#B8A8E5",
      "#D5E5C0",
      "#E8D8E8",
      "#9BBBD9"
    )[seq_len(n_bulk)],
    unique(mat_data$bulk)
  )

  row_anno <- ComplexHeatmap::rowAnnotation(
    # "Tumor" = stringr::str_extract(dataset_label, "BRCA_HER2|TNBC|LUNG|OV"),
    "Bulk Data" = stringr::str_extract(dataset_label, "TCGA_[^_]*|GSE[0-9]*"),
    annotation_name_gp = grid::gpar(fontface = "bold"),
    annotation_name_rot = 60,
    col = list(
      #   "Tumor" = palette_tumor,
      "Bulk Data" = palette_bulk
    ),
    gap = grid::unit(0.6, "mm")
  )

  # Cap values > 3 (corresponding to padj < 0.001, -log10(0.001) = 3)
  wide_mat_data[wide_mat_data > 3] <- 3
  wide_mat_data[wide_mat_data < -3] <- -3

  # significance legend: bi-directional discrete bins
  # Negative direction (down-regulated): blue palette
  # Positive direction (up-regulated): red palette
  lgd <- ComplexHeatmap::Legend(
    labels = c(
      "< 0.001 (down)",
      "< 0.01 (down)",
      "< 0.05 (down)",
      ">= 0.05",
      "< 0.05 (up)",
      "< 0.01 (up)",
      "< 0.001 (up)"
    ),
    title = ComplexHeatmap::gt_render(
      "<span> log10(p-adjusted) <br> with direction</span>"
    ),
    legend_gp = grid::gpar(
      fill = c(
        "#2166AC",
        "#67A9CF",
        "#D1E5F0",
        "#E8E8E8",
        "#FDDBC7",
        "#EF8A62",
        "#B2182B"
      )
    )
  )

  htmap <- ComplexHeatmap::Heatmap(
    matrix = wide_mat_data,
    name = "-log10(padj)",
    col = col_fun,
    na_col = "#BBBBBB",

    cluster_rows = TRUE,
    cluster_columns = FALSE,
    column_dend_height = grid::unit(3, "cm"),
    row_dend_width = grid::unit(3, "cm"),

    # top_annotation = ComplexHeatmap::columnAnnotation(
    #   Distribution = ComplexHeatmap::anno_boxplot(
    #     wide_mat_data,
    #     gp = grid::gpar(
    #       fill = c(
    #         "#B8D9A0",
    #         "#E8E0DB",
    #         "#A8D9D8",
    #         "#F5C0A8",
    #         "#D5C9E8",
    #         "#9BBBD9",
    #         "#B8A8E5",
    #         "#D5E5C0",
    #         "#E8D8E8",
    #         "#88D5B0",
    #         "#E5B8E5",
    #         "#88E0D5",
    #         "#A8D5A8",
    #         "#E5A8D0",
    #         "#F5C8E5",
    #         "#F0E890",
    #         "#F5E8B8",
    #         "#B8BCA0",
    #         "#C8E5D0",
    #         "#B888D5",
    #         "#C5B8E8",
    #         "#E5C090",
    #         "#D5A0A8",
    #         "#88A8D5",
    #         "#F5D890",
    #         "#D8E5B0",
    #         "#B8D5D0",
    #         "#E0E8D8",
    #         "#E0A8E5",
    #         "#A888D5",
    #         "#C5B0B8",
    #         "#E5A8E8",
    #         "#E8A8D5",
    #         "#D0C5B8",
    #         "#F5A098",
    #         "#B0E5D0",
    #         "#C0D8E5",
    #         "#90C8D5",
    #         "#E8A0B0",
    #         "#F5E0B8"
    #       )
    #     ),
    #   ),
    #   height = grid::unit(5, "cm")
    # ),
    left_annotation = c(row_anno),

    # column_title = "Hallmarks Gene Set",
    # column_title_side = "bottom",
    column_names_gp = grid::gpar(fontsize = 14),
    row_names_gp = grid::gpar(fontsize = 14),

    row_labels = gsub(".*_", "", dataset_label) %>%
      stringr::str_replace("SGL", "LP_SGL"),
    row_split = stringr::str_replace(rownames(wide_mat_data), "_SGL", "SGL") %>%
      gsub("_.*", "", .) %>%
      stringr::str_replace("BRCA", "BRCA HER2"),
    show_parent_dend_line = FALSE,

    column_labels = wrap_every_n_chars(
      colnames(wide_mat_data),
      width = chr_width
    ),
    column_names_rot = 60,
    column_order = col_order,

    heatmap_height = grid::unit(0.1, "npc"),
    height = grid::unit(0.7, "npc"),
    width = grid::unit(0.85, "npc"),

    show_heatmap_legend = FALSE
  )

  Cairo::CairoPNG(
    filename = filename,
    dpi = 400,
    width = 10000,
    height = 10000
  )
  ComplexHeatmap::draw(htmap, annotation_legend_list = list(lgd))
  dev.off()
}

# rpkgkit::make_func_call_explicit(use_packages = c("ComplexHeatmap", "grid"))
