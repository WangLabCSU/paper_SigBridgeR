setwd(file.path(usethis::proj_path(), "5_GSEA/all_viz"))

library(dplyr)

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

combined_stats_file <- data.table::fread("gsea_res_all_stats.csv")


# * 缩写基因集名称映射函数
shorten_gene_sets <- function(x) {
  # x: character vector，原始列名/基因集名称
  short <- x

  # 1. 保留方向性修饰词，用箭头/缩写替代
  short <- gsub("Genes up-regulated", "UP:", short, ignore.case = TRUE)
  short <- gsub("Genes down-regulated", "DOWN:", short, ignore.case = TRUE)
  short <- gsub("up-regulated", "UP:", short, ignore.case = TRUE)
  short <- gsub("down-regulated", "DOWN:", short, ignore.case = TRUE)

  # 2. UP:/DOWN: 标签后去除多余的前置词
  short <- gsub(
    "(UP:|DOWN:)\\s*(by |during |in response to |through |of )",
    "\\1",
    short
  )

  # 3. 去除 "Genes " 前缀（普通描述性前缀）
  short <- gsub(
    "Genes (encoding |involved in |defining |defined |important for |regulating |mediating |involve in )",
    "",
    short,
    ignore.case = TRUE
  )
  short <- gsub("Genes specifically ", "", short, ignore.case = TRUE)
  short <- gsub("^Genes ", "", short)

  # 4. 去除其他冗余前缀
  short <- gsub(
    "(A subgroup of |processing of |production of |formation of |development of |metabolism of |components of )",
    "",
    short,
    ignore.case = TRUE
  )
  short <- gsub("^genes ", "", short)

  # 5. 去除长注释后缀
  short <- gsub(
    "(, as in.*|, e\\.g\\..*|, a cellular.*|, which.*|; also.*)",
    "",
    short
  )

  # 6. 去除尾部不完整的标点/空格
  short <- gsub("[,\\s]+$", "", short)
  short <- gsub("\\.$", "", short)

  # 7. 截断过长的名称
  short <- ifelse(
    nchar(short) > 45,
    paste0(substr(short, 1, 42), "..."),
    short
  )

  short
}

plot_heatmap2 <- function(
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
    as.matrix()

  # 简化基因集列名
  colnames(wide_mat_data) <- shorten_gene_sets(colnames(wide_mat_data))

  # 唯一区分标签
  dataset_label <- rownames(wide_mat_data)

  palette_tumor <- c(
    "BRCA_HER2" = "#773232",
    "TNBC" = "#12298f",
    "LUNG" = "#23964f",
    "OV" = "#5a229b"
  )
  palette_pheno <- c(
    "survival" = "#ddcbaf",
    "binary" = "#b3cad4"
  )
  palette_bulk <- setNames(
    c(
      "#B8D9A0",
      "#E8E0DB",
      "#A8D9D8",
      "#F5C0A8",
      "#D5C9E8",
      "#9BBBD9",
      "#B8A8E5",
      "#D5E5C0",
      "#E8D8E8",
      "#88D5B0"
    ),
    unique(mat_data$bulk)
  )

  row_anno_tissue <- ComplexHeatmap::rowAnnotation(
    "Phenotype" = stringr::str_extract(dataset_label, "binary|survival"),
    "Bulk Data" = stringr::str_extract(dataset_label, "TCGA_[^_]*|GSE[0-9]*"),
    annotation_name_gp = grid::gpar(fontface = "bold"),
    col = list(
      "Phenotype" = palette_pheno,
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
    },
    col = col_fun,
    na_col = "#FFFFFF",

    cluster_rows = FALSE,
    cluster_columns = TRUE,

    top_annotation = ComplexHeatmap::columnAnnotation(
      Dist = ComplexHeatmap::anno_boxplot(
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
    left_annotation = c(row_anno_tissue),

    # column_title = "Hallmarks Gene Set",
    # column_title_side = "bottom",
    column_names_gp = grid::gpar(fontsize = 10),
    row_names_gp = grid::gpar(fontsize = 10),

    row_labels = gsub(".*_", "", dataset_label) %>%
      stringr::str_replace("SGL", "LP_SGL"),
    row_split = stringr::str_replace(rownames(wide_mat_data), "_SGL", "SGL") %>%
      gsub("_.*", "", .),

    heatmap_height = grid::unit(0.1, "npc"),
    height = grid::unit(0.75, "npc"),
  )

  CairoPNG("gsea_res_all_stats.png", dpi = 400, width = 4000, height = 10000)
  ComplexHeatmap::draw(htmap)
  dev.off()
}
