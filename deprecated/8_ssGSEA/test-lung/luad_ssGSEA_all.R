# ! 运行此脚本前，先运行/home/yyx/R/Project/R_code/SigBridgeR/Tmp/ssGSEA/test-luad/luad_ssGSEA_survival.R
# ! 把框架打好

library(ComplexHeatmap)

# ? 这个脚本的作用是将survival和bianry phenotype列注释结合起来
es.mat = qs::qread("luad_ssgsea_score.qs")

CreateHtmapAnno.survival <- function(
  seurat_obj,
  es.mat,
  color = "#123123"
) {
  # 提取元数据
  meta <- seurat_obj@meta.data %>%
    tibble::rownames_to_column("cell")

  # 处理表达矩阵并合并元数据
  es.df <- t(es.mat) %>%
    as.data.frame() %>%
    tibble::rownames_to_column("cell") %>%
    left_join(meta, by = "cell")

  # 定义可能的注释列
  possible_annotations <- c(
    "scissor",
    "scPAS",
    "scAB",
    "scPP",
    "DEGAS",
    "LP_SGL",
    "PIPET"
  )

  # 检查哪些注释列存在
  existing_annotations <- possible_annotations %>%
    purrr::keep(~ .x %in% colnames(es.df))

  # 如果没有找到任何注释列，返回空的 HeatmapAnnotation
  if (length(existing_annotations) == 0) {
    return(HeatmapAnnotation())
  }

  # 创建注释列表
  annotation_list <- list()
  color_list <- list()

  # 为每个存在的注释列配置颜色
  if ("scissor" %in% existing_annotations) {
    annotation_list$scissor_survival <- es.df$scissor
    color_list$scissor_survival <- c(
      "Neutral" = "#CECECE",
      "Positive" = "#ff3333",
      "Negative" = "#386c9b"
    )
  }

  if ("scPAS" %in% existing_annotations) {
    annotation_list$scPAS_survival <- es.df$scPAS
    color_list$scPAS_survival <- c(
      "Neutral" = "#CECECE",
      "Positive" = "#ff3333",
      "Negative" = "#386c9b"
    )
  }

  if ("scAB" %in% existing_annotations) {
    annotation_list$scAB_survival <- es.df$scAB
    color_list$scAB_survival <- c(
      "Other" = "#CECECE",
      "Positive" = "#ff3333"
    )
  }

  if ("scPP" %in% existing_annotations) {
    annotation_list$scPP_survival <- es.df$scPP
    color_list$scPP_survival <- c(
      "Neutral" = "#CECECE",
      "Positive" = "#ff3333",
      "Negative" = "#386c9b"
    )
  }
  if ("DEGAS" %in% existing_annotations) {
    annotation_list$DEGAS_survival <- es.df$DEGAS
    color_list$DEGAS_survival <- c(
      "Other" = "#CECECE",
      "Positive" = "#ff3333"
    )
  }

  if ("LP_SGL" %in% existing_annotations) {
    annotation_list$LP_SGL_survival <- es.df$LP_SGL
    color_list$LP_SGL_survival <- c(
      "Neutral" = "#CECECE",
      "Positive" = "#ff3333",
      "Negative" = "#386c9b"
    )
  }

  if ("PIPET" %in% existing_annotations) {
    annotation_list$PIPET_survival <- es.df$PIPET
    color_list$PIPET_survival <- c(
      "Neutral" = "#CECECE",
      "Positive" = "#ff3333",
      "Negative" = "#386c9b"
    )
  }

  # 创建 HeatmapAnnotation
  HeatmapAnnotation(
    df = annotation_list,
    col = color_list,
    annotation_label = list(
      scissor_survival = "scissor survival",
      scPAS_survival = "scPAS survival",
      scAB_survival = "scAB survival",
      scPP_survival = "scPP survival",
      DEGAS_survival = "DEGAS survival",
      LP_SGL_survival = "LP_SGL survival",
      PIPET_survival = "PIPET survival"
    ),
    annotation_name_gp = gpar(col = color, fontface = "bold"),
    show_legend = FALSE,
    height = unit(1, "cm"),
    gap = unit(0.6, "mm")
  )
}

CreateHtmapAnno.binary <- function(
  seurat_obj,
  es.mat,
  color = "#123123"
) {
  # 提取元数据
  meta <- seurat_obj@meta.data %>%
    tibble::rownames_to_column("cell")

  # 处理表达矩阵并合并元数据
  es.df <- t(es.mat) %>%
    as.data.frame() %>%
    tibble::rownames_to_column("cell") %>%
    left_join(meta, by = "cell")

  # 定义可能的注释列
  possible_annotations <- c(
    "scissor",
    "scPAS",
    "scAB",
    "scPP",
    "DEGAS",
    "LP_SGL",
    "PIPET"
  )

  # 检查哪些注释列存在
  existing_annotations <- possible_annotations %>%
    purrr::keep(~ .x %in% colnames(es.df))

  # 如果没有找到任何注释列，返回空的 HeatmapAnnotation
  if (length(existing_annotations) == 0) {
    return(HeatmapAnnotation())
  }

  # 创建注释列表
  annotation_list <- list()
  color_list <- list()

  # 为每个存在的注释列配置颜色
  if ("scissor" %in% existing_annotations) {
    annotation_list$scissor_binary <- es.df$scissor
    color_list$scissor_binary <- c(
      "Neutral" = "#CECECE",
      "Positive" = "#ff3333",
      "Negative" = "#386c9b"
    )
  }

  if ("scPAS" %in% existing_annotations) {
    annotation_list$scPAS_binary <- es.df$scPAS
    color_list$scPAS_binary <- c(
      "Neutral" = "#CECECE",
      "Positive" = "#ff3333",
      "Negative" = "#386c9b"
    )
  }

  if ("scAB" %in% existing_annotations) {
    annotation_list$scAB_binary <- es.df$scAB
    color_list$scAB_binary <- c(
      "Other" = "#CECECE",
      "Positive" = "#ff3333"
    )
  }

  if ("scPP" %in% existing_annotations) {
    annotation_list$scPP_binary <- es.df$scPP
    color_list$scPP_binary <- c(
      "Neutral" = "#CECECE",
      "Positive" = "#ff3333",
      "Negative" = "#386c9b"
    )
  }
  if ("DEGAS" %in% existing_annotations) {
    annotation_list$DEGAS_binary <- es.df$DEGAS
    color_list$DEGAS_binary <- c(
      "Other" = "#CECECE",
      "Positive" = "#ff3333"
    )
  }

  if ("LP_SGL" %in% existing_annotations) {
    annotation_list$LP_SGL_binary <- es.df$LP_SGL
    color_list$LP_SGL_binary <- c(
      "Neutral" = "#CECECE",
      "Positive" = "#ff3333",
      "Negative" = "#386c9b"
    )
  }

  if ("PIPET" %in% existing_annotations) {
    annotation_list$PIPET_binary <- es.df$PIPET
    color_list$PIPET_binary <- c(
      "Neutral" = "#CECECE",
      "Positive" = "#ff3333",
      "Negative" = "#386c9b"
    )
  }

  # 创建 HeatmapAnnotation
  HeatmapAnnotation(
    df = annotation_list,
    col = color_list,
    annotation_label = list(
      scissor_binary = "scissor binary",
      scPAS_binary = "scPAS binary",
      scAB_binary = "scAB binary",
      scPP_binary = "scPP binary",
      DEGAS_binary = "DEGAS binary",
      LP_SGL_binary = "LP_SGL binary",
      PIPET_binary = "PIPET binary"
    ),
    annotation_name_gp = gpar(col = color, fontface = "bold"),
    show_legend = FALSE,
    height = unit(1, "cm"),
    gap = unit(0.6, "mm")
  )
}

GSE3141_survival = qs::qread(
  "/home/data/sigbridger/benchmark_data/lung/GSE3141/GSE3141_luad_merged_seurat.qs"
)

GSE8894_survival = qs::qread(
  "/home/data/sigbridger/benchmark_data/lung/GSE8894/GSE8894_luad_merged_seurat.qs"
)
GSE31210_survival = qs::qread(
  "/home/data/sigbridger/benchmark_data/lung/GSE31210/GSE31210_luad_merged_seurat.qs"
)
tcga_survival = qs::qread(
  "/home/data/sigbridger/benchmark_data/lung/TCGA-LUAD/tcga_luad_merged_seurat.qs"
)

col_anno_GSE3141_sur = CreateHtmapAnno.survival(
  GSE3141_survival,
  es.mat = es.mat,
  color = "#00913fff"
)

col_anno_GSE8894_sur = CreateHtmapAnno.survival(
  GSE8894_survival,
  es.mat = es.mat,
  color = "#cf7815ff"
)
col_anno_GSE31210_sur = CreateHtmapAnno.survival(
  GSE31210_survival,
  es.mat = es.mat,
  color = "#bf009cff"
)
col_anno_tcga_sur = CreateHtmapAnno.survival(
  tcga_survival,
  es.mat = es.mat,
  color = "#00c7baff"
)

tcga_binary <- qs::qread(
  "/home/data/sigbridger/benchmark_binary/lung/TCGA-LUAD/tcga_luad_merged_seurat.qs",
  nthreads = 4
)

col_anno_tcga_binary = CreateHtmapAnno.binary(
  tcga_binary,
  es.mat = es.mat,
  color = "#00c7baff"
)

meta <- GSE3141_survival[[]]
meta$cell <- colnames(GSE3141_survival)
es.df <- t(es.mat) %>% # 行为细胞，列为 cluster
  as.data.frame() %>%
  tibble::rownames_to_column("cell") %>%
  left_join(meta, by = "cell")


order <- order(
  es.df$seurat_clusters,
  es.df$lung_cancer,
  es.df$Celltype,
  es.df$cnv_status,
  es.df$Sample
)
es.mat_ordered <- es.mat[, order, drop = FALSE]

bulk_lgd = Legend(
  at = c("GSE3141", "GSE8894", "GSE31210", "TCGA"),
  title = "Bulk Sample for Methods",
  legend_gp = gpar(
    fill = c("#00913fff", "#cf7815ff", "#bf009cff", "#00c7baff")
  )
)


# ! cell type + cnv status

pallete_celltype <- randomcoloR::distinctColorPalette(
  length(unique(es.df$Celltype)),
  runTsne = T
)
col_anno_types <- HeatmapAnnotation(
  cell_type = es.df$Celltype[order],
  cnv_status = es.df$cnv_status[order],
  annotation_name_gp = gpar(fontface = "bold"),
  col = list(
    cell_type = setNames(pallete_celltype, unique(es.df$Celltype)),
    cnv_status = setNames(pallete_celltype[1:2], unique(es.df$cnv_status))
  ),
  height = unit(1, "cm"),
  gap = unit(0.6, "mm")
)

# ! cluster

pallete_cluster <- randomcoloR::distinctColorPalette(
  length(unique(es.df$seurat_clusters)),
  runTsne = T
)
col_anno_cluster <- HeatmapAnnotation(
  cluster = es.df$seurat_clusters[order],
  annotation_name_gp = gpar(fontface = "bold"),
  col = list(
    cluster = setNames(pallete_cluster, unique(es.df$seurat_clusters))
  ),
  height = unit(1, "cm"),
  gap = unit(0.6, "mm")
)


# ! sample

pallete_sample <- randomcoloR::distinctColorPalette(
  length(unique(es.df$Sample)),
  runTsne = T
)
col_anno_sample <- HeatmapAnnotation(
  Sample = es.df$Sample[order],
  annotation_name_gp = gpar(fontface = "bold"),
  col = list(
    Sample = setNames(pallete_sample, unique(es.df$Sample))
  ),
  height = unit(1, "cm"),
  gap = unit(0.6, "mm")
)


# ! ssGSEA score distribution

col_anno_distribution <- HeatmapAnnotation(
  ssGSEA = anno_points(
    as.vector(es.mat_ordered),
    height = unit(4, "cm"),
    size = unit(0.5, "mm"),
    ylim = c(-1, 1),
    axis_param = list(
      at = c(-0.9, -0.6, -0.3, 0, 0.3, 0.6, 0.9)
    ),
    gp = gpar(col = "#9f9f9fff", alpha = 0.5),
  ),
  gap = unit(0.6, "mm")
)

screen_lgd = Legend(
  at = c("Positive", "Negative", "Neutral", "Other"),
  title = "Screen Result Group",
  legend_gp = gpar(
    fill = c(
      "#ff3333",
      "#386c9b",
      "#CECECE",
      "#CECECE"
    )
  ),
)
# ! 间隔
row_anno_gap = HeatmapAnnotation(
  empty = anno_empty(border = FALSE), # 让热图扁一点
  height = unit(0.6, "mm"),
  gap = unit(0.6, "mm")
)

# ! 依据cluster 分组
column_split = es.df$seurat_clusters |> sort()

# # ! mean score of each cluster ,按照分数大小排序了，abandoned
cluster_mean_score = purrr::map_dbl(
  sort(unique(es.df$seurat_clusters)),
  function(x) {
    mean(es.df$lung_cancer[es.df$seurat_clusters == x])
  }
)
cluster_mean_score = round(cluster_mean_score, 2)


htmap_cmb = ComplexHeatmap::Heatmap(
  es.mat_ordered,
  name = "ssGSEA\nScore",
  top_annotation = c(
    col_anno_cluster,
    row_anno_gap,
    col_anno_types,
    row_anno_gap,
    col_anno_sample,
    row_anno_gap,
    col_anno_GSE3141_sur,
    row_anno_gap,
    col_anno_GSE8894_sur,
    row_anno_gap,
    col_anno_GSE31210_sur,
    row_anno_gap,
    col_anno_tcga_sur,
    row_anno_gap,
    col_anno_tcga_binary,
    row_anno_gap,
    col_anno_distribution
  ),
  right_annotation = rowAnnotation(
    ssGSEA = anno_empty(border = FALSE)
  ),
  show_column_names = FALSE,
  show_row_names = FALSE,
  cluster_columns = FALSE,
  cluster_rows = FALSE,
  use_raster = TRUE,
  column_split = column_split,
  column_title_rot = 90,
  column_title_side = "bottom",
  column_title_gp = gpar(fontface = "bold", fontsize = 8),
  column_title = cluster_mean_score,
  heatmap_width = unit(3, "npc")
)
# draw(
#     htmap_cmb,
#     annotation_legend_list = list(bulk_lgd, screen_lgd)
# )

ragg::agg_png(
  file.path(plot_dir, "../lung_ssGSEA_heatmap_all.png"),
  width = 1600,
  height = 1000
)
draw(
  htmap_cmb,
  annotation_legend_list = list(bulk_lgd, screen_lgd)
)
dev.off()
