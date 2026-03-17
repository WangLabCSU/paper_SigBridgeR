# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd("/home/yyx/R/Project/R_code/SigBridgeR/Tmp/ssGSEA/test-lung")

set.seed(123)

library(Seurat)
library(GSVA)
library(ComplexHeatmap)
library(tidyr)
library(dplyr)
library(ggplot2)

# ! luad-sc-123902, bulk-GSE3141, binary phenotype

seurat <- qs::qread(
  "/home/data/sigbridger/benchmark_binary/lung/TCGA-LUAD/tcga_luad_merged_seurat.qs",
  nthreads = 4
)
es.mat = qs::qread("luad_ssgsea_score.qs")

meta <- seurat[[]]
meta$cell <- colnames(seurat)
es.df <- t(es.mat) %>% # 行为细胞，列为 cluster
  as.data.frame() %>%
  tibble::rownames_to_column("cell") %>%
  left_join(meta, by = "cell")

es.df$scissor <- factor(
  es.df$scissor,
  levels = c("Neutral", "Negative", "Positive")
)
es.df$scPAS <- factor(
  es.df$scPAS,
  levels = c("Neutral", "Negative", "Positive")
)
es.df$scAB <- factor(
  es.df$scAB,
  levels = c("Other", "Positive")
)

table(es.df$scissor)
#  Neutral Negative Positive
#    30687      230     1261

table(es.df$scPAS)
#  Neutral Negative Positive
#    31628       70      480

table(es.df$scAB)
#    Other Positive
#  105    32073

# 7. limma 差异检验（以 cluster0 为例）
library(limma)
scissor_design <- model.matrix(~scissor, data = es.df)
scpas_design <- model.matrix(~scPAS, data = es.df)
scab_design <- model.matrix(~scAB, data = es.df)


scissor_fit <- lmFit(es.mat, scissor_design)
scpas_fit <- lmFit(es.mat, scpas_design)
scab_fit <- lmFit(es.mat, scab_design)

scissor_fit <- eBayes(scissor_fit)
scpas_fit <- eBayes(scpas_fit)
scab_fit <- eBayes(scab_fit)

scissor_toptable_pos = topTable(
  scissor_fit,
  coef = "scissorPositive",
  number = 100
) # 比较Positive vs Neutral
scpas_toptable_pos = topTable(scpas_fit, coef = "scPASPositive", number = 100) # 比较Positive vs Neutral
scab_toptable_pos = topTable(scab_fit, coef = "scABPositive", number = 100) # 比较Positive vs Other
# # 或者
# toptable_neg = topTable(fit, coef = "scissorNegative", number = 100) # 比较Negative vs Neutral

plot_dir = "/home/yyx/R/Project/R_code/SigBridgeR/Tmp/ssGSEA/test-lung/binary_plot"
if (!dir.exists(plot_dir)) {
  dir.create(plot_dir, recursive = TRUE)
}

# # ! 可视化-基础热图

order <- order(
  es.df$seurat_clusters,
  es.df$lung_cancer,
  es.df$Celltype,
  es.df$cnv_status,
  es.df$Sample
)
es.mat_ordered <- es.mat[, order, drop = FALSE]

# ! 筛选

# # 列注释（细胞分组）
col_anno_screen <- HeatmapAnnotation(
  scissor = es.df$scissor[order],
  scPAS = es.df$scPAS[order],
  scAB = es.df$scAB[order],
  # scPP = es.df$scPP[order],
  annotation_name_gp = gpar(col = "#00913fff", fontface = "bold"),
  show_legend = FALSE,
  col = list(
    scissor = c(
      "Neutral" = "#CECECE",
      "Positive" = "#ff3333",
      "Negative" = "#386c9b"
    ),
    scPAS = c(
      "Neutral" = "#CECECE",
      "Positive" = "#ff3333",
      "Negative" = "#386c9b"
    ),
    scAB = c(
      "Other" = "#CECECE",
      "Positive" = "#ff3333"
    )
    # ,
    # scPP = c(
    #     "Neutral" = "#CECECE",
    #     "Positive" = "#ff3333",
    #     "Negative" = "#386c9b"
    # )
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


# ! 细胞类型 + cnv_status
set.seed(123)
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

CreateHtmapAnno <- function(seurat_obj, es.mat, color = "#123123") {
  # 提取元数据
  meta <- seurat_obj@meta.data %>%
    tibble::rownames_to_column("cell")

  # 处理表达矩阵并合并元数据
  es.df <- t(es.mat) %>%
    as.data.frame() %>%
    tibble::rownames_to_column("cell") %>%
    left_join(meta, by = "cell")

  # 定义可能的注释列
  possible_annotations <- c("scissor", "scPAS", "scAB", "scPP")

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
    annotation_list$scissor <- es.df$scissor
    color_list$scissor <- c(
      "Neutral" = "#CECECE",
      "Positive" = "#ff3333",
      "Negative" = "#386c9b"
    )
  }

  if ("scPAS" %in% existing_annotations) {
    annotation_list$scPAS <- es.df$scPAS
    color_list$scPAS <- c(
      "Neutral" = "#CECECE",
      "Positive" = "#ff3333",
      "Negative" = "#386c9b"
    )
  }

  if ("scAB" %in% existing_annotations) {
    annotation_list$scAB <- es.df$scAB
    color_list$scAB <- c(
      "Other" = "#CECECE",
      "Positive" = "#ff3333"
    )
  }

  if ("scPP" %in% existing_annotations) {
    annotation_list$scPP <- es.df$scPP
    color_list$scPP <- c(
      "Neutral" = "#CECECE",
      "Positive" = "#ff3333",
      "Negative" = "#386c9b"
    )
  }

  # 创建 HeatmapAnnotation
  HeatmapAnnotation(
    df = annotation_list,
    col = color_list,
    annotation_name_gp = gpar(col = color, fontface = "bold"),
    show_legend = FALSE,
    height = unit(1, "cm"),
    gap = unit(0.6, "mm")
  )
}

# ? 继续组合heatmap
bulk_lgd = Legend(
  at = c(
    # "GSE3141",
    # "GSE8894",
    # "GSE31210",
    "TCGA"
  ),
  title = "Bulk Sample for Methods",
  legend_gp = gpar(
    fill = c(
      "#00913fff"
      # , "#cf7815ff", "#bf009cff", "#00c7baff"
    )
  )
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
  )
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
    col_anno_screen,
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
  heatmap_width = unit(3, "npc"),
  row_title = "Binary phenotype"
)
draw(
  htmap_cmb,
  annotation_legend_list = list(bulk_lgd, screen_lgd)
)

ragg::agg_png(
  file.path(plot_dir, "lung_ssGSEA_heatmap_binary.png"),
  width = 1600,
  height = 1000
)
draw(
  htmap_cmb,
  annotation_legend_list = list(bulk_lgd, screen_lgd)
)
dev.off()

# -----------------------------
# 平均得分
MeanScore = function(seurat_obj, es.mat, sample = "GSE") {
  meta <- seurat_obj@meta.data %>%
    tibble::rownames_to_column("cell")

  # 处理表达矩阵并合并元数据
  es.df <- t(es.mat) %>%
    as.data.frame() %>%
    tibble::rownames_to_column("cell") %>%
    left_join(meta, by = "cell")

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

  mean_scores_scissor <- if ("scissor" %in% existing_annotations) {
    es.df %>%
      group_by(scissor) %>%
      summarise(across("lung_cancer", mean, na.rm = TRUE)) %>%
      pivot_longer(
        cols = "lung_cancer",
        names_to = "cluster",
        values_to = "mean_score"
      ) %>%
      pivot_longer(
        cols = "scissor",
        names_to = "screen_method",
        values_to = "screen_ressult"
      )
  } else {
    NULL
  }
  mean_scores_scpas <- if ("scPAS" %in% existing_annotations) {
    es.df %>%
      group_by(scPAS) %>%
      summarise(across("lung_cancer", mean, na.rm = TRUE)) %>%
      pivot_longer(
        cols = "lung_cancer",
        names_to = "cluster",
        values_to = "mean_score"
      ) %>%
      pivot_longer(
        cols = "scPAS",
        names_to = "screen_method",
        values_to = "screen_ressult"
      )
  } else {
    NULL
  }
  mean_scores_scab <- if ("scAB" %in% existing_annotations) {
    es.df %>%
      group_by(scAB) %>%
      summarise(across("lung_cancer", mean, na.rm = TRUE)) %>%
      pivot_longer(
        cols = "lung_cancer",
        names_to = "cluster",
        values_to = "mean_score"
      ) %>%
      pivot_longer(
        cols = "scAB",
        names_to = "screen_method",
        values_to = "screen_ressult"
      )
  } else {
    NULL
  }
  mean_scores_scpp <- if ("scPP" %in% existing_annotations) {
    es.df %>%
      group_by(scPP) %>%
      summarise(across("lung_cancer", mean, na.rm = TRUE)) %>%
      pivot_longer(
        cols = "lung_cancer",
        names_to = "cluster",
        values_to = "mean_score"
      ) %>%
      pivot_longer(
        cols = "scPP",
        names_to = "screen_method",
        values_to = "screen_ressult"
      )
  } else {
    NULL
  }
  mean_scores_degas <- if ("DEGAS" %in% existing_annotations) {
    es.df %>%
      group_by(DEGAS) %>%
      summarise(across("lung_cancer", mean, na.rm = TRUE)) %>%
      pivot_longer(
        cols = "lung_cancer",
        names_to = "cluster",
        values_to = "mean_score"
      ) %>%
      pivot_longer(
        cols = "DEGAS",
        names_to = "screen_method",
        values_to = "screen_ressult"
      )
  } else {
    NULL
  }
  mean_scores_lp_sgl <- if ("LP_SGL" %in% existing_annotations) {
    es.df %>%
      group_by(LP_SGL) %>%
      summarise(across("lung_cancer", mean, na.rm = TRUE)) %>%
      pivot_longer(
        cols = "lung_cancer",
        names_to = "cluster",
        values_to = "mean_score"
      ) %>%
      pivot_longer(
        cols = "LP_SGL",
        names_to = "screen_method",
        values_to = "screen_ressult"
      )
  } else {
    NULL
  }
  mean_scores_pipet <- if ("PIPET" %in% existing_annotations) {
    es.df %>%
      group_by(PIPET) %>%
      summarise(across("lung_cancer", mean, na.rm = TRUE)) %>%
      pivot_longer(
        cols = "lung_cancer",
        names_to = "cluster",
        values_to = "mean_score"
      ) %>%
      pivot_longer(
        cols = "PIPET",
        names_to = "screen_method",
        values_to = "screen_ressult"
      )
  } else {
    NULL
  }

  mean_scores <- rbind(
    mean_scores_scissor,
    mean_scores_scpas,
    mean_scores_scab,
    mean_scores_scpp,
    mean_scores_degas,
    mean_scores_lp_sgl,
    mean_scores_pipet
  ) %>%
    unite("screen", screen_method, screen_ressult, sep = "_") %>%
    mutate("Sample" = sample)
}


mean_score_ctga = MeanScore(seurat, es.mat, "TCGA")

p <- ggplot(
  mean_score_ctga,
  aes(x = Sample, y = screen, fill = mean_score)
) +
  geom_point(size = 6, alpha = 0.9, shape = 21, color = "black") +
  geom_text(
    aes(label = round(mean_score, 3)),
    color = "black",
    size = 3,
    hjust = -0.5
  ) +
  scale_fill_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    name = "Mean Score"
  ) +
  labs(
    title = "Mean ssGSEA Scores by Bulk Sample Group",
    subtitle = "sc-GSE123902, survival as phenotype",
    y = "Mean ssGSEA Score",
    x = "Bulk Sample"
  ) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(
  file.path(plot_dir, "lung_mean_scores_binary.png"),
  plot = p,
  width = 4,
  height = 7,
  dpi = 400
)
