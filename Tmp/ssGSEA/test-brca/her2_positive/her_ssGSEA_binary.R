library(data.table)
library(Seurat)
library(GSVA)
library(ComplexHeatmap)
library(dplyr)
library(tidyr)
library(ggplot2)

setwd(
  "/home/yyx/R/Project/R_code/SigBridgeR/Tmp/ssGSEA/test-brca/her2_positive"
)
set.seed(123)

es.mat <- qs::qread("her2_ssgsea_score.qs")

# ! brca-GSE161529_her2, bulk - GSE42568, pheno-binary-normal vs. cancer

seurat <- qs::qread(
  "/home/data/sigbridger/benchmark_binary/brca/HER2/GSE42568_brca_merged_seurat.qs",
  nthreads = 4
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

meta <- seurat@meta.data
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
# es.df$scPP <- factor(
#     es.df$scPP,
#     levels = c("Neutral", "Negative", "Positive")
# )

table(es.df$scissor)
#  Neutral Negative Positive
#    22103     3654     6115
table(es.df$scPAS)
#  Neutral Negative Positive
#    31272      600        0

table(es.df$scAB)
#    Other Positive
#    13631    18241
# table(es.df$scPP)

library(limma)
scissor_design <- model.matrix(~scissor, data = es.df)
scpas_design <- model.matrix(~scPAS, data = es.df)
scab_design <- model.matrix(~scAB, data = es.df)
# scpp_design <- model.matrix(~scPP, data = es.df)

scissor_fit <- lmFit(es.mat, scissor_design)
scpas_fit <- lmFit(es.mat, scpas_design)
scab_fit <- lmFit(es.mat, scab_design)
# scpp_fit <- lmFit(es.mat, scpp_design)

scissor_fit <- eBayes(scissor_fit)
scpas_fit <- eBayes(scpas_fit)
scab_fit <- eBayes(scab_fit)
# scpp_fit <- eBayes(scpp_fit)

scissor_toptable_pos = topTable(
  scissor_fit,
  coef = "scissorPositive",
  number = 100
) # 比较Positive vs Neutral
scpas_toptable_pos = topTable(scpas_fit, coef = "scPASPositive", number = 100) # ! NA 存在，因为没有Positive细胞
scab_toptable_pos = topTable(scab_fit, coef = "scABPositive", number = 100) # 比较Positive vs Other
# scpp_toptable_pos = topTable(scpp_fit, coef = "scPPPositive", number = 100) # 比较Positive vs Neutral
# # 或者
# toptable_neg = topTable(fit, coef = "scissorNegative", number = 100) # 比较Negative vs Neutral

plot_dir = "/home/yyx/R/Project/R_code/SigBridgeR/Tmp/ssGSEA/test-brca/her2_positive/binary_plot"
if (!dir.exists(plot_dir)) {
  dir.create(plot_dir, recursive = TRUE)
}

# # ! 可视化-基础热图
set.seed(123)
pallete <- randomcoloR::distinctColorPalette(
  length(unique(es.df$seurat_clusters)),
  runTsne = T
)
order <- order(
  es.df$seurat_clusters,
  es.df$Her2,
  es.df$scissor,
  # es.df$scPAS,
  es.df$scAB
  # ,es.df$scPP
)
es.mat_ordered <- es.mat[, order, drop = FALSE]

col_anno_GSE42568 <- CreateHtmapAnno(
  seurat,
  es.mat = es.mat_ordered,
  color = "#00913fff"
)

col_anno_cluster <- HeatmapAnnotation(
  cluster = es.df$seurat_clusters[order],
  annotation_name_gp = gpar(fontface = "bold"),
  col = list(
    cluster = setNames(pallete, unique(es.df$seurat_clusters))
  ),
  height = unit(1, "cm"),
  gap = unit(0.6, "mm")
)
col_anno_distribution <- HeatmapAnnotation(
  ssGSEA = anno_points(
    as.vector(es.mat_ordered),
    height = unit(4, "cm"),
    size = unit(0.7, "mm"),
    ylim = c(-1, 1),
    axis_param = list(
      at = c(-0.9, -0.6, -0.3, 0, 0.3, 0.6, 0.9)
    ),
    gp = gpar(col = "#9f9f9fff", alpha = 0.5)
  )
)


# ! bulk-GSE162228
seurat2 <- qs::qread(
  "/home/data/sigbridger/benchmark_binary/brca/HER2/GSE162228_brca_merged_seurat.qs"
)

col_anno_GSE162228 <- CreateHtmapAnno(
  seurat2,
  es.mat = es.mat_ordered,
  color = "#cf7815ff"
)

# ! HER2 Group
pallete_her2 <- randomcoloR::distinctColorPalette(
  length(unique(es.df$group)),
  runTsne = T
)
col_anno_her2 <- HeatmapAnnotation(
  Her2_group = es.df$group[order],
  annotation_name_gp = gpar(fontface = "bold"),
  annotation_label = "Her2 Group",
  col = list(
    HER2_group = setNames(pallete_her2, unique(es.df$group))
  ),
  height = unit(1, "cm"),
  gap = unit(0.6, "mm")
)


bulk_lgd = Legend(
  at = c("GSE42568", "GSE162228"),
  title = "Bulk Sample for Methods",
  legend_gp = gpar(
    fill = c("#00913fff", "#cf7815ff")
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
  ),
)

# ! 依据cluster 分组
column_split = es.df$seurat_clusters |> sort()

# # ! mean score of each cluster ,按照分数大小排序了，abandoned
cluster_mean_score = purrr::map_dbl(
  sort(unique(es.df$seurat_clusters)),
  function(x) {
    mean(es.df$Her2[es.df$seurat_clusters == x])
  }
)
cluster_mean_score = round(cluster_mean_score, 2)


# ! 每个anno之间的间距
row_anno_gap = HeatmapAnnotation(
  empty = anno_empty(border = FALSE), # 让热图扁一点
  height = unit(0.6, "mm")
)


htmap_cmb = ComplexHeatmap::Heatmap(
  es.mat_ordered,
  name = "ssGSEA\nScore",
  top_annotation = c(
    col_anno_cluster,
    row_anno_gap,
    col_anno_her2,
    row_anno_gap,
    col_anno_GSE42568,
    row_anno_gap,
    col_anno_GSE162228,
    row_anno_gap,
    col_anno_distribution,
    row_anno_gap
  ),
  right_annotation = rowAnnotation(
    ssGSEA = anno_empty(border = FALSE)
  ),
  # bottom_annotation = HeatmapAnnotation(
  #     empty = anno_empty(border = FALSE), # 让热图扁一点
  #     height = unit(6, "cm")
  # ),
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
draw(
  htmap_cmb,
  annotation_legend_list = list(bulk_lgd, screen_lgd)
)

# ----------------------------------------------------------------------
# ! 可视化-平均得分

MeanScore = function(seurat_obj, es.mat, sample = "GSE") {
  meta <- seurat_obj@meta.data %>%
    tibble::rownames_to_column("cell")

  # 处理表达矩阵并合并元数据
  es.df <- t(es.mat) %>%
    as.data.frame() %>%
    tibble::rownames_to_column("cell") %>%
    left_join(meta, by = "cell")
  mean_scores_scissor <- if ("scissor" %in% colnames(es.df)) {
    es.df %>%
      group_by(scissor) %>%
      summarise(across("Her2", mean, na.rm = TRUE)) %>%
      pivot_longer(
        cols = "Her2",
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
  mean_scores_scpas <- if ("scPAS" %in% colnames(es.df)) {
    es.df %>%
      group_by(scPAS) %>%
      summarise(across("Her2", mean, na.rm = TRUE)) %>%
      pivot_longer(
        cols = "Her2",
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
  mean_scores_scab <- if ("scAB" %in% colnames(es.df)) {
    es.df %>%
      group_by(scAB) %>%
      summarise(across("Her2", mean, na.rm = TRUE)) %>%
      pivot_longer(
        cols = "Her2",
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
  mean_scores_scpp <- if ("scPP" %in% colnames(es.df)) {
    es.df %>%
      group_by(scPP) %>%
      summarise(across("Her2", mean, na.rm = TRUE)) %>%
      pivot_longer(
        cols = "Her2",
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

  mean_scores <- rbind(
    mean_scores_scissor,
    mean_scores_scpas,
    mean_scores_scab,
    mean_scores_scpp
  ) %>%
    unite("screen", screen_method, screen_ressult, sep = "_") %>%
    mutate("Sample" = sample)
}

mean_scores_GSE42568 <- MeanScore(seurat, es.mat, sample = "GSE42568")
mean_scores_GSE162228 <- MeanScore(seurat2, es.mat, sample = "GSE162228")

mean_scores = grep("mean_scores_", ls(), value = T) %>%
  purrr::map_df(
    ~ {
      get(.x)
    }
  ) %>%
  rbind()

p <- ggplot(mean_scores, aes(x = Sample, y = screen, fill = mean_score)) +
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
    subtitle = "sc data-GSE161528 HER2, binary phenotype",
    y = "Mean ssGSEA Score",
    x = "Bulk Sample"
  ) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(
  file.path(plot_dir, "her2_mean_scores_binary.png"),
  plot = p,
  width = 6,
  height = 7,
  dpi = 400
)
