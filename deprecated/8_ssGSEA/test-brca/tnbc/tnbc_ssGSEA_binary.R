library(data.table)
library(Seurat)
library(GSVA)
library(ComplexHeatmap)
library(dplyr)
library(tidyr)
library(ggplot2)

setwd("/home/yyx/R/Project/R_code/SigBridgeR/Tmp/ssGSEA/test-brca/tnbc")

data_path = "/home/yyx/R/Project/R_code/SigBridgeR/Tmp/ssGSEA/BRCA marker"

set.seed(123)


# ! brca-GSE161529_tnbc, bulk-GSE42569, binary phenotype

seurat <- qs::qread(
  "/home/data/sigbridger/benchmark_binary/brca/TNBC/GSE42568_tnbc_merged_seurat.qs",
  nthreads = 4
)
# seurat <- SCTransform(seurat)
# DefaultAssay(seurat) <- "SCT" # 推荐用 SCT 或 logNorm 计数

# 2. 计算各 cluster 的 marker（对比其余所有细胞）
# ! from CellMarker2.0, ovary cancer cell
markers <- data.table::fread(file.path(data_path, "markers.csv"))
markers = list("TNBC" = markers$TNBC_marker)


# 3. 每个 cluster 取 top 100 基因（可改）
# Idents(seurat) <- "scissor"

# markers <- FindAllMarkers(
#     seurat,
#     only.pos = TRUE, # 只要上调
#     min.pct = 0.2,
#     logfc.threshold = 0.2
# )
# genesets <- markers %>%
#     group_by(cluster) %>%
#     top_n(n = 100, wt = avg_log2FC)

# # %>% # 按 logFC 排序
# genesets = split(x = genesets$gene, f = genesets$cluster) # 变成 list，名 = cluster

# 4. 提取单细胞表达矩阵（基因 × 细胞）
expr <- as.matrix(SeuratObject::LayerData(seurat, layer = "data"))


# 5. 构造 GSVA 参数对象（≥1.50 版写法）
ssgsea_param <- gsvaParam(
  exprData = expr,
  geneSets = markers,
  # kcdf = auto # *自动选择了gaussian
)

set.seed(123)
# ! 决速步
es.mat <- gsva(ssgsea_param) # 行为 cluster-通路，列为细胞

qs::qsave(es.mat, "tnbc_ssgsea_score.qs") # PID=960788
es.mat <- qs::qread("tnbc_ssgsea_score.qs")

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
# es.df$scPP <- factor(
#     es.df$scPP,
#     levels = c("Neutral", "Negative", "Positive")
# )

table(es.df$scissor)
#  Neutral Negative Positive
#    42529     6652     5469
table(es.df$scPAS)
#  Neutral Negative Positive
#    54164      486        0

table(es.df$scAB)
#    Other Positive
#    22899    31751
# table(es.df$scPP)

# 7. limma 差异检验（以 cluster0 为例）
library(limma)
scissor_design <- model.matrix(~scissor, data = es.df)
scpas_design <- model.matrix(~scPAS, data = es.df)
scab_design <- model.matrix(~scAB, data = es.df)
scpp_design <- model.matrix(~scPP, data = es.df)

scissor_fit <- lmFit(es.mat, scissor_design)
scpas_fit <- lmFit(es.mat, scpas_design)
scab_fit <- lmFit(es.mat, scab_design)
scpp_fit <- lmFit(es.mat, scpp_design)

scissor_fit <- eBayes(scissor_fit)
scpas_fit <- eBayes(scpas_fit)
scab_fit <- eBayes(scab_fit)
scpp_fit <- eBayes(scpp_fit)

scissor_toptable_pos = topTable(
  scissor_fit,
  coef = "scissorPositive",
  number = 100
) # 比较Positive vs Neutral
scpas_toptable_pos = topTable(scpas_fit, coef = "scPASPositive", number = 100) # 比较Positive vs Neutral
scab_toptable_pos = topTable(scab_fit, coef = "scABPositive", number = 100) # 比较Positive vs Other
scpp_toptable_pos = topTable(scpp_fit, coef = "scPPPositive", number = 100) # 比较Positive vs Neutral
# # 或者
# toptable_neg = topTable(fit, coef = "scissorNegative", number = 100) # 比较Negative vs Neutral

plot_dir = "/home/yyx/R/Project/R_code/SigBridgeR/Tmp/ssGSEA/test-brca/tnbc/binary_plot"
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
  es.df$TNBC,
  es.df$scissor,
  es.df$scPAS,
  es.df$scAB
  # ,    es.df$scPP
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
    size = unit(0.5, "mm"),
    ylim = c(-1, 1),
    axis_param = list(
      at = c(-0.9, -0.6, -0.3, 0, 0.3, 0.6, 0.9)
    ),
    gp = gpar(col = "#9f9f9fff", alpha = 0.5),
  ),
  gap = unit(0.6, "mm")
)

# ! bulk-GSE162228
seurat2 <- qs::qread(
  "/home/data/sigbridger/benchmark_binary/brca/TNBC/GSE162228_tnbc_merged_seurat.qs",
  nthreads = 4
)

col_anno_GSE162228 <- CreateHtmapAnno(
  seurat2,
  es.mat = es.mat_ordered,
  color = "#cf7815ff"
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

set.seed(334)
pallete_tnbc <- randomcoloR::distinctColorPalette(
  length(unique(es.df$group))
  # ,    runTsne = T
)
col_anno_tnbc <- HeatmapAnnotation(
  TNBC_group = es.df$group[order],
  annotation_name_gp = gpar(fontface = "bold"),
  annotation_label = "TNBC Group",
  col = list(
    HER2_group = setNames(pallete_tnbc, unique(es.df$group))
  ),
  height = unit(1, "cm"),
  gap = unit(0.6, "mm")
)


# ! 间隔
row_anno_gap = HeatmapAnnotation(
  empty = anno_empty(border = FALSE), # 让热图扁一点
  height = unit(0.6, "mm")
)

# ! 依据cluster 分组
column_split = es.df$seurat_clusters |> sort()

# # ! mean score of each cluster ,按照分数大小排序了，abandoned
cluster_mean_score = purrr::map_dbl(
  sort(unique(es.df$seurat_clusters)),
  function(x) {
    mean(es.df$TNBC[es.df$seurat_clusters == x])
  }
)
cluster_mean_score = round(cluster_mean_score, 2)

htmap_cmb = ComplexHeatmap::Heatmap(
  es.mat_ordered,
  name = "ssGSEA\nScore",
  top_annotation = c(
    col_anno_cluster,
    row_anno_gap,
    col_anno_tnbc,
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
  filename = file.path(plot_dir, "ssGSEA_binary_tnbc_heatmap.png"),
  width = 1600,
  height = 1000
)
draw(
  htmap_cmb,
  annotation_legend_list = list(bulk_lgd, screen_lgd)
)
dev.off()
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
      summarise(across("TNBC", mean, na.rm = TRUE)) %>%
      pivot_longer(
        cols = "TNBC",
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
      summarise(across("TNBC", mean, na.rm = TRUE)) %>%
      pivot_longer(
        cols = "TNBC",
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
      summarise(across("TNBC", mean, na.rm = TRUE)) %>%
      pivot_longer(
        cols = "TNBC",
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
      summarise(across("TNBC", mean, na.rm = TRUE)) %>%
      pivot_longer(
        cols = "TNBC",
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
    subtitle = "sc-GSE161528 TNBC, binary phenotype",
    y = "Mean ssGSEA Score",
    x = "Bulk Sample"
  ) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(
  file.path(plot_dir, "TNBC_mean_scores_binary.png"),
  plot = p,
  width = 6,
  height = 7,
  dpi = 400
)
