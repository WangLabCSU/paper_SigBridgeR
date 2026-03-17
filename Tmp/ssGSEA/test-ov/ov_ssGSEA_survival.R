# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd("/home/yyx/R/Project/R_code/SigBridgeR/Tmp/ssGSEA/test-ov")

set.seed(123)

library(Seurat)
library(GSVA)
library(ComplexHeatmap)
# library(clusterProfiler)
library(dplyr)
library(ggplot2)

# ! hgsoc-GSE165897, bulk-GSE140082

seurat <- qs::qread(
  "/home/data/sigbridger/benchmark_data/ov/GSE165897/GSE140082_ov_merged_seurat.qs",
  nthreads = 4
)
# seurat <- SCTransform(seurat)
# DefaultAssay(seurat) <- "SCT" # 推荐用 SCT 或 logNorm 计数

# 2. 计算各 cluster 的 marker（对比其余所有细胞）
# ! from CellMarker2.0, ovary cancer cell
markers <- data.table::fread("../cellmarker2/ovary_cancer.csv")
markers <- markers %>% arrange(desc(count))
markers2 <- data.table::fread("../cell_taxonomy/ovary_cancer.csv")

markers_all <- c(markers$marker, markers2$ovary_cancer) |> unique()

genesets = list("ovary_cancer" = markers_all)


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
  geneSets = genesets,
  # kcdf = auto # *自动选择了gaussian
)

set.seed(123)
# ! 决速步
es.mat <- gsva(ssgsea_param) # 行为 cluster-通路，列为细胞

if (!dir.exists('ssGSEA/test-ov')) {
  dir.create('ssGSEA/test-ov', recursive = TRUE)
}

qs::qsave(es.mat, "ov_ssgsea_score.qs") # PID=1933839
es.mat = qs::qread("ov_ssgsea_score.qs")

# ! 这里应该scPP运行时有bug，导致Negative和Neutral的细胞数都是0, 建议直接读取原始Screen文件
# scpp_result <- qs::qread(
#   "/home/data/sigbridger/benchmark_data/ov/GSE165897/GSE165897_GSE140082_scpp_result.qs"
# )
# boxplot(t(es.mat), las = 2)

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
es.df$scPP <- factor(
  es.df$scPP,
  levels = c("Neutral", "Negative", "Positive")
)
es.df$DEGAS <- factor(
  es.df$DEGAS,
  levels = c("Other", "Positive")
)
es.df$scPP <- factor(
  es.df$LP_SGL,
  levels = c("Neutral", "Negative", "Positive")
)


# # # 7. limma 差异检验（以 cluster0 为例）
# library(limma)
# scissor_design <- model.matrix(~scissor, data = es.df)
# scpas_design <- model.matrix(~scPAS, data = es.df)
# scab_design <- model.matrix(~scAB, data = es.df)
# scpp_design <- model.matrix(~scPP, data = es.df)

# scissor_fit <- lmFit(es.mat, scissor_design)
# scpas_fit <- lmFit(es.mat, scpas_design)
# scab_fit <- lmFit(es.mat, scab_design)
# scpp_fit <- lmFit(es.mat, scpp_design)

# scissor_fit <- eBayes(scissor_fit)
# scpas_fit <- eBayes(scpas_fit)
# scab_fit <- eBayes(scab_fit)
# scpp_fit <- eBayes(scpp_fit)

# scissor_toptable_pos = topTable(
#   scissor_fit,
#   coef = "scissorPositive",
#   number = 100
# ) # 比较Positive vs Neutral
# scpas_toptable_pos = topTable(
#   scpas_fit,
#   coef = "scPASPositive",
#   number = 100
# ) # 比较Positive vs Neutral
# scab_toptable_pos = topTable(
#   scab_fit,
#   coef = "scABPositive",
#   number = 100
# ) # 比较Positive vs Other
# scpp_toptable_pos = topTable(
#   scpp_fit,
#   coef = "scPPPositive",
#   number = 100
# ) # 比较Positive vs Negative

# # 或者
# toptable_neg = topTable(fit, coef = "scissorNegative", number = 100) # 比较Negative vs Neutral

plot_dir = "/home/yyx/R/Project/R_code/SigBridgeR/Tmp/ssGSEA/test-ov/plot"
if (!dir.exists(plot_dir)) {
  dir.create(plot_dir, recursive = TRUE)
}

order <- order(
  es.df$seurat_clusters,
  es.df$ovary_cancer,
  es.df$cell_type,
  es.df$scissor,
  es.df$scPAS,
  es.df$scAB,
  es.df$DEGAS,
  es.df$scPP,
  es.df$LP_SGL,
  es.df$cell_subtype
)
es.mat_ordered <- es.mat[, order, drop = FALSE]


# ! screen

# # 列注释（细胞分组）
col_anno_screen <- HeatmapAnnotation(
  scissor = es.df$scissor[order],
  scPAS = es.df$scPAS[order],
  scAB = es.df$scAB[order],
  scPP = es.df$scPP[order],
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
    ),
    scPP = c(
      "Neutral" = "#CECECE",
      "Positive" = "#ff3333",
      "Negative" = "#386c9b"
    )
  ),
  height = unit(1, "cm"),
  gap = unit(0.6, "mm")
)

# ! cell type + subtype
set.seed(123)
pallete_subtype <- randomcoloR::distinctColorPalette(
  length(unique(es.df$cell_subtype)),
  runTsne = T
)
col_anno_types <- HeatmapAnnotation(
  cell_type = es.df$cell_type[order],
  cell_subtype = es.df$cell_subtype[order],
  annotation_label = list(
    cell_type = "Cell Type",
    cell_subtype = "Cell Subtype"
  ),
  annotation_name_gp = gpar(fontface = "bold"),
  col = list(
    cell_type = setNames(pallete_subtype[1:3], unique(es.df$cell_type)),
    cell_subtype = setNames(pallete_subtype, unique(es.df$cell_subtype))
  ),
  height = unit(1, "cm"),
  gap = unit(0.6, "mm")
)

# ! distribution

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

# htmap <- ComplexHeatmap::Heatmap(
#     es.mat_ordered,
#     name = "ssGSEA\nScore",
#     top_annotation = c(
#         col_anno_types,
#         col_anno_screen,
#         col_anno_distribution
#     ),
#     show_column_names = FALSE,
#     show_row_names = FALSE,
#     cluster_columns = FALSE,
#     cluster_rows = FALSE,
#     use_raster = TRUE
# )

# png(
#     filename = file.path(plot_dir, "GSE140082_basic_heatmap.png"),
#     width = 1080
# )
# draw(htmap)
# dev.off()

# # ! 可视化-平均得分

library(ggplot2)
library(tidyr)
# 计算各组的平均得分
mean_scores_scissor <- es.df %>%
  group_by(scissor) %>%
  summarise(across("ovary_cancer", mean, na.rm = TRUE)) %>%
  pivot_longer(
    cols = "ovary_cancer",
    names_to = "cluster",
    values_to = "mean_score"
  ) %>%
  pivot_longer(
    cols = "scissor",
    names_to = "screen_method",
    values_to = "screen_ressult"
  )
mean_scores_scpas <- es.df %>%
  group_by(scPAS) %>%
  summarise(across("ovary_cancer", mean, na.rm = TRUE)) %>%
  pivot_longer(
    cols = "ovary_cancer",
    names_to = "cluster",
    values_to = "mean_score"
  ) %>%
  pivot_longer(
    cols = "scPAS",
    names_to = "screen_method",
    values_to = "screen_ressult"
  )
mean_scores_scab <- es.df %>%
  group_by(scAB) %>%
  summarise(across("ovary_cancer", mean, na.rm = TRUE)) %>%
  pivot_longer(
    cols = "ovary_cancer",
    names_to = "cluster",
    values_to = "mean_score"
  ) %>%
  pivot_longer(
    cols = "scAB",
    names_to = "screen_method",
    values_to = "screen_ressult"
  )
mean_scores_scpp <- es.df %>%
  group_by(scPP) %>%
  summarise(across("ovary_cancer", mean, na.rm = TRUE)) %>%
  pivot_longer(
    cols = "ovary_cancer",
    names_to = "cluster",
    values_to = "mean_score"
  ) %>%
  pivot_longer(
    cols = "scPP",
    names_to = "screen_method",
    values_to = "screen_ressult"
  )

mean_scores <- rbind(
  mean_scores_scissor,
  mean_scores_scpas,
  mean_scores_scab,
  mean_scores_scpp
) %>%
  unite("screen", screen_method, screen_ressult, sep = "_")
mean_scores$bulk_sample <- "GSE140082"
mean_scores$sc_sample <- "GSE165897"

plot <- ggplot(mean_scores, aes(x = cluster, y = screen, fill = mean_score)) +
  geom_point(size = 6, alpha = 0.9, shape = 21, color = "black") +
  scale_fill_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    name = "Mean Score"
  ) +
  labs(
    title = "Mean ssGSEA Scores by Group",
    y = "Mean ssGSEA Score",
    x = "Cluster"
  ) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(
  file.path(plot_dir, "GSE140082_mean_scores.png"),
  plot = plot,
  width = 4,
  height = 10,
  dpi = 400
)

# -------------------------------------------------------------------------

# ! 另外一个Bulk-GSE32062，只有筛选不同

seurat2 <- qs::qread(
  "/home/data/sigbridger/benchmark_data/ov/GSE165897/GSE32062_soc_merged_seurat.qs",
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

meta2 <- seurat2@meta.data
meta2$cell <- colnames(seurat2)
es.df2 <- t(es.mat) %>% # 行为细胞，列为 cluster
  as.data.frame() %>%
  tibble::rownames_to_column("cell") %>%
  left_join(meta2, by = "cell")

es.df2$scissor <- factor(
  es.df2$scissor,
  levels = c("Neutral", "Negative", "Positive")
)
es.df2$scPAS <- factor(
  es.df2$scPAS,
  levels = c("Neutral", "Negative", "Positive")
)
es.df2$scAB <- factor(
  es.df2$scAB,
  levels = c("Other", "Positive")
)

table(es.df2$scissor)
#  Neutral Negative Positive
#    47670     1598     1597

table(es.df2$scPAS)
#  Neutral Negative Positive
#    45270     3933     1662

table(es.df2$scAB)
#    Other Positive
#    47389     3476

# # 7. limma 差异检验（以 cluster0 为例）
library(limma)
scissor_design2 <- model.matrix(~scissor, data = es.df2)
scpas_design2 <- model.matrix(~scPAS, data = es.df2)
scab_design2 <- model.matrix(~scAB, data = es.df2)


scissor_fit2 <- lmFit(es.mat, scissor_design2)
scpas_fit2 <- lmFit(es.mat, scpas_design2)
scab_fit2 <- lmFit(es.mat, scab_design2)

scissor_fit2 <- eBayes(scissor_fit2)
scpas_fit2 <- eBayes(scpas_fit2)
scab_fit2 <- eBayes(scab_fit2)

scissor_toptable_pos2 = topTable(
  scissor_fit2,
  coef = "scissorPositive",
  number = 100
) # 比较Positive vs Neutral
scpas_toptable_pos2 = topTable(
  scpas_fit2,
  coef = "scPASPositive",
  number = 100
) # 比较Positive vs Neutral
scab_toptable_pos2 = topTable(
  scab_fit2,
  coef = "scABPositive",
  number = 100
) # 比较Positive vs Other

es.mat_ordered <- es.mat[, order, drop = FALSE]

# # 列注释（细胞分组）
col_anno_screen2 <- CreateHtmapAnno(
  seurat2,
  es.mat,
  color = "#cf7815ff"
)

# htmap2 <- Heatmap(
#     es.mat_ordered,
#     name = "ssGSEA\nScore",
#     top_annotation = col_anno_screen2,
#     show_column_names = FALSE,
#     show_row_names = FALSE,
#     cluster_columns = FALSE,
#     use_raster = TRUE
# )

# png(
#     filename = file.path(plot_dir, "GSE140082_basic_heatmap.png"),
#     width = 1080
# )
# draw(htmap)
# dev.off()

# # ! 可视化-平均得分

library(ggplot2)
library(tidyr)
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
      summarise(across("ovary_cancer", mean, na.rm = TRUE)) %>%
      pivot_longer(
        cols = "ovary_cancer",
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
      summarise(across("ovary_cancer", mean, na.rm = TRUE)) %>%
      pivot_longer(
        cols = "ovary_cancer",
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
      summarise(across("ovary_cancer", mean, na.rm = TRUE)) %>%
      pivot_longer(
        cols = "ovary_cancer",
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
      summarise(across("ovary_cancer", mean, na.rm = TRUE)) %>%
      pivot_longer(
        cols = "ovary_cancer",
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


# 计算各组的平均得分
mean_scores_scissor2 <- es.df2 %>%
  group_by(scissor) %>%
  summarise(across("ovary_cancer", mean, na.rm = TRUE)) %>%
  pivot_longer(
    cols = "ovary_cancer",
    names_to = "cluster",
    values_to = "mean_score"
  ) %>%
  pivot_longer(
    cols = "scissor",
    names_to = "screen_method",
    values_to = "screen_ressult"
  )
mean_scores_scpas2 <- es.df2 %>%
  group_by(scPAS) %>%
  summarise(across("ovary_cancer", mean, na.rm = TRUE)) %>%
  pivot_longer(
    cols = "ovary_cancer",
    names_to = "cluster",
    values_to = "mean_score"
  ) %>%
  pivot_longer(
    cols = "scPAS",
    names_to = "screen_method",
    values_to = "screen_ressult"
  )
mean_scores_scab2 <- es.df2 %>%
  group_by(scAB) %>%
  summarise(across("ovary_cancer", mean, na.rm = TRUE)) %>%
  pivot_longer(
    cols = "ovary_cancer",
    names_to = "cluster",
    values_to = "mean_score"
  ) %>%
  pivot_longer(
    cols = "scAB",
    names_to = "screen_method",
    values_to = "screen_ressult"
  )

mean_scores2 <- rbind(
  mean_scores_scissor2,
  mean_scores_scpas2,
  mean_scores_scab2
) %>%
  unite("screen", screen_method, screen_ressult, sep = "_")
mean_scores2$bulk_sample <- "GSE32062"
mean_scores2$sc_sample <- "GSE165897"

plot <- ggplot(mean_scores2, aes(x = cluster, y = screen, fill = mean_score)) +
  geom_point(size = 6, alpha = 0.9, shape = 21, color = "black") +
  scale_fill_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    name = "Mean Score"
  ) +
  labs(
    title = "Mean ssGSEA Scores by Group",
    y = "Mean ssGSEA Score",
    x = "Cluster"
  ) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(
  file.path(plot_dir, "GSE140082_mean_scores.png"),
  plot = plot,
  width = 4,
  height = 10,
  dpi = 400
)


# -----------------------------------------------------------------
# ! 拼图
# ? heatmap
bulk_lgd = Legend(
  at = c("GSE140082", "GSE32062"),
  title = "Bulk Sample for Methods",
  legend_gp = gpar(fill = c("#00913fff", "#cf7815ff")),
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
  height = unit(0.6, "mm")
)

# ! 依据cluster 分组
column_split = es.df$seurat_clusters |> sort()

# # ! mean score of each cluster ,按照分数大小排序了，abandoned
cluster_mean_score = purrr::map_dbl(
  sort(unique(es.df$seurat_clusters)),
  function(x) {
    mean(es.df$ovary_cancer[es.df$seurat_clusters == x])
  }
)
cluster_mean_score = round(cluster_mean_score, 2)

# ! treatment phase

pallete_tp <- randomcoloR::distinctColorPalette(
  length(unique(es.df$treatment_phase))
)
col_anno_treatment <- HeatmapAnnotation(
  treatment_phase = es.df$treatment_phase[order],
  annotation_label = "Treatment Phase",
  annotation_name_gp = gpar(fontface = "bold"),
  col = list(
    treatment_phase = setNames(pallete_tp, unique(es.df$treatment_phase))
  ),
  height = unit(1, "cm"),
  gap = unit(0.6, "mm")
)

# ! anatomical location

pallete_al <- randomcoloR::distinctColorPalette(
  length(unique(es.df$anatomical_location))
)
col_anno_al <- HeatmapAnnotation(
  anatomical_location = es.df$anatomical_location[order],
  annotation_label = "Anatomical Location",
  annotation_name_gp = gpar(fontface = "bold"),
  col = list(
    anatomical_location = setNames(
      pallete_al,
      unique(es.df$anatomical_location)
    )
  ),
  height = unit(1, "cm"),
  gap = unit(0.6, "mm")
)

# ! seurat cluster

pallete_cluster <- randomcoloR::distinctColorPalette(
  length(unique(es.df$seurat_clusters)),
  runTsne = T
)
col_anno_cluster = HeatmapAnnotation(
  Cluster = es.df$seurat_clusters[order],
  annotation_label = "Cluster",
  annotation_name_gp = gpar(fontface = "bold"),
  col = list(
    Cluster = setNames(
      pallete_cluster,
      unique(es.df$seurat_clusters)
    )
  ),
  height = unit(1, "cm"),
  gap = unit(0.6, "mm")
)

# ! 拼图

htmap_cmb = ComplexHeatmap::Heatmap(
  es.mat_ordered,
  name = "ssGSEA\nScore",
  top_annotation = c(
    col_anno_cluster,
    row_anno_gap,
    col_anno_types,
    row_anno_gap,
    col_anno_treatment,
    row_anno_gap,
    col_anno_al,
    row_anno_gap,
    col_anno_screen,
    row_anno_gap,
    col_anno_screen2,
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
  row_title = "Survival as phenotype"
)
# draw(
#     htmap_cmb,
#     annotation_legend_list = list(bulk_lgd, screen_lgd)
# )

ragg::agg_png(
  file.path(plot_dir, "ov_ssGSEA_heatmap_survival.png"),
  width = 1600,
  height = 1000
)
draw(
  htmap_cmb,
  annotation_legend_list = list(bulk_lgd, screen_lgd)
)
dev.off()

# ! 拼图
# ? 平均得分
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
      summarise(across("ovary_cancer", mean, na.rm = TRUE)) %>%
      pivot_longer(
        cols = "ovary_cancer",
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
      summarise(across("ovary_cancer", mean, na.rm = TRUE)) %>%
      pivot_longer(
        cols = "ovary_cancer",
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
      summarise(across("ovary_cancer", mean, na.rm = TRUE)) %>%
      pivot_longer(
        cols = "ovary_cancer",
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
      summarise(across("ovary_cancer", mean, na.rm = TRUE)) %>%
      pivot_longer(
        cols = "ovary_cancer",
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
      summarise(across("ovary_cancer", mean, na.rm = TRUE)) %>%
      pivot_longer(
        cols = "ovary_cancer",
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
      summarise(across("ovary_cancer", mean, na.rm = TRUE)) %>%
      pivot_longer(
        cols = "ovary_cancer",
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
      summarise(across("ovary_cancer", mean, na.rm = TRUE)) %>%
      pivot_longer(
        cols = "ovary_cancer",
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
mean_score_GSE140082 = MeanScore(seurat, es.mat, "GSE140082")
mean_score_GSE32062 = MeanScore(seurat2, es.mat, "GSE32062")


mean_scores = grep("mean_score_", ls(), value = T) %>%
  purrr::map_df(
    ~ {
      get(.x)
    }
  ) %>%
  rbind()

library(ggplot2)
p <- ggplot(
  mean_scores,
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
    title = "Mean ssGSEA Scores by Group",
    subtitle = "Ovary cancer, sc-GSE165897, survival as phenotype",
    y = "Mean ssGSEA Score",
    x = "Under different bulk sample"
  ) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(
  file.path(plot_dir, "GSE140082_GSE32062_mean_scores_survival.png"),
  width = 5,
  height = 7,
  plot = p,
  dpi = 400
)
