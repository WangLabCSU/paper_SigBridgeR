# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd("/home/yyx/R/Project/R_code/SigBridgeR/Tmp/ssGSEA/test-lung")

set.seed(123)

library(Seurat)
library(GSVA)
library(ComplexHeatmap)
# library(clusterProfiler)
library(dplyr)
library(ggplot2)

# ! luad-sc-123902, bulk-GSE3141

seurat <- qs::qread(
  "/home/data/sigbridger/benchmark_data/lung/GSE3141/GSE3141_luad_merged_seurat.qs",
  nthreads = 4
)
# seurat <- SCTransform(seurat)
# DefaultAssay(seurat) <- "SCT" # 推荐用 SCT 或 logNorm 计数

# 2. 计算各 cluster 的 marker（对比其余所有细胞）
# ! from CellMarker2.0, ovary cancer cell
markers <- data.table::fread("../cellmarker2/lung_cacner.csv")
# ! FROM cell taxonomy https://ngdc.cncb.ac.cn/celltaxonomy/celltype/CT:00001067
markers2 <- data.table::fread("../cell_taxonomy/lung_cancer.csv")
marker_all = c(markers$marker, markers2$lung_cancer)
genesets = list("lung_cancer" = unique(marker_all))


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
  geneSets = genesets
  # kcdf = auto # *自动选择了gaussian
)

set.seed(123)
# ! 决速步
es.mat <- gsva(ssgsea_param) # 行为 cluster-通路，列为细胞

plot_dir = "/home/yyx/R/Project/R_code/SigBridgeR/Tmp/ssGSEA/test-lung/plot"
if (!dir.exists('ssGSEA/test-lung')) {
  dir.create('ssGSEA/test-lung', recursive = TRUE)
}

qs::qsave(es.mat, "luad_ssgsea_score.qs") # PID=1233279


es.mat = qs::qread("luad_ssgsea_score.qs")

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

table(es.df$scissor)
# Neutral Negative Positive
#    31398      133      647

table(es.df$scPAS)
#  Neutral Negative Positive
#    31393      412      373

table(es.df$scAB)
#    Other Positive
#    29140     3038

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

plot_dir = "/home/yyx/R/Project/R_code/SigBridgeR/Tmp/ssGSEA/test-lung/plot"
if (!dir.exists(plot_dir)) {
  dir.create(plot_dir, recursive = TRUE)
}

# # ! 可视化-基础热图
set.seed(123)

order <- order(
  es.df$seurat_clusters,
  es.df$lung_cancer,
  es.df$Celltype,
  es.df$cnv_status,
  es.df$Sample
)
es.mat_ordered <- es.mat[, order, drop = FALSE]

# ! screen

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


# ---------------------------------------------------------

# ! 可视化-差异显著性的热图
# ! 只有一行，没有意义
# 提取limma结果中的显著差异cluster
# scissor_pos_sig_clusters <- scissor_toptable_pos %>%
#     filter(adj.P.Val < 0.05) %>%
#     rownames()
# scpas_pos_sig_clusters <- scpas_toptable_pos %>%
#     filter(adj.P.Val < 0.05) %>%
#     rownames()
# scab_pos_sig_clusters <- scab_toptable_pos %>%
#     filter(adj.P.Val < 0.05) %>%
#     rownames()

# if (length(pos_sig_clusters) > 0) {
#     sig_es.mat <- es.mat[pos_sig_clusters, scissor_order]

#     pos_sig_htmap = Heatmap(
#         sig_es.mat,
#         name = "ssGSEA\nScore",
#         top_annotation = col_anno,
#         show_column_names = FALSE,
#         cluster_rows = TRUE,
#         cluster_columns = FALSE,
#         row_names_gp = gpar(fontsize = 10),
#         column_title = "Significantly Different Clusters"
#     )

#     png(
#         filename = file.path(plot_dir, "pos_sig_heatmap.png"),
#         width = 1080
#     )
#     draw(pos_sig_htmap)
#     dev.off()
# }

# neg_sig_clusters <- toptable_neg %>%
#     filter(adj.P.Val < 0.05) %>%
#     rownames()

# if (length(neg_sig_clusters) > 0) {
#     sig_es.mat <- es.mat[neg_sig_clusters, scissor_order]

#     neg_sig_htmap = Heatmap(
#         sig_es.mat,
#         name = "ssGSEA\nScore",
#         top_annotation = col_anno,
#         show_column_names = FALSE,
#         cluster_rows = TRUE,
#         cluster_columns = FALSE,
#         row_names_gp = gpar(fontsize = 10),
#         column_title = "Significantly Different Clusters"
#     )

#     png(
#         filename = file.path(plot_dir, "neg_sig_heatmap.png"),
#         width = 1080
#     )
#     draw(neg_sig_htmap)
#     dev.off()
# }

# # ! 可视化-小提琴图
# ! 只有一行，没有意义

# library(tidyr)
# es.long <- es.df %>%
#     select(cell, scissor, all_of(scissor_pos_sig_clusters)) %>%
#     pivot_longer(
#         cols = all_of(sig_clusters),
#         names_to = "cluster",
#         values_to = "score"
#     )

# vio = ggplot(es.long, aes(x = scissor, y = score, fill = scissor)) +
#     geom_violin(alpha = 0.7) +
#     geom_boxplot(width = 0.2, alpha = 0.8) +
#     facet_wrap(~cluster, scales = "free_y") +
#     scale_fill_manual(
#         values = c("Neutral" = "grey", "Negative" = "blue", "Positive" = "red")
#     ) +
#     labs(
#         title = "ssGSEA Scores of Significant Clusters",
#         y = "ssGSEA Score",
#         x = "Scissor Group"
#     ) +
#     theme_bw()
# ggsave(
#     filename = file.path(plot_dir, "vio.png"),
#     plot = vio,
#     width = 10,
#     height = 10
# )

# # ! 可视化-在降维图上显示ssGSEA得分
# ! 只有一行，没有意义
# 提取UMAP坐标
# umap_df <- as.data.frame(seurat@reductions$umap@cell.embeddings)
# umap_df$cell <- rownames(umap_df)
# umap_df <- left_join(umap_df, es.df, by = "cell")

# 绘制每个显著cluster的UMAP
# plot_umap_cluster <- function(cluster_name) {
#     ggplot(
#         umap_df,
#         aes(x = umap_1, y = umap_2, color = .data[[cluster_name]])
#     ) +
#         geom_point(size = 0.5, alpha = 0.7) +
#         scale_color_gradient2(
#             low = "blue",
#             mid = "white",
#             high = "red",
#             midpoint = median(umap_df[[cluster_name]])
#         ) +
#         labs(title = paste("UMAP:", cluster_name), color = "ssGSEA\nScore") +
#         theme_classic() +
#         facet_wrap(~scissor) # 按组分别显示
# }

# for (cluster in sig_clusters) {
#     png(
#         filename = file.path(plot_dir, glue::glue("umap_{cluster}.png")),
#         width = 600,
#         height = 600
#     )
#     print(plot_umap_cluster(cluster))
#     dev.off()
# }

# # ! 可视化-相关性热图
# ! 只有一行，没有意义, 只能算出一个格子

# # 计算cluster间的相关性
# cluster_cor <- cor(t(es.mat))

# # 绘制相关性热图
# cor_htmap = Heatmap(
#     cluster_cor,
#     name = "Correlation",
#     col = circlize::colorRamp2(c(-1, 0, 1), c("blue", "white", "red")),
#     row_names_gp = gpar(fontsize = 8),
#     column_names_gp = gpar(fontsize = 8),
#     column_title = "Correlation between Clusters"
# )

# ---------------------------------------------------------
# ? 继续组合heatmap

seurat_GSE8894 = qs::qread(
  "/home/data/sigbridger/benchmark_data/lung/GSE8894/GSE8894_luad_merged_seurat.qs"
)
seurat_GSE31210 = qs::qread(
  "/home/data/sigbridger/benchmark_data/lung/GSE31210/GSE31210_luad_merged_seurat.qs"
)
seurat_tcga = qs::qread(
  "/home/data/sigbridger/benchmark_data/lung/TCGA-LUAD/tcga_luad_merged_seurat.qs"
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

col_anno_screen_GSE8894 <- CreateHtmapAnno(
  seurat_GSE8894,
  es.mat,
  color = "#cf7815ff"
)
col_anno_screen_GSE31210 <- CreateHtmapAnno(
  seurat_GSE31210,
  es.mat,
  color = "#bf009cff"
)
col_anno_screen_tcga <- CreateHtmapAnno(
  seurat_tcga,
  es.mat,
  color = "#00c7baff"
)

bulk_lgd = Legend(
  at = c("GSE3141", "GSE8894", "GSE31210", "TCGA"),
  title = "Bulk Sample for Methods",
  legend_gp = gpar(
    fill = c("#00913fff", "#cf7815ff", "#bf009cff", "#00c7baff")
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
    col_anno_screen_GSE8894,
    row_anno_gap,
    col_anno_screen_GSE31210,
    row_anno_gap,
    col_anno_screen_tcga,
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
  file.path(plot_dir, "lung_ssGSEA_heatmap_survival.png"),
  width = 1600,
  height = 1000
)
draw(
  htmap_cmb,
  annotation_legend_list = list(bulk_lgd, screen_lgd)
)
dev.off()


# ---------------------------------------------------------

# # ! 可视化-平均得分
library(tidyr)
# 计算各组的平均得分
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


mean_score_GSE3141 = MeanScore(seurat, es.mat, "GSE3141")
mean_score_GSE8894 = MeanScore(seurat_GSE8894, es.mat, "GSE8894")
mean_score_GSE31210 = MeanScore(seurat_GSE31210, es.mat, "GSE31210")
mean_score_tcga = MeanScore(seurat_tcga, es.mat, "TCGA")

mean_scores = grep("mean_score_", ls(), value = T) %>%
  purrr::map_df(
    ~ {
      get(.x)
    }
  ) %>%
  rbind()

library(ggplot2)
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
    subtitle = "sc-GSE123902, survival as phenotype",
    y = "Mean ssGSEA Score",
    x = "Bulk Sample"
  ) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(
  file.path(plot_dir, "lung_mean_scores_survival.png"),
  plot = p,
  width = 6,
  height = 7,
  dpi = 400
)
