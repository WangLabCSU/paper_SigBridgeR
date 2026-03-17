# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd("/home/yyx/R/Project/R_code/SigBridgeR/Tmp")

set.seed(123)

library(Seurat)
library(GSVA)
library(ComplexHeatmap)
# library(clusterProfiler)
library(dplyr)
library(ggplot2)

# ! brca-her2, bulk-GSE42568, sc-GSE161529, pheno=survival, 4 methods

seurat <- qs::qread(
  "/home/yyx/R/Project/R_code/SigBridgeR/Tmp/benchmark_data/brca/HER2/GSE42568_her2_merged_seurat.qs",
  nthreads = 4
)
# seurat <- SCTransform(seurat)
# DefaultAssay(seurat) <- "SCT" # 推荐用 SCT 或 logNorm 计数

# 2. 计算各 cluster 的 marker（对比其余所有细胞）
# ! from CellMarker2.0, her2 breast cancer cell
markers <- data.table::fread("ssGSEA/cellmarker2/her2_breast_cancer.csv")
genesets = list("HER2_Breast" = markers$marker)


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

if (dir.exists('ssGSEA/test-ov')) {
  dir.create('ssGSEA/test-ov', recursive = TRUE)
}

qs::qsave(es.mat, "ssGSEA/test-ov/cluster_ssgsea_score.qs") # PID=1507167
# es.mat = qs::qread("ssGSEA/test-ov/cluster_ssgsea_score.qs")

# # boxplot(t(es.mat), las = 2)

# meta <- seurat@meta.data
# meta$cell <- colnames(seurat)
# es.df <- t(es.mat) %>% # 行为细胞，列为 cluster
#     as.data.frame() %>%
#     tibble::rownames_to_column("cell") %>%
#     left_join(meta, by = "cell")

# es.df$scissor <- factor(
#     es.df$scissor,
#     levels = c("Neutral", "Negative", "Positive")
# )

# # table(es.df$scissor)

# # 7. limma 差异检验（以 cluster0 为例）
# library(limma)
# design <- model.matrix(~scissor, data = es.df)

# # head(design)
# #   (Intercept) scissorNegative scissorPositive
# # 1           1               0               0
# # 2           1               0               0
# # 3           1               0               1
# # 4           1               0               1
# # 5           1               0               1
# # 6           1               0               1

# fit <- lmFit(es.mat, design)
# fit <- eBayes(fit)

# toptable_pos = topTable(fit, coef = "scissorPositive", number = 100) # 比较Positive vs Neutral
# # 或者
# toptable_neg = topTable(fit, coef = "scissorNegative", number = 100) # 比较Negative vs Neutral

# plot_dir = "/home/yyx/R/Project/R_code/SigBridgeR/Tmp/ssGSEA/test-ov/plot"
# # ! 可视化-基础热图
# scissor_order <- order(es.df$scissor)
# es.mat_ordered <- es.mat[, scissor_order]

# # 列注释（细胞分组）
# col_anno <- HeatmapAnnotation(
#     scissor = es.df$scissor[scissor_order],
#     col = list(
#         scissor = c(
#             "Neutral" = "#CECECE",
#             "Other" = "#CECECE",
#             "Positive" = "#ff3333",
#             "Negative" = "#386c9b"
#         )
#     )
# )
# htmap <- Heatmap(
#     es.mat_ordered,
#     name = "ssGSEA\nScore",
#     top_annotation = col_anno,
#     show_column_names = FALSE,
#     cluster_columns = FALSE,
#     column_split = es.df$scissor[scissor_order],
#     row_names_gp = gpar(fontsize = 8)
# )
# png(
#     filename = file.path(plot_dir, "basic_heatmap.png"),
#     width = 1080
# )
# draw(htmap)
# dev.off()

# # ! 可视化-差异显著性的热图
# # 提取limma结果中的显著差异cluster
# pos_sig_clusters <- toptable_pos %>%
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

# library(tidyr)
# es.long <- es.df %>%
#     select(cell, scissor, all_of(sig_clusters)) %>%
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
# # 提取UMAP坐标
# umap_df <- as.data.frame(seurat@reductions$umap@cell.embeddings)
# umap_df$cell <- rownames(umap_df)
# umap_df <- left_join(umap_df, es.df, by = "cell")

# # 绘制每个显著cluster的UMAP
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

# # ! 可视化-平均得分

# # 计算各组的平均得分
# mean_scores <- es.df %>%
#     group_by(scissor) %>%
#     summarise(across(starts_with("cluster"), mean, na.rm = TRUE)) %>%
#     pivot_longer(
#         cols = starts_with("cluster"),
#         names_to = "cluster",
#         values_to = "mean_score"
#     )

# # 只显示显著cluster
# mean_scores_sig <- mean_scores %>%
#     filter(cluster %in% sig_clusters)

# ggplot(mean_scores_sig, aes(x = cluster, y = mean_score, fill = scissor)) +
#     geom_bar(stat = "identity", position = "dodge") +
#     scale_fill_manual(
#         values = c("Neutral" = "grey", "Negative" = "blue", "Positive" = "red")
#     ) +
#     labs(
#         title = "Mean ssGSEA Scores by Group",
#         y = "Mean ssGSEA Score",
#         x = "Cluster"
#     ) +
#     theme_classic() +
#     theme(axis.text.x = element_text(angle = 45, hjust = 1))
