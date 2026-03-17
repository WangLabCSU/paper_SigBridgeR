library(data.table)
library(Seurat)
library(GSVA)
library(ComplexHeatmap)
library(dplyr)
library(tidyr)
library(ggplot2)

setwd(
  "/home/yyx/R/Project/R_code/SigBridgeR/Tmp/ssGSEA/test-brca/only_her2"
)

data_path = "/home/yyx/R/Project/R_code/SigBridgeR/Tmp/ssGSEA/BRCA marker"

set.seed(123)


# ! brca-GSE161529_her2, bulk - GSE42568, pheno-survival

seurat <- qs::qread(
  "/home/data/sigbridger/benchmark_data/brca/HER2/GSE42568_her2_merged_seurat.qs",
  nthreads = 4
)
seurat2 <- qs::qread(
  "/home/data/sigbridger/benchmark_data/brca/HER2/GSE162228_her2_merged_seurat.qs",
  nthreads = 4
)

data = SeuratObject::LayerData(seurat, assay = "RNA")
genes = rownames(data)
"ERBB2" %chin% genes
# seurat <- SCTransform(seurat)
# DefaultAssay(seurat) <- "SCT" # 推荐用 SCT 或 logNorm 计数

# 2. 计算各 cluster 的 marker（对比其余所有细胞）
# ! from CellMarker2.0, her2+ brca
markers = list("Her2" = "ERBB2")

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

qs::qsave(es.mat, "only_her2_ssgsea_score.qs") # PID=965591


# ----------------------------------------------------------------------
# ! 可视化-平均得分
MeanScore = function(seurat_obj, es.mat, sample = "GSE") {
  meta <- seurat_obj@meta.data %>%
    tibble::rownames_to_column("cell")

  # 处理表达矩阵并合并元数据
  es.df <- Matrix::t(es.mat) %>%
    as.data.frame() %>%
    tibble::rownames_to_column("cell") %>%
    dplyr::left_join(meta, by = "cell")
  mean_scores_scissor <- if ("scissor" %in% colnames(es.df)) {
    es.df %>%
      dplyr::group_by(scissor) %>%
      dplyr::summarise(dplyr::across("Her2", mean, na.rm = TRUE)) %>%
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
      dplyr::group_by(scPAS) %>%
      dplyr::summarise(dplyr::across("Her2", mean, na.rm = TRUE)) %>%
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
      dplyr::group_by(scAB) %>%
      dplyr::summarise(dplyr::across("Her2", mean, na.rm = TRUE)) %>%
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
      dplyr::group_by(scPP) %>%
      dplyr::summarise(dplyr::across("Her2", mean, na.rm = TRUE)) %>%
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

  mean_scores <- AUCell::rbind(
    mean_scores_scissor,
    mean_scores_scpas,
    mean_scores_scab,
    mean_scores_scpp
  ) %>%
    tidyr::unite("screen", screen_method, screen_ressult, sep = "_") %>%
    dplyr::mutate("Sample" = sample)
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
    subtitle = "sc-GSE161528 HER2, survival as phenotype",
    y = "Mean ssGSEA Score",
    x = "Bulk Sample"
  ) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(
  file.path(plot_dir, "her2_mean_scores_survival.png"),
  plot = p,
  width = 6,
  height = 7,
  dpi = 400
)
