library(data.table)
library(Seurat)
library(GSVA)
library(ComplexHeatmap)
library(dplyr)
library(tidyr)
library(ggplot2)

setwd(
  "/home/yyx/R/Project/R_code/SigBridgeR/Tmp/ssGSEA/test-brca/tumor_tnbc"
)

data_path = "/home/yyx/R/Project/R_code/SigBridgeR/Tmp/ssGSEA/cellmarker2"

set.seed(123)


# ! brca-GSE161529_tnbc, bulk - GSE42568, pheno-survival

seurat <- qs::qread(
  "/home/data/sigbridger/benchmark_data/brca/seurat_tnbc.qs",
  nthreads = 4
)
# seurat <- SCTransform(seurat)
# DefaultAssay(seurat) <- "SCT" # 推荐用 SCT 或 logNorm 计数

# 2. 计算各 cluster 的 marker（对比其余所有细胞）
# ! from CellMarker2.0, breast cancer marker
markers <- data.table::fread(file.path(data_path, "breast_cancer.csv"))
markers = markers$`Cell Marker`
markers = strsplit(markers, split = ",")
markers = unlist(markers)
markers = gsub(pattern = " ", replacement = "", x = markers)
markers = toupper(markers)

expand_slash <- function(x) {
  unlist(lapply(x, function(s) {
    if (!grepl("/", s)) {
      return(s)
    }

    parts <- strsplit(s, "/")[[1]]
    first <- parts[1]

    # 检查除 first 外的部分是否全为数字
    suffixes <- parts[-1]
    if (length(suffixes) > 0 && all(grepl("^\\d+$", suffixes))) {
      # 尝试从 first 中提取 prefix + 首数字
      m <- regexec("^(.*?)(\\d+)$", first)
      match <- regmatches(first, m)[[1]]
      if (length(match) == 3) {
        prefix <- match[2]
        first_num <- match[3]
        nums <- c(first_num, suffixes)
        return(paste0(prefix, nums))
      }
    }
    # fallback: 直接返回所有 parts（按 / 拆分）
    parts
  }))
}

markers = expand_slash(markers)


# ! FROM Cell Taxonomy https://ngdc.cncb.ac.cn/celltaxonomy/celltype/CT:00001067
markers2 <- data.table::fread(
  "/home/yyx/R/Project/R_code/SigBridgeR/Tmp/ssGSEA/cell_taxonomy/breast_cancer.csv"
)
markers2 = markers2$`Breast Cancer`

markers = list("BRCA" = unique(c(markers, markers2)))

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

qs::qsave(es.mat, "tum_tnbc_ssgsea_score.qs") # PID=40521
# es.mat <- qs::qread("tum_tnbc_ssgsea_score.qs")

cli::cli_alert_success("ssGSEA on TNBC done!")

es.mat <- qs::qread("tum_tnbc_ssgsea_score.qs")


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

data_path = "/home/data/sigbridger/benchmark_data/brca/TNBC"
seurat_GSE42568 = qs::qread(
  file.path(data_path, "GSE42568_tnbc_merged_seurat.qs"),
  nthreads = 4L
)

meta <- seurat_GSE42568@meta.data
meta$cell <- colnames(seurat_GSE42568)
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

table(es.df$scissor)
#  Neutral Negative Positive
#    24743     3342     3787
table(es.df$scPAS)
#  Neutral Negative Positive
#    31540      195      137

table(es.df$scAB)
#   Other Positive
#    50848     3802
table(es.df$scPP)
#  Neutral Negative Positive
#    28410     1936     1526

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

plot_dir = "/home/yyx/R/Project/R_code/SigBridgeR/Tmp/ssGSEA/test-brca/tumor_tnbc/plot"
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
  es.df$BRCA,
  es.df$scissor,
  es.df$scPAS,
  es.df$scAB,
  es.df$scPP
)
es.mat_ordered <- es.mat[, order, drop = FALSE]

col_anno_GSE42568 <- CreateHtmapAnno(
  seurat_GSE42568,
  es.mat = es.mat_ordered,
  color = "#00913fff"
)

col_anno_cluster <- HeatmapAnnotation(
  Cluster = es.df$seurat_clusters[order],
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
seurat_GSE162228 <- qs::qread(
  file.path(data_path, "GSE162228_tnbc_merged_seurat.qs"),
  nthreads = 4L
)

col_anno_GSE162228 <- CreateHtmapAnno(
  seurat_GSE162228,
  es.mat = es.mat_ordered,
  color = "#cf7815ff"
)

# ! bulk TCGA-BRCA
seurat_TCGA_BRCA <- qs::qread(
  file.path(data_path, "tcga_tnbc_merged_seurat.qs"),
  nthreads = 4L
)
col_anno_TCGA_BRCA <- CreateHtmapAnno(
  seurat_TCGA_BRCA,
  es.mat = es.mat_ordered,
  color = "#bf009cff"
)

bulk_lgd = Legend(
  at = c("GSE42568", "GSE162228", "TCGA-BRCA"),
  title = "Bulk Sample for Methods",
  legend_gp = gpar(
    fill = c("#00913fff", "#cf7815ff", "#bf009cff")
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

# ! TNBC Group
set.seed(333)
pallete_tnbc <- randomcoloR::distinctColorPalette(
  length(unique(es.df$group))
  # ,    runTsne = T
)
col_anno_tnbc <- HeatmapAnnotation(
  Her2_group = es.df$group[order],
  annotation_name_gp = gpar(fontface = "bold"),
  annotation_label = "Her2 Group",
  col = list(
    TNBC_group = setNames(pallete_tnbc, unique(es.df$group))
  ),
  height = unit(1, "cm"),
  gap = unit(0.6, "mm")
)

# ! 每个anno之间的间距
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
    mean(es.df$BRCA[es.df$seurat_clusters == x])
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
    col_anno_TCGA_BRCA,
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
  heatmap_width = unit(3, "npc")
)
draw(
  htmap_cmb,
  annotation_legend_list = list(bulk_lgd, screen_lgd)
)

ragg::ang_png(
  file = file.path(plot_dir, "tum_ssGSEA_tnbc_survival_heatmap.png"),
  width = 16,
  height = 10
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
      summarise(across("BRCA", mean, na.rm = TRUE)) %>%
      pivot_longer(
        cols = "BRCA",
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
      summarise(across("BRCA", mean, na.rm = TRUE)) %>%
      pivot_longer(
        cols = "BRCA",
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
      summarise(across("BRCA", mean, na.rm = TRUE)) %>%
      pivot_longer(
        cols = "BRCA",
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
      summarise(across("BRCA", mean, na.rm = TRUE)) %>%
      pivot_longer(
        cols = "BRCA",
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
      summarise(across("BRCA", mean, na.rm = TRUE)) %>%
      pivot_longer(
        cols = "BRCA",
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
      summarise(across("BRCA", mean, na.rm = TRUE)) %>%
      pivot_longer(
        cols = "BRCA",
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
      summarise(across("BRCA", mean, na.rm = TRUE)) %>%
      pivot_longer(
        cols = "BRCA",
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


mean_scores_GSE42568 <- MeanScore(seurat_GSE42568, es.mat, sample = "GSE42568")
mean_scores_GSE162228 <- MeanScore(
  seurat_GSE162228,
  es.mat,
  sample = "GSE162228"
)
mean_scores_tcga <- MeanScore(seurat_TCGA_BRCA, es.mat, sample = "TCGA-BRCA")

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
    subtitle = "sc-GSE161528 TNBC, survival as phenotype",
    y = "Mean ssGSEA Score",
    x = "Bulk Sample"
  ) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(
  file.path(plot_dir, "tnbc_mean_scores_survival.png"),
  plot = p,
  width = 5.5,
  height = 7,
  dpi = 400
)
