# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd(usethis::proj_path())

set.seed(123)

# library(Seurat)
library(GSVA)
# library(ComplexHeatmap)
# library(clusterProfiler)
# library(dplyr)
# library(ggplot2)

# ! hgsoc-GSE165897, bulk-GSE140082

seurat <- qs::qread(
  "/home/data/sigbridger/benchmark_data/brca/seurat_her2.qs",
  nthreads = 4
)
# seurat <- SCTransform(seurat)
# DefaultAssay(seurat) <- "SCT" # 推荐用 SCT 或 logNorm 计数

expr <- as.matrix(SeuratObject::LayerData(seurat, layer = "data"))

# 2. 计算各 cluster 的 marker（对比其余所有细胞）
# ! from Negative control markes
# ! random selected
markers <- data.table::fread(
  "Tmp/ssGSEA_negative_compare/brca/her2_random20_markers_100rep.csv"
)

mat_out <- NULL


sample <- colnames(markers)
for (i in seq_len(ncol(markers))) {
  genesets <- list(markers[[i]])
  names(genesets) <- paste0('Neg_ctrl_', sample[i])

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
  # 5. 构造 GSVA 参数对象（≥1.50 版写法）
  ssgsea_param <- gsvaParam(
    exprData = expr,
    geneSets = genesets,
    # kcdf = auto # *自动选择了gaussian
  )

  set.seed(123)
  # ! 决速步
  es.mat <- gsva(ssgsea_param) # 行为 cluster-通路，列为细胞

  # 拼接到结果矩阵
  mat_out <- cbind(mat_out, es.mat)
  gc(verbose = FALSE)
  cli::cli_h2("Finish {i}")
}

# 保存完整的结果矩阵
qs::qsave(
  mat_out,
  "Tmp/ssGSEA_negative_compare/brca/her2_Sample_100_ssgsea_score.qs"
)
# PID=1000308
