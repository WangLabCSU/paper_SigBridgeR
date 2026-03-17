# ! ----------------可视化----------------
setwd(
  "/home/yyx/R/Project/R_code/SigBridgeR/Tmp/method_compare/survival/brca/her2/GSE42568"
)
# * cnv_status标记为肿瘤细胞，作为判断依据
scpas_random_search <- data.table::data.table()
for (i in 1:6) {
  for (j in 1:10) {
    file_path <- paste0("scpas/", i, j, "_scpas_random_search.csv")
    if (file.exists(file_path)) {
      scpas_random_search <- cbind(
        scpas_random_search,
        data.table::fread(file_path)
      )
    }
  }
}

keep_first_V1 <- function(df) {
  nms <- colnames(df)
  # 找出所有 "V1" 的位置
  v1_idx <- which(nms == "V1")
  # 保留第一个 V1，其余 V1 列剔除
  keep_idx <- setdiff(seq_along(nms), v1_idx[-1])
  df[, ..keep_idx]
}
scpas_random_search <- keep_first_V1(scpas_random_search)

scpas_random_search[, c(
  "process_12",
  "process_18",
  "process_13",
  "process_19",
  "process_22",
  "process_26",
  "process_35",
  "process_45",
  "process_410",
  "process_59"
)] <- FALSE

data.table::fwrite(
  scpas_random_search,
  file = "scpas_random_search.csv",
  row.names = TRUE
)

seurat <- qs::qread(file.path(
  "/home/data/sigbridger/benchmark_data/brca",
  "seurat_her2.qs"
))
seurat_tum = readRDS(
  '/home/data/data-resource/single-cell/BRCA/GSE161529_Seurat/SeuratObject_HER2Tum.rds'
)
tumor_cells = rownames(seurat_tum[[]])
scpas_random_search$benchmark = ifelse(
  colnames(seurat) %chin% tumor_cells,
  TRUE,
  FALSE
)

# * 计算指标
source(
  "/home/yyx/R/Project/R_code/SigBridgeR/Tmp/method_compare/compute_metrics.R"
)
metrics = compute_metrics(scpas_random_search)
tmetrics = data.table::transpose(metrics)
colnames(tmetrics) = rownames(metrics)
rownames(tmetrics) = grep('process', colnames(scpas_random_search), value = T)

set.seed(42568)
arg_samples <- data.frame(
  "nfeature" = sample(seq(500, 5000, by = 100), 50, replace = TRUE), # 第1维
  "imputation" = sample(c("None", "KNN", "ALRA"), 50, replace = TRUE), # 第2维
  "independent" = sample(c(FALSE, TRUE), 50, replace = TRUE) # 第3维
) %>%
  add_row(nfeature = 2000, imputation = "None", independent = TRUE)

library(data.table)
arg_samples = data.table::as.data.table(arg_samples)
arg_samples[, names(tmetrics) := tmetrics]

data.table::fwrite(
  arg_samples,
  file = "scpas_arg_samples.csv",
  row.names = TRUE
)

# * 图
arg_samples <- arg_samples %>%
  tidyr::unite("y_axis", imputation, independent, sep = " | ")

library(ggplot2)
p <- ggplot(arg_samples, aes(x = nfeature, y = y_axis, fill = Accuracy)) +
  geom_point(size = 6, alpha = .9, shape = 21, color = "black") +
  scale_fill_gradient(
    low = "white",
    high = "red",
    name = "Accuracy"
  ) +
  # 1. 保证横轴按真实数值排序，并强制每 500 一个刻度
  scale_x_continuous(
    breaks = seq(500, 5000, by = 250),
    labels = seq(500, 5000, by = 250)
  ) +
  labs(
    title = "Validation of the Screening Efficiency of scPAS under Random Parameters",
    subtitle = "x = nfeature, y = imputation | independent",
    x = "Number of features",
    y = NULL
  ) +
  theme_minimal(base_size = 14) + # 全局字体基准
  theme(
    # 2. 轴文字放大
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 13),
    # 3. x 轴 45° 倾斜
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    # 4. 图例文字放大
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 13)
  )

ggsave(
  filename = "scpas_acc.png",
  plot = p,
  width = 10,
  height = 8,
  dpi = 400
)
