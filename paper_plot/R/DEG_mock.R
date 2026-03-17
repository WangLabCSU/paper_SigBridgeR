# ============================================================
# DEG 模拟数据分析：火山图 + ComplexHeatmap热图
# ============================================================

library(ggplot2)
library(ggrepel)
library(dplyr)
library(tibble)
library(ComplexHeatmap)
library(circlize)

# ============================================================
# 1. 生成DEG模拟数据
# ============================================================
set.seed(12345)

# 参数设置
n_genes <- 20000 # 总基因数
n_case <- 10 # 病例组样本数
n_ctrl <- 10 # 对照组样本数
n_up <- 50 # 上调基因数
n_down <- 60 # 下调基因数
n_ns <- n_genes - n_up - n_down

# 组内样本相关性（模拟生物学重复）
rho <- 0.75 # 组内相关系数
Sigma_ctrl <- matrix(rho, n_ctrl, n_ctrl) + diag(1 - rho, n_ctrl)
Sigma_case <- matrix(rho, n_case, n_case) + diag(1 - rho, n_case)
L_ctrl <- chol(Sigma_ctrl)
L_case <- chol(Sigma_case)

# ---------------- 1. 生成基础表达量（模拟真实RNA-seq特性） ----------------
# 步骤1: 生成基础均值（模拟基因表达丰度分布：长尾分布）
base_mu <- 2^runif(n_genes, 3, 14) # 8 ~ 16384 counts，符合真实分布

# 步骤2: 为差异基因分配效应大小（严格控制分布）
log2fc_true <- numeric(n_genes)
log2fc_true[1:n_up] <- rnorm(n_up, mean = 1.5, sd = 0.2) # 上调集中分布
log2fc_true[(n_up + 1):(n_up + n_down)] <- rnorm(n_down, mean = -1.5, sd = 0.2) # 下调集中分布
# NS基因log2FC接近0（标准差0.08，非常集中）
log2fc_true[(n_up + n_down + 1):n_genes] <- rnorm(n_ns, mean = 0, sd = 0.08)

# 步骤3: 生成表达矩阵（负二项分布 + 组内相关性）
expr_matrix <- matrix(0, nrow = n_genes, ncol = n_ctrl + n_case)
gene_ids <- paste0("Gene_", sprintf("%05d", 1:n_genes))

for (i in 1:n_genes) {
  # 根据基础表达量设置离散度（高表达→低离散度）
  disp <- 1 / sqrt(base_mu[i]) * 2 + 0.1 # 均值-方差负相关

  # 生成对照组（多元正态→指数变换→负二项）
  z_ctrl <- rnorm(n_ctrl) %*% L_ctrl
  mu_ctrl <- base_mu[i] * exp(z_ctrl * 0.2) # 添加组内变异
  expr_matrix[i, 1:n_ctrl] <- rnbinom(n_ctrl, mu = mu_ctrl, size = 1 / disp)

  # 生成病例组（添加log2FC效应）
  fc_effect <- 2^log2fc_true[i]
  z_case <- rnorm(n_case) %*% L_case
  mu_case <- base_mu[i] * fc_effect * exp(z_case * 0.2)
  expr_matrix[i, (n_ctrl + 1):(n_ctrl + n_case)] <- rnbinom(
    n_case,
    mu = mu_case,
    size = 1 / disp
  )
}

# 添加少量技术零值（模拟dropout）
dropout_idx <- which(matrix(
  runif(n_genes * (n_ctrl + n_case)) < 0.02,
  nrow = n_genes,
  ncol = n_ctrl + n_case
))
expr_matrix[dropout_idx] <- 0

rownames(expr_matrix) <- gene_ids
colnames(expr_matrix) <- c(paste0("Ctrl_", 1:n_ctrl), paste0("Case_", 1:n_case))

# ---------------- 2. 真实计算DEG（t-test + BH校正） ----------------
deg_results <- tibble(
  gene = gene_ids,
  log2FC = numeric(n_genes),
  pvalue = numeric(n_genes),
  baseMean = rowMeans(expr_matrix)
)

# 逐基因t检验（更真实，避免limma的过度平滑）
for (i in 1:n_genes) {
  ctrl_vals <- expr_matrix[i, 1:n_ctrl]
  case_vals <- expr_matrix[i, (n_ctrl + 1):(n_ctrl + n_case)]

  # 避免除零：添加伪计数
  ctrl_log <- log2(ctrl_vals + 1)
  case_log <- log2(case_vals + 1)

  deg_results$log2FC[i] <- mean(case_log) - mean(ctrl_log)

  # 低表达基因直接设为NS（避免假阳性）
  if (deg_results$baseMean[i] < 10) {
    deg_results$pvalue[i] <- runif(1, 0.3, 0.99) # 高p值
  } else {
    tt <- t.test(case_log, ctrl_log, var.equal = FALSE)
    deg_results$pvalue[i] <- tt$p.value
  }
}

# 多重检验校正
deg_results <- deg_results %>%
  mutate(
    padj = p.adjust(pvalue, method = "BH"),
    sig = case_when(
      padj < 0.05 & log2FC > 0.8 ~ "Up",
      padj < 0.05 & log2FC < -0.8 ~ "Down",
      TRUE ~ "NS"
    ),
    # 选择最显著的基因用于label（避免过度标注）
    label = ifelse(padj < 1e-10 & abs(log2FC) > 1.8, gene, NA_character_)
  ) %>%
  # 按|log2FC|排序，确保热图展示最显著基因
  arrange(desc(abs(log2FC)))

deg_results$label[deg_results$sig %in% c("Up", "Down")] = deg_results$gene[
  deg_results$sig %in% c("Up", "Down")
]

# ============================================================
# 2. 绘制火山图
# ============================================================
volcano_plot <- ggplot(
  deg_results,
  aes(x = log2FC, y = -log10(padj), color = sig)
) +
  # NS基因：浅灰色小点
  geom_point(
    data = filter(deg_results, sig == "NS"),
    size = 0.8,
    alpha = 0.3,
    color = "gray70"
  ) +
  # 显著基因：大点+饱和色
  geom_point(data = filter(deg_results, sig != "NS"), size = 1.6, alpha = 0.9) +
  # 阈值线：灰色虚线
  geom_vline(
    xintercept = c(-0.8, 0.8),
    linetype = "dashed",
    color = "gray50",
    linewidth = 0.7
  ) +
  geom_hline(
    yintercept = -log10(0.05),
    linetype = "dashed",
    color = "gray50",
    linewidth = 0.7
  ) +
  # 智能标注（仅标注极端显著基因，避免拥挤）
  geom_text_repel(
    data = filter(deg_results, !is.na(label)),
    aes(label = label),
    size = 2.5,
    box.padding = 0.4,
    point.padding = 0.3,
    segment.color = "gray60",
    segment.size = 0.25,
    max.overlaps = 10,
    force = 3,
    seed = 123
  ) +
  scale_color_manual(
    values = c("Up" = "#D62728", "Down" = "#1F77B4", "NS" = "gray70"),
    labels = c(
      "Up" = "Up-regulated",
      "Down" = "Down-regulated",
      "NS" = "Not significant"
    )
  ) +
  labs(
    title = "Volcano Plot",
    x = "log2(Fold Change)",
    y = "-log10(Adjusted p-value)",
    color = "Significance"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    legend.position = c(0.92, 0.88),
    legend.background = element_rect(
      fill = "white",
      color = "black",
      size = 0.3
    ),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 0.5)
  ) +
  coord_cartesian(xlim = c(-3, 3), ylim = c(0, 4))

print(volcano_plot)
ggsave(volcano_plot, filename = "plot/volcano_plot.png", width = 8, height = 10)

# ============================================================
# 3. 绘制ComplexHeatmap热图
# ============================================================
# 3.1 准备列注释（样本分组）
sig_genes <- deg_results %>%
  filter(abs(log2FC) > 1 & padj < 0.01) %>%
  pull(gene)

top_genes <- head(sig_genes, 60)
heatmap_data_raw <- expr_matrix[top_genes, ]
heatmap_zscore <- t(apply(log2(heatmap_data_raw + 1), 1, scale))

col_anno_df <- data.frame(
  Group = factor(
    rep(c("Non_Positive", "Positive"), c(n_ctrl, n_case)),
    levels = c("Non_Positive", "Positive")
  )
)
rownames(col_anno_df) <- colnames(heatmap_data_raw)

col_ha <- HeatmapAnnotation(
  Group = col_anno_df$Group,
  col = list(Group = c("Non_Positive" = "#8A8A8A", "Positive" = "#FF3333")),
  annotation_name_gp = gpar(fontsize = 10, fontface = "bold")
)

row_anno_df <- deg_results %>%
  filter(gene %in% top_genes) %>%
  arrange(match(gene, top_genes)) %>%
  mutate(Direction = ifelse(log2FC > 0, "Up", "Down")) %>%
  select(Direction)

row_ha <- rowAnnotation(
  Direction = row_anno_df$Direction,
  col = list(Direction = c("Up" = "#D62728", "Down" = "#1F77B4")),
  annotation_name_gp = gpar(fontsize = 10, fontface = "bold"),
  width = unit(4, "mm")
)

ht_opt(
  heatmap_column_names_gp = gpar(fontsize = 8),
  heatmap_row_names_gp = gpar(fontsize = 8)
)

heatmap_plot <- Heatmap(
  matrix = heatmap_zscore,
  name = "Z-score",
  col = colorRamp2(c(-2.5, 0, 2.5), c("#2C3E50", "#ECF0F1", "#E74C3C")),
  top_annotation = col_ha,
  right_annotation = row_ha,
  row_names_side = "left",
  row_names_gp = gpar(fontsize = 7.5),
  column_names_gp = gpar(fontsize = 9, angle = 45),
  column_split = col_anno_df$Group,
  row_split = factor(row_anno_df$Direction, levels = c("Down", "Up")),
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  show_row_dend = TRUE,
  show_column_dend = FALSE,
  row_dend_reorder = TRUE,
  width = unit(9, "cm"),
  height = unit(11, "cm"),
  heatmap_legend_param = list(
    title = "Expression\n(Z-score)",
    direction = "horizontal",
    title_position = "topcenter"
  )
) +
  Heatmap(
    matrix = deg_results %>%
      filter(gene %in% top_genes) %>%
      arrange(match(gene, top_genes)) %>%
      pull(log2FC),
    name = "log2FC",
    width = unit(5, "mm"),
    col = colorRamp2(c(-2.2, 0, 2.2), c("#1F77B4", "white", "#D62728")),
    show_row_names = FALSE,
    border = "white"
  )

pdf("plot/DEG_heatmap_plot.pdf")
draw(
  heatmap_plot,
  heatmap_legend_side = "bottom",
  padding = unit(c(2, 5, 2, 2), "mm")
)
dev.off()
