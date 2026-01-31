library(ggplot2)
library(dplyr)
library(tidyr)

# ========== 1. 生成模拟轨迹数据 ==========
set.seed(123)

# 参数设置
n_cells_main <- 300 # 主干细胞数
n_cells_branch1 <- 150 # 分支1细胞数
n_cells_branch2 <- 120 # 分支2细胞数

# 主干轨迹：沿曲线生成
pseudotime_main <- sort(runif(n_cells_main, 0, 1))
x_main <- pseudotime_main * 5
y_main <- 0.5 * sin(2 * pi * pseudotime_main) # 正弦波动模拟生物噪声

# 分支1：从主干中点分出
branch_point <- 0.5
idx_branch <- which.min(abs(pseudotime_main - branch_point))
x_bp <- x_main[idx_branch]
y_bp <- y_main[idx_branch]

pseudotime_b1 <- sort(runif(n_cells_branch1, 0, 0.5))
x_branch1 <- x_bp + pseudotime_b1 * 3 * cos(pi / 4)
y_branch1 <- y_bp + pseudotime_b1 * 3 * sin(pi / 4)

# 分支2：从同一点向另一方向分出
pseudotime_b2 <- sort(runif(n_cells_branch2, 0, 0.5))
x_branch2 <- x_bp + pseudotime_b2 * 3 * cos(-pi / 4)
y_branch2 <- y_bp + pseudotime_b2 * 3 * sin(-pi / 4)

# 合并为数据框
trajectory_data <- bind_rows(
  tibble(
    cell_id = paste0("cell_", 1:n_cells_main),
    x = x_main + rnorm(n_cells_main, sd = 0.15), # 添加空间噪声
    y = y_main + rnorm(n_cells_main, sd = 0.15),
    pseudotime = pseudotime_main,
    branch = "Main",
    lineage = "Lineage_A"
  ),
  tibble(
    cell_id = paste0(
      "cell_",
      (n_cells_main + 1):(n_cells_main + n_cells_branch1)
    ),
    x = x_branch1 + rnorm(n_cells_branch1, sd = 0.12),
    y = y_branch1 + rnorm(n_cells_branch1, sd = 0.12),
    pseudotime = branch_point + pseudotime_b1,
    branch = "Branch_1",
    lineage = "Lineage_B"
  ),
  tibble(
    cell_id = paste0(
      "cell_",
      (n_cells_main + n_cells_branch1 + 1):(n_cells_main +
        n_cells_branch1 +
        n_cells_branch2)
    ),
    x = x_branch2 + rnorm(n_cells_branch2, sd = 0.12),
    y = y_branch2 + rnorm(n_cells_branch2, sd = 0.12),
    pseudotime = branch_point + pseudotime_b2,
    branch = "Branch_2",
    lineage = "Lineage_C"
  )
)

# ========== 2. 模拟基因表达（沿轨迹动态变化） ==========
# 基因1：在分支1高表达
trajectory_data <- trajectory_data %>%
  mutate(
    gene_marker1 = ifelse(
      branch == "Branch_1",
      5 + 3 * (pseudotime - branch_point),
      1 + rnorm(n(), sd = 0.5)
    ),
    gene_marker1 = pmax(0, gene_marker1) # 避免负值
  )

# ========== 3. 可视化轨迹 ==========
# 3.1 轨迹结构 + 伪时间着色
p1 <- ggplot(trajectory_data, aes(x = x, y = y, color = pseudotime)) +
  geom_point(alpha = 0.7, size = 1.5) +
  scale_color_viridis_c(name = "Pseudotime", option = "plasma") +
  labs(
    title = "Simulated Cell Trajectory (Pseudotime)",
    x = "Dimension 1",
    y = "Dimension 2"
  ) +
  cowplot::theme_cowplot() +
  theme(legend.position = "right")

# 3.2 按分支/谱系着色
p2 <- ggplot(trajectory_data, aes(x = x, y = y, color = lineage)) +
  geom_point(alpha = 0.8, size = 1.8) +
  scale_color_manual(
    values = c(
      "Lineage_A" = "#2E86AB",
      "Lineage_B" = "#A23B72",
      "Lineage_C" = "#F18F01"
    ),
    name = "Lineage"
  ) +
  labs(
    title = "Trajectory with Lineage Annotation",
    x = "Dimension 1",
    y = "Dimension 2"
  ) +
  cowplot::theme_cowplot() +
  theme(legend.position = "right")

# 3.3 基因表达沿轨迹变化
p3 <- ggplot(
  trajectory_data,
  aes(x = pseudotime, y = gene_marker1, color = lineage)
) +
  geom_point(alpha = 0.6, size = 1) +
  geom_smooth(method = "loess", se = TRUE, span = 0.3) +
  scale_color_manual(
    values = c(
      "Lineage_A" = "#2E86AB",
      "Lineage_B" = "#A23B72",
      "Lineage_C" = "#F18F01"
    )
  ) +
  labs(
    title = "Gene Expression Along Pseudotime",
    x = "Pseudotime",
    y = "Expression (log scale)"
  ) +
  cowplot::theme_cowplot() +
  theme(legend.position = "none")

# ========== 4. 展示结果 ==========
ggsave("plot/cell_trajectory_pseudotime.png", p1, height = 6, width = 10)
ggsave("plot/cell_trajectory_lineage.png", p2, height = 6, width = 10)
ggsave("plot/cell_trajectory_gene_expression.png", p3, height = 6, width = 10)
