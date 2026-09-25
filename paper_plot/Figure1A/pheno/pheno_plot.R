setwd(file.path(
  usethis::proj_path(),
  "paper_plot/Figure1A/pheno"
))

library(ggplot2)
library(dplyr)
library(tidyr)

data_dir <- "/home/data/sigbridger/benchmark_data/"

pheno <- qs::qread(file.path(data_dir, "lung/TCGA_LUAD_pheno.qs"))

# ----------------------------------------------------------------------------------------------------------
# * Viz: 表型元数据可用性热图（迷你元素）
# 行 = 表型字段，列 = 样本；每个字段一种颜色（调色盘），缺失值留白。
# 只保留轴标题与图表题，其余元素（图例、网格、x 轴刻度标签）全部去掉。

# 仅保留有信息量的字段：剔除样本/患者 ID、常量列，以及整列为空的列
fields <- pheno %>%
  select(-all_of(c("sample", "patient"))) %>%
  select(where(~ length(unique(.x[!is.na(.x)])) > 1L)) %>%
  colnames()

# 调色盘：按 fields 顺序，每个字段取一种颜色
field_palette <- c(
  "#B8D9A0",
  "#E8E0DB",
  "#A8D9D8",
  "#F5C0A8",
  "#D5C9E8",
  "#9BBBD9",
  "#B8A8E5",
  "#D5E5C0",
  "#E8D8E8",
  "#88D5B0",
  "#E5B8E5",
  "#88E0D5",
  "#A8D5A8",
  "#E5A8D0",
  "#F5C8E5",
  "#F0E890",
  "#F5E8B8",
  "#B8BCA0",
  "#C8E5D0",
  "#B888D5",
  "#C5B8E8",
  "#E5C090",
  "#D5A0A8",
  "#88A8D5",
  "#F5D890",
  "#D8E5B0",
  "#B8D5D0",
  "#E0E8D8",
  "#E0A8E5",
  "#A888D5",
  "#C5B0B8",
  "#E5A8E8",
  "#E8A8D5",
  "#D0C5B8",
  "#F5A098",
  "#B0E5D0",
  "#C0D8E5",
  "#90C8D5",
  "#E8A0B0",
  "#F5E0B8"
)
field_cols <- setNames(field_palette[seq_along(fields)], fields)

# 长表：(样本, 字段) 一行；有值用该字段的颜色，缺失留白
# 字段类型混合（数值/字符），统一转成字符后再合并为一列
avail <- pheno %>%
  select(sample, all_of(fields)) %>%
  mutate(across(-sample, as.character)) %>%
  pivot_longer(-sample, names_to = "field", values_to = "value") %>%
  mutate(fill_key = ifelse(is.na(value), "Missing", as.character(field)))

# 样本按 tissue_type 分组：肿瘤（tumor == 1）在前、正常在后，组内保持原文件顺序
sample_levels <- pheno$sample[order(-pheno$tumor, seq_len(nrow(pheno)))]

avail <- avail %>%
  mutate(
    sample = factor(sample, levels = sample_levels),
    field = factor(field, levels = fields)
  )

p_pheno <- avail %>%
  ggplot(aes(x = sample, y = field, fill = fill_key)) +
  geom_tile(color = NA) +
  scale_fill_manual(values = c(field_cols, "Missing" = "white")) +
  scale_x_discrete(expand = c(0, 0), breaks = NULL) +
  scale_y_discrete(expand = c(0, 0)) +
  labs(
    title = "Phenotype metadata",
    x = "Sample",
    y = "Phenotype field"
  ) +
  theme_minimal(base_size = 8) +
  theme(
    legend.position = "none",
    panel.grid = element_blank(),
    axis.ticks = element_blank(),
    axis.text.y = element_text(size = 6.5, color = "black"),
    axis.title = element_text(size = 9, face = "bold"),
    plot.title = element_text(size = 10, face = "bold", hjust = 0),
    plot.title.position = "plot"
  )

ggplot2::ggsave("pheno_heatmap.png", p_pheno, width = 6, height = 3, dpi = 400)

p_pheno
