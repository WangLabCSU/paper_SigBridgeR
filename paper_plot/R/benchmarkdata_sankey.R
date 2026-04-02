library(ggplot2)
library(ggalluvial)
library(dplyr)

setwd(usethis::proj_path())

# 示例数据创建
set.seed(123)
data_survival <- data.frame(
  single_cell = c(
    rep("GSE161529\nHER2", 3),
    rep("GSE161529\nTNBC", 3),
    rep("GSE123902", 4),
    rep("GSE1678897", 2)
  ),
  bulk = c(
    rep(c("GSE42568", "GSE162228", "TCGA-BRCA"), 2),
    c("GSE3141", "GSE8894", "GSE31210", "TCGA-LUAD"),
    c("GSE32062", "GSE140082")
  ),
  phenotype = c("Survival")
)

data_binary <- data.frame(
  single_cell = c(
    rep("GSE161529\nHER2", 3),
    rep("GSE161529\nTNBC", 3),
    "GSE123902",
    rep("GSE1678897", 2)
  ),
  bulk = c(
    rep(c("GSE42568", "GSE162228", "TCGA-BRCA"), 2),
    c("TCGA-LUAD"),
    c("GSE9891", "GSE140082")
  ),
  phenotype = c("Binary")
)

data <- rbind(data_survival, data_binary)
# 按搭配关系排序：根据single_cell的顺序来排序bulk
data_custom_ordered <- data %>%
  mutate(
    # 保持single_cell的原始顺序或自定义顺序
    single_cell = factor(
      single_cell,
      levels = c(
        "GSE161529\nHER2",
        "GSE161529\nTNBC",
        "GSE123902",
        "GSE1678897"
      )
    ),

    # 根据single_cell的搭配来排序bulk
    bulk = factor(bulk, levels = unique(bulk[order(single_cell)])),

    phenotype = factor(phenotype)
  )

# 桑基图
p <- ggplot(
  data_custom_ordered,
  aes(
    axis1 = single_cell,
    axis2 = bulk,
    axis3 = phenotype
  )
) +
  geom_alluvium(aes(fill = single_cell), alpha = 0.25) +
  scale_fill_manual(
    values = c(
      "GSE161529\nHER2" = "#aa3824",
      "GSE161529\nTNBC" = "#a76723",
      "GSE1678897" = "#408e22",
      "GSE123902" = "#2261a8"
    )
  ) +
  geom_stratum(
    alpha = 0.25,
    color = "black"
  ) +
  geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  scale_x_discrete(
    limits = c("Single Cell Datasets", "Bulk Datasets", "Phenotype Type"),
    expand = c(0.15, 0.05)
  ) +
  theme_minimal() +
  labs(title = "Data Benchmarking") +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 12),
    axis.text.y = element_blank()
  )

ggsave("paper_plot/plot/benchmarkdata_sankey.png", p, width = 11, height = 8)
