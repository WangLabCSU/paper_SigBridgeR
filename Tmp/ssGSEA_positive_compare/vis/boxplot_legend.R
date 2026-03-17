setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(ggplot2)

temp <- data.frame(
  cluster = c("Positive", "Negative", "Neutral"),
  ssgsea_score = c(1, 2, 3)
)

group_colors <- c(
  "Positive" = "#ff3333",
  "Negative" = "#386c9b",
  "Neutral" = "#CECECE"
)

p <- ggplot(
  temp,
  aes(
    x = cluster,
    y = ssgsea_score,
    fill = cluster, # box color
  )
) +
  geom_boxplot() +
  scale_fill_manual(
    name = "Group", # 图例标题
    values = group_colors # 颜色映射
  )

ggsave(
  filename = file.path("boxplot_legend.png"),
  plot = p,
  width = 10,
  height = 8,
  dpi = 400
)
