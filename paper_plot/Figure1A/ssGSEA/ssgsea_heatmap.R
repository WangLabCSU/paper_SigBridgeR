library(dplyr)
library(tidyr)
library(ggplot2)

setwd(file.path(
  usethis::proj_path(),
  "paper_plot/Figure1A/ssGSEA"
))

score <- qs::qread(
  "../../../4_positive_ctrl/esmat/binary/luad/TCGA_LUAD/ssGSEA_score_TCGA_LUAD.qs"
)

method_cols <- c(
  "scissor",
  "scPAS",
  "scPP",
  "DEGAS",
  "LP_SGL",
  "PIPET",
  "SCIPAC",
  "scAB"
)

plot_data <- score |>
  select(
    pos_ssGSEA = binary_TCGA_LUAD_pos_ssGSEA,
    neg_ssGSEA = binary_TCGA_LUAD_neg_ssGSEA,
    all_of(method_cols)
  ) |>
  pivot_longer(
    cols = all_of(method_cols),
    names_to = "Method",
    values_to = "Classification"
  ) |>
  mutate(
    Classification = recode(
      Classification,
      Positive = "Positive",
      Negative = "Negative",
      Neutral = "NS",
      Other = "NS"
    ),
    Method = factor(Method, levels = method_cols),
    Classification = factor(
      Classification,
      levels = c("Positive", "Negative", "NS")
    )
  ) |>
  pivot_longer(
    cols = c(pos_ssGSEA, neg_ssGSEA),
    names_to = "ssGSEA_type",
    values_to = "ssGSEA_score"
  ) |>
  group_by(ssGSEA_type, Method, Classification) |>
  summarise(
    mean_ssGSEA_score = mean(ssGSEA_score, na.rm = TRUE),
    .groups = "drop"
  ) |>
  mutate(
    ssGSEA_type = recode(
      ssGSEA_type,
      pos_ssGSEA = "pos_ssGSEA",
      neg_ssGSEA = "neg_ssGSEA"
    )
  )

ggplot(
  plot_data,
  aes(x = Classification, y = Method, fill = mean_ssGSEA_score)
) +
  geom_tile(color = "white", linewidth = 0.3) +
  facet_wrap(~ssGSEA_type) +
  scale_fill_gradient(
    low = "white",
    high = "#a02020",
    name = "Mean ssGSEA\nscore"
  ) +
  labs(x = NULL, y = NULL) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.background = element_blank(),
    panel.grid = element_blank()
  )

ggsave("ssgse_htmap.png")
