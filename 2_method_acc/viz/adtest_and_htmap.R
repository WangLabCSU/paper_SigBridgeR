setwd(file.path(usethis::proj_path(), "2_method_acc/viz"))

library(dplyr)
library(data.table)
source("pairwise_ad_test.R")
source("ad_p_htmap.R")


combined2 <- data.table::fread("combined.csv")

# stats <- combined2[,
#   list(
#     N = .N,
#     mean_f1 = mean(F1, na.rm = TRUE),
#     sd_f1 = sd(F1, na.rm = TRUE),
#     mean_acc = mean(Accuracy, na.rm = TRUE),
#     sd_acc = sd(Accuracy, na.rm = TRUE),
#     median_f1 = median(F1, na.rm = TRUE),
#     median_acc = median(Accuracy, na.rm = TRUE),
#     min_f1 = min(F1, na.rm = TRUE),
#     max_f1 = max(F1, na.rm = TRUE),
#     min_acc = min(Accuracy, na.rm = TRUE),
#     max_acc = max(Accuracy, na.rm = TRUE)
#   ),
#   by = method_name
# ]

# stats2 <- combined2[,
#   list(
#     N = .N,
#     mean_f1 = mean(F1, na.rm = TRUE),
#     sd_f1 = sd(F1, na.rm = TRUE),
#     mean_acc = mean(Accuracy, na.rm = TRUE),
#     sd_acc = sd(Accuracy, na.rm = TRUE),
#     median_f1 = median(F1, na.rm = TRUE),
#     median_acc = median(Accuracy, na.rm = TRUE),
#     min_f1 = min(F1, na.rm = TRUE),
#     max_f1 = max(F1, na.rm = TRUE),
#     min_acc = min(Accuracy, na.rm = TRUE),
#     max_acc = max(Accuracy, na.rm = TRUE)
#   ),
#   by = list(method_name, tissue, screen_method)
# ]

ad_test_f1 <- pairwise_dist_test(
  dt = combined2,
  group_col = "method_name",
  value_col = "F1"
)
ad_test_acc <- pairwise_dist_test(
  dt = combined2,
  group_col = "method_name",
  value_col = "Accuracy"
)

data.table::fwrite(ad_test_f1, "ad_test_f1.csv")
data.table::fwrite(ad_test_acc, "ad_test_acc.csv")

ad_test_f1 <- data.table::fread("ad_test_f1.csv")
ad_test_acc <- data.table::fread("ad_test_acc.csv")

convert_to_p_mat <- function(wide) {
  wide <- tibble::as_tibble(wide)
  row_names <- unlist(wide$group1)
  total <- union(colnames(wide[, -1]), row_names)

  missing <- setdiff(total, row_names)
  missing_tbl <- data.frame(rep(NA, length(total)))
  missing_tbl <- t(missing_tbl)
  colnames(missing_tbl) <- total
  rownames(missing_tbl) <- missing
  res <- dplyr::bind_rows(
    tibble::column_to_rownames(wide, "group1"),
    as.data.frame(missing_tbl)
  )

  res <- res[sort(rownames(res)), sort(colnames(res))]
  res <- as.matrix(res)
  diag(res) <- 0L
  res <- Matrix::forceSymmetric(res, uplo = "U")
  as.matrix(res)
}


wide_f1 <- ad_test_f1[, .(group1, group2, `AD.aympt. P-value`)]
wide_f1 <- data.table::dcast(
  wide_f1,
  formula = group1 ~ group2,
  value.var = "AD.aympt. P-value"
) %>%
  convert_to_p_mat()

ad_p_htmap(
  corr = wide_f1,
  p_mat = wide_f1,
  filename = "ad_p_htmap_f1.png",
  width = 3500,
  height = 3000,
  tolerance = 1e-240
)

wide_acc <- ad_test_acc[, .(group1, group2, `AD.aympt. P-value`)]
wide_acc <- data.table::dcast(
  wide_acc,
  formula = group1 ~ group2,
  value.var = "AD.aympt. P-value"
) %>%
  convert_to_p_mat()
ad_p_htmap(
  corr = wide_acc,
  p_mat = wide_acc,
  filename = "ad_p_htmap_acc.png",
  width = 3500,
  height = 3000,
  tolerance = 1e-240
)
