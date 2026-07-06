# * 缩写基因集名称映射函数
shorten_gene_sets <- function(x) {
  # x: character vector，原始列名/基因集名称
  short <- x

  # 1. 保留方向性修饰词，用箭头/缩写替代
  short <- gsub("Genes up-regulated", "UP:", short, ignore.case = TRUE)
  short <- gsub("Genes down-regulated", "DOWN:", short, ignore.case = TRUE)
  short <- gsub("up-regulated", "UP:", short, ignore.case = TRUE)
  short <- gsub("down-regulated", "DOWN:", short, ignore.case = TRUE)

  # 2. UP:/DOWN: 标签后去除多余的前置词
  short <- gsub(
    "(UP:|DOWN:)\\s*(by |during |in response to |through |of )",
    "\\1",
    short
  )

  # 3. 去除 "Genes " 前缀（普通描述性前缀）
  short <- gsub(
    "Genes (encoding |involved in |defining |defined |important for |regulating |mediating |involve in )",
    "",
    short,
    ignore.case = TRUE
  )
  short <- gsub("Genes specifically ", "", short, ignore.case = TRUE)
  short <- gsub("^Genes ", "", short)

  # 4. 去除其他冗余前缀
  short <- gsub(
    "(A subgroup of |processing of |production of |formation of |development of |metabolism of |components of )",
    "",
    short,
    ignore.case = TRUE
  )
  short <- gsub("^genes ", "", short)

  # 5. 去除长注释后缀
  short <- gsub(
    "(, as in.*|, e\\.g\\..*|, a cellular.*|, which.*|; also.*)",
    "",
    short
  )

  # 6. 去除尾部不完整的标点/空格
  short <- gsub("[,\\s]+$", "", short)
  short <- gsub("\\.$", "", short)

  # 7. 截断过长的名称
  short <- ifelse(
    nchar(short) > 45,
    paste0(substr(short, 1, 42), "..."),
    short
  )

  short
}
