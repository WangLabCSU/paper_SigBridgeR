log2_transform_bulk <- function(
  x,
  method = c("cpm", "tpm", "rle", "none"),
  pseudocount = "adaptive",
  lib_size = NULL,
  gene_length = NULL,
  return_sparse = NULL,
  round_digits = NULL
) {
  # 🔧 依赖检查（轻量，仅 base + Matrix）
  if (!requireNamespace("Matrix", quietly = TRUE)) {
    stop("Package 'Matrix' required for sparse support.")
  }
  library(Matrix, warn.conflicts = FALSE)

  # 🧪 输入检查
  method <- match.arg(method)
  #   stopifnot(
  #     is.numeric(x),
  #     is.matrix(x) || inherits(x, "Matrix") || is.data.frame(x)
  #   )

  # 标准化为 matrix / dgCMatrix（data.frame → matrix）
  if (is.data.frame(x)) {
    rn <- rownames(x)
    cn <- colnames(x)
    x <- as.matrix(x)
    if (!is.null(rn)) {
      rownames(x) <- rn
    }
    if (!is.null(cn)) colnames(x) <- cn
  }

  # 记录原始类型（用于返回）
  is_sparse_input <- inherits(x, "sparseMatrix")
  return_sparse <- if (is.null(return_sparse)) {
    is_sparse_input
  } else {
    return_sparse
  }

  # 确保为 dgCMatrix（高效稀疏操作）
  if (is_sparse_input && !inherits(x, "dgCMatrix")) {
    x <- as(x, "dgCMatrix")
  } else if (!is_sparse_input && is.matrix(x)) {
    # 可选：若密度 < 0.1，自动转稀疏（可注释掉）
    # if (mean(x != 0) < 0.1) x <- as(x, "dgCMatrix")
  }

  # 🔍 检查是否已是 log-scale（防重复转换）
  # 策略：若最小值 > 0 且非整数 → 警告
  if (!is.null(dim(x)) && length(x) > 0) {
    x_min <- if (inherits(x, "dgCMatrix")) {
      if (x@x[1] <= 0) min(x@x) else min(c(0, x@x))
    } else {
      min(x)
    }
    if (x_min > 0 && !all.equal(x, round(x), check.attributes = FALSE)) {
      warning(
        "Input appears already log-transformed (min > 0 & non-integer). Consider method='none'."
      )
    }
  }

  # 📏 计算文库大小（若未提供）
  if (is.null(lib_size)) {
    if (inherits(x, "dgCMatrix")) {
      lib_size <- colSums(x)
    } else {
      lib_size <- colSums(x, na.rm = TRUE)
    }
    # 替换 0 → 中位数（避免 Inf）
    if (any(lib_size == 0)) {
      lib_size[lib_size == 0] <- median(lib_size[lib_size > 0])
    }
  } else {
    stopifnot(length(lib_size) == ncol(x))
  }

  # 🔁 归一化
  x_norm <- switch(
    method,
    cpm = {
      # CPM = counts / lib_size * 1e6
      if (inherits(x, "dgCMatrix")) {
        sweep2_sp(x, lib_size / 1e6, margin = 2, FUN = "/")
      } else {
        sweep(x, 2, lib_size / 1e6, "/")
      }
    },
    tpm = {
      stopifnot(!is.null(gene_length))
      stopifnot(length(gene_length) == nrow(x))
      # RPKM = (counts × 1e9) / (gene_len × lib_size)
      # TPM = RPKM / sum(RPKM) × 1e6
      if (inherits(x, "dgCMatrix")) {
        # Step 1: RPKM-like
        x_rpkm <- sweep2_sp(
          x,
          gene_length * lib_size,
          margin = c(1, 2),
          FUN = "/"
        ) *
          1e9
        # Step 2: per-sample scaling
        colsums <- colSums(x_rpkm)
        sweep2_sp(x_rpkm, colsums / 1e6, margin = 2, FUN = "/")
      } else {
        rpkm <- (x * 1e9) / (gene_length %o% lib_size)
        tpm <- sweep(rpkm, 2, colSums(rpkm) / 1e6, "/")
        tpm
      }
    },
    rle = {
      # RLE (like DESeq2 size factors)
      geo_means <- apply(x + 1, 1, function(row) exp(mean(log(row + 1)))) # +1 for zeros
      size_factors <- apply(x + 1, 2, function(col) {
        median((col + 1) / geo_means, na.rm = TRUE)
      })
      if (inherits(x, "dgCMatrix")) {
        sweep2_sp(x, size_factors, margin = 2, FUN = "/")
      } else {
        sweep(x, 2, size_factors, "/")
      }
    },
    none = x
  )

  # 🔢 确定 pseudocount
  pc <- switch(
    pseudocount,
    fixed = if (is.numeric(pseudocount) && length(pseudocount) == 1) {
      pseudocount
    } else {
      1
    },
    adaptive = {
      # 自适应：1% 非零表达的中位数，但 ≥ 0.5
      if (inherits(x_norm, "dgCMatrix")) {
        nonzero_vals <- x_norm@x[x_norm@x > 0]
      } else {
        nonzero_vals <- x_norm[x_norm > 0]
      }
      if (length(nonzero_vals) == 0) {
        1
      } else {
        max(0.5, quantile(nonzero_vals, 0.01, na.rm = TRUE))
      }
    },
    median_nonzero = {
      if (inherits(x_norm, "dgCMatrix")) {
        median(x_norm@x[x_norm@x > 0], na.rm = TRUE)
      } else {
        median(x_norm[x_norm > 0], na.rm = TRUE)
      }
    },
    {
      if (!is.numeric(pseudocount) || length(pseudocount) != 1) {
        stop("Invalid pseudocount")
      }
      pseudocount
    }
  )

  # 📈 log2 转换
  if (inherits(x_norm, "dgCMatrix")) {
    x_log <- x_norm
    x_log@x <- log2(x_norm@x + pc)
  } else {
    x_log <- log2(x_norm + pc)
  }

  # ✂️ 可选：四舍五入（节省存储，如 float32）
  if (!is.null(round_digits)) {
    if (inherits(x_log, "dgCMatrix")) {
      x_log@x <- round(x_log@x, round_digits)
    } else {
      x_log <- round(x_log, round_digits)
    }
  }

  # 🔄 返回同类型
  if (return_sparse && !inherits(x_log, "sparseMatrix")) {
    x_log <- as(x_log, "dgCMatrix")
  } else if (!return_sparse && inherits(x_log, "sparseMatrix")) {
    x_log <- as.matrix(x_log)
  }

  # 🏷️ 保留行列名
  if (!is.null(rownames(x))) {
    rownames(x_log) <- rownames(x)
  }
  if (!is.null(colnames(x))) {
    colnames(x_log) <- colnames(x)
  }

  return(x_log)
}

# 🔧 辅助函数：稀疏矩阵 sweep（高效列/行操作）
sweep2_sp <- function(x, STATS, margin, FUN) {
  # x: dgCMatrix, STATS: vector, margin: 1=row / 2=col, FUN: "/" or "*"
  stopifnot(inherits(x, "dgCMatrix"))
  if (margin == 2) {
    # 列操作：x[, j] / STATS[j]
    for (j in seq_along(STATS)) {
      idx <- (x@p[j] + 1):x@p[j + 1]
      if (length(idx) > 0 && x@p[j] < x@p[j + 1]) {
        x@x[idx] <- match.fun(FUN)(x@x[idx], STATS[j])
      }
    }
  } else if (margin == 1) {
    # 行操作：需转置或循环
    stop("Row-wise sweep for sparse not implemented (use t() if needed)")
  }
  x
}
