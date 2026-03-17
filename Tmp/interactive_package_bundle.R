interactive_package_bundle <- function(pkg_path = ".", output_name = NULL) {
  # --- 1. 初始化路径与环境 ---
  pkg_path <- normalizePath(pkg_path, mustWork = TRUE)

  # 确定输出文件名
  if (is.null(output_name)) {
    desc_file <- file.path(pkg_path, "DESCRIPTION")
    if (file.exists(desc_file)) {
      desc <- read.dcf(desc_file)
      pkg_name <- desc[1, "Package"]
      pkg_ver <- desc[1, "Version"]
      output_name <- sprintf("%s_%s_minimal.tar.gz", pkg_name, pkg_ver)
    } else {
      output_name <- "package_minimal.tar.gz"
    }
  }

  # 输出文件的绝对路径 (非常重要，因为后面会切换工作目录)
  output_path <- file.path(pkg_path, output_name)

  # --- 2. 扫描根目录文件 ---
  all_items <- list.files(pkg_path, full.names = FALSE, all.files = TRUE)
  ignore_items <- c(
    ".",
    "..",
    ".git",
    ".gitignore",
    ".Rproj",
    ".Rproj.user",
    ".Rbuildignore",
    ".DS_Store",
    output_name
  )
  all_items <- all_items[!all_items %in% ignore_items]

  if (length(all_items) == 0) {
    stop("未找到任何可打包的文件。")
  }

  # --- 3. 设置默认选中状态 ---
  essentials <- c("DESCRIPTION", "NAMESPACE", "R", "man", "src", "LICENSE")
  exclusions <- c("tests", "vignettes", "data")

  selection_status <- setNames(rep(FALSE, length(all_items)), all_items)

  for (item in all_items) {
    if (item %in% essentials) {
      selection_status[item] <- TRUE
    } else if (item %in% exclusions) {
      selection_status[item] <- FALSE
    } else {
      selection_status[item] <- TRUE
    }
  }

  # --- 4. 交互式选择循环 ---
  cat("\n=== R 包交互式打包工具 ===\n")
  cat("当前路径:", pkg_path, "\n")
  cat("操作指南：输入数字 toggle 选中状态，输入 0 开始打包。\n\n")

  repeat {
    cat(sprintf("%-4s %-20s %s\n", "ID", "Name", "Status"))
    cat(paste(rep("-", 40), collapse = ""), "\n")

    for (i in seq_along(all_items)) {
      item <- all_items[i]
      status <- if (selection_status[item]) "[x] 保留" else "[ ] 排除"
      note <- ""
      if (item == "inst" && selection_status[item]) {
        note <- " (自动排除 inst/pkgdown)"
      }
      cat(sprintf("%-4d %-20s %s%s\n", i, item, status, note))
    }
    cat("\n")

    choice <- readline(prompt = "输入文件 ID 切换状态 (0 确认打包): ")
    choice <- suppressWarnings(as.integer(choice))

    if (is.na(choice)) {
      cat("无效输入，请输入数字。\n\n")
      next
    }

    if (choice == 0) {
      break
    }

    if (choice > 0 && choice <= length(all_items)) {
      target_item <- all_items[choice]
      selection_status[target_item] <- !selection_status[target_item]
      cat(sprintf(
        "已切换 '%s' 状态为：%s\n\n",
        target_item,
        ifelse(selection_status[target_item], "保留", "排除")
      ))
    } else {
      cat("ID 超出范围。\n\n")
    }
  }

  # --- 5. 执行打包 ---
  cat("\n正在准备文件...\n")

  temp_dir <- tempfile(pattern = "pkg_build_")
  dir.create(temp_dir, showWarnings = FALSE, recursive = TRUE)

  selected_items <- names(selection_status)[selection_status == TRUE]

  if (length(selected_items) == 0) {
    unlink(temp_dir, recursive = TRUE)
    stop("未选择任何文件，打包取消。")
  }

  # 复制文件到临时目录
  copy_success <- TRUE
  for (item in selected_items) {
    src <- file.path(pkg_path, item)
    dst <- file.path(temp_dir, item)

    # 特殊处理 inst 目录
    if (item == "inst" && dir.exists(src)) {
      cat("处理 inst 目录 (排除 pkgdown)...\n")
      dir.create(dst, showWarnings = FALSE)
      inst_contents <- list.files(src, full.names = FALSE)
      for (sub_item in inst_contents) {
        if (sub_item != "pkgdown") {
          sub_src <- file.path(src, sub_item)
          # 修复：移除 warn 参数，增加 overwrite
          res <- file.copy(sub_src, dst, recursive = TRUE, overwrite = TRUE)
          if (!all(res)) copy_success <- FALSE
        }
      }
    } else {
      # 普通文件/目录复制
      # 修复：移除 warn 参数
      res <- file.copy(src, dst, recursive = TRUE, overwrite = TRUE)
      if (!all(res)) {
        warning(sprintf("复制失败：%s -> %s", src, dst))
        copy_success <- FALSE
      }
    }
  }

  if (!copy_success) {
    unlink(temp_dir, recursive = TRUE)
    stop("文件复制过程中出现错误，请检查权限或路径。")
  }

  # 验证临时目录内容 (调试用)
  temp_contents <- list.files(temp_dir)
  cat("临时目录内容检查:", paste(temp_contents, collapse = ", "), "\n")

  # 确保选中的文件都在临时目录中
  missing_files <- setdiff(selected_items, temp_contents)
  if (length(missing_files) > 0) {
    unlink(temp_dir, recursive = TRUE)
    stop(sprintf(
      "致命错误：以下文件未成功复制到临时目录：%s",
      paste(missing_files, collapse = ", ")
    ))
  }

  # 压缩 (使用 tar.gz 格式)
  cat("正在压缩...\n")

  # 关键修复：保存当前工作目录，切换到临时目录，执行 tar，再恢复
  old_wd <- getwd()
  setwd(temp_dir)

  tryCatch(
    {
      # 确保 output_path 是绝对路径，否则切换目录后会存错地方
      utils::tar(output_path, files = selected_items, compression = "gzip")

      if (file.exists(output_path)) {
        cat("\n=== 打包完成 ===\n")
        cat("输出文件:", output_path, "\n")
        cat("文件大小:", file.info(output_path)$size, "bytes\n")
      } else {
        stop("tar 命令执行完毕但未生成文件。")
      }
    },
    error = function(e) {
      cat("打包失败:", conditionMessage(e), "\n")
    },
    finally = {
      # 恢复工作目录
      setwd(old_wd)
      # 清理临时目录
      unlink(temp_dir, recursive = TRUE)
    }
  )

  return(invisible(output_path))
}
