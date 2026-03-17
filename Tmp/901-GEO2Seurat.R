GEO2Seurat <- function(
  geo_accession,
  min.cells = 3,
  min.features = 200,
  mt.pattern = "^MT-",
  qc.filter = TRUE,
  ncores = 1,
  save.rds = NULL
) {
  # 加载必要包（静默）
  suppressPackageStartupMessages({
    require(GEOquery)
    require(Seurat)
    require(cli)
    require(Matrix)
    require(future)
    require(doParallel)
  })

  # 设置并行
  if (ncores > 1) {
    future::plan("multicore", workers = ncores)
    cli::cli_alert_info("使用 {ncores} 个核心进行并行处理")
  }

  # 1. 下载GEO数据 ------------------------------------------------------------
  cli::cli_h1("步骤1: 下载GEO数据 {.field {geo_accession}}")
  cli::cli_progress_step("下载GEO数据集...")

  gse <- tryCatch(
    {
      GEOquery::getGEO(geo_accession, GSEMatrix = TRUE)
    },
    error = function(e) {
      cli::cli_abort("无法下载GEO数据: {e$message}")
    }
  )

  if (length(gse) == 0) {
    cli::cli_abort("未找到GEO数据，请检查accession number是否正确")
  }

  cli::cli_progress_done()
  cli::cli_alert_success("成功下载数据集，包含 {length(gse)} 个平台")

  # 2. 提取表达矩阵 ----------------------------------------------------------
  cli::cli_h2("步骤2: 提取表达数据")

  expr_data <- tryCatch(
    {
      exprs(gse[[1]])
    },
    error = function(e) {
      cli::cli_abort("无法提取表达矩阵: {e$message}")
    }
  )

  # 转换为稀疏矩阵节省内存
  if (!inherits(expr_data, "sparseMatrix")) {
    cli::cli_progress_step("将稠密矩阵转换为稀疏矩阵...")
    expr_data <- as(expr_data, "sparseMatrix")
  }

  cli::cli_alert_info(
    "表达矩阵维度: {nrow(expr_data)} 个基因 x {ncol(expr_data)} 个样本"
  )

  # 3. 数据预处理 ------------------------------------------------------------
  cli::cli_h2("步骤3: 数据预处理")

  # 并行过滤低表达基因
  cli::cli_progress_step("过滤低表达基因...")
  keep_genes <- Matrix::rowSums(expr_data > 0) >= min.cells
  expr_data <- expr_data[keep_genes, ]
  cli::cli_alert_success(
    "过滤后保留 {sum(keep_genes)}/{length(keep_genes)} 个基因"
  )

  # 转置矩阵
  cli::cli_progress_step("转置矩阵...")
  expr_data <- Matrix::t(expr_data)

  # 4. 创建Seurat对象 -------------------------------------------------------
  cli::cli_h2("步骤4: 创建Seurat对象")

  cli::cli_progress_step("初始化Seurat对象...")
  seurat_obj <- tryCatch(
    {
      CreateSeuratObject(
        counts = expr_data,
        project = geo_accession,
        min.cells = min.cells,
        min.features = min.features
      )
    },
    error = function(e) {
      cli::cli_abort("创建Seurat对象失败: {e$message}")
    }
  )

  # 添加元数据
  pdata <- pData(gse[[1]])
  if (ncol(pdata) > 0) {
    cli::cli_progress_step("添加样本元数据...")
    seurat_obj <- AddMetaData(seurat_obj, metadata = pdata)
  }

  # 添加基因注释
  fdata <- fData(gse[[1]])
  if (!is.null(fdata)) {
    cli::cli_progress_step("添加基因注释...")
    seurat_obj[["RNA"]]@meta.features <- cbind(
      seurat_obj[["RNA"]]@meta.features,
      fdata[match(rownames(seurat_obj), rownames(fdata)), ]
    )
  }

  # 5. 质量控制 -------------------------------------------------------------
  if (qc.filter) {
    cli::cli_h2("步骤5: 质量控制")

    # 计算线粒体比例
    if (!is.null(mt.pattern)) {
      cli::cli_progress_step("计算线粒体基因比例...")
      seurat_obj[["percent.mt"]] <- PercentageFeatureSet(
        seurat_obj,
        pattern = mt.pattern
      )
    }

    # 可视化QC指标
    qc_plot <- VlnPlot(
      seurat_obj,
      features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
      ncol = 3,
      pt.size = 0.1
    )
    print(qc_plot)

    # 交互式QC过滤
    cli::cli_alert("当前细胞数: {ncol(seurat_obj)}")
    filter_params <- list(
      nFeature_RNA = c(min.features, Inf),
      percent.mt = c(0, 100)
    )

    cli::cli_text("建议QC过滤参数:")
    cli::cli_ul()
    cli::cli_li("nFeature_RNA: 200-2500")
    cli::cli_li("percent.mt: <5-10%")
    cli::cli_end()

    # 应用过滤
    cli::cli_progress_step("应用QC过滤...")
    seurat_obj <- subset(
      seurat_obj,
      subset = nFeature_RNA > filter_params$nFeature_RNA[1] &
        nFeature_RNA < filter_params$nFeature_RNA[2] &
        percent.mt > filter_params$percent.mt[1] &
        percent.mt < filter_params$percent.mt[2]
    )
    cli::cli_alert_success("过滤后保留 {ncol(seurat_obj)} 个细胞")
  }

  # 6. 数据标准化 -----------------------------------------------------------
  cli::cli_h2("步骤6: 数据标准化")

  cli::cli_progress_step("执行LogNormalize标准化...")
  seurat_obj <- NormalizeData(
    seurat_obj,
    normalization.method = "LogNormalize",
    scale.factor = 10000,
    verbose = FALSE
  )

  # 7. 保存结果 -------------------------------------------------------------
  if (!is.null(save.rds)) {
    cli::cli_h2("步骤7: 保存结果")

    cli::cli_progress_step("保存Seurat对象到 {.file {save.rds}}...")
    saveRDS(seurat_obj, file = save.rds)
    cli::cli_alert_success("保存成功")
  }

  # 返回结果
  cli::cli_alert_success("转换完成!")
  cli::cli_alert("最终对象信息:")
  print(seurat_obj)

  return(seurat_obj)
}
