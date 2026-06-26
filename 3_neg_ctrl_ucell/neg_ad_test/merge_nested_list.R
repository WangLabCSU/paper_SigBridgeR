#' 合并稀疏嵌套列表
#'
#' @param lst 包含子列表的列表 (List of lists)
#' @param is_empty 判断元素是否为空的函数，默认兼容 NULL、空向量和 0行 tbl
#' @return 合并后的单一列表
merge_nested_lists <- function(
  lst,
  is_empty = function(x) {
    is.null(x) ||
      length(x) == 0L ||
      (is.data.frame(x) && nrow(x) == 0L)
  }
) {
  # 1. 边界条件校验
  if (!is.list(lst) || length(lst) == 0L) {
    return(list())
  }

  # 2. 校验子列表是否等长
  lengths <- lengths(lst)
  if (length(unique(lengths)) > 1L) {
    stop("All sub-lists must have the same length.")
  }

  n_slots <- lengths[1L]
  if (n_slots == 0L) {
    return(list())
  }

  # 3. 预分配结果列表 (避免动态扩容带来的内存开销)
  res <- vector("list", n_slots)

  # 4. 保留原始槽位名称 (如果存在)
  slot_names <- names(lst[[1L]])
  if (!is.null(slot_names)) {
    names(res) <- slot_names
  }

  # 5. 遍历槽位与子列表
  for (i in seq_len(n_slots)) {
    for (j in seq_along(lst)) {
      elem <- lst[[j]][[i]]

      # 一旦发现非空元素，立即填入并利用先验知识 break 剪枝
      if (!is_empty(elem)) {
        res[[i]] <- elem
        break
      }
    }
  }

  return(res)
}
