FindDuplicates <- function(data, row_names) {
  dup_flags <- duplicated(row_names) | duplicated(row_names, fromLast = TRUE)
  dup_indices <- which(dup_flags)
  dup_names <- unique(row_names[dup_flags])

  if (length(dup_names) == 0) {
    message("No duplicates found.")
    return(invisible(NULL))
  }

  for (name in dup_names) {
    idx <- which(row_names == name)
    message(
      "\nDuplicate: '",
      name,
      "' | Positions: ",
      paste(idx, collapse = ", ")
    )
    print(data[idx, 1:2], row.names = FALSE) # 只打印前两列，隐藏行号
  }
}
