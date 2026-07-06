#' 在字符向量每个元素中按固定宽度插入换行符
#' @param x character vector
#' @param width 每隔 width 个字符插入一个 \n
#' @param mode 换行模式："hyphen" — 单词中断时加 -；"word" — 尽量按单词边界换行
#' @return 与 x 等长的字符向量，保留了 names 属性
wrap_every_n_chars <- function(x, width = 20, mode = c("word", "hyphen")) {
  mode <- match.arg(mode)

  vapply(
    x,
    \(s) {
      if (is.na(s) || nchar(s) <= width) {
        return(s)
      }

      if (mode == "hyphen") {
        # 固定每 width 个字符断开，加 -
        chars <- strsplit(s, "")[[1L]]
        n <- length(chars)
        g <- ceiling(seq_len(n) / width)
        pieces <- trimws(tapply(chars, g, paste, collapse = ""))
        paste(pieces, collapse = "-\n")
      } else {
        # word 模式：按单词边界换行
        words <- strsplit(s, " ")[[1L]]
        lines <- character()
        current <- ""

        # 辅助函数：把超长单词切成若干段，段尾加 -，并去掉续段前导的 -
        split_long_word <- function(word) {
          chars <- strsplit(word, "")[[1L]]
          n <- length(chars)
          g <- ceiling(seq_len(n) / width)
          pieces <- tapply(chars, g, paste, collapse = "")
          # 续段如果以 - 开头则去掉（原词的连字符已表达断词）
          for (j in seq_along(pieces)[-1L]) {
            pieces[j] <- sub("^-", "", pieces[j])
          }
          list(
            head = paste0(pieces[-length(pieces)], "-"),
            tail = pieces[length(pieces)]
          )
        }

        for (word in words) {
          if (nchar(word) > width) {
            if (nchar(current) > 0) {
              lines <- c(lines, current)
              current <- ""
            }
            spl <- split_long_word(word)
            lines <- c(lines, spl$head)
            current <- spl$tail
            next
          }

          if (nchar(current) == 0) {
            current <- word
          } else if (nchar(current) + 1L + nchar(word) <= width) {
            current <- paste(current, word, sep = " ")
          } else {
            lines <- c(lines, current)
            current <- word
          }
        }
        if (nchar(current) > 0) {
          lines <- c(lines, current)
        }
        paste(lines, collapse = "\n")
      }
    },
    character(1),
    USE.NAMES = FALSE
  )
}
