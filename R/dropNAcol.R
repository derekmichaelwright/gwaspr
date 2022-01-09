#' dropNAcol
#'
#' Removes columns without any values.
#' @param x data.
#' @return data with empty columns removed.
#' @export

dropNAcol <- function(x) {
  x[,colSums(is.na(x)) < nrow(x)]
}
