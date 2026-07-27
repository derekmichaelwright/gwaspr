#' mv
#'
#' A copy of the plyr::mapvalues function.
#' @param x The factor or vector to modify.
#' @param from A vector of the items to replace.
#' @param to A vector of replacement values.
#' @param warn_missing print a message if any of the old values are not actually present in x.
#' @return plyr::mapvalues output
#' @export

mv <- function (x, from, to, warn_missing = TRUE) {
  plyr::mapvalues(x = x, from = from, to = to, warn_missing = warn_missing)
}
