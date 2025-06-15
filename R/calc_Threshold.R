#' calc_Threshold.
#'
#' calculate thresholds.
#' @param myG Genotype file
#' @param p desired p value
#' @return GWAS threshold
#' @export

calc_Threshold <- function(myG, p = 0.05) {
  -log10( p / (nrow(myG)) )
}
