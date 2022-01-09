#' list_Threshold.
#'
#' calculate thresholds.
#' @param num_markers Number of markers
#' @param p desired p value
#' @return GWAS threshold
#' @export

list_Threshold <- function(num_markers, p = 0.05) {
  -log10( p / (nrow(num_markers)-1) )
}
