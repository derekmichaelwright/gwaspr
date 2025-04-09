#' list_Traits
#'
#' Get a list of traits with GWAS results.
#' @param folder Folder containing GWAS results.
#' @return List of traits with GWAS results.
#' @export

list_Traits <- function(folder = "GWAS_Results/") {
  #
  fnames <- grep(".GWAS.Results", list.files(folder))
  fnames <- list.files(folder)[fnames]
  fnames <- gsub("GAPIT.|GLM.|MLM.|CMLM.|MLMM.|SUPER.|FarmCPU.|Blink.|BLINK.|.GWAS.Results|Association.GWAS_Results.|.csv","",fnames)
  fnames
  unique(fnames)
}
