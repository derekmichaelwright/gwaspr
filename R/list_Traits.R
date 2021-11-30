#' list_Traits
#'
#' Get a list of traits with GWAS results.
#' @param folder Folder containing GWAS results.
#' @return List of traits with GWAS results.
#' @export

list_Traits <- function(folder = "C:/gitfolder/gwas_tutorial/Results", fullfilename = F) {
  #
  fnames <- grep(".GWAS.Results", list.files(folder))
  fnames <- list.files(folder)[fnames]
  fnames <- gsub("GAPIT.|GLM.|MLM.|CMLM.|MLMM.|SUPER.|FarmCPU.|Blink.|.GWAS.Results.csv","",fnames)
  fnames
  unique(fnames)
}
