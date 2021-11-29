#' list_Results
#'
#' Get list of GWAS results.
#' @param folder Folder containing GWAS results.
#' @return List of GWAS results.
#' @export

list_Results <- function(folder = "C:/gitfolder/gwas_tutorial/Results", fullfilename = F) {
  #
  fnames <- grep(".GWAS.Results", list.files(folder))
  list.files(folder)[fnames]
}
