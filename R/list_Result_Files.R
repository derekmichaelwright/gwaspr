#' list_Result_Files
#'
#' Get list of GWAS results.
#' @param folder Folder containing GWAS results.
#' @return List of GWAS results.
#' @export

list_Result_Files <- function(folder) {
  fnames <- grep(".GWAS.Results", list.files(folder))
  list.files(folder)[fnames]
}
