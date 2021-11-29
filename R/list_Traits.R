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
  myT <- NULL
  # i <- fnames[3]
  for(i in fnames) {
    myTi <- substr(i, gregexpr("\\.", i)[[1]][2]+1, gregexpr("\\.", i)[[1]][3]-1)
    myT <- c(myT, myTi)
  }
  #
  unique(myT)
}
