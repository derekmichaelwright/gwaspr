#' is_Ran
#'
#' Create a table of which models have been run on each trait.
#' @param folder Folder containing GWAS results.
#' @param myY Phenotype file with traits as column names which will be checked for GWAS results.
#' @return A table of which models have been run on each trait.
#' @export

is_Ran <- function(folder = "GWAS_Results/", myY) {
  myYs <- colnames(myY)[-1]
  myRs <- list_Traits(folder)
  xx <- data.frame(Traits = myYs, GWAS = ifelse(myYs %in% myRs, "X", ""))
  xx
}
