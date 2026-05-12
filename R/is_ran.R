#' is_Ran
#'
#' Create a table of which models have been run on each trait.
#' @param myY Phenotype file with traits as column names which will be checked for GWAS results.
#' @param folder Folder containing GWAS results.
#' @return A table of which models have been run on each trait.
#' @export

is_Ran <- function(myY, folder = "GWAS_Results/") {
  myYs <- colnames(myY)[-1]
  myRs <- list_Traits(folder)
  xx <- data.frame(Traits = myYs, GWAS = ifelse(myYs %in% myRs, "X", ""))
  xx
}
