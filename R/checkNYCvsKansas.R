#' checkTaxaNames
#'
#' Check to see if the names in your myY and myG file match.
#' @param xG GWAS genotype object. Note:  needs to be in hapmap format and read in with header = F.
#' @param xY GWAS phenotype object.
#' @return A list of names not present in both datasets.
#' @export

checkTaxaNames <- function(xG = myG, xY = myY) {
  #
  print("The following taxa are not present in your genotype data")
  #
  print( setdiff(xY$Name, t(xG[1,])) ) # should be character(0)
  #
  print("The following taxa are not present in your genotype data")
  #
  print( setdiff(t(xG[1,12:ncol(xG)]), xY$Name) ) # should show just the first 11 columns of myG
}

#i <- myTraits[1]
#folder = "vignettes/GWAS_Results/"; deleteKansas = F
