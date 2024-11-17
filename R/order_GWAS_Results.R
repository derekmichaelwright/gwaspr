#' order_GWAS_Results
#'
#' Read in GWAS results, then arrange by P.value and write as csv..
#' @param folder Folder containing GWAS results.
#' @param files The files to read.
#' @return Rearranged files of GWAS results.
#' @export

order_GWAS_Results <- function(folder, files) {
  #
  for(i in files){
    xx <- read.csv(paste0(folder, i)) %>%
      arrange(P.value)
    write.csv(xx, paste0(folder, i), row.names = F)
  }
}
