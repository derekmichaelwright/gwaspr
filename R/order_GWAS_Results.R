#' order_GWAS_Results
#'
#' Read in GWAS results, then arrange by P.value and write as csv..
#' @param folder Folder containing GWAS results.
#' @param files The files to read.
#' @return Rearranged files of GWAS results.
#' @export

order_GWAS_Results <- function(
    folder = "GWAS_Results/",
    files = list_Result_Files(folder)) {
  #
  x1 <- is_ran(folder = folder) %>% dropNAcol()
  x2 <- is_ordered(folder = folder)
  x2 <- x2[,colnames(x1)]
  x2 <- x2[rowSums(is.na(x2)) > 0,]
  #
  if(nrow(x2) > 0) {
    files <- files[grepl(x2$Trait, files)]
    #
    for(i in files){
      xx <- read.csv(paste0(folder, i)) %>%
        arrange(P.value)
      write.csv(xx, paste0(folder, i), row.names = F)
    }
  }
  message("GWAS Results files are ordered")
}
