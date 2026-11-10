#' order_GWAS_Results
#'
#' Read in GWAS results, then arrange by P.value and write as csv..
#' @param folder Folder containing GWAS results.
#' @param files The files to read. Default will read all files.
#' @param nrows Number of rows to read.
#' @return Rearranged files of GWAS results.
#' @export

order_GWAS_Results <- function(
    folder = "GWAS_Results/",
    files = list_Result_Files(folder),
    nrows = 10000) {
  #
  xx <- is_Ordered(folder = folder, nrows = nrows)
  xx[xx==""] <- NA
  xx <- xx %>% dropNAcol()
  xx[is.na(xx)] <- ""
  xx <- xx[rowSums(xx=="") > 0,]
  #
  if(nrow(xx) > 0) {
    files <- files[grepl(paste(xx$Trait, collapse="|"), files)]
    #i<-files[1]
    for(i in files){
      xi <- read.csv(paste0(folder, i)) %>%
        arrange(P.value) %>% mutate(pval = -log10(P.value))
      xi1 <- xi %>% filter(P.value != 0)
      xi2 <- xi %>% filter(P.value == 0)
      xi <- bind_rows(xi1, xi2)
      write.csv(xi, paste0(folder, i), row.names = F)
    }
  }
  message("GWAS Results files are ordered")
}
#i<- "GAPIT.Association.GWAS_Results.FarmCPU.TP.Dotted_LDP(Kansas).csv"
