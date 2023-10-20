#' table_Results_Summary
#'
#' Create a table of which models have been run on each trait.
#' @param folder Folder containing GWAS results.
#' @param isOrdered set to True to see if results files are ordered based on P.value
#' @return A table of which models have been run on each trait.
#' @export

table_Results_Summary <- function(folder = "GWAS_Results/", isOrdered = F) {
  yy <- data.frame(Trait = list_Traits(folder),
                   MLM = NA, MLMM = NA, FarmCPU = NA, BLINK = NA)
  xx <- list_Result_Files(folder)
  #i<-xx[303]
  for(i in xx) {
    myTrait <- gsub("GAPIT.Association.GWAS_Results.|MLM.|MLMM.|FarmCPU.|BLINK.|.csv", "", i)
    myMod <- gsub("GAPIT.Association.GWAS_Results.|.csv", "", i)
    myMod <- substr(myMod, 1, regexpr("\\.", myMod)-1)
    yy[yy$Trait==myTrait,colnames(yy)[colnames(yy)==myMod]] <- "X"
  }
  if(isOrdered) {
    yy <- yy %>% mutate(MLM = NA, MLMM = NA, FarmCPU = NA, BLINK = NA)
    for(i in xx) {
      myTrait <- gsub("GAPIT.Association.GWAS_Results.|MLM.|MLMM.|FarmCPU.|BLINK.|.csv", "", i)
      myMod <- gsub("GAPIT.Association.GWAS_Results.|.csv", "", i)
      myMod <- substr(myMod, 1, regexpr("\\.", myMod)-1)
      #
      xi <- read.csv(paste0(folder, i), nrows = 100)
      xy <- xi %>% arrange(P.value)
      myLog <- sum(!xi$SNP == xy$SNP) == 0
      if(myLog) {
        yy[yy$Trait==myTrait,colnames(yy)[colnames(yy)==myMod]] <- "X"
      }
    }
  }
  yy
}
