#' is_ran
#'
#' Create a table of which models have been run on each trait.
#' @param folder Folder containing GWAS results.
#' @return A table of which models have been run on each trait.
#' @export

is_ran <- function(folder = "GWAS_Results/") {
  #
  yy <- data.frame(Trait = list_Traits(folder),
                   GLM = NA, MLM = NA, MLMM = NA, 
                   FarmCPU = NA, BLINK = NA, CMLM = NA, SUPER = NA)
  xx <- list_Result_Files(folder)
  #
  for(i in xx) {
    myTrait <- gsub("GAPIT.Association.GWAS_Results.|MLM.|MLMM.|FarmCPU.|BLINK.|GLM.|CMLM.|SUPER.|\\(NYC\\)|\\(Kansas\\)|.csv", "", i)
    myMod <- gsub("GAPIT.Association.GWAS_Results.|.csv", "", i)
    myMod <- substr(myMod, 1, regexpr("\\.", myMod)-1)
    #
    yy[yy$Trait==myTrait,colnames(yy)[colnames(yy)==myMod]] <- "X"
  }
  #
  yy
}
