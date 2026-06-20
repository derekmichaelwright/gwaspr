#' run_Summary
#'
#' Create a table of which models have been run on each trait.
#' @param folder Folder containing GWAS results.
#' @return A table of which models have been run on each trait.
#' @export

run_Summary <- function(folder = "GWAS_Results/") {
  #
  yy <- data.frame(Trait = list_Traits(folder),
                   GLM = "", MLM = "", MLMM = "",
                   FarmCPU = "", BLINK = "", CMLM = "", SUPER = "")
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
