#' is_ordered
#'
#' Create a table of which GWAS results files are ordered by p.value.
#' @param folder Folder containing GWAS results.
#' @return A table of which models have been run on each trait.
#' @export

is_ordered <- function(folder = "GWAS_Results/") {
  #
  yy <- data.frame(Trait = list_Traits(folder),
                   GLM = NA, MLM = NA, MLMM = NA, 
                   FarmCPU = NA, BLINK = NA, CMLM = NA, SUPER = NA, 
                   FarmCPU_Kansas = NA, BLINK_Kansas = NA)
  xx <- list_Result_Files(folder)
  #
  for(i in xx) {
    isK <- grepl("\\(Kansas\\)", i)
    myTrait <- gsub("GAPIT.Association.GWAS_Results.|MLM.|MLMM.|FarmCPU.|BLINK.|GLM.|CMLM.|SUPER.|\\(NYC\\)|\\(Kansas\\)|.csv", "", i)
    myMod <- gsub("GAPIT.Association.GWAS_Results.|.csv", "", i)
    myMod <- substr(myMod, 1, regexpr("\\.", myMod)-1)
    if(isK == T) { myMod <- paste0(myMod, "_Kansas") }
    #
    xi <- read.csv(paste0(folder, i), nrows = 100)
    xy <- xi %>% arrange(P.value)
    myLog <- sum(!xi$SNP == xy$SNP) == 0
    if(myLog) { yy[yy$Trait==myTrait,colnames(yy)[colnames(yy)==myMod]] <- "X" }
  }
  #
  yy
}
