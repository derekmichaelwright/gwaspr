#' is_Ordered
#'
#' Create a table of which GWAS results files are ordered by p.value.
#' @param folder Folder containing GWAS results.
#' @param nrows Number of rows to read.
#' @return A table of which models have been run on each trait.
#' @export

is_Ordered <- function(folder = "GWAS_Results/", nrows = 10000) {
  #
  yy <- data.frame(Trait = list_Traits(folder),
                   GLM = "", MLM = "", MLMM = "",
                   FarmCPU = "", BLINK = "", CMLM = "", SUPER = "",
                   FarmCPU_Kansas = "", BLINK_Kansas = "")
  xx <- list_Result_Files(folder)
  #
  for(i in xx) {
    isK <- grepl("\\(Kansas\\)", i)
    myTrait <- gsub("GAPIT.Association.GWAS_Results.|MLM.|MLMM.|FarmCPU.|BLINK.|GLM.|CMLM.|SUPER.|\\(NYC\\)|\\(Kansas\\)|.csv", "", i)
    myMod <- gsub("GAPIT.Association.GWAS_Results.|.csv", "", i)
    myMod <- substr(myMod, 1, regexpr("\\.", myMod)-1)
    if(isK == T) { myMod <- paste0(myMod, "_Kansas") }
    #
    xi <- read.csv(paste0(folder, i), nrows = nrows)
    xy <- xi %>% arrange(P.value)
    xy1 <- xy %>% filter(P.value != 0)
    xy2 <- xy %>% filter(P.value == 0)
    xy <- bind_rows(xy1, xy2)
    myLog <- sum(!xi$SNP[1:100] == xy$SNP[1:100]) == 0
    if(myLog) { yy[yy$Trait==myTrait,colnames(yy)[colnames(yy)==myMod]] <- "X" }
  }
  #
  yy
}
#i<-"GAPIT.Association.GWAS_Results.MLMM.TP.DM.Mar_LDP(NYC).csv"
