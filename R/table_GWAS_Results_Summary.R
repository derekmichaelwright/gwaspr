#' table_GWAS_Results_Summary
#'
#' Create a summary using the output from `table_GWAS_Results()`.
#' @param xx Object from `table_GWAS_Results()`.
#' @param onlySig Logical, if TRUE, any suggested associations will be removed.
#' @param binMarkers Logical, if TRUE, markers will be bined based on `binSize`.
#' @param binSize range on each side of marker to bin. default = 1,000,000.
#' @return A table of significant GWAS results.
#' @export

table_GWAS_Results_Summary <- function(xx, onlySig = F, binMarkers = F, binSize = 1000000) {
  #
  if(onlySig == T) { xx <- xx %>% filter(Threshold == "Significant") }
  xx <- xx %>% arrange(desc(negLog10_P))
  #
  output <- NULL
  #
  if(binMarkers == T) {
    while(nrow(xx) > 0) {
      x1 <- xx %>% slice(1)
      xi <- xx %>% filter(Chr == x1$Chr, Pos >= x1$Pos - binSize, Pos <= x1$Pos + binSize)
      x1 <- x1 %>%
        mutate(Hits = nrow(xi),
               min_P = min(xi$P.value),
               max_P = max(xi$P.value),
               min_HBP = min(xi$H.B.P.Value),
               max_HBP = max(xi$H.B.P.Value),
               MAF = mean(MAF),
               Models = paste(unique(xi$Model), collapse=";"),
               Traits = paste(unique(xi$Trait), collapse=";"),
               min_negLog10_P = -log10(max_P),
               max_negLog10_P = -log10(min_P),
               min_negLog10_HBP = -log10(max_HBP),
               max_negLog10_HBP = -log10(min_HBP) )
      #
      x1 <- x1 %>%
        dplyr::select(SNP, Chr, Pos, Hits, MAF,
               max_negLog10_P, min_negLog10_P,
               Models, Traits,
               max_negLog10_HBP, min_negLog10_HBP,
               min_P, max_P, min_HBP, max_HBP)
      #
      output <- rbind(output, x1)
      #
      xx <- xx %>% filter(!SNP %in% xi$SNP)
    }
  }
  #
  if(binMarkers == F) {
    while(nrow(xx) > 0) {
      x1 <- xx %>% slice(1)
      xi <- xx %>% filter(SNP %in% x1$SNP)
      x1 <- x1 %>%
        mutate(Hits = nrow(xi),
               min_P = min(xi$P.value),
               max_P = max(xi$P.value),
               min_HBP = min(xi$H.B.P.Value),
               max_HBP = max(xi$H.B.P.Value),
               MAF = mean(MAF),
               Models = paste(unique(xi$Model), collapse=";"),
               Traits = paste(unique(xi$Trait), collapse=";"),
               min_negLog10_P = -log10(max_P),
               max_negLog10_P = -log10(min_P),
               min_negLog10_HBP = -log10(max_HBP),
               max_negLog10_HBP = -log10(min_HBP) )
      #
      x1 <- x1 %>%
        dplyr::select(SNP, Chr, Pos, Hits, MAF,
               max_negLog10_P, min_negLog10_P,
               Models, Traits,
               max_negLog10_HBP, min_negLog10_HBP,
               min_P, max_P, min_HBP, max_HBP)
      #
      output <- rbind(output, x1)
      #
      xx <- xx %>% filter(SNP != xi$SNP)
    }
  }
  #
  output %>% arrange(desc(max_negLog10_P))
}

#folder = "GWAS_Results/"; fnames = list_Result_Files(folder); nrowstoread = 1000
#threshold = 6; sug.threshold = NULL
#xx <- table_GWAS_Results("GWAS_Results/", threshold = 6.8, sug.threshold = 5)
#yy<-xx
#xx<-yy
