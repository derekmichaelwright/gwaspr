#' table_GWAS_Results
#'
#' Create a table of significant GWAS results.
#' @param folder Folder containing GWAS results.
#' @param fnames The files to read.
#' @param threshold Significant threshold.
#' @param sug.threshold Suggestive threshold.
#' @param nrowstoread Number of rows to read.
#' @param useHBPvalues Logical, should H.B.P.Values be uses.
#' @return A table of significant GWAS results.
#' @export

table_GWAS_Results <- function(
    folder = "GWAS_Results/",
    fnames = list_Result_Files(folder),
    nrowstoread = 1000,
    threshold = 6,
    sug.threshold = NULL) {
  #
  output <- NULL
  for(i in fnames){
    trait <- substr(i, gregexpr("GWAS_Results", i)[[1]][1]+13,
                    gregexpr(".csv", i)[[1]][1]-1 )
    model <- substr(trait, 1, gregexpr("\\.", trait)[[1]][1]-1 )
    trait <- substr(trait, gregexpr("\\.", trait)[[1]][1]+1, nchar(trait) )
    #
    oi <- read.csv(paste0(folder, i), nrows = nrowstoread) %>%
      mutate(Model = model, Trait = trait,
             negLog10_P = -log10(P.value),
             negLog10_HBP = -log10(H.B.P.Value))
    #
    oi <- oi %>% mutate(Threshold = ifelse(negLog10_P < threshold, "Suggestive", "Significant"))
    oi <- oi %>% mutate(H.B.Threshold = ifelse(negLog10_HBP < 1.3, "Not Significant", "Significant"))
    #
    if (!is.null(sug.threshold)) {
      oi <- oi %>% filter(negLog10_P >= sug.threshold)
    }
    else {
      oi <- oi %>% filter(negLog10_P > threshold)
    }
    #
    output <- bind_rows(output, oi)
  }
  #
  if(sum(colnames(output)=="nobs")>0) { output <- select(output, -nobs) }
  #
  output %>% arrange(desc(negLog10_P))
}
