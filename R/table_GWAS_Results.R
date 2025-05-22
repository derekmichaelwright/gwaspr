#' table_GWAS_Results
#'
#' Create a table of significant GWAS results.
#' @param folder Folder containing GWAS results.
#' @param files The files to read.
#' @param threshold Significant threshold.
#' @param sug.threshold Suggestive threshold.
#' @param nrowstoread Number of rows to read.
#' @param useHBPvalues Logical, should H.B.P.Values be uses.
#' @return A table of significant GWAS results.
#' @export

table_GWAS_Results <- function(
    folder = "GWAS_Results/",
    files = list_Result_Files(folder),
    nrowstoread = 1000,
    threshold = 6,
    sug.threshold = NULL,
    useHBPvalues = F) {
  #
  output <- NULL
  for(i in files){
    trait <- substr(i, gregexpr("GWAS_Results", i)[[1]][1]+13,
                    gregexpr(".csv", i)[[1]][1]-1 )
    model <- substr(trait, 1, gregexpr("\\.", trait)[[1]][1]-1 )
    trait <- substr(trait, gregexpr("\\.", trait)[[1]][1]+1, nchar(trait) )
    #
    oi <- read.csv(paste0(folder, i), nrows = nrowstoread) %>%
      mutate(Model = model, Trait = trait,
             negLog10_P = -log10(P.value),
             negLog10_HBP = -log10(H.B.P.Value))
    if(useHBPvalues == T) {
      oi <- oi %>% mutate(pvals = negLog10_HBP)
    } else {
      oi <- oi %>% mutate(pvals = negLog10_P)
    }
    ####
    oi <- oi %>%
      mutate(Threshold = ifelse(pvals <
                                  threshold, "Suggestive", "Significant"))
    if (!is.null(sug.threshold)) {
      oi <- oi %>% filter(pvals >= sug.threshold)
    }
    else {
      oi <- oi %>% filter(pvals > threshold)
    }
    #
    output <- bind_rows(output, oi)
  }
  #
  if(sum(colnames(output)=="nobs")>0) { output <- select(output, -nobs) }
  #
  #
  output %>% arrange(desc(negLog10_P)) %>% select(-pvals) #%>%
    #select(SNP, Chromosome, Position, Model, Trait, P.value, `-log10(p)`, MAF, Threshold,
     #      FDR_Adjusted_P.values, effect, Rsquare.of.Model.without.SNP, Rsquare.of.Model.with.SNP)
}
