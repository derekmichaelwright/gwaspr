#' table_GWAS_Results
#'
#' Create a table of significant GWAS results.
#' @param folder Folder containing GWAS results.
#' @param files The files to read.
#' @param threshold Significant threshold.
#' @param sug.threshold Suggestive threshold.
#' @param nrowstoread Number of rows to read.
#' @return A table of significant GWAS results.
#' @export

table_GWAS_Results <- function(folder, files = list.files(folder), nrowstoread = 1000,
                               threshold = 6, sug.threshold = NULL) {
  #
  output <- NULL
  for(i in files){
    trait <- substr(i, gregexpr("GWAS_Results", i)[[1]][1]+13,
                    gregexpr(".csv", i)[[1]][1]-1 )
    model <- substr(trait, 1, gregexpr("\\.", trait)[[1]][1]-1 )
    trait <- substr(trait, gregexpr("\\.", trait)[[1]][1]+1, nchar(trait) )
    #
    oi <- read.csv(paste0(folder, i), nrows = nrowstoread) %>%
      mutate(`-log10(p)` = -log10(P.value),
             Model = model,
             Trait = trait,
             Threshold = ifelse(`-log10(p)` < threshold, "Suggestive", "Significant"))
    #
    if(!is.null(sug.threshold)) {
      oi <- oi %>% filter(`-log10(p)` >= sug.threshold)
    } else{ oi <- oi %>% filter(`-log10(p)` > threshold) }
    #
    output <- bind_rows(output, oi)
  }
  #
  if(sum(colnames(output)=="nobs")>0) { output <- select(output, -nobs) }
  #
  #
  output %>% arrange(desc(`-log10(p)`)) #%>%
    #select(SNP, Chromosome, Position, Model, Trait, P.value, `-log10(p)`, MAF, Threshold,
     #      FDR_Adjusted_P.values, effect, Rsquare.of.Model.without.SNP, Rsquare.of.Model.with.SNP)
}
