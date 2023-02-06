#' old_table_GWAS_Results
#'
#' Create a table of significant GWAS results. Note: this function is meant for the old version of GAPIT.
#' @param folder Folder containing GWAS results.
#' @param files The files to read.
#' @param threshold Significant threshold.
#' @param sug.threshold Suggestive threshold.
#' @param rowread Number of rows to read.
#' @return A table of significant GWAS results.
#' @export

old_table_GWAS_Results <- function(folder, files,
                               threshold = 6, sug.threshold = NULL,
                               rowread = 2000) {
  #
  output <- NULL
  for(i in files){
    trait <- substr(i, gregexpr("GAPIT.", i)[[1]][1]+6,
                    gregexpr(".GWAS.Results.csv", i)[[1]][1]-1 )
    model <- substr(trait, 1, gregexpr("\\.", trait)[[1]][1]-1 )
    trait <- substr(trait, gregexpr("\\.", trait)[[1]][1]+1, nchar(trait) )
    #
    oi <- read.csv(paste0(folder, i), nrows = rowread) %>%
      mutate(`-log10(p)` = -log10(P.value),
             Model = model,
             Trait = trait,
             Threshold = ifelse(`-log10(p)` < threshold, "Suggestive", "Significant"))
    output <- bind_rows(output, oi)
  }
  #
  if(sum(colnames(output)=="nobs")>0) { output <- select(output, -nobs) }
  #
  if(!is.null(sug.threshold)) {
    output <- output %>% filter(`-log10(p)` >= sug.threshold)
  } else{ output <- output %>% filter(`-log10(p)` > threshold) }
  #
  #colnames(output)[c(2:3,5)] <- c("Chromosome","Position","MAF")
  output %>% arrange(desc(`-log10(p)`)) %>%
    select(SNP, Chromosome, Position, Model, Trait, P.value, `-log10(p)`, maf, Threshold,
           FDR_Adjusted_P.values, effect, Rsquare.of.Model.without.SNP, Rsquare.of.Model.with.SNP)
}
