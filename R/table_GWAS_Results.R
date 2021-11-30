#' table_GWAS_Results
#'
#' Create a table of significant GWAS results.
#' @param folder Folder containing GWAS results.
#' @param file The file to read.
#' @param threshold Significant threshold.
#' @param threshold2 Suggestive threshold.
#' @param rowread Number of rows to read.
#' @return A table of significant GWAS results.
#' @export

table_GWAS_Results <- function(folder, file,
                              threshold = 6, threshold2 = NULL,
                              rowread = 2000) {
  #
  trait <- substr(file, gregexpr("GAPIT.", file)[[1]][1]+6,
                  gregexpr(".GWAS.Results.csv", file)[[1]][1]-1 )
  model <- substr(trait, 1, gregexpr("\\.", trait)[[1]][1]-1 )
  trait <- substr(trait, gregexpr("\\.", trait)[[1]][1]+1, nchar(trait) )
  #
  output <- read.csv(paste0(folder, file), nrows = rowread) %>%
    mutate(`-log10(p)` = -log10(P.value))
  #
  if(sum(colnames(output)=="nobs")>0) { output <- select(output, -nobs) }
  #
  if(!is.null(threshold2)) {
    output <- output %>% filter(-log10(P.value) >= threshold2)
  } else{ output <- output %>% filter(-log10(P.value) > threshold) }
  #
  output <- output %>%
    mutate(`-log10(p)` = -log10(P.value),
           Model = model,
           Trait = trait,
           Threshold = ifelse(`-log10(p)` < threshold, "Suggestive", "Significant"))
  output
}
