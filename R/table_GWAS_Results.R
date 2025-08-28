#' table_GWAS_Results
#'
#' Create a table of significant GWAS results.
#' @param folder Folder containing GWAS results.
#' @param fnames The files to read.
#' @param threshold Significant threshold.
#' @param sug.threshold Suggestive threshold.
#' @param nrowstoread Number of rows to read.
#' @param useHBPvalues Logical, if TRUE, H.B.P.Values will be uses.
#' @param skyline Which skyline type to use. Can be "NYC" or "Kansas". If left NULL, it will use the highest P.value.
#' @return A table of significant GWAS results.
#' @export

table_GWAS_Results <- function(
    folder = "GWAS_Results/",
    fnames = list_Result_Files(folder),
    nrowstoread = 1000,
    threshold = 6,
    sug.threshold = NULL,
    skyline = NULL
    ) {
  #
  #if(removeKansas == T) { fnames <- fnames[!grepl("Kansas", fnames)] }
  #
  output <- NULL
  #
  for(i in fnames){
    trait <- substr(i, gregexpr("GWAS_Results", i)[[1]][1]+13,
                    gregexpr(".csv", i)[[1]][1]-1 )
    trait <- gsub("\\(Kansas\\)|\\(NYC\\)", "", trait)
    model <- substr(trait, 1, gregexpr("\\.", trait)[[1]][1]-1 )
    trait <- substr(trait, gregexpr("\\.", trait)[[1]][1]+1, nchar(trait) )
    sky <- substr(i, gregexpr("\\(", i)[[1]][1] + 1, gregexpr("\\)", i)[[1]][1] - 1)
    #
    oi <- read.csv(paste0(folder, i), nrows = nrowstoread) %>%
      mutate(Model = model, Type = sky, Trait = trait,
             negLog10_P = -log10(P.value),
             negLog10_HBP = -log10(H.B.P.Value))
    #
    oi <- oi %>% mutate(Threshold = ifelse(negLog10_P < threshold, "Suggestive", "Significant"))
    #
    if (!is.null(sug.threshold)) {
      oi <- oi %>% filter(negLog10_P >= sug.threshold)
    } else {
      oi <- oi %>% filter(negLog10_P > threshold)
    }
    #
    output <- bind_rows(output, oi)
  }
  #
  if(sum(colnames(output)=="nobs")>0) { output <- dplyr::select(output, -nobs) }
  #
  if(!is.null(skyline)) {
    if(skyline == "NYC")    { output <- output %>% filter(!paste(Model, Type) %in% c("FarmCPU Kansas", "BLINK Kansas")) }
    if(skyline == "Kansas") { output <- output %>% filter(!paste(Model, Type) %in% c("FarmCPU NYC", "BLINK NYC")) }
  }
  #
  output  <- output %>% arrange(desc(P.value)) %>% filter(!duplicated(paste(SNP, Model, P.value)))
  #
  output %>% arrange(desc(negLog10_P))
}

#folder = "GWAS_Results/"; fnames = list_Result_Files(folder); nrowstoread = 1000
#threshold = 6; sug.threshold = NULL; removeKansas = T
