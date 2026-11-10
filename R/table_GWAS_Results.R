#' table_GWAS_Results
#'
#' Create a table of significant GWAS results.
#' @param folder Folder containing GWAS results.
#' @param threshold Significant threshold.
#' @param sug.threshold Suggestive threshold.
#' @param nrows Number of rows to read.
#' @param skyline Which skyline type to use. Can be "NYC" or "Kansas". If NULL, it will use the highest P.value.
#' @param models GWAS models to use.
#' @return A table of significant GWAS results.
#' @export

table_GWAS_Results <- function(
    folder = "GWAS_Results/",
    traits = list_Traits(folder),
    threshold = 6,
    sug.threshold = NULL,
    nrows = 10000,
    skyline = "NYC",
    models = c("MLM", "MLMM", "FarmCPU", "BLINK", "GLM", "CMLM", "SUPER")
    ) {
  #
  xx <- NULL
  #
  fnames <- list_Result_Files(folder)
  fnames <- fnames[grepl(paste(traits,collapse="|"), fnames)]
  for(i in fnames){
    trait <- substr(i, gregexpr("GWAS_Results", i)[[1]][1]+13,
                    gregexpr(".csv", i)[[1]][1]-1 )
    trait <- gsub("\\(Kansas\\)|\\(NYC\\)", "", trait)
    model <- substr(trait, 1, gregexpr("\\.", trait)[[1]][1]-1 )
    trait <- substr(trait, gregexpr("\\.", trait)[[1]][1]+1, nchar(trait) )
    sky <- substr(i, gregexpr("\\(", i)[[1]][1] + 1, gregexpr("\\)", i)[[1]][1] - 1)
    #
    xi <- read.csv(paste0(folder, i), nrows = nrows) %>%
      mutate(Model = model, Type = sky, Trait = trait,
             negLog10_P = -log10(P.value),
             negLog10_HBP = -log10(H.B.P.Value))
    #
    xi <- xi %>% mutate(Threshold = ifelse(negLog10_P < threshold, "Suggestive", "Significant"))
    #
    if (!is.null(sug.threshold)) {
      xi <- xi %>% filter(negLog10_P >= sug.threshold)
    } else {
      xi <- xi %>% filter(negLog10_P > threshold)
    }
    #
    xx <- bind_rows(xx, xi)
  }
  #
  if(sum(colnames(xx)=="nobs")>0) { xx <- dplyr::select(xx, -nobs) }
  #
  if(!is.null(skyline)) {
    if(skyline == "NYC")    { xx <- xx %>% filter(!paste(Model, Type) %in% c("FarmCPU Kansas", "BLINK Kansas")) }
    if(skyline == "Kansas") { xx <- xx %>% filter(!paste(Model, Type) %in% c("FarmCPU NYC", "BLINK NYC")) }
  }
  #
  xx  <- xx %>%
    arrange(desc(negLog10_P)) %>%
    filter(!duplicated(paste(SNP, Model, P.value))) %>%
    filter(Model %in% models) %>%
    select(SNP, Chr, Pos, negLog10_P, Trait, Model, Threshold, everything()) %>%
    arrange(desc(negLog10_P))
  #
  xx
}

#folder = "GWAS_Results/"; fnames = list_Result_Files(folder); nrowstoread = 1000
#threshold = 6; sug.threshold = NULL; removeKansas = T
