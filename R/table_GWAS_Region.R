#' table_GWAS_Region
#'
#' Create a table of significant GWAS results.
#' @param folder Folder containing GWAS results.
#' @param chr Chromosome to plot.
#' @param pos1 Start position on chromosome.
#' @param pos2 End position on chromosome.
#' @param fnames The files to read.
#' @param threshold Significant threshold.
#' @param sug.threshold Suggestive threshold.
#' @param nrowstoread Number of rows to read.
#' @param useHBPvalues Logical, if TRUE, H.B.P.Values will be uses.
#' @param skyline Which skyline type to use. Can be "NYC" or "Kansas". If NULL, it will use the highest P.value.
#' @param models GWAS models to use.
#' @return A table of significant GWAS results.
#' @export

table_GWAS_Region <- function(
    folder = "GWAS_Results/",
    chr = 1,
    pos1 = 350000000,
    pos2 = 370000000,
    fnames = list_Result_Files(folder),
    threshold = 6,
    sug.threshold = NULL,
    nrowstoread = 1000,
    skyline = "NYC",
    models = c("MLM", "MLMM", "FarmCPU", "BLINK", "GLM", "CMLM", "SUPER")
    ) {
  #
  xx <- NULL
  #
  for(i in fnames){
    trait <- substr(i, gregexpr("GWAS_Results", i)[[1]][1]+13,
                    gregexpr(".csv", i)[[1]][1]-1 )
    trait <- gsub("\\(Kansas\\)|\\(NYC\\)", "", trait)
    model <- substr(trait, 1, gregexpr("\\.", trait)[[1]][1]-1 )
    trait <- substr(trait, gregexpr("\\.", trait)[[1]][1]+1, nchar(trait) )
    sky <- substr(i, gregexpr("\\(", i)[[1]][1] + 1, gregexpr("\\)", i)[[1]][1] - 1)
    #
    xi <- read.csv(paste0(folder, i), nrows = nrowstoread) %>%
      filter(Chr == chr, Pos > pos1, Pos < pos2) %>%
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
    if(nrow(xi) > 0) { xx <- bind_rows(xx, xi) }
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
    filter(Model %in% models)
  #
  xx <- xx %>% group_by(Trait) %>%
    mutate(Scaled_P = scales::rescale(negLog10_P))
  xx %>% arrange(desc(negLog10_P))
  #
  xi <- read.csv(paste0(folder, i)) %>%
    filter(Chr == chr, Pos > pos1, Pos < pos2)
  #
  ggplot(xi, aes(x = Pos / 100000000)) +
    geom_blank() +
    geom_point(data = xx, aes(y = Trait, size = Scaled_P, alpha = negLog10_P)) +
    facet_grid(. ~ Chr) +
    theme_gwaspr() +
    labs(y = "log10", x = "100 Mbp")
}


#folder = "GWAS_Results/"; fnames = list_Result_Files(folder); nrowstoread = 1000
#threshold = 6; sug.threshold = NULL; removeKansas = T
