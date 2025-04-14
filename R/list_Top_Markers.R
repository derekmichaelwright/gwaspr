#' list_Top_Markers
#'
#' Finds the markers with the highest association on each chromosome.
#' @param folder Folder containing GWAS results.
#' @param traits GWAS trait.
#' @param models GWAS models to include.
#' @param chroms Chromosomes to include.
#' @param n Number per trait and model.
#' @param threshold filters results with -log10(p) below threshold.
#' @return Table of top results.
#' @export

list_Top_Markers <- function(
    folder = "GWAS_Results/",
    traits = list_Traits(folder),
    models = c("MLM", "FarmCPU", "BLINK", "MLMM", "GLM", "CMLM", "SUPER"),
    threshold = 5,
    chroms = 1:50,
    n = 3) {
  #
  fnames <- list_Result_Files(folder)
  #
  xx <- NULL
  for(i in fnames) {
    trait <- substr(i, gregexpr("GWAS_Results", i)[[1]][1]+13,
                    gregexpr(".csv", i)[[1]][1]-1 )
    model <- substr(trait, 1, gregexpr("\\.", trait)[[1]][1]-1 )
    trait <- substr(trait, gregexpr("\\.", trait)[[1]][1]+1, nchar(trait) )
    #
    if(model %in% models & trait %in% traits) {
      xi <- read.csv(paste0(folder, i)) %>%
        filter(Chr %in% chroms) %>%
        mutate(`-log10(p)` = -log10(P.value)) %>%
        group_by(Chr) %>%
        filter(`-log10(p)` >= threshold) %>%
        slice(1:n) %>%
        ungroup() %>%
        mutate(Trait = trait, Model = model)
      #
      xx <- bind_rows(xx, xi)
    }
  }
  xx <- xx %>% mutate(Traits = NA, Models = NA)
  for(i in 1:nrow(xx)) {
    xi <- xx %>% filter(SNP == xx$SNP[i])
    xx$Traits[i] <- paste(unique(xi$Trait), collapse = "; ")
    xx$Models[i] <- paste(unique(xi$Model), collapse = "; ")
  }
  #
  xx %>% arrange(desc(`-log10(p)`)) %>%
    filter(!duplicated(SNP)) %>%
    select(SNP, Chr, Pos, Traits, Models, Max_LogP = `-log10(p)`)
}
