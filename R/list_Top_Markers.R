#' list_Top_Markers
#'
#' Finds the markers with the highest association on each chromosome.
#' @param folder Folder containing GWAS results.
#' @param traits GWAS trait.
#' @param chr Chromosomes to include.
#' @param models GWAS models to include.
#' @param threshold filters results with -log10(p) below threshold.
#' @param n Number per trait and model.
#' @return Table of top results.
#' @export

list_Top_Markers <- function(
    folder = "GWAS_Results/",
    traits = list_Traits(folder),
    chr = 1:50,
    models = c("MLM", "FarmCPU", "BLINK", "MLMM", "GLM", "CMLM", "SUPER"),
    threshold = 5,
    n = 5) {
  #
  fnames <- list_Result_Files(folder)
  #
  xx <- NULL
  for(i in fnames) {
    trait <- substr(i, gregexpr("GWAS_Results", i)[[1]][1]+13,
                    gregexpr(".csv", i)[[1]][1]-1 )
    trait <- gsub("\\(Kansas\\)|\\(NYC\\)", "", trait)
    model <- substr(trait, 1, gregexpr("\\.", trait)[[1]][1]-1 )
    trait <- substr(trait, gregexpr("\\.", trait)[[1]][1]+1, nchar(trait) )
    sky <- substr(i, gregexpr("\\(", i)[[1]][1] + 1, gregexpr("\\)", i)[[1]][1] - 1)
    #
    if(model %in% models & trait %in% traits) {
      xi <- read.csv(paste0(folder, i)) %>%
        filter(Chr %in% chr) %>%
        mutate(`-log10(p)` = -log10(P.value)) %>%
        group_by(Chr) %>%
        filter(`-log10(p)` >= threshold) %>%
        slice(1:n) %>%
        ungroup() %>%
        mutate(Trait = trait, Model = model, Type = sky)
      #
      xx <- bind_rows(xx, xi)
    }
  }
  xx <- xx %>% mutate(Traits = NA, Models = NA, Type = NA)
  for(i in 1:nrow(xx)) {
    xi <- xx %>% filter(SNP == xx$SNP[i])
    xx$Traits[i] <- paste(unique(xi$Trait), collapse = "; ")
    xx$Models[i] <- paste(unique(xi$Model), collapse = "; ")
  }
  #
  xx %>% filter(P.value > 0) %>%
    arrange(desc(`-log10(p)`)) %>%
    filter(!duplicated(SNP)) %>%
    select(SNP, Chr, Pos, Traits, Models, Max_LogP = `-log10(p)`)
}
