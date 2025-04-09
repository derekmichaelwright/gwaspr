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
    n = 1) {
  #
  fname <- grepl("GWAS_Results", list.files(folder)) &
    grepl(traits, list.files(folder)) &
    grepl(paste0("\\.",models,"\\."), list.files(folder))
  fname <- list.files(folder)[fname]
  fname
  x <- read.csv(paste0(folder, fname))
  #
  x <- x %>% filter(Chr %in% chroms) %>%
    group_by(Chr) %>%
    top_n(., n = n, -log10(P.value)) %>%
    mutate(`-log10(p)` = round(-log10(P.value),2)) %>%
    arrange(Chr, rev(`-log10(p)`)) %>%
    as.data.frame() %>%
    select(SNP, Chr, Pos, `-log10(p)`)
  #
  if(!is.null(threshold)) { x <- x %>% filter(`-log10(p)` > threshold) }
  #
  x
}
