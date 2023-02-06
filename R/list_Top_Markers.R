#' list_Top_Markers
#'
#' Finds the markers with the highest association on each chromosome.
#' @param trait GWAS trait.
#' @param method GWAS method.
#' @param chroms Chromosomes to include.
#' @param n Number per chromosome.
#' @param threshold filters results with -log10(p) below threshold.
#' @return Table of top results.
#' @export

list_Top_Markers <- function(trait, model, folder = NULL, threshold = NULL, chroms = 1:50, n = 1) {
  #
  fname <- grepl("GWAS_Results", list.files(folder)) &
    grepl(trait, list.files(folder)) &
    grepl(paste0("\\.",model,"\\."), list.files(folder))
  fname <- list.files(folder)[fname]
  fname
  x <- read.csv(paste0(folder, fname))
  #colnames(x)[2:3] <- c("Chromosome", "Position")
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
