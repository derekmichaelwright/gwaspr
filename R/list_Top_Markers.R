#' list_Top_Markers
#'
#' Finds the markers with the highest association on each chromosome.
#' @param x data.
#' @param chroms Chromosomes to include.
#' @param n Number per chromosome.
#' @param threshold_filter Logical, filters non-significant results.
#' @return Table of top results.
#' @export

list_Top_Markers <- function(x, chroms = 1:7, n = 1, threshold = NULL) {
  #
  threshold <- -log10(0.05 / (nrow(x)) )
  #
  x <- x %>% filter(Chromosome %in% chroms) %>%
    group_by(Chromosome) %>%
    top_n(., n = n, -log10(P.value)) %>%
    mutate(`-log10(p)` = round(-log10(P.value),2)) %>%
    arrange(Chromosome, rev(`-log10(p)`)) %>%
    as.data.frame() %>%
    select(SNP, CHR=Chromosome, POS=Position, `-log10(p)`)
  #
  if(!is.null(threshold)) { x <- x %>% filter(`-log10(p)` > threshold) }
  #
  x
}
