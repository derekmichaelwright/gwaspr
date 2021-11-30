#' list_Top_Markers
#'
#' Finds the markers with the highest association on each chromosome.
#' @param x data.
#' @param chroms Chromosomes to include.
#' @param n Number per chromosome.
#' @param threshold_filter Logical, filters non-significant results.
#' @return Table of top results.
#' @export

list_Top_Markers <- function(x, chroms = 1:7, n = 1, threshold_filter = T) {
  #
  threshold <- -log10(0.05 / (nrow(x)) )
  #
  x <- x %>% filter(Chromosome %in% chroms) %>%
    group_by(Chromosome) %>%
    top_n(., n = n, -log10(P.value)) %>%
    mutate(NegLog10 = round(-log10(P.value),2)) %>%
    arrange(Chromosome, rev(NegLog10)) %>%
    as.data.frame() %>%
    select(SNP, CHR=Chromosome, POS=Position, NegLog10)
  #
  if(threshold_filter == T) { x <- x %>% filter(NegLog10 > threshold) }
  #
  x
}
