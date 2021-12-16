#' list_Top_Markers
#'
#' Finds the markers with the highest association on each chromosome.
#' @param trait GWAS trait.
#' @param method GWAS method.
#' @param chroms Chromosomes to include.
#' @param n Number per chromosome.
#' @param threshold_filter Logical, filters non-significant results.
#' @return Table of top results.
#' @export

list_Top_Markers <- function(trait, model, folder = NULL, threshold = NULL, chroms = 1:50, n = 1) {
  #
  fname <- grepl(".GWAS.Results", list.files(folder)) &
    grepl(trait, list.files(folder)) &
    grepl(paste0("\\.",model,"\\."), list.files(folder))
  fname <- list.files(folder)[fname]
  fname
  x <- read.csv(paste0(folder, fname))
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

#setwd("C:/gitfolder/gwas_tutorial")
#library(gwaspr)
#trait <-"DTF_Nepal_2017"
#method <- "MLM"
#folder <- "GWAS_Results/"
#folder <- "Results/"
#files <- list_Result_Files(folder)
#files <- files[!grepl("GLM", files)]
#traits <- unique(gsub("GAPIT.|GLM.|MLM.|CMLM.|MLMM.|SUPER.|FarmCPU.|Blink.|.GWAS.Results.csv", "", files))[c(1,3)]

#models <- c("MLM","MLMM","FarmCPU","Blink")
