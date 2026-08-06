#' calc_LD_Decay
#'
#' Calculates pairwise LD for markers within your myG genotype file.
#' @param xG GWAS genotype object. Note: needs to be in hapmap format.
#' @param outputFolder Folder to place RData files into.
#' @param markerNum Number of markers to randomly select for each chromosome.
#' @return RData files in the outputFolder location of LD data.
#' @export

calc_LD_Decay <- function(xG = myG, outputFolder, markerNum = 200 ) {
  #
  dna <- data.frame(stringsAsFactors = F,
                    Symbol = c("A", "C", "G", "T", "U",
                               "R", "Y", "S", "W", "K", "M", "N"),
                    Value  = c("A/A","C/C","G/G","T/T","U/U",
                               "A/G","C/T","G/C","A/T","G/T","A/C","N/N") )
  #
  xx <- xG[,c(-2,-5,-6,-7,-8,-9,-10,-11)]
  #
  for(i in 4:ncol(xx)) {
    xx[xx[,i]=="N", i] <- NA
    xx[,i] <- mv(xx[,i], dna$Symbol, dna$Value)
  }
  #
  LD_Decay <- function(x=xx, folder = outputFolder, Chr = 1, Num = 200) {
    xc <- x %>% filter(chrom == Chr)
    xr <- xc %>% column_to_rownames("rs")
    xr <- xr[,c(-1,-2)]
    # get random row numbers for random marker selection
    rr <- round(runif(Num, 1, nrow(xr)))
    # Create loop to make sure there are no duplicates
    while(sum(duplicated(rr))>0) {
      ra <- round(runif(sum(duplicated(rr)), 1, nrow(xr)))
      rr <- rr[!duplicated(rr)]
      rr <- c(rr, ra)
    }
    rr <- rr[order(rr)]
    # Subset of random markers
    xr <- xr[rr,]
    xi <- xr %>% t() %>% as.data.frame()
    # Calculate LD
    myLD <- genetics::LD(genetics::makeGenotypes(xi))
    save(myLD, file = paste0(folder, "LD_Chrom_", Chr, "_", Num, ".Rdata"))
  }
  # Create a loop to do all chromosomes
  for(i in unique(xG$chrom)) {
    LD_Decay(x = xx, folder = outputFolder, Chr = i, Num = markerNum)
  }
  #
}

#myG <- read.csv("vignettes/gwaspr_myG_hmp.csv", header = T)
#xG = myG; outputFolder = "vignettes/LD_Decay/"; markerNum = 200
