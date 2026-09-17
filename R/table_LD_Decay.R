#' table_LD_Decay
#'
#' Calculates pairwise LD for markers within your myG genotype file.
#' @param xG GWAS genotype object. Note: needs to be in hapmap format.
#' @param outputFolder Folder to place RData files into.
#' @return RData files in the outputFolder location of LD data.
#' @export

table_LD_Decay <- function(xG = myG, outputFolder ) {
  #
  yy <- NULL
  myFiles <- list.files(outputFolder)
  #i<-"LD_Chrom_1_2000.Rdata"
  for(i in myFiles) {
    myChr <- gsub("LD_Chrom_", "", i)
    myChr <- substr(myChr, 1, regexpr("", myChr)  )
    myNum <- gsub(paste0("LD_Chrom_", myChr,"_"), "", i)
    myNum <- gsub(".Rdata","", myNum)
    xc <- xG %>% filter(chrom == myChr)
    #
    load(paste0(outputFolder, i))
    xi <- myLD$`R^2` %>% as.data.frame() %>%
      rownames_to_column("SNP1") %>%
      gather(SNP2, LD, 2:ncol(.)) %>%
      filter(!is.na(LD)) %>%
      mutate(Chr = myChr,
             SNP1_d = plyr::mapvalues(SNP1, xc$rs, xc$pos, warn_missing = F),
             SNP2_d = plyr::mapvalues(SNP2, xc$rs, xc$pos, warn_missing = F),
             Distance = as.numeric(SNP2_d) - as.numeric(SNP1_d)) %>%
      arrange(Distance, rev(LD))
    xii <- xi %>% filter(Distance < 1000000)
    myloess <- stats::loess(LD ~ Distance, data = xii, span=0.50)
    xi <- xi %>%
      mutate(Moving_Avg = movingAverage(LD, n = 100),
             Loess_Chr = ifelse(Distance < 1000000, predict(myloess), NA) )
    #
    myThresh <- xi %>% filter(Loess_Chr < 0.2)
    myThresh <- min(myThresh$Distance, na.rm = T)
    yi <- data.frame(Chr = myChr, Num = myNum, Threshold = myThresh)
    #
    yy <- bind_rows(yy, yi)
    #
  }
  #
  yy <- yy %>% spread(Num, Threshold)
  yy
}

#outputFolder = "vignettes/LD_Decay/"
