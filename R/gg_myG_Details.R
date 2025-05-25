#' gg_myG_Details
#'
#' Analyses your genotype data and outputs 2 summary files with marker details and 3 plots.
#' @param filename GWAS genotype object. Note: needs to be in hapmap format.
#' @param myPrefix Prefix for file names.
#' @return 2 .csv files with marker details & 3 Marker plots.
#' @export

gg_myG_Details <- function(filname, myPrefix = "") {
  #
  xx <- read.csv(filename, header = T)
  #
  xNames <- colnames(xx)[12:ncol(xx)]
  #
  numAlleles <- (max(nchar(xx$alleles)) + 1) / 2
  xx <- xx %>% separate(alleles, into = paste0("M",1:numAlleles), sep = "/", remove = F)
  #
  xx_AF <- function(x) { sum(x %in% x[1], na.rm = T)-1 }
  #
  for(i in 1:numAlleles) { xx[,paste0("numM",i)] <- apply(xx[,c(paste0("M",i), xNames)], 1, xx_AF) }
  #
  xx$numMajor <- apply(xx[,paste0("numM",1:numAlleles)], 1, max)
  #
  zz <- xx[,paste0("numM",1:numAlleles)]
  zz <- as.data.frame(zz == xx[,"numMajor"])
  #
  for(i in 1:nrow(zz)) {
    for(j in numAlleles:2) {
      if( zz[i,j] == T & sum(zz[i,1:(j-1)]) > 0 ) { zz[i,j] <- F }
    } }
  zz <- as.matrix(zz)
  #
  for(i in 1:nrow(xx)) {
    xx$Major[i] <- xx[i,paste0("M",1:numAlleles)][zz[i,]]
    xx$Minors[i] <- paste(dropNAcol(xx[i,paste0("M",1:numAlleles)][!zz[i,]]), collapse = "/")
  }
  xx$Major <- as.character(xx$Major)
  xx$Minors <- as.character(xx$Minors)
  #
  xx$numTotal <- apply(xx[,paste0("numM",1:numAlleles)], 1, sum)
  #
  xx_Het <- function(x) { sum(x %in% c("R","Y","S","W","K","M"), na.rm = T) }
  xx_Homo <- function(x) { sum(x %in% c("A", "C", "G", "T"), na.rm = T) }
  xx_N <- function(x) { sum(x %in% "N", na.rm = T) }
  #
  xx <- xx %>%
    mutate(numN     = apply(xx[,xNames], 1, xx_N),
           numHet   = apply(xx[,xNames], 1, xx_Het),
           numHomo  = apply(xx[,xNames], 1, xx_Homo),
           numMinors = numTotal - numMajor,
           #
           Het = numHet / (numHet + numHomo),
           MAF = (numMinors*2 + numHet) / (numMinors*2 + numMajor*2 + numHet*2),
           MCF = numN / (numHet + numHomo + numN) )
  xg <- xx %>%
    select(-xNames, -strand, -assembly, -center, -protLSID, -assayLSID, -panel, -QCcode) %>%
    select(rs, chrom, pos, alleles, Major, Minors, everything())
  #
  write.csv(xg, paste0(myPrefix, "_01_myG_Markers.csv"), row.names = F)
  # Plot
  mp1 <- ggplot(xx, aes(x = Het * 100)) +
    geom_histogram(color = "black", fill = "darkgreen", alpha = 0.7) +
    theme_gwaspr() +
    labs(title = "Marker Details (myG)",
         subtitle = "Heterozygosity", x = NULL, y = NULL)
  mp2 <- ggplot(xx, aes(x = MAF * 100)) +
    geom_histogram(color = "black", fill = "darkgreen", alpha = 0.7) +
    scale_x_continuous(breaks = seq(5,100, by = 5),
                       minor_breaks = seq(5,100, by = 5)) +
    theme_gwaspr() +
    labs(subtitle = "Minor Allele Frequency", x = "Percent", y = NULL)
  mp3 <- ggplot(xx, aes(x = MCF * 100)) +
    geom_histogram(color = "black", fill = "darkgreen", alpha = 0.7) +
    theme_gwaspr() +
    labs(subtitle = "Missing Call Frequency", x = NULL, y = NULL)
  mpA <- ggarrange(mp1, mp2, mp3, ncol = 3, nrow = 1, align = "h")
  ggsave(paste0(myPrefix, "_04_myG_Markers.png"), mpA, width = 12, height = 4)
  #
  #
  #
  yy <- data.frame(Name = xNames)
  #
  for(i in xNames) {
    xi <- xx[,i]
    numHet <- sum(xi %in% c("R","Y","S","W","K","M"), na.rm = T)
    numHomo <- sum(xi %in% c("A", "C", "G", "T"), na.rm = T)
    numMiss <- sum(xi %in% "N", na.rm = T)
    numMajor <- sum(xi == xx$Major, na.rm = T)
    numMinors <- sum(grepl(xi, xx$Minors), na.rm = T)
    #
    yy$numMajor[yy$Name == i] <- numMajor
    yy$numMinors[yy$Name == i] <- numMinors
    yy$numHet[yy$Name == i] <- numHet
    yy$numHomo[yy$Name == i] <- numHomo
    yy$numMiss[yy$Name == i] <- numMiss
    #
    yy$Het[yy$Name == i] <- numHet / (numHet + numHomo)
    yy$MAF[yy$Name == i] <- (numMinors*2 + numHet) / (numMinors*2 + numMajor*2 + numHet*2)
    yy$MCF[yy$Name == i] <- numMiss / (numHet + numHomo + numMiss)
  }
  write.csv(yy, paste0(myPrefix, "_02_myG_Genotypes.csv"), row.names = F)
  # Plot
  mp4 <- ggplot(yy, aes(x = Het * 100)) +
    geom_histogram(color = "black", fill = "darkgreen", alpha = 0.7) +
    theme_gwaspr() +
    labs(title = "Genotype Details (myY)",
         subtitle = "Heterozygosity", x = NULL, y = NULL)
  mp5 <- ggplot(yy, aes(x = MAF * 100)) +
    geom_histogram(color = "black", fill = "darkgreen", alpha = 0.7) +
    theme_gwaspr() +
    labs(subtitle = "Minor Allele Frequency", x = "Percent", y = NULL)
  mp6 <- ggplot(yy, aes(x = MCF * 100)) +
    geom_histogram(color = "black", fill = "darkgreen", alpha = 0.7) +
    theme_gwaspr() +
    labs(subtitle = "Missing Call Frequency", x = NULL, y = NULL)
  mpB <- ggarrange(mp4, mp5, mp6, ncol = 3, nrow = 1, align = "h")
  ggsave(paste0(myPrefix,"_05_myG_Genotypes.png"), mpB, width = 12, height = 4)
  #
  mp <- ggarrange(mpA, mpB, ncol = 1, nrow = 2, align = "v")
  ggsave(paste0(myPrefix,"_03_myG_Details.png"), mp, width = 12, height = 8)
}
