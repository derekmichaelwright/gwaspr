#' gg_LD_Decay
#'
#' Creates a number of LD decay plots.
#' @param xG GWAS genotype object. Note: needs to be in hapmap format.
#' @param outputFolder Folder to place RData files into.
#' @param markerNum Number of markers to randomly select for each chromosome.
#' @return LD decay plots.
#' @export

gg_LD_Decay <- function(xG = myG, outputFolder, markerNum = 1000) {
  #
  # Create function
  movingAverage <- function(x, n = 5) {
    stats::filter(x, rep(1 / n, n), sides = 2)
  }
  # Prep data
  xx <- NULL
  for(i in unique(xG$chrom)) {
    xc <- xG %>% filter(chrom == i)
    load(paste0(outputFolder,"LD_Chrom_",i,"_",markerNum,".Rdata"))
    xi <- myLD$`R^2` %>% as.data.frame() %>%
      rownames_to_column("SNP1") %>%
      gather(SNP2, LD, 2:ncol(.)) %>%
      filter(!is.na(LD)) %>%
      mutate(Chr = i,
             SNP1_d = plyr::mapvalues(SNP1, xc$rs, xc$pos, warn_missing = F),
             SNP2_d = plyr::mapvalues(SNP2, xc$rs, xc$pos, warn_missing = F),
             Distance = as.numeric(SNP2_d) - as.numeric(SNP1_d)) %>%
      arrange(Distance, rev(LD))
    #
    xii <- xi %>% filter(Distance < 1000000)
    myloess <- stats::loess(LD ~ Distance, data = xii, span=0.50)
    #
    xi <- xi %>%
      mutate(Moving_Avg = movingAverage(LD, n = 100),
             Loess_Chr = ifelse(Distance < 1000000, predict(myloess), NA))
    #
    xx <- bind_rows(xx, xi)
  }
  #
  x1 <- xx %>% group_by(Chr) %>%
    summarise(Mean_LD = mean(Moving_Avg, na.rm = T))
  x2 <- xx %>% filter(Loess_Chr < 0.2) %>% group_by(Chr) %>%
    summarise(Threshold_0.2 = min(Distance, na.rm = T))
  #x3 <- xx %>% left_join(x1, by = "Chr") %>%
  #  group_by(Chr) %>% filter(Moving_Avg < Mean_LD) %>%
  #  summarise(Threshold_0.1 = min(Distance, na.rm = T))
  myChr <-  x1 %>% left_join(x2, by = "Chr") #%>%
  #left_join(x3, by = "Chr")
  #
  xx <- left_join(xx, myChr, by = "Chr") %>%
    mutate(Chr = as.factor(Chr))
  yy <- xx %>% filter(Distance < 1000000) %>%
    arrange(Distance, rev(LD))
  myloess <- stats::loess(LD ~ Distance, data = yy, span=0.50)
  yy <- yy %>%
    mutate(Moving_Avg = movingAverage(LD, n = 100),
           Loess_Geno = ifelse(Distance < 1000000, predict(myloess), NA))
  #
  x1 <- yy %>% summarise(Mean_LD = mean(Moving_Avg, na.rm = T))
  x2 <- yy %>% filter(Loess_Chr < 0.2) %>%
    summarise(Threshold_0.2 = min(Distance, na.rm = T))
  myGeno <- cross_join(x1, x2)
  # Plot first 1 Mbp
  mp1 <- ggplot(yy, aes(x = Distance/1000)) +
    geom_point(aes(y = LD), size = 0.3, pch = 15, alpha = 0.3) +
    geom_line(aes(y = Moving_Avg), alpha = 0.7) +
    geom_line(aes(y = Loess_Geno)) +
    geom_hline(yintercept = 0.2, color = "blue", lty = 2) +
    geom_hline(data = myGeno, aes(yintercept = Mean_LD), color = "red", lty = 2) +
    geom_vline(data = myGeno, lty = 2, linewidth = 0.3, alpha = 0.8,
               aes(xintercept = Threshold_0.2/1000)) +
    scale_x_continuous(seq(0, 1000, by = 100), expand = c(0,10)) +
    facet_wrap(paste("Threshold = ", myGeno$Threshold_0.2)) +
    theme_gwaspr(legend.position = "none",
                 axis.title.y = ggtext::element_markdown()) +
    labs(title = "LD of all markers within 1 Mbp", y = "r^2", x = "Kbp")
  # Plot first 1 Mbp
  mp2 <- ggplot(yy, aes(x = Distance/1000, y = Loess_Chr)) +
    geom_line(aes(y = Moving_Avg), alpha = 0.5) +
    geom_line() +
    geom_hline(data = myChr, aes(yintercept = Mean_LD), color = "red", lty = 2) +
    geom_hline(yintercept = 0.2, color = "blue", lty = 2) +
    scale_y_continuous(breaks = seq(0, 0.5, by = 0.1), limits = c(0,0.5)) +
    facet_wrap(paste("Chr =", Chr) + paste("Threshold =", Threshold_0.2) ~ .,
               ncol = 7, scales = "free_x") +
    geom_vline(data = myChr, lty = 2, size = 0.3, alpha = 0.8,
               aes(xintercept = Threshold_0.2/1000)) +
    theme_gwaspr(legend.position = "none",
                 axis.title.y = ggtext::element_markdown()) +
    labs(title = "By Chromosome", y = "r^2", x = "Kbp")
  # Plot full chromsomes
  mp3 <- ggplot(xx, aes(x = Distance/1000000, y = Moving_Avg)) +
    geom_line(size = 0.5, alpha = 0.5) +
    geom_hline(data = myChr, aes(yintercept = Mean_LD), color = "red", lty = 2) +
    geom_hline(yintercept = 0.2, color = "blue", lty = 2) +
    scale_y_continuous(breaks = seq(0, 0.5, by = 0.1), limits = c(0,0.5)) +
    facet_wrap(paste("Chr =", Chr) ~ ., ncol = 7, scales = "free_x") +
    theme_gwaspr(legend.position = "none",
                 axis.title.y = ggtext::element_markdown()) +
    labs(title = "All Markers", y = "r^2", x = "Mbp")
  #yy <- yy %>% filter(Distance < 10000)
  # Append
  mp <- ggarrange(mp1, mp2, mp3, ncol = 1, nrow = 3)
  mp
}

#myG <- read.csv("vignettes/gwaspr_myG_hmp.csv", header = T)
#xG = myG; outputFolder = "vignettes/LD_Decay/"; markerNum = 1000
