#' gg_LDHeatmap
#'
#' Creates a manhattan plot.
#' @param xG GWAS genotype object. Note: needs to be in hapmap format.
#' @param chr Chromosome to plot.
#' @param pos1 Start position within the selected chromosome.
#' @param pos2 End position within the selected chromosome.
#' @param metric Which LD calculation to use. Default is "D'".
#' @param threshold Value for selecting linked markers. default is "0.9".
#' @param myMs Markers to highlight within the plot.
#' @param myTitle Title for the plot.
#' @param axisTextSize Text size for the axis labels (genotype names).
#' @param nameTrim String used to trim marker names.
#' @param color.low Color for gradient low.
#' @param color.mid Color for gradient mid.
#' @param color.high Color for gradient high.
#' @return A LD Heatmap plot.
#' @export

gg_LDheatmap <- function(xG = myG,
                         chr, pos1, pos2,
                         metric = "R^2",
                         threshold = 0.9,
                         myMs = NULL,
                         myTitle = NULL,
                         axisTextSize = NULL,
                         nameTrim = NULL,
                         color.low = "white",
                         color.mid = "goldenrod1",
                         color.high = "darkred") {
  #
  print("starting")
  myCaption <- paste(myMs, collapse = "|")
  xG <- xG[,c(-2,-5,-6,-7,-8,-9,-10,-11)]
  xG <- xG %>% filter(chrom == chr, pos >= pos1, pos <= pos2)
  #
  if(!is.null(nameTrim)) { xG$rs <- gsub(nameTrim, "", xG$rs) }
  if(!is.null(nameTrim) & !is.null(myMs)) { myMs <- gsub(nameTrim, "", myMs) }
  #
  xC <- xG %>% select(rs, chrom, pos)
  xC <- xC %>% mutate(scaled_pos = scales::rescale(1:nrow(xC), to = c(min(xC$pos),max(xC$pos))),
                      rank = 1:nrow(xC))
  myScale <- scales::rescale(c(pos1, xC$pos, pos2), to = c(1,nrow(xC)))
  xC$scaled_rank <- myScale[c(-1,-length(myScale))]
  #
  dna <- data.frame(stringsAsFactors = F,
                    Symbol = c("A", "C", "G", "T", "U", "R", "Y", "S", "W", "K", "M"),
                    Value  = c("A/A","C/C","G/G","T/T","U/U", "A/G","C/T","G/C","A/T","G/T","A/C") )
  #
  for(i in 4:ncol(xG)) {
    xG[xG[,i]=="N", i] <- NA
    xG[,i] <-  plyr::mapvalues(xG[,i], dna$Symbol, dna$Value, warn_missing = FALSE)
  }
  #
  xG <- xG %>% column_to_rownames("rs")
  xG <- xG[,c(-1,-2)]
  #
  xG <- xG %>% t() %>% as.data.frame()
  print("prepping genotypes for LD")
  xG <- genetics::makeGenotypes(xG)
  print("calcualting LD")
  myLD <- genetics::LD(xG)
  print("LD finsished")
  #
  xx <- myLD[metric][[1]] %>% as.data.frame()
  xx[1,1] <- 0
  xx[nrow(xx),ncol(xx)] <- 0
  xx <- xx %>%
    rownames_to_column("SNP1") %>%
    gather(SNP2, LD, 2:ncol(.)) %>%
    filter(!is.na(LD), SNP1 %in% xC$rs, SNP2 %in% xC$rs) %>%
    mutate(Chr = chr,
           SNP1_d = plyr::mapvalues(SNP1, xC$rs, xC$pos, warn_missing = F),
           SNP2_d = plyr::mapvalues(SNP2, xC$rs, xC$pos, warn_missing = F) ) %>%
    arrange(SNP1_d)
  #
  myLength <- pos2 - pos1
  xm <- data.frame(SNP1 = myMs, SNP2 = myMs)
  xl <- xx %>% filter(SNP1 %in% myMs | SNP2 %in% myMs) %>%
    filter(LD > threshold) %>%
    mutate(SNP1 = ifelse(SNP1 %in% myMs, SNP2, SNP1),
           SNP2 = ifelse(SNP2 %in% myMs, SNP1, SNP2))

  #xl <- data.frame(SNP1 = xC$rs, SNP2 = xC$rs)
  #
  myPmin <- xx$SNP1[1]
  myPmax <- xx$SNP1[nrow(xx)]
  print("ploting results")
  mp2 <- ggplot(xx, aes(x = SNP1, y = SNP2)) +
    geom_tile(aes(fill = LD)) +
    geom_rect(data = xm, xmin = myPmin, aes(xmax = SNP1, ymin = SNP2, ymax = SNP2),  color = alpha("black",0.1), fill = NULL) +
    geom_rect(data = xm, aes(xmin = SNP1, xmax = SNP2, ymin = SNP2), ymax = myPmax, color = alpha("black",0.1), fill = NULL) +
    scale_fill_gradient2(name = metric, low = color.low, mid = color.mid, high = color.high, midpoint = 0.5) +
    scale_x_discrete(position = "top") +
    theme_gwaspr(legend.position = "bottom",
                 axis.text.x = element_text(angle = 90, vjust = 1),
                 axis.text = element_text(size = axisTextSize)) +
    theme(panel.grid = element_blank()) +
    labs(x = NULL, y = NULL, caption = myCaption)
  if(!is.null(myMs)) {
    mp2 <- mp2 +
      geom_point(data = xm, pch = 8, size = 2) +
      geom_point(data = xl, pch = 16, size = 1)
    xm <- xC %>% filter(rs %in% myMs)
    xl <- xC %>% filter(rs %in% xl$SNP1)
  }
  mySubtitle <-paste("Chr", chr, "| pos", format(pos1, scientific = F), "-", format(pos2, scientific = F), "| Length", format(myLength,scientific = F))
  mp1 <- ggplot(xC) +
    geom_hline(yintercept = 0) +
    geom_segment(aes(x = rank, xend = scaled_rank, y = 1, yend = 0), linewidth = 0.1) +
    scale_x_continuous(breaks = xC$scaled_pos, labels = xC$rs, expand = c(0,0.5),
                       position = "top") +
    theme_void() +
    theme(panel.grid = element_blank(),
          plot.background = element_rect(fill = "white")) +
    scale_y_reverse() +
    labs(y = NULL, x = NULL, title = myTitle,
         subtitle = mySubtitle)
  if(!is.null(myMs)) {
    mp1 <- mp1 +
      geom_point(data = xl, aes(x = scaled_rank, y = 0), pch = 16, alpha = 0.7, color = "darkred") +
      geom_point(data = xl, aes(x = rank, y = 1), pch = 16, alpha = 0.7, color = "darkred") +
      geom_point(data = xm, aes(x = scaled_rank, y = 0), pch = 8, color = "darkred") +
      geom_point(data = xm, aes(x = rank, y = 1), pch = 8, color = "darkred")
  }
  #
  mp <- ggarrange(mp1, mp2, nrow = 2, ncol = 1, align = "v", heights = c(0.2,1))
  colnames(xx)[3] <- metric
  list(mp, xx)
  #mp
}

#
#xG = myG; chr = 6; pos1 = 2500000; pos2 = 4000000; myMs = "Lcu.1GRN.Chr6p3269280";
#xG = myG; chr = 2; pos1 = 7500000; pos2 = 8000000; myMs = "Lcu.1GRN.Chr2p7836680";
#metric = "D'";  axisTextSize = 3; threshold = 0.9;
#nameTrim = "Lcu.1GRN.Chr6"; myTitle = "Lentil Diversity Panel"
#color.low = "white"; color.mid = "goldenrod1"; color.high = "darkred"

#metric = "R^2"
