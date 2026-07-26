#' gg_LDHeatmap
#'
#' Creates a manhattan plot.
#' @param xG GWAS genotype object. Note: needs to be in hapmap format.
#' @param chr Chromosome to plot.
#' @param pos1 Start position within the selected chromosome.
#' @param pos2 End position within the selected chromosome.
#' @param myMs Markers to highlight within the plot.
#' @param myTitle Title for the plot.
#' @param axisTextSize Text size for the axis labels (genotype names).
#' @param nameTrim String used to trim marker names.
#' @return A LD Heatmap plot.
#' @export

gg_LDheatmap <- function(xG = myG, chr = 6, pos1 = 0, pos2 = 6000000,
                         myMs = NULL,
                         myTitle = NULL,
                         axisTextSize = NULL,
                         nameTrim = NULL ) {
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
  xx <- myLD$`R^2` %>% as.data.frame()
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
  xl <- data.frame(SNP1 = xC$rs, SNP2 = xC$rs)
  #
  print("ploting results")
  mp2 <- ggplot(xx, aes(x = SNP1, y = SNP2)) +
    geom_tile(aes(fill = LD)) +
    scale_fill_gradient2(low = "grey90", mid = "goldenrod1", high = "darkred", midpoint = 0.5) + #"grey90"
    scale_x_discrete(position = "top") +
    #scale_y_discrete(limits = rev) +
    theme_gwaspr(legend.position = "bottom",
                 axis.text.x = element_text(angle = 90, vjust = 1),
                 axis.text = element_text(size = axisTextSize)) +
    theme(panel.grid = element_blank()) +
    labs(x = NULL, y = NULL, caption = myCaption)
  if(!is.null(myMs)) {
    mp2 <- mp2 + geom_point(data = xm, pch = 8)
    xm <- xC %>% filter(rs %in% myMs)
  }
  mySubtitle <-paste("Chr", chr, "| pos", format(pos1, scientific = F), "-", format(pos2, scientific = F), "| Length", format(myLength,scientific = F))
  mp1 <- ggplot(xC) +
    geom_hline(yintercept = 0) +
    #geom_hline(yintercept = 1) +
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
      geom_point(data = xm, aes(x = scaled_rank, y = 0), pch = 8, color = "darkred") +
      geom_point(data = xm, aes(x = rank, y = 1), pch = 8, color = "darkred")
  }
  #
  mp <- ggarrange(mp1, mp2, nrow = 2, ncol = 1, align = "v", heights = c(0.2,1))
  mp
}

#
#xG = myG_LDP; chr = 1; pos1 = 365986872-500000; pos2 = 365986872+500000
#myMs = "Lcu.1GRN.Chr1p365986872"; axistextsize = 4; nameTrim = "Lcu.1GRN.Chr1"

