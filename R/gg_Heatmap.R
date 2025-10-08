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

gg_LDheatmap <- function(xg = myG, chr = 6, pos1 = 0, pos2 = 6000000,
                         myMs = NULL,
                         myTitle = NULL,
                         axisTextSize = 0.5,
                         nameTrim = NULL) {
  #
  print("starting")
  myCaption <- paste(myMs, collapse = "|")
  xg <- xg[,c(-2,-5,-6,-7,-8,-9,-10,-11)]
  xg <- xg %>% filter(chrom == chr, pos >= pos1, pos <= pos2)
  #
  if(!is.null(nameTrim)) { xg$rs <- gsub(nameTrim, "", xg$rs) }
  if(!is.null(nameTrim) & !is.null(myMs)) { myMs <- gsub(nameTrim, "", myMs) }
  #
  xc <- xg %>% select(rs, chrom, pos)
  xc <- xc %>% mutate(scaled_pos = scales::rescale(1:nrow(xc), to = c(min(xc$pos),max(xc$pos))) )
  #
  dna <- data.frame(stringsAsFactors = F,
                    Symbol = c("A", "C", "G", "T", "U",
                               "R", "Y", "S", "W", "K", "M"),
                    Value  = c("A/A","C/C","G/G","T/T","U/U",
                               "A/G","C/T","G/C","A/T","G/T","A/C") )
  for(i in 4:ncol(xg)) {
    xg[xg[,i]=="N", i] <- NA
    xg[,i] <-  plyr::mapvalues(xg[,i], dna$Symbol, dna$Value, warn_missing = FALSE)
  }
  #
  xg <- xg %>% column_to_rownames("rs")
  xg <- xg[,c(-1,-2)]
  #
  xg <- xg %>% t() %>% as.data.frame()
  print("prepping genotypes for LD")
  xg <- genetics::makeGenotypes(xg)
  print("calcualting LD")
  myLD <- genetics::LD(xg)
  print("LD finsished")
  #
  xx <- myLD$`R^2` %>% as.data.frame()
  xx[1,1] <- 0
  xx[nrow(xx),ncol(xx)] <- 0
  xx <- xx %>%
    rownames_to_column("SNP1") %>%
    gather(SNP2, LD, 2:ncol(.)) %>%
    filter(!is.na(LD), SNP1 %in% xc$rs, SNP2 %in% xc$rs) %>%
    mutate(Chr = chr,
           SNP1_d = plyr::mapvalues(SNP1, xc$rs, xc$pos, warn_missing = F),
           SNP2_d = plyr::mapvalues(SNP2, xc$rs, xc$pos, warn_missing = F) ) %>%
    arrange(SNP1_d)
  #
  myExpand <- nrow(xx) * 0.0000035
  xm <- data.frame(SNP1 = myMs, SNP2 = myMs)
  xl <- data.frame(SNP1 = xc$rs, SNP2 = xc$rs)
  #
  print("ploting results")
  mp2 <- ggplot(xx, aes(x = SNP1, y = SNP2)) +
    geom_tile(aes(fill = LD)) +
    geom_point(data = xm, pch = 8) +
    scale_fill_gradient(low = "grey90", high = "darkred") +
    scale_x_discrete(position = "top") +
    theme_gwaspr(legend.position = "none",
                 axis.text.x = element_text(angle = 90, vjust = 1),
                 axis.text = element_text(size = axisTextSize)) +
    theme(panel.grid = element_blank()) +
    labs(x = NULL, y = NULL, caption = myCaption)
  #mp2
  xm <- xc %>% filter(rs %in% myMs)
  mp1 <- ggplot(xc) +
    geom_hline(yintercept = 0) +
    #geom_hline(yintercept = 1) +
    geom_segment(aes(x = scaled_pos, xend = pos, y = 1, yend = 0), linewidth = 0.1) +
    geom_point(data = xm, aes(x = pos, y = 0), pch = 8, color = "darkred") +
    geom_point(data = xm, aes(x = scaled_pos, y = 1), pch = 8, color = "darkred") +
    scale_x_continuous(breaks = xc$scaled_pos, labels = xc$rs, expand = c(myExpand,0),
                       position = "top") +
    theme_void() +
    theme(panel.grid = element_blank(),
          plot.background = element_rect(fill = "white")) +
    scale_y_reverse() +
    labs(y = NULL, x = NULL, title = myTitle,
         subtitle = paste("Chr", chr, "| pos", format(pos1, scientific = F), "-", format(pos2, scientific = F)))
  #mp1
  ggarrange(mp1, mp2, nrow = 2, ncol = 1, align = "v", heights = c(0.4,1))
}

#
#xg = myG_LDP; chr = 1; pos1 = 365986872-500000; pos2 = 365986872+500000
#myMs = "Lcu.1GRN.Chr1p365986872"; axistextsize = 4; nameTrim = "Lcu.1GRN.Chr1"

