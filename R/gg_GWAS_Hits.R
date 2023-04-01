#' gg_GWAS_Hits
#'
#' Creates a summary GWAS plot of significant associations.
#' @param xx Table of significant GWAS results. See ?table_GWAS_Results().
#' @param myG Genotype data.
#' @param myTs List of traits to use.
#' @param myR Range for binning GWAS hits.
#' @param myTitle Title for horizontal facet.
#' @param sigMin Minimum number of hits to plot.
#' @param myCV (optional) for filtering if you have a "CV" column.
#' @param vlines Markers to be labelled with a vertical red line.
#' @param vline.colors colors for each vertical line.
#' @return A GWAS Hits plot.
#' @export

gg_GWAS_Hits <- function(xx, myG, myTs, myR = 2000000, myTitle = "",
                         sigMin = 0, myCV = NULL,
                         vlines = NULL, vline.colors = rep("red", length(vlines)) ) {
  #
  myCols <- c("darkgreen", "darkorange3", "darkslategray4","darkblue")
  xx <- xx %>% filter(Trait %in% myTs) %>% mutate(Hits = NA)
  myG <- myG %>% select(SNP=1, Chr=3, Pos=4)
  #
  i<-1
  for(i in 1:nrow(xx)) {
    myChr <- xx$Chr[i]
    myPos <- xx$Pos[i]
    myMod <- xx$Model[i]
    xi <- xx %>%
      filter(Chr == myChr, Model == myMod,
             Pos > myPos-myR, Pos < myPos+myR)
    xx$Hits[i] <- nrow(xi)
    xi <- xi %>% filter(!(SNP == xx$SNP[i] & Trait == xx$Trait[i]))
    xx$Pos[xx$SNP %in% xi$SNP & xx$Trait %in% xi$Trait] <- NA
  }
  xx <- xx %>% filter(!is.na(Pos), Hits > sigMin) %>%
    mutate(Model = factor(Model, levels = c("MLM", "FarmCPU", "BLINK")))
  #
  # myTitle <- paste0(myTitle, " (", length(myTs), ")")
  if(!is.null(myCV)) { xx <- xx %>% filter(CV %in% myCV) }
  ymax <- max(xx$Hits)
  #
  # group_by(Model) %>% summarise(Count = n())
  mp <- ggplot(xx, aes(x = Pos / 100000000) ) +
    geom_blank(data = myG) #+
  #geom_hline(yintercept = sigThresh, alpha = 0.5, color = "darkred")
  if(!is.null(vlines)) {
    myGM <- myG %>%
      filter(SNP %in% vlines) %>%
      mutate(SNP = factor(SNP, levels = vlines)) %>%
      arrange(SNP)
    mp <- mp +
      geom_vline(data = myGM, alpha = 0.7,
                 aes(xintercept = Pos/1e+08, color = SNP)) +
      scale_color_manual(name = NULL, values = vline.colors)
  }
  mp <- mp +
    geom_point(aes(y = Hits, size = Hits, key1 = SNP,
                   fill = Model, shape = Model), alpha = 0.7) +
    facet_grid(paste(myTitle) ~ paste("Chr", Chr),
               space = "free_x", scales = "free_x") +
    scale_shape_manual(values = c(21,24:25), breaks = c("MLM","FarmCPU","BLINK")) +
    scale_fill_manual(values = myCols, breaks = c("MLM","FarmCPU","BLINK")) +
    scale_size_continuous(range = c(0.5,3)) +
    scale_y_continuous(labels = scales::label_number(accuracy = 1),
                       limits = c(0,ymax)) +
    theme_AGL +
    theme(legend.position = "bottom") +
    guides(color = F, size = F) +
    labs(x = "100 Mbp", y = "Significant Associations")
  mp
}
