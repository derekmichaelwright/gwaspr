#' gg_GWAS_Hits
#'
#' Creates a summary GWAS plot of significant associations.
#' @param xx Table of significant GWAS results. See ?table_GWAS_Results().
#' @param xG Genotype data in hapmap format.
#' @param xCV Optional, for filtering if you have a "CV" column.
#' @param traits List of traits to use.
#' @param range Range for binning GWAS hits.
#' @param title Title for horizontal facet.
#' @param sigMin Minimum number of hits to plot.
#' @param models Models to read.
#' @param model.colors Colors for each model.
#' @param vlines Markers to be labelled with a vertical red line.
#' @param vline.colors colors for each vertical line.
#' @param vline.types lty for each vertical line.
#' @param legend.rows Number of rows for the legend.
#' @param facet Logical. facet by GWAS model.
#' @param caption Logical. should a caption be added.
#' @return A GWAS Hits plot.
#' @export

gg_GWAS_Hits <- function(
    folder = "GWAS_Results/",
    xx, xG, xCV = NULL,
    traits = list_Traits(folder),
    groups = NULL,
    chr = NULL, pos1 = NULL, pos2 = NULL,
    range = 2000000,
    title = "",
    sigMin = 0,
    models =  c("BLINK", "FarmCPU", "MLMM", "MLM", "GLM", "CMLM", "SUPER"),
    model.colors = gwaspr_ModelColors,
    model.shapes = c(21:25,21:22),
    vlines = NULL,
    vline.colors = rep("red", length(vlines)),
    vline.types = rep(1, length(vlines)),
    vline.legend = T,
    legend.rows = 1,
    chr.unit = "100 Mbp",
    facet = F,
    caption = T,
    skyline = "NYC"
    ) {
  #
  check1 <- run_Summary(folder = folder) %>% dropNAcol()
  check2 <- is_Ordered(folder = folder) %>%
    select(colnames(check1))
  if(sum(is.na(check2)) > 0) {
    warning("Some of your GWAS results files might not be ordered by pvalue")
    warning("use order_GWAS_Results() before making these summary plots")
    warning("use is_ordered() to check if GWAS results are properly ordered for use with this function")
  }
  #
  fnames <- list_Result_Files(folder)
  fnames <- fnames[grepl(paste(models,collapse="|"),fnames)]
  fnames <- fnames[grepl(paste(c(paste0(traits, ".csv"), paste0(traits, "\\(")), collapse="|"), fnames)]
  #
  if(!is.null(skyline)) {
    if(skyline == "NYC") { fnames <- fnames[!grepl("\\(Kansas\\)", fnames)] }
    if(skyline == "Kansas") {
      fnames <- fnames[!(grepl("\\(NYC\\)", fnames) & grepl("FarmCPU", fnames))]
      fnames <- fnames[!(grepl("\\(NYC\\)", fnames) & grepl("BLINK",   fnames))]
    }
  }
  #
  #
  xP <- NULL
  # i<- fnames[1]
  for(i in fnames) {
    xPi <- table_GWAS_Results(folder = folder, fnames = i, skyline = skyline)
    #
    if(nrow(xPi)>0) { xP <- bind_rows(xP, xPi) }
  }
  #
  xP <- xP %>% mutate(Pvalue = negLog10_P)
  #
  if(chr.unit == "1 kbp")   { x.unit = 1000 }
  if(chr.unit == "10 kbp")  { x.unit = 10000 }
  if(chr.unit == "50 kbp")  { x.unit = 50000 }
  if(chr.unit == "100 kbp") { x.unit = 100000 }
  if(chr.unit == "1 Mbp")   { x.unit = 1000000 }
  if(chr.unit == "10 Mbp")  { x.unit = 10000000 }
  if(chr.unit == "50 Mbp")  { x.unit = 50000000 }
  if(chr.unit == "100 Mbp") { x.unit = 100000000 }
  if(chr.unit == "1 Gbp")   { x.unit = 1000000000 }
  if(!chr.unit %in% c("1 kbp", "10 kbp", "50 kbp", "100 kbp",
                      "1 Mbp", "10 Mbp", "50 Mbp", "100 Mbp", "1 Gbp")) { print("error in chr.unit") }
  #
  xP <- xP %>%
    filter(!is.na(SNP)) %>%
    arrange(Chr, Pos, P.value, Trait) %>%
    mutate(Model = factor(Model, levels = models),
           Trait = factor(Trait, levels = traits),
           Group = "") %>%
    filter(!is.na(Trait)) %>%
    arrange(desc(Model))
  #
  xG <- read.csv(paste0(folder, fnames[grepl(traits[1],fnames)][1])) %>%
    mutate(Trait = xP$Trait[1],
           Trait = factor(Trait, levels = traits),
           Group = "", Model = NULL)
  #
  if(!is.null(chr)) {
    xG <- xG %>% filter(Chr %in% chr)
    xP <- xP %>% filter(Chr %in% chr)
    if(!is.null(genes)) { genes <- genes %>% filter(Chr %in% chr) }

    if(!is.null(pos1)) {
      xG <- xG %>% filter(Pos > pos1, Pos < pos2)
      xP <- xP %>% filter(Pos > pos1, Pos < pos2)
      if(!is.null(genes)) { genes <- genes %>% filter(Start > pos1, End < pos2) }
    }
  }
  #
  if(!is.null(groups)) {
    xP <- xP %>%
      mutate(Group = plyr::mapvalues(Trait, traits, groups),
             Group = factor(Group, levels = unique(groups)))
    xG <- xG %>%
      mutate(Group = plyr::mapvalues(Trait, traits, groups),
             Group = factor(Group, levels = unique(groups)))
  }
  #
  #
  #
  xx <- xP %>% filter(Trait %in% traits, Model %in% models) %>% mutate(Hits = NA)
  xG <- xG %>% dplyr::select(SNP=1, Chr=3, Pos=4)
  # i<-1
  for(i in 1:nrow(xx)) {
    myChr <- xx$Chr[i]
    myPos <- xx$Pos[i]
    myMod <- xx$Model[i]
    myGrp <- xx$Group[i]
    xi <- xx %>%
      filter(Chr == myChr, Model == myMod, Group == myGrp,
             Pos > myPos-range, Pos < myPos+range)
    xx$Hits[i] <- nrow(xi)
    xi <- xi %>% filter(!(SNP == xx$SNP[i] & Trait == xx$Trait[i]))
    xx$Pos[xx$SNP %in% xi$SNP & xx$Trait %in% xi$Trait] <- NA
  }
  xx <- xx %>% filter(!is.na(Pos), Hits > sigMin) %>%
    mutate(Model = factor(Model, levels = models))
  #
  if(!is.null(xCV)) { xx <- xx %>% filter(CV %in% xCV) }
  ymax <- max(xx$Hits)
  #
  myCaption <- paste("Bin size =", range, "\n", "Traits =",
                     paste(unique(xx$Trait), collapse = ", "))
  #
  mp <- ggplot(xx, aes(x = Pos / 100000000) ) + geom_blank(data = xG)
  #
  if(!is.null(vlines)) {
    xGM <- xG %>%
      filter(SNP %in% vlines) %>%
      mutate(SNP = factor(SNP, levels = vlines)) %>%
      arrange(SNP)
    mp <- mp +
      geom_vline(data = xGM, alpha = 0.7,
                 aes(xintercept = Pos/1e+08, color = SNP, lty = SNP)) +
      scale_color_manual(name = NULL, values = vline.colors) +
      scale_linetype_manual(name = NULL, values = vline.types)
  }
  mp <- mp +
    geom_point(aes(y = Hits, size = Hits, key1 = SNP,
                   fill = Model, shape = Model), alpha = 0.7) +
    facet_grid(. ~ paste("Chr", Chr), space = "free_x", scales = "free") +
    scale_shape_manual(values = model.shapes)+
    scale_fill_manual(values = model.colors)+
    scale_size_continuous(range = c(0.5,3)) +
    theme_gwaspr() +
    theme(legend.position = "bottom") +
    guides(fill = guide_legend(nrow = legend.rows, override.aes = list(size = 1.5)),
           color = guide_legend(nrow = legend.rows, byrow = T),
           size = "none",
           shape = guide_legend(nrow = legend.rows, byrow = T) ) +
    labs(title = title, x = "100 Mbp", y = "Significant Associations", caption = myCaption)
  #
  if(facet == T) {
    mp <- mp + facet_grid(Model ~ paste("Chr", Chr), space = "free_x", scales = "free")
  }
  mp
}

