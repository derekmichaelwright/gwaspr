#' gg_Marker_Box
#'
#' [Creates a marker plot with myG and myY objects.](https://derekmichaelwright.github.io/gwaspr/articles/gg_Marker.html)
#' @param xG GWAS genotype object. Note: needs to be in hapmap format.
#' @param xY GWAS phenotype object.
#' @param traits Traits to plot.
#' @param markers Markers to plot.
#' @param marker.colors Colors to fill in the violin and boxplots.
#' @param remove.hets Logical, Whether to remove hets or not. advisded if plotting multiple markers.
#' @param plot.violin Logical, whether or not to plot violins.
#' @param plot.box Logical, whether or not to plot the boxplots.
#' @param plot.points Logical, whether or not to plot points.
#' @param box.width Width for the boxplot.
#' @param point.size Size for the points.
#' @param point.beeswarm Logical. If False (the default), will plot points with `geom_quasirandom`. If TRUE, will plot points with `geom_beeswarm`.
#' @param myncol Number of columns for facetting when plotting multiple traits.
#' @param title Title for the plot.
#' @param legend.rows Number of rows for the legend.
#' @param subtitle Subtitle for the plot. Defaults to the list of markers.
#' @param yLab Label for the y-axis.
#' @param cv.source Where to get your `cv.name` from. Default is "xG", while "xY" is the other option.
#' @param cv.name Covariable data for points.
#' @param cv.colors Covariable colors for filling points.
#' @param cv.label Label for the covariate.
#' @return Marker plot.
#' @export

gg_Marker_Box <- function (
    xG,
    xY,
    traits,
    markers,
    marker.colors = gwaspr_Colors,
    remove.hets = T,
    plot.violin = T,
    plot.box = T,
    plot.points = T,
    box.width = 0.1,
    point.size = 1,
    point.beeswarm = F,
    myncol = NULL,
    title = NULL,
    legend.rows = 1,
    subtitle = paste(markers, collapse = "\n"),
    yLab = traits,
    cv.source = "xG",
    cv.name = NULL,
    cv.colors = NULL,
    cv.label = NULL
    ) {
  #
  myLab <- paste(markers, collapse = "\n")
  #
  xT <- xY %>% select(Name=1, traits) %>%
    gather(Trait, Value, traits)
  #
  xx <- xG %>% rename(SNP=1) %>%
    filter(SNP %in% markers) %>%
    dplyr::select(-2,-3,-4,-5,-6,-7,-8,-9,-10,-11) %>%
    column_to_rownames("SNP") %>%
    t() %>% as.data.frame() %>%
    select(markers) %>%
    mutate(Alleles = NA)
  #
  if(remove.hets == T) {
    for(i in 1:length(markers)) { xx <- xx[xx[,i] %in% c("A","T","G","C","AA","TT","GG","CC"),] }
  }
  #
  for(i in 1:nrow(xx)) { xx$Alleles[i] <- paste(xx[i,1:length(markers)], collapse = "-") }
  #
  xx <- xx %>% rownames_to_column("Name") %>%
    left_join(xT, by = "Name") %>%
    filter(!is.na(Value))
  #
  yy <- xx %>% filter(Trait == traits[1]) %>%
    group_by(Alleles) %>%
    summarise(Value = mean(Value, na.rm = T)) %>%
    arrange(Value)
  xx <- xx %>%
    mutate(Alleles = factor(Alleles, levels = rev(yy$Alleles)),
           Trait = factor(Trait, levels = traits))
  # Prep covariable data
  if(!is.null(cv.name)) {
    if(cv.source == "xY") {
      xx <- xx %>%
        left_join(xY %>% select(Name=1, cv.name), by = "Name") %>%
        filter(!is.na(get(cv.name)))
    }
    if(cv.source == "xG") {
      xCV <- xG %>% rename(SNP=1) %>%
        filter(SNP == cv.name) %>%
        dplyr::select(-2,-3,-4,-5,-6,-7,-8,-9,-10,-11) %>%
        column_to_rownames("SNP") %>%
        t() %>% as.data.frame() %>%
        rownames_to_column("Name")
      xx <- xx %>%
        left_join(xCV, by = "Name") %>%
        filter(!is.na(get(cv.name)))
      if(remove.hets == T) { xx <- xx %>% filter(get(cv.name) %in% c("A","T","G","C","AA","TT","GG","CC")) }
    }
  }
  #
  # Plot
  mp <- ggplot(xx, aes(x = Alleles, y = Value))
  if(plot.violin == T & plot.box == T) {
    mp <- mp + geom_violin(aes(fill = Alleles), alpha = 0.3) +
      geom_boxplot(fill = "white", width = box.width, outlier.shape = NA)
  }
  if(plot.violin == T & plot.box == F) {
    mp <- mp + geom_violin(aes(fill = Alleles), alpha = 0.3)
  }
  if(plot.violin == F & plot.box == T) {
    mp <- mp + geom_boxplot(aes(fill = Alleles), alpha = 0.5, width = box.width, outlier.shape = NA)
  }
  mp <- mp +
    facet_wrap(Trait ~ ., scales = "free_y", ncol = myncol) +
    scale_fill_manual(name = NULL, values = marker.colors, guide = F) +
    theme_gwaspr(legend.position = "none",
                 axis.text.x = element_text(angle = 45, hjust = 1) ) +
    labs(title = title, subtitle = subtitle, x = NULL, y = yLab)
  if (plot.points == T) {
    if(is.null(cv.name)) {
      if(point.beeswarm == T) {
        mp <- mp + geom_beeswarm(size = point.size, alpha = 0.8, pch = 16, method = "center")
      } else {
        mp <- mp + geom_quasirandom(size = point.size, alpha = 0.8, pch = 16)
      }
    } else {
      if(point.beeswarm == T) {
        mp <- mp + geom_beeswarm(aes(color = get(cv.name)), size = point.size, alpha = 0.8, pch = 16, method = "center")
      } else {
        mp <- mp + geom_quasirandom(aes(color = get(cv.name)), size = point.size, alpha = 0.8, pch = 16)
      }
      mp <- mp +
        scale_color_manual(name = cv.label, values = cv.colors) +
        theme(legend.position = "bottom") +
        guides(color = guide_legend(nrow = legend.rows)) #, override.aes = list(size = 2)
    }
  }
  mp
}

#xG = myG; xY = myY;
#traits = c("Canopy.Height_Ba16", "Canopy.Width_Ba16", "Canopy.Height_Ba17", "Canopy.Width_Ba17")
#markers = myMarkers
#marker.colors = c("darkgreen", "darkgoldenrod3", "darkred", "steelblue4", "darkslategray", "maroon4", "purple4", "darkblue")
#box.width = 0.1; plot.points = T; plot.violin = F; plot.box = T; myncol = 4
#yLab = ""; subtitle = ""; point.size = 2; title = ""; cv.colors = gwaspr_Colors; cv.label = "CV"
#legend.rows = 1; point.beeswarm = F; remove.hets = T;
