#' gg_Marker_Bubble
#'
#' [Creates a marker scatter plot with myG and myY objects.](https://derekmichaelwright.github.io/gwaspr/articles/gg_Marker_Bubble.html)
#' @param xG GWAS genotype object. Note: needs to be in hapmap format.
#' @param xY GWAS phenotype object.
#' @param trait.x Trait for the x-axis.
#' @param trait.y Trait for the y-axis.
#' @param facet.rows Trait for the row facetting.
#' @param facet.cols Trait for the col facetting.
#' @param markers Markers to plot.
#' @param marker.colors Colors to fill the marker bubbles.
#' @param remove.hets Logical, Whether to remove hets or not. advised if plotting multiple markers.
#' @param remove.N Logical, Whether to remove N or NN.
#' @param legend.pos Position for the legend.
#' @param legend.rows Number of rows for the legend.
#' @param title Title for the plot.
#' @param subtitle Subtitle for the plot. Defaults to the list of markers.
#' @return Marker plot.
#' @export

gg_Marker_Bubble <- function (
    xG,
    xY,
    trait.x = NULL,
    trait.y = "Alleles",
    facet.rows = NULL,
    facet.cols = NULL,
    markers,
    marker.colors = gwaspr_Colors,
    remove.hets = F,
    remove.N = T,
    legend.pos = "bottom",
    legend.rows = 1,
    title = NULL,
    subtitle = paste(markers, collapse = "\n")
    ) {
  # Prep data
  xx <- xG %>% rename(SNP=1) %>%
    filter(SNP %in% markers) %>%
    dplyr::select(-2,-3,-4,-5,-6,-7,-8,-9,-10,-11) %>%
    column_to_rownames("SNP") %>%
    t() %>% as.data.frame() %>%
    #select(markers) %>%
    mutate(Alleles = NA)
  #
  if(remove.hets == T) { for(i in 1:length(markers)) { xx <- xx[xx[,i] %in% c("A","T","G","C","AA","TT","GG","CC"),] } }
  if(remove.hets == T) { for(i in 1:length(markers)) { xx <- xx[!xx[,i] %in% c("N","NN"),] } }
  #
  for(k in 1:nrow(xx)) { xx$Alleles[k] <- paste(xx[k,1:length(markers)], collapse = "-") }
  #
  xx <- xx %>%
    select(Alleles) %>%
    mutate(Alleles2=Alleles) %>%
    rownames_to_column("Name") %>%
    left_join(xY, by = "Name")
  #
  ifelse(is.null(trait.x), xx <- xx %>% mutate(Tx = ""), xx <- xx %>% rename(Tx=trait.x))
  ifelse(is.null(trait.y), xx <- xx %>% mutate(Ty = ""), xx <- xx %>% rename(Ty=trait.y))
  ifelse(is.null(facet.rows), xx <- xx %>% mutate(Fr = ""), xx <- xx %>% rename(Fr=facet.rows))
  ifelse(is.null(facet.cols), xx <- xx %>% mutate(Fc = ""), xx <- xx %>% rename(Fc=facet.cols))
  #
  xx <- xx %>%
    rename(Alleles=Alleles2) %>%
    group_by(Alleles, Tx, Ty, Fr, Fc) %>%
    summarise(Count = n()) %>%
    ungroup()
  #
  ggplot(xx, aes(x = Tx, y = Ty)) +
    geom_point(aes(fill = Alleles, size = Count), pch = 21, alpha = 0.7) +
    geom_text(aes(label = Count)) +
    facet_grid(Fr ~ Fc, scales= "free_x", space = "free_x", drop = T) +
    scale_y_discrete(sec.axis = dup_axis(name = facet.rows)) +
    scale_x_discrete(sec.axis = dup_axis(name = facet.cols)) +
    scale_fill_manual(values = gwaspr_Colors) +
    scale_size_continuous(range = c(1,15), guide = "none") +
    theme_gwaspr(legend.position = legend.pos,
                 axis.text.x.bottom = element_text(angle = 45, hjust = 1),
                 axis.text.x.top = element_blank(),
                 axis.ticks.x.top = element_blank(),
                 axis.text.y.right = element_blank(),
                 axis.ticks.y.right = element_blank()) +
    labs(x = trait.x, y = trait.y,
         title = title, subtitle = subtitle)
}

#xG = myG; xY = myY
#trait.x = "SubRegion"
#trait.y = "Alleles"
#facet.rows = "Region"
#facet.cols = "Cotyledon_Color"
#markers = "Lcu.1GRN.Chr1p365986872"
#marker.axis = ""
#marker.colors = gwaspr_Colors; remove.hets = T; remove.N = T; legend.rows = 1;
#title = NULL; subtitle = paste(markers, collapse = "\n")

#trait.x = "Cotyledon_Color"
#trait.y = "Alleles"
#facet.rows = NULL
#facet.cols = NULL
