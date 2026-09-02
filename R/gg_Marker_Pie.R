#' gg_Marker_Pie
#'
#' [Creates a marker plot with myG and myY objects.](https://derekmichaelwright.github.io/gwaspr/articles/gg_Marker.html)
#' @param xG GWAS genotype object. Note:  needs to be in hapmap format.
#' @param xY GWAS phenotype object.
#' @param trait Trait to plot.
#' @param trait.label Label for the Trait.
#' @param markers Markers to plot.
#' @param marker.colors Color palette.
#' @param title Title for the plot.
#' @param subtitle Subtitle for the plot. Defaults to the list of markers.
#' @param ncol number of columns for facetting.
#' @return Marker plot.
#' @export

gg_Marker_Pie <- function (
    xG,
    xY,
    trait,
    trait.label = trait,
    markers,
    marker.colors = gwaspr_Colors,
    title = NULL,
    subtitle = paste(markers, collapse = "\n"),
    ncol = NULL
    ) {
 #
 xY <- xY %>% dplyr::select(Name=1, myTrait=trait)
 xx <- xG %>% rename(SNP=1) %>%
   filter(SNP %in% markers) %>%
   dplyr::select(-2,-3,-4,-5,-6,-7,-8,-9,-10,-11) %>%
   column_to_rownames("SNP") %>%
   t() %>% as.data.frame() %>%
   select(markers) %>%
   mutate(Alleles = NA)
 for(i in 1:length(markers)) { xx <- xx[xx[,i] %in% c("A","T","G","C","AA","TT","GG","CC"),] }
 #
 for(i in 1:nrow(xx)) { xx$Alleles[i] <- paste(xx[i,1:length(markers)], collapse = "-") }
 #
 xx <- xx %>%
   rownames_to_column("Name") %>%
   left_join(xY, by = "Name") %>%
   filter(!is.na(myTrait)) %>%
   group_by(Alleles) %>%
   mutate(AlleleCount = n()) %>%
   group_by(Alleles, myTrait) %>%
   mutate(TraitCount = n(),
          myTrait = factor(myTrait)) %>%
   ungroup() %>%
   mutate(Percent = 100* TraitCount / AlleleCount) %>%
   filter(!duplicated(paste(Alleles, myTrait, TraitCount, AlleleCount, Percent)))
 #
 xx <- xx %>%
   group_by(Alleles) %>%
   mutate(TraitPos = cumsum(TraitCount),
          myTrait = factor(myTrait))
 #
 # Plot
 mp <- ggplot(xx, aes(x = "x", y = Percent, fill = myTrait)) +
   geom_col(alpha = 0.7) +
   geom_text(aes(label = TraitCount), position = position_stack(vjust = 0.5)) +
   geom_text(aes(label = paste("n =",AlleleCount), x = "x_empty"), y = 50) +
   coord_polar("y", start = 0) +
   facet_wrap(. ~ Alleles, scales = "free", ncol = ncol) +
   scale_fill_manual(name = trait.label, values = marker.colors) +
   scale_x_discrete(limits = c("x_empty", "x")) +
   theme_gwaspr_pie(legend.position = "bottom") +
   labs(title = title, subtitle = subtitle, y = NULL)
 mp
}

#xG = myG; xY = myY;
#trait = "Disease.Score_Ba16"
#markers = myMarkers[1]
#marker.colors = c("darkgreen", "darkgoldenrod3", "darkred", "steelblue4", "darkslategray", "maroon4", "purple4", "darkblue")
