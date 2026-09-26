#' gg_Marker_Pie
#'
#' [Creates a marker plot with myG and myY objects.](https://derekmichaelwright.github.io/gwaspr/articles/gg_Marker.html)
#' @param xG GWAS genotype object. Note:  needs to be in hapmap format.
#' @param xY GWAS phenotype object.
#' @param trait Trait to plot.
#' @param trait.label Label for the trait.
#' @param trait.levels Factor levels for the trait.
#' @param markers Markers to plot.
#' @param marker.colors Color palette.
#' @param remove.hets Logical, Whether to remove hets or not. advised if plotting multiple markers.
#' @param title Title for the plot.
#' @param subtitle Subtitle for the plot. Defaults to the list of markers.
#' @param ncol number of columns for facetting.
#' @param legend.rows number of rows in legend.
#' @param groupByTrait Logical, if TRUE, will make pies of each trait instead of each marker.
#' @return Marker plot.
#' @export

gg_Marker_Pie <- function (
    xG,
    xY,
    trait,
    trait.label = trait,
    trait.levels = NULL,
    markers,
    marker.colors = gwaspr_Colors,
    remove.hets = T,
    title = NULL,
    subtitle = paste(markers, collapse = "\n"),
    ncol = NULL,
    legend.rows = 1,
    removeHets = T,
    groupByTrait = F
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
 #
 if(remove.hets == T) { for(i in 1:length(markers)) { xx <- xx[xx[,i] %in% c("A","T","G","C","AA","TT","GG","CC"),] } }
 if(remove.hets == F) { for(i in 1:length(markers)) { xx <- xx[!xx[,i] %in% c("N","NN"),] } }
 #
 for(i in 1:nrow(xx)) { xx$Alleles[i] <- paste(xx[i,1:length(markers)], collapse = "-") }
 #
 xx <- xx %>%
   rownames_to_column("Name") %>%
   left_join(xY, by = "Name") %>%
   filter(!is.na(myTrait)) %>%
   group_by(Alleles) %>%
   mutate(AlleleCount = n()) %>%
   ungroup() %>% group_by(Alleles, myTrait) %>%
   mutate(TraitAlleleCount = n(),
          myTrait = factor(myTrait)) %>%
   ungroup() %>% group_by(myTrait) %>%
   mutate(TraitCount = n()) %>%
   ungroup() %>%
   mutate(Percent = 100* TraitCount / AlleleCount) %>%
   filter(!duplicated(paste(Alleles, myTrait, TraitAlleleCount, AlleleCount, Percent)))
 #
 if(is.null(trait.levels)) {
   xx <- xx %>%
     group_by(Alleles) %>%
     mutate(TraitPos = cumsum(TraitCount),
            myTrait = factor(myTrait))
 }
 if(!is.null(trait.levels)) {
   xx <- xx %>%
     group_by(Alleles) %>%
     mutate(TraitPos = cumsum(TraitCount),
            myTrait = factor(myTrait, levels = trait.levels))
 }
 #
 # Plot
 if(groupByTrait == F) {
   mp <- ggplot(xx, aes(x = "x", y = Percent, fill = myTrait)) +
     geom_col(alpha = 0.7) +
     geom_text(aes(label = TraitAlleleCount), position = position_stack(vjust = 0.5)) +
     geom_text(aes(label = paste("n =",AlleleCount), x = "x_empty"), y = 50) +
     coord_polar("y", start = 0) +
     facet_wrap(. ~ Alleles, scales = "free", ncol = ncol) +
     scale_fill_manual(name = trait.label, values = marker.colors) +
     scale_x_discrete(limits = c("x_empty", "x")) +
     theme_gwaspr_pie(legend.position = "bottom") +
     guides(fill = guide_legend(nrow = legend.rows)) +
     labs(title = title, subtitle = subtitle, y = NULL)
 }
 if(groupByTrait == T) {
   mp <- ggplot(xx, aes(x = "x", y = Percent, fill = Alleles)) +
     geom_col(alpha = 0.7) +
     geom_text(aes(label = TraitAlleleCount), position = position_stack(vjust = 0.5)) +
     geom_text(aes(label = paste("n =", TraitCount), x = "x_empty"), y = 50) +
     coord_polar("y", start = 0) +
     facet_wrap(. ~ myTrait, scales = "free", ncol = ncol) +
     scale_fill_manual(name = trait.label, values = marker.colors) +
     scale_x_discrete(limits = c("x_empty", "x")) +
     theme_gwaspr_pie(legend.position = "bottom") +
     guides(fill = guide_legend(nrow = legend.rows)) +
     labs(title = title, subtitle = subtitle, y = NULL)
 }
 mp
}

#xG = myG; xY = myY;
#trait = "Disease.Score_Ba16"
#markers = myMarkers[1]
#marker.colors = c("darkgreen", "darkgoldenrod3", "darkred", "steelblue4", "darkslategray", "maroon4", "purple4", "darkblue")
