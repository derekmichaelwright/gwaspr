#' gg_Marker
#'
#' Creates a marker plot with myG and myY objects.
#' @param xG GWAS genotype object. Note:  needs to be in hapmap format.
#' @param xY GWAS phenotype object.
#' @param trait Trait to plot.
#' @param markers Markers to plot.
#' @param marker.colors Color palette.
#' @return Marker plot.
#' @export

gg_Marker_Pie <- function (
    xG,
    xY,
    trait,
    myncol = NULL,
    markers,
    marker.colors = c("darkgreen", "darkgoldenrod3", "darkred", "steelblue4",
                      "darkslategray", "maroon4", "purple4", "darkblue") 
    ) {
 #
 title <- paste(markers, collapse = "\n")
 #
 xY <- xY %>% select(1, myTrait=trait)
 xx <- xG %>% rename(SNP=1) %>%
   filter(SNP %in% markers) %>%
   select(-2,-3,-4,-5,-6,-7,-8,-9,-10,-11) %>%
   column_to_rownames("SNP") %>%
   t() %>% as.data.frame() %>% 
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
 # Plot
 mp <- ggplot(xx, aes(x = "", y = Percent)) +
   geom_col(aes(fill = myTrait), stat = "identity",  color = "black", alpha = 0.5) +
   coord_polar("y", start = 0) +
   facet_grid(. ~ paste0(Alleles,"\nn = ", AlleleCount), scales = "free") +
   scale_fill_manual(name = NULL, values = marker.colors) +
   theme_gwaspr_pie(legend.position = "bottom") +
   labs(title = title, y = NULL, x = NULL)
 mp
}

#xG = myG; xY = myY;  markers = myMarkers[1]
#trait = "Disease.Score_Ba16"
#marker.colors = c("darkgreen", "darkgoldenrod3", "darkred", "steelblue4", "darkslategray", "maroon4", "purple4", "darkblue") 
#line.color=F; myncol = NULL; plot.histogram = T; plot.density = T
