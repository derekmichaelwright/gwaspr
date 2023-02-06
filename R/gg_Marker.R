#' gg_Marker
#'
#' Creates a marker plot with myG and myY objects.
#' @param myG GWAS genotype object. Note:
#' @param myY GWAS phenotype object.
#' @param trait Trait to plot.
#' @param marker Marker to plot.
#' @param marker2 Second marker to plot.
#' @param type Type of plot to make.
#' @param points Logical, whether or not to plot points
#' @param colors Color palette.
#' @return Marker plot.
#' @export

gg_Marker <- function(myG, myY, trait, marker, marker2 = NULL,
                      type = "bar", points = T,
                      colors = c("darkgreen", "darkgoldenrod3", "darkred", "steelblue4",
                                 "darkslategray", "maroon4", "purple4", "darkblue")) {
  gwaspr_dna <- data.frame(stringsAsFactors = F,
                           Symbol = c("A", "C", "G", "T", "U", "R", "Y", "S", "W", "K", "M", "N"),
                           Value  = c("AA","CC","GG","TT","UU","AG","CT","GC","AT","GT","AC","NN") )
  #
  x1 <- myG %>% filter(rs == marker)
  if(!is.null(marker2)) {
    x2 <- myG %>% filter(rs == marker2) %>%
      gather(Name, Allele2, 12:ncol(.)) %>% select(Name, Allele2)
  }
  fname <- paste(trait, marker, marker2, "png", sep = ".")
  title <- paste(trait, marker, marker2, sep = " | ")
  xx <- x1 %>% gather(Name, Allele, 12:ncol(.)) %>%
    select(Name, Allele) %>%
    mutate(Allele = plyr::mapvalues(Allele, gwaspr_dna$Symbol, gwaspr_dna$Value),
           Allele = factor(Allele, levels = gwaspr_dna$Value)) %>%
    filter(Allele != "NN")
  if(!is.null(marker2)) { xx <- xx %>% left_join(x2, by = "Name") }
  xx <- xx %>% left_join(myY, by = "Name") %>%
    filter(!is.na(get(trait))) %>%
    mutate(Allele = plyr::mapvalues(Allele, gwaspr_dna$Symbol, gwaspr_dna$Value))
  if(type == "bar") {
    mp <- ggplot(xx, aes(x = get(trait), fill = get(trait))) +
      geom_bar(color = "black", alpha = 0.7) +
      facet_grid(. ~ Allele) +
      scale_fill_manual(name = NULL, values = colors) +
      theme_gwaspr(legend.position = "none") +
      labs(title = title, x = NULL)
  }
  if(type == "box") {
    if(!is.null(marker2)) {
      mp <- ggplot(xx, aes(x = Allele, y = get(trait))) +
        geom_violin(fill = "grey90") +
        geom_boxplot(width = 0.1) +
        geom_beeswarm(aes(color = Allele2), alpha = 0.8) +
        scale_color_manual(name = NULL, values = colors) +
        theme_gwaspr() +
        labs(title = title, y = trait, x = NULL)

    }else {
      mp <- ggplot(xx, aes(x = Allele, y = get(trait))) +
        geom_violin(fill = "grey90") +
        geom_boxplot(width = 0.1) +
        theme_gwaspr() +
        labs(title = title, x = NULL)
      if(points == T) {
        mp <- mp + geom_beeswarm(aes(color = Allele), alpha = 0.5) +
          scale_color_manual(name = NULL, values = colors)
      }
    }
  }
  if(type == "double bar") {
    xx <- xx %>% mutate(Allele2 = plyr::mapvalues(Allele2, gwaspr_dna$Symbol, gwaspr_dna$Value))
    mp <- ggplot(xx, aes(x = get(trait), fill = Allele2)) +
      geom_bar(position = position_dodge(preserve = "single"), color = "black") +
      facet_grid(. ~ Allele) +
      scale_fill_manual(name = marker2, values = colors) +
      theme_gwaspr(legend.position = "bottom") +
      labs(title = title, x = NULL)
  }
  if(type == "double box") {
    xx <- xx %>% mutate(Allele2 = plyr::mapvalues(Allele2, gwaspr_dna$Symbol, gwaspr_dna$Value))
    mp <- ggplot(xx, aes(x = Allele, y = get(trait), color = Allele2)) +
      geom_violin(fill = "grey70") +
      geom_boxplot(width = 0.1, position = "dodge") +
      geom_beeswarm(aes(color = Allele2), alpha = 0.8) +
      scale_color_manual(name = NULL, values = colors) +
      labs(title = title, x = NULL)
  }
  mp
}
