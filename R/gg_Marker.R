#' gg_Manh
#'
#' Removes columns without any values.
#' @param trait .
#' @param marker
#' @param marker2
#' @param type
#' @param points
#' @param colors
#' @return Plot.
#' @export

gg_Marker <- function(trait, marker, marker2 = NULL, 
                      type = "bar", points = T,
                      colors = c("darkgreen", "darkgoldenrod3", "darkred", "steelblue4",
                                 "darkslategray", "maroon4", "purple4", "darkblue")) {
  x1 <- myG %>% filter(rs == marker)
  if(!is.null(marker2)) { 
    x2 <- myG %>% filter(rs == marker2) %>% 
      gather(Name, Allele2, 12:ncol(.)) %>% select(Name, Allele2)
  }
  fname <- paste(trait, marker, marker2, "png", sep = ".")
  title <- paste(trait, marker, marker2, sep = " | ")
  xx <- x1 %>% gather(Name, Allele, 12:ncol(.)) %>%
    select(Name, Allele) %>%
    mutate(Allele = plyr::mapvalues(Allele, dna$Symbol, dna$Value),
           Allele = factor(Allele, levels = dna$Value)) %>%
    filter(Allele != "NN")
  if(!is.null(marker2)) { xx <- xx %>% left_join(x2, by = "Name") }
  xx <- xx %>% left_join(myY, by = "Name") %>%
    filter(!is.na(get(trait))) %>%
    mutate(Allele = plyr::mapvalues(Allele, dna$Symbol, dna$Value))
  if(type == "bar") {
    mp <- ggplot(xx, aes(x = get(trait), fill = get(trait))) + 
      geom_bar(color = "black") + 
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
        labs(title = title, x = NULL)
      
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
    xx <- xx %>% mutate(Allele2 = plyr::mapvalues(Allele2, dna$Symbol, dna$Value))
    mp <- ggplot(xx, aes(x = get(trait), fill = Allele2)) + 
      geom_bar(position = position_dodge(preserve = "single"), color = "black") + 
      facet_grid(. ~ Allele) +
      scale_fill_manual(name = marker2, values = colors) +
      theme_gwaspr(legend.position = "bottom") +
      labs(title = title, x = NULL)
  }
  if(type == "double box") {
    xx <- xx %>% mutate(Allele2 = plyr::mapvalues(Allele2, dna$Symbol, dna$Value))
    mp <- ggplot(xx, aes(x = Allele, y = get(trait), color = Allele2)) + 
      geom_violin(fill = "grey70") + 
      geom_boxplot(width = 0.1, position = "dodge") +
      geom_beeswarm(aes(color = Allele2), alpha = 0.8) +
      scale_color_manual(name = NULL, values = colors) +
      labs(title = title, x = NULL)
  }
  mp
}