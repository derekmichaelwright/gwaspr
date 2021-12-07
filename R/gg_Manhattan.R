#' gg_Manhattan
#'
#' Creates a manhattan plot.
#' @param folder Folder containing GWAS results.
#' @param trait The trait to read.
#' @param title A title for the plot.
#' @param markers Markers to be labelled.
#' @param labels Labels to be used for markers.
#' @param lines Logical, value of whether or not to include vertical lines with markers.
#' @param facet Logical, whether or not to produce a facetted or multi-model plot.
#' @param qq Logical, whether or not to add a QQ plot
#' @param pmax A max value for the y-axis.
#' @param models Models to read.
#' @param colors1 Colors for each chromosome.
#' @param colors2 Colors for each model.
#' @return A manhattan plot.
#' @export

gg_Manhattan <- function(folder, trait, title = trait, threshold = NULL, threshold2 = NULL,
                         markers = NULL, labels = markers,
                         lines = F, facet = T, qq = T, pmax = NULL,
                         models = c("GLM","MLM","CMLM","MLMM","SUPER","FarmCPU","Blink"),
                         colors1 = c("darkgreen","darkgoldenrod3","darkgreen","darkgoldenrod3",
                                     "darkgreen","darkgoldenrod3","darkgreen"),
                         colors2 = c("darkgreen", "darkred", "darkorange3",
                                     "steelblue", "darkorchid4", "darkgoldenrod2")) {
  #
  fnames <- grep(paste0(trait,".GWAS.Results"), list.files(folder))
  fnames <- list.files(folder)[fnames]
  xx <- NULL
  # i <- fnames[3]
  for(i in fnames) {
    mod <- substr(i, gregexpr("\\.", i)[[1]][1]+1, gregexpr("\\.", i)[[1]][2]-1)
    xi <- read.csv(paste0(folder, i))
    if(sum(colnames(xi)=="nobs")>0) { xi <- select(xi, -nobs) }
    xi <- xi %>%
      mutate(Model = mod,
             `-log10(p)`     = -log10(P.value),
             `-log10(p)_exp` = -log10((rank(P.value, ties.method="first")-.5)/nrow(.)))
    xx <- bind_rows(xx, xi)
  }
  xx <- xx %>% filter(Model %in% models) %>%
    mutate(Model = factor(Model, levels = models)) %>%
    arrange(desc(Model))
  #
  if(is.null(threshold)) { threshold <- -log10(0.05 / nrow(xi)) }
  if(!is.null(pmax)) {
    xx <- xx %>% mutate(`-log10(p)` = ifelse(`-log10(p)` > pmax, pmax, `-log10(p)`))
  }
  x1 <- xx %>% filter(`-log10(p)` < threshold)
  x2 <- xx %>% filter(`-log10(p)` > threshold)
  # Man plot
  mp1 <- ggplot(x1, aes(x = Position / 1000000, y = `-log10(p)`))
  if(!is.null(markers) & lines == T) {
    mp1 <- mp1 +
      geom_vline(data = xx %>% filter(SNP %in% markers),
                 aes(xintercept = Position / 1000000), alpha = 0.5)
  }
  if(facet == T) {
    mp1 <- mp1 +
      geom_hline(yintercept = threshold, color = "red") +
      geom_hline(yintercept = threshold2, color = "blue") +
      geom_point(aes(color = factor(Chromosome)), pch = 1) +
      geom_point(data = x2, pch = 21, size = 1.5, color = "black", fill = "darkred") +
      facet_grid(Model ~ Chromosome, scales = "free", space = "free_x") +
      scale_x_continuous(breaks = seq(0, 2000, by = 100)) +
      scale_color_manual(values = colors1) +
      theme_gwaspr(legend.position = "none",
                   axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5),
                   axis.title.y = element_markdown()) +
      labs(title = title, y = "-log<sub>10</sub>(*p*)", x = "Mbp")
    if(!is.null(markers)) {
      xx <- xx %>% mutate(Label = ifelse(SNP %in% markers, SNP, NA),
                          Label = plyr::mapvalues(Label, markers, labels))
      mp1 <- mp1 + geom_text_repel(data = xx %>% filter(SNP %in% markers),
                                   aes(label = Label), size = 2)
    }
    # QQ plot
    mp2 <- ggplot(x1, aes(y = `-log10(p)`, x = `-log10(p)_exp`)) +
      geom_hline(yintercept = threshold, color = "red") +
      geom_hline(yintercept = threshold2, color = "blue") +
      geom_point(pch = 1, color = colors1[1], alpha = 0.8) +
      geom_point(data = x2, pch = 21, color = "black", fill = "darkred", alpha = 0.8) +
      geom_abline() +
      facet_grid(Model ~ "QQ", scales = "free_y") +
      theme_gwaspr() +
      labs(title = "", y = NULL, x = "Expected")
    # Append plots
    if(qq == T) { mp <- ggpubr::ggarrange(mp1, mp2, ncol = 2, widths = c(4,1))
    } else { mp <- mp1 }
  } else {
    mp1 <- mp1 +
      geom_hline(yintercept = threshold, color = "red") +
      geom_hline(yintercept = threshold2, color = "blue") +
      geom_point(size = 0.1, aes(color = Model), pch = 1) +
      geom_point(data = x2, aes(color = Model), alpha = 0.8) +
      facet_grid(. ~ Chromosome, scales = "free", space = "free_x") +
      scale_x_continuous(breaks = seq(0, 2000, by = 100)) +
      scale_color_manual(values = colors2) +
      theme_gwaspr(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5),
                   axis.title.y = element_markdown()) +
      labs(title = title, y = "-log<sub>10</sub>(*p*)", x = "Mbp")
    if(!is.null(markers)) {
      xx <- xx %>% mutate(Label = ifelse(SNP %in% markers, SNP, NA),
                          Label = plyr::mapvalues(Label, markers, labels))
      mp1 <- mp1 + geom_text_repel(data = xx %>% filter(SNP %in% markers),
                                   aes(label = Label), size = 2)
    }
    # QQ plot
    mp2 <- ggplot(x1, aes(y = `-log10(p)`, x = `-log10(p)_exp`)) +
      geom_hline(yintercept = threshold, color = "red") +
      geom_hline(yintercept = threshold2, color = "blue") +
      geom_point(pch = 1, aes(color = Model)) +
      geom_point(data = x2, aes(color = Model)) +
      geom_abline() +
      facet_grid(. ~ "QQ", scales = "free_y") +
      scale_color_manual(values = colors2) +
      theme_gwaspr() +
      labs(title = "", y = NULL, x = "Expected")
    # Append plots
    if(qq == T) { mp <- ggpubr::ggarrange(mp1, mp2, ncol = 2, widths = c(4,1), align = "h",
                                          legend = "bottom", common.legend = T)
    } else { mp <- mp1 }
  }
  mp
}

