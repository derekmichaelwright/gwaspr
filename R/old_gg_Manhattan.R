#' old_gg_Manhattan
#'
#' Creates a manhattan plot. Note: this function is meant for the old version of GAPIT.
#' @param folder Folder containing GWAS results.
#' @param trait The trait to read.
#' @param title A title for the plot.
#' @param markers Markers to be labelled.
#' @param labels Labels to be used for markers.
#' @param vlines Markers which will be used as a location for a vertical lines.
#' @param vline.colors colors for each vertical line.
#' @param vline.legend Logical, wheterh or not to add a legend for the vlines.
#' @param facet Logical, whether or not to produce a facetted or multi-model plot.
#' @param qq Logical, whether or not to add a QQ plot
#' @param pmax A max value for the y-axis.
#' @param models Models to read.
#' @param chrom.colors Colors for each chromosome. Used if `facet = T`.
#' @param model.colors Colors for each model. Used if `facet = F`.
#' @return A manhattan plot.
#' @export

old_gg_Manhattan <- function(folder, trait, title = trait, threshold = NULL, sug.threshold = NULL,
                         markers = NULL, labels = markers,
                         vlines = markers, vline.colors = rep("red",length(vlines)), vline.legend = F,
                         facet = T, qq = T, pmax = NULL,
                         models = c("MLM","MLMM","FarmCPU","Blink","GLM"),
                         chrom.colors = c("darkgreen","darkgoldenrod3","darkgreen","darkgoldenrod3",
                                          "darkgreen","darkgoldenrod3","darkgreen"),
                         model.colors = c("darkgreen", "darkred", "darkorange3",
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
  mp1 <- ggplot(x1, aes(x = Pos / 100000000, y = `-log10(p)`))
  if(!is.null(vlines)) {
    mp1 <- mp1 +
      geom_vline(data = xx %>% filter(SNP %in% vlines), alpha = 0.4,
                 aes(xintercept = Pos / 100000000, color = SNP))
    if(vline.legend == T) {
      mp1 <- mp1 + scale_color_manual(name = NULL, values = vline.colors)
    } else {
      mp1 <- mp1 + scale_color_manual(name = NULL, values = vline.colors, guide = F)
    }

  }
  if(facet == T) {
    mp1 <- mp1 +
      geom_hline(yintercept = threshold, color = "red", alpha = 0.8, size = 0.5) +
      geom_hline(yintercept = sug.threshold, color = "blue", alpha = 0.8, size = 0.5) +
      geom_point(aes(fill = factor(Chr)), pch = 21, size = 1, color = alpha("white", 0)) +
      geom_point(data = x2, pch = 21, size = 1.5, color = "black", fill = "darkred", alpha = 0.8) +
      facet_grid(Model ~ Chr, scales = "free", space = "free_x") +
      scale_fill_manual(values = chrom.colors) +
      scale_x_continuous(breaks = 0:20) +
      theme_gwaspr(legend.position = "none",
                   axis.title.y = element_markdown()) +
      labs(title = title, y = "-log<sub>10</sub>(*p*)", x = "100 Mbp")
    if(!is.null(markers)) {
      xx <- xx %>% mutate(Label = ifelse(SNP %in% markers, SNP, NA),
                          Label = plyr::mapvalues(Label, markers, labels)) %>%
        filter(`-log10(p)` > min(threshold, sug.threshold))
      mp1 <- mp1 + geom_text_repel(data = xx %>% filter(SNP %in% markers),
                                   aes(label = Label), size = 2)
    }
    # QQ plot
    mp2 <- ggplot(x1, aes(y = `-log10(p)`, x = `-log10(p)_exp`)) +
      geom_hline(yintercept = threshold, color = "red", alpha = 0.8, size = 0.5) +
      geom_hline(yintercept = sug.threshold, color = "blue", alpha = 0.8, size = 0.5) +
      geom_point(pch = 1, color = chrom.colors[1], alpha = 0.8) +
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
      geom_hline(yintercept = threshold, color = "red", alpha = 0.8, size = 0.5) +
      geom_hline(yintercept = sug.threshold, color = "blue", alpha = 0.8, size = 0.5) +
      geom_point(size = 0.1, aes(fill = Model), pch = 21, color = alpha("white", 0)) +
      geom_point(data = x2, aes(fill = Model), pch = 21, size = 1.25, alpha = 0.8) +
      facet_grid(. ~ Chr, scales = "free", space = "free_x") +
      scale_fill_manual(values = model.colors) +
      scale_x_continuous(breaks = 0:20) +
      theme_gwaspr(axis.title.y = element_markdown()) +
      labs(title = title, y = "-log<sub>10</sub>(*p*)", x = "100 Mbp")
    if(!is.null(markers)) {
      xx <- xx %>% mutate(Label = ifelse(SNP %in% markers, SNP, NA),
                          Label = plyr::mapvalues(Label, markers, labels)) %>%
        filter(`-log10(p)` > min(threshold, sug.threshold))
      mp1 <- mp1 + geom_text_repel(data = xx %>% filter(SNP %in% markers),
                                   aes(label = Label), size = 2)
    }
    # QQ plot
    mp2 <- ggplot(x1, aes(y = `-log10(p)`, x = `-log10(p)_exp`)) +
      geom_hline(yintercept = threshold, color = "red", alpha = 0.8, size = 0.5) +
      geom_hline(yintercept = sug.threshold, color = "blue", alpha = 0.8, size = 0.5) +
      geom_point(pch = 1, aes(color = Model)) +
      geom_point(data = x2, aes(color = Model)) +
      geom_abline() +
      facet_grid(. ~ "QQ", scales = "free_y") +
      scale_color_manual(values = model.colors) +
      theme_gwaspr() +
      labs(title = "", y = NULL, x = "Expected")
    # Append plots
    if(qq == T) {
      mp <- ggpubr::ggarrange(mp1, mp2, ncol = 2, widths = c(4,1), align = "h",
                              legend = "bottom", common.legend = T)
    } else { mp <- mp1 }
  }
  mp
}
