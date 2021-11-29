#' gg_Manh
#'
#' Creates a manhattan plot.
#' @param folder Folder containing GWAS results.
#' @param trait The trait to read.
#' @param subtitle A subtitle for the plot.
#' @param markers Markers to be labelled.
#' @param labels Labels to be used for markers.
#' @param lines Logical value of whether or not to include vertical lines with markers.
#' @param models Models to read.
#' @param colors Colors for each chromosome
#' @return A manhattan plot.
#' @export

gg_Manh <- function(folder = NULL,
                   trait = NULL,
                   subtitle = NULL,
                   markers = NULL,
                   labels = markers,
                   lines = F,
                   facet = T,
                   models = c("GLM","MLM","CMLM","MLMM","SUPER","FarmCPU","Blink"),
                   colors = c("darkgreen","darkgoldenrod3","darkgreen","darkgoldenrod3",
                              "darkgreen","darkgoldenrod3","darkgreen")) {
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
             negLog10     = -log10(P.value),
             negLog10_exp = -log10((rank(P.value, ties.method="first")-.5)/nrow(.)))
    xx <- bind_rows(xx, xi)
  }
  xx <- xx %>% filter(Model %in% models) %>%
    mutate(Model = factor(Model, levels = models)) %>%
    arrange(desc(Model))
  threshold <- -log10(0.05 / nrow(xi))
  x1 <- xx %>% filter(negLog10 < threshold)
  x2 <- xx %>% filter(negLog10 > threshold)
  # Man plot
  mp1 <- ggplot(x1, aes(x = Position / 1000000, y = -log10(P.value)))
  if(!is.null(markers) & lines == T) {
    mp1 <- mp1 +
      geom_vline(data = xx %>% filter(SNP %in% markers),
                 aes(xintercept = Position / 1000000), alpha = 0.5)
  }
  if(facet == T) {
    mp1 <- mp1 +
      geom_hline(yintercept = threshold) +
      geom_point(aes(color = factor(Chromosome)), pch = 1, alpha = 0.8) +
      geom_point(data = x2, pch = 21, size = 1.5, color = "black", fill = "darkred", alpha = 0.8) +
      facet_grid(Model ~ Chromosome, scales = "free") +
      scale_x_continuous(breaks = seq(100, 700, by = 100)) +
      scale_color_manual(values = colors) +
      theme_gwaspr(legend.position = "none",
                   axis.text.x = element_text(angle = 90, hjust = 0.5)) +
      labs(title = paste(trait, subtitle), y = "-log10(p)", x = "Mbp")
    if(!is.null(markers)) {
      xx <- xx %>% mutate(Label = ifelse(SNP %in% markers, SNP, NA),
                          Label = plyr::mapvalues(Label, markers, labels))
      mp1 <- mp1 + geom_text_repel(data = xx %>% filter(SNP %in% markers),
                                   aes(label = Label), size = 2)
    }
    # QQ plot
    mp2 <- ggplot(x1, aes(y = negLog10, x = negLog10_exp)) +
      geom_point(pch = 1, color = colors[1], alpha = 0.8) +
      geom_point(data = x2, pch = 21, color = "black", fill = "darkred", alpha = 0.8) +
      geom_hline(yintercept = threshold) + geom_abline() +
      facet_grid(Model ~ "QQ", scales = "free_y") +
      theme_gwaspr() +
      labs(title = "", y = NULL, x = "Expected")
    # Append and save plots
    mp <- ggpubr::ggarrange(mp1, mp2, ncol = 2, widths = c(4,1))
  } else {
    colors2 <- c("darkred", "darkorange3", "steelblue", "darkgreen", "darkorchid4", "darkgoldenrod2")
    mp1 <- mp1 +
      geom_hline(yintercept = threshold) +
      geom_point(size = 0.1, aes(color = Model), pch = 1) +
      geom_point(data = x2, aes(color = Model)) +
      facet_grid(. ~ Chromosome, scales = "free") +
      scale_x_continuous(breaks = seq(100, 700, by = 100)) +
      scale_color_manual(values = colors2) +
      theme_gwaspr(axis.text.x = element_text(angle = 90, hjust = 0.5)) +
      labs(title = paste(trait, subtitle), y = "-log10(p)", x = "Mbp")
    if(!is.null(markers)) {
      xx <- xx %>% mutate(Label = ifelse(SNP %in% markers, SNP, NA),
                          Label = plyr::mapvalues(Label, markers, labels))
      mp1 <- mp1 + geom_text_repel(data = xx %>% filter(SNP %in% markers),
                                   aes(label = Label), size = 2)
    }
    # QQ plot
    mp2 <- ggplot(x1, aes(y = negLog10, x = negLog10_exp)) +
      geom_point(pch = 1, aes(color = Model)) +
      geom_point(data = x2, aes(color = Model)) +
      geom_hline(yintercept = threshold) + geom_abline() +
      facet_grid(. ~ "QQ", scales = "free_y") +
      scale_color_manual(values = colors2) +
      theme_gwaspr() +
      labs(title = "", y = NULL, x = "Expected")
    # Append and save plots
    mp <- ggpubr::ggarrange(mp1, mp2, ncol = 2, widths = c(4,1), align = "h",
                            legend = "bottom", common.legend = T)
  }
  mp
}

