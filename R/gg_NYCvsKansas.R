#' gg_NYCvsKansas
#'
#' Creates a manhattan plot comparing the NYC and Kansas results.
#' @param folder Folder containing GWAS results.
#' @param trait The trait to read.
#' @param title A title for the plot.
#' @param threshold Significant Threshold.
#' @param chr.colors Colors for each chromosome.
#' @param chr.unit Unit for the x-axis. Can be one of c("kbp","100 kbp","Mbp","100 Mbp","Gbp").
#' @return A manhattan plot.
#' @export

gg_NYCvsKansas <- function (
    folder = "GWAS_Results/",
    trait = list_Traits(folder)[1],
    title = trait,
    threshold = NULL,
    chr.colors = rep(c("darkgreen", "darkgoldenrod3"), 30),
    chr.unit = "100 Mbp"
    ) {
  #
  # Read in files
  #
  fnames <- list_Result_Files(folder)
  fnames <- fnames[grepl("BLINK|FarmCPU", fnames)]
  fnames <- fnames[grepl(paste0(trait, c(".csv","\\("), collapse="|"), fnames)]
  #
  xx <- NULL
  for (i in fnames) {
    mod <- substr(i, gregexpr("GWAS_Results", i)[[1]][1] + 13, gregexpr(".csv", i)[[1]][1] - 1)
    mod <- substr(mod, 1, gregexpr("\\.", mod)[[1]][1] - 1)
    sky <- substr(i, gregexpr("\\(", i)[[1]][1] + 1, gregexpr("\\)", i)[[1]][1] - 1)
    #
    xi <- read.csv(paste0(folder, i))
    if (sum(colnames(xi) == "nobs") > 0) { xi <- dplyr::select(xi, -nobs) }
    if (sum(colnames(xi) == "effect") > 0) { xi <- rename(xi, Effect=effect) }
    xi <- xi %>%
      mutate(Model = mod,
             negLog10_P = -log10(P.value),
             negLog10_Exp = -log10((rank(P.value, ties.method = "first") - 0.5)/nrow(.)),
             Type = sky )
    xx <- bind_rows(xx, xi)
  }
  #
  xx  <- xx %>% filter(!duplicated(paste(SNP, Model, P.value)))
  #
  # Prep data
  #
  #
  if (is.null(threshold)) {
    threshold <- -log10(0.05/nrow(xi))
  }
  #
  xx <- xx %>% mutate(Sig.level = ifelse(negLog10_P >= threshold, "Sig","Not Sig"))
  x2 <- xx %>% filter(negLog10_P > threshold)
  #
  if(chr.unit == "kbp")     { x.unit = 1000 }
  if(chr.unit == "100 kbp") { x.unit = 100000 }
  if(chr.unit == "Mbp")     { x.unit = 1000000 }
  if(chr.unit == "100 Mbp") { x.unit = 100000000 }
  if(chr.unit == "Gbp")     { x.unit = 1000000000 }
  if(!chr.unit %in% c("kbp", "100 kbp", "Mbp", "100 Mbp", "Gbp")) { print("error in chr.unit") }
  #
  myBreaks <- 0:(round(max(xx$Pos)/x.unit))
  #
  # Start Plots
  #
  mp1 <- ggplot(xx, aes(x = Pos/x.unit, y = negLog10_P)) +
    theme_gwaspr(axis.title.y = element_markdown()) +
    labs(title = title, y = "-log<sub>10</sub>(*p*)", x = chr.unit)
  #
  mp2 <- ggplot(xx, aes(y = negLog10_P, x = negLog10_Exp)) +
    theme_gwaspr() +
    labs(title = "", y = NULL, x = "Expected")
  #
  # Add threshold lines
  #
  mp1 <- mp1 +
    geom_hline(yintercept = threshold, color = "red", alpha = 0.8, linewidth = 0.5)
  mp2 <- mp2 +
    geom_hline(yintercept = threshold, color = "red", alpha = 0.8, linewidth = 0.5)
  #
  # Plot facetted by model
  #
  mp1 <- mp1 +
    geom_point(aes(fill = factor(Chr), size = Sig.level), pch = 21, color = alpha("white", 0)) +
    geom_point(data = x2, pch = 21, size = 1.5, color = "black", fill = "darkred", alpha = 0.8) +
    facet_grid(Model+Type ~ Chr, scales = "free", space = "free_x") +
    scale_fill_manual(name = NULL, values = alpha(chr.colors, 0.8), guide = "none") +
    scale_size_manual(name = NULL, values = c(0.4,1.25,0.75), guide = "none") +
    scale_x_continuous(breaks = myBreaks, minor_breaks = myBreaks) +
    #guides(color = guide_legend(nrow = legend.rows, byrow = T, override.aes = list(alpha = 1))) +
    theme(legend.position = "bottom")
  #
  mp2 <- mp2 +
    geom_point(pch = 1, color = chr.colors[1], alpha = 0.8) +
    geom_point(data = x2, pch = 21, color = "black", fill = "darkred", alpha = 0.8) +
    geom_abline() +
    facet_grid(Model+Type ~ "QQ", scales = "free_y")
  #
  mp <- ggarrange(mp1, mp2, ncol = 2, widths = c(4,1), align = "h",
                  legend = "bottom", common.legend = T)
  #
  # Output Plot
  #
  mp
}

#folder = "vignettes/GWAS_Results/"; trait = "DTF_Nepal_2017";
#title = trait; threshold = NULL;
#chr.colors = rep(c("darkgreen", "darkgoldenrod3"), 30); chr.unit = "100 Mbp"
