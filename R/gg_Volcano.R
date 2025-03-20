#' gg_Volcano
#'
#' Creates a volcano plot.
#' @param folder Folder containing GWAS results.
#' @param trait The trait to read.
#' @param title A title for the plot.
#' @param markers Markers to be labelled.
#' @param labels Labels to be used for markers.
#' @param models Models to read.
#' @return A volcano plot.
#' @export

gg_Volcano <- function(folder, trait, title = trait,
                       markers = NULL, labels = markers,
                       models = c("MLM", "FarmCPU", "BLINK", "MLMM", "GLM")) {
  fnames <- list.files(folder)[grepl("GWAS_Results", list.files(folder))]
  fnames <- fnames[grepl(paste0(trait,".csv"), fnames)]
  xx <- NULL
  # i <- fnames[1]
  for(i in fnames) {
    mod <- substr(i, gregexpr("GWAS_Results", i)[[1]][1]+13, gregexpr(".csv", i)[[1]][1]-1)
    mod <- substr(mod, 1, gregexpr("\\.", mod)[[1]][1]-1)
    xi <- read.csv(paste0(folder, i))
    if(sum(colnames(xi)=="nobs")>0) { xi <- select(xi, -nobs) }
    xi <- xi %>%
      mutate(Model = mod,
             `-log10(p)`     = -log10(P.value),
             `-log10(p)_exp` = -log10((rank(P.value, ties.method="first")-.5)/nrow(.)),
             Sig = ifelse(`-log10(p)` > -log10(0.05/nrow(.)), "Significant", "Non-Significant"))
    xx <- bind_rows(xx, xi)
  }
  xx <- xx %>% filter(Model %in% models) %>%
    mutate(Model = factor(Model, levels = models)) %>%
    filter(!is.na(Model), !is.na(Effect))
  xm <- xx %>% filter(SNP %in% markers) %>%
    mutate(SNP = plyr::mapvalues(SNP, markers, labels))
  threshold <- -log10(0.05/nrow(xi))
  # Volcano plot
  ggplot(xx, aes(x = Effect, y = `-log10(p)`)) +
    geom_hline(yintercept = threshold, color = "red") +
    geom_point(aes(color = Sig, shape = factor(Chr)), alpha = 0.7) +
    geom_text_repel(data = xm, size = 1, aes(label = SNP)) +
    facet_wrap(Model ~ ., scales = "free", ncol = length(unique(xx$Model))) +
    scale_color_manual(values = c("darkgreen", "black")) +
    scale_shape_manual(name = "Chr", values = c(1:7)) +
    guides(color = F) +
    theme_gwaspr(axis.title.y = element_markdown()) +
    labs(title = title, y = "-log<sub>10</sub>(*p*)")
}
