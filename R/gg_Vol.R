#' gg_Vol
#'
#' Creates a volcano plot.
#' @param folder Folder containing GWAS results.
#' @param trait The trait to read.
#' @param subtitle A subtitle for the plot.
#' @param markers Markers to be labelled.
#' @param labels Labels to be used for markers.
#' @param models Models to read.
#' @return A volcano plot.
#' @export

gg_Vol <- function(folder, 
                   trait, 
                   subtitle = NULL,
                   markers = NULL, 
                   labels = markers,
                   models = c("GLM","MLM","MLMM","CMLM","SUPER","FarmCPU","Blink")) {
  fnames <- grep(paste0(trait,".GWAS.Results"), list.files(folder))
  fnames <- list.files(folder)[fnames]
  xx <- NULL
  # i <- fnames[1]
  for(i in fnames) {
    mod <- substr(i, gregexpr("\\.", i)[[1]][1]+1, gregexpr("\\.", i)[[1]][2]-1)
    xi <- read.csv(paste0(folder, i))
    if(sum(colnames(xi)=="nobs")>0) { xi <- select(xi, -nobs) }
    xi <- xi %>% 
      mutate(Model = mod,
             negLog10     = -log10(P.value),
             negLog10_exp = -log10((rank(P.value, ties.method="first")-.5)/nrow(.)),
             Sig = ifelse(negLog10 > -log10(0.05/nrow(.)), "Significant", "Non-Significant"))
    xx <- bind_rows(xx, xi)
  }
  xx <- xx %>% filter(Model %in% models) %>%
    mutate(Model = factor(Model, levels = models)) %>% 
    filter(!is.na(Model), !is.na(effect))
  xm <- xx %>% filter(SNP %in% markers) %>% 
    mutate(SNP = plyr::mapvalues(SNP, markers, labels))
  threshold <- -log10(0.05/nrow(xi))
  # Volcano plot
  ggplot(xx, aes(x = effect, y = -log10(P.value))) + 
    geom_hline(yintercept = threshold) +
    geom_point(aes(color = Sig, shape = factor(Chromosome))) +
    geom_point(data = xm, aes(shape = factor(Chromosome)), color = "darkred") +
    geom_text_repel(data = xm, size = 1, aes(label = SNP)) +
    facet_wrap(Model ~ ., scales = "free", ncol = length(unique(xx$Model))) +
    scale_color_manual(values = c("darkgreen", "black"), guide = F) +
    scale_fill_manual(values = c(alpha("darkgreen", 0), "darkred")) +
    scale_shape_manual(name = "Chr", values = c(1:7)) +
    guides(color = F) +
    theme_gwaspr() +
    labs(title = paste(trait, subtitle, sep = " | "), y = "-log10(p)")
}