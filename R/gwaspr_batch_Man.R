#' gwaspr_batch_Man
#'
#' [Creates Manhattan plots for all results in your specified folder.](https://derekmichaelwright.github.io/gwaspr/articles/gwaspr_batch_Man.html)
#' @param folder Folder containing GWAS results.
#' @param output Output folder for Manhattan plots.
#' @param traits The traits to read.
#' @param models Models to read.
#' @return Manhattan plots.
#' @export

gwaspr_batch_Man <- function(
    folder = "vignettes/GWAS_Results/",
    output = "vignettes/figures/",
    traits = list_Traits(folder),
    models =  c("MLM", "MLMM", "FarmCPU", "BLINK", "GLM", "CMLM", "SUPER")
    ) {
  #
  xx <- myResults <- table_GWAS_Results(folder = folder)
  # Create for loop i<-traits[1]
  for(i in traits) {
    # Prep data
    xi <- xx %>% filter(Trait == i) %>% #filter(!duplicated(SNP)) %>%
      group_by(Chr) %>%
      top_n(negLog10_P, n = 1) %>%
      arrange(Chr)
    # Plot
    mp <- gg_Manhattan(folder = folder, traits = i, models = models)
    # Save
    ggsave(paste0(output, "gwaspr_bm_", i, ".png"),
           mp, width = 12, height = 3.5, bg = "white" )
  }
}
