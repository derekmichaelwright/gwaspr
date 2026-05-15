#' clean_GWAS_Results
#'
#' Delete extra GWAS files in the folder containing GWAS result files from GAPIT.
#' @param folder Folder containing GWAS results.
#' @param remove_nonResults Delelets all non .csv files (all .pdf files)
#' @param remove_Kansas Delete Kansas files.
#' @param remove_NYC Delete NYC files.
#' @return A table of which models have been run on each trait.
#' @export

clean_GWAS_Results <- function(folder = "GWAS_Results/",
                               remove_nonResults = F,
                               remove_Kansas = F,
                               remove_NYC = F) {
  #
  if(remove_nonResults == T) {
    to_be_deleted <- list.files(folder)[grepl(".pdf", list.files(folder))]
    if(length(to_be_deleted) > 0) { file.remove(paste0(folder, to_be_deleted)) }
  }
  #
  if(remove_Kansas == T) {
    to_be_deleted <- list.files(folder)[grepl("\\(Kansas\\)", list.files(folder))]
    if(length(to_be_deleted) > 0) { file.remove(paste0(folder, to_be_deleted)) }
  }
  #
  if(remove_NYC == T) {
    to_be_deleted <- list.files(folder)[grepl("\\(NYC\\)", list.files(folder))]
    if(length(to_be_deleted) > 0) { file.remove(paste0(folder, to_be_deleted)) }
  }
}
