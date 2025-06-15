#' checkNYCvsKansas
#'
#' Check if NYC and Kansas results are identical, if so, the Kansas file will be deleted.
#' @param folder Folder containing GWAS results.
#' @param deteleKansas Logical, if TRUE, will delete any `Kansas` files with no difference between the `NYC` files.
#' @return A table of which runs were identical and deleted.
#' @export

checkNYCvsKansas <- function(folder = "GWAS_Results/", deleteKansas = F) {
  #
  myTraits <- list_Traits(folder)
  fnames <- list_Result_Files(folder)
  output <- data.frame(Trait = myTraits, FarmCPU = NA, BLINK = NA)
  #
  for(i in myTraits) {
    #
    # FarmCPU
    #
    f1 <- fnames[grepl(i,fnames) & grepl("\\(NYC\\)",fnames) & grepl("FarmCPU",fnames)]
    f2 <- fnames[grepl(i,fnames) & grepl("\\(Kansas\\)",fnames) & grepl("FarmCPU",fnames)]
    x1 <- read.csv(paste0(folder, f1)) %>% select(SNP, P.value) %>% arrange(SNP)
    x2 <- read.csv(paste0(folder, f2)) %>% select(SNP, P.value) %>% arrange(SNP)
    diffs <- sum(x1 != x2)
    output$FarmCPU[output$Trait == i] <- diffs
    #if(diffs == 0) { if(deleteKansas == T) { file.remove(paste0(folder, f2)) } }
    #
    # BLINK
    #
    f1 <- fnames[grepl(i,fnames) & grepl("\\(NYC\\)",fnames) & grepl("BLINK",fnames)]
    f2 <- fnames[grepl(i,fnames) & grepl("\\(Kansas\\)",fnames) & grepl("BLINK",fnames)]
    x1 <- read.csv(paste0(folder, f1)) %>% select(SNP, P.value) %>% arrange(SNP)
    x2 <- read.csv(paste0(folder, f2)) %>% select(SNP, P.value) %>% arrange(SNP)
    diffs <- sum(x1 != x2)
    output$BLINK[output$Trait == i] <- diffs
    #if(diffs == 0) { if(deleteKansas == T) { file.remove(paste0(folder, f2)) } }
  }
  #
  output
}

#i <- myTraits[1]
