#' convert_IUPAC
#'
#' Converts hapmap formatted genotype files to and from IUPAC codes.
#' @param xG GWAS genotype object. Note: needs to be in hapmap format.
#' @param toIUPAC Logical, convert to IUPAC (TRUE) or from IUPAC (FALSE).
#' @return A LD Heatmap plot.
#' @export

convert_IUPAC <- function(xg, toIUPAC = T) {
  #
  dna <- data.frame(stringsAsFactors = F,
                    Symbol = c("A", "C", "G", "T", "R", "Y", "S", "W", "K", "M", "N", "R", "Y", "S", "W", "K", "M"),
                    Value  = c("AA","CC","GG","TT","AG","CT","GC","AT","GT","AC","NN", "GA","TC","CG","TA","TG","CA") )
  #
  if(toIUPAC == T) {
    # Convert to single letter IUPAC
    for(i in 12:ncol(xg)) { xg[,i] <- plyr::mapvalues(xg[,i], dna$Value, dna$Symbol) }
  }
  if(toIUPAC == F) {
    # Convert from IUPAC to double letter
    for(i in 12:ncol(xg)) { xg[,i] <- plyr::mapvalues(xg[,i], dna$Symbol, dna$Value) }
  }
  xg
}

#xg = myG_LDP

