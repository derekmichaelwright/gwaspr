% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gg_GWAS_Summary.R
\name{gg_GWAS_Summary}
\alias{gg_GWAS_Summary}
\title{gg_GWAS_Summary}
\usage{
gg_GWAS_Summary(
  folder = NULL,
  traits = list_Traits(),
  threshold = -log10(5e-08),
  sug.threshold = -log10(5e-06),
  models = c("MLM", "MLMM", "FarmCPU", "BLINK", "GLM"),
  colors = c("darkgreen", "darkred", "darkorange3", "steelblue", "darkgoldenrod2"),
  shapes = 21:25,
  hlines = NULL,
  vlines = NULL,
  vline.colors = rep("red", length(vlines)),
  vline.legend = T,
  title = NULL,
  caption = paste0("Sig Threshold = ", threshold, " = Large\\nSuggestive = ",
    sug.threshold, " = Small"),
  rowread = 2000,
  legend.position = "bottom",
  lrows = 1
)
}
\arguments{
\item{folder}{Folder containing GWAS results.}

\item{traits}{The traits to read.}

\item{threshold}{Significant threshold.}

\item{sug.threshold}{Suggestive threshold.}

\item{models}{Models to read.}

\item{colors}{Colors for each model.}

\item{shapes}{The shape values to use for the different models. e.g., 21:25}

\item{hlines}{Locations for horizontal lines. e.g., hlines = c(1.5,2.5).}

\item{vlines}{Markers to be labelled with a vertical red line.}

\item{vline.colors}{colors for each vertical line.}

\item{vline.legend}{Logical, display of vline color legend.}

\item{title}{A title for the plot.}

\item{caption}{A caption for the plot.}

\item{rowread}{Number of rows to read for each GWAS results file.}

\item{lrows}{Number of rows for the legend.}
}
\value{
A GWAS summary plot.
}
\description{
Creates a summary GWAS plot of significant associations.
}
