% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gg_Manhattan_Zoom.R
\name{gg_Manhattan_Zoom}
\alias{gg_Manhattan_Zoom}
\title{gg_Manhattan_Zoom}
\usage{
gg_Manhattan_Zoom(
  folder = "GWAS_Results/",
  trait = list_Traits(folder)[1],
  chrom,
  pos1,
  pos2,
  title = trait,
  threshold = NULL,
  sug.threshold = NULL,
  markers = NULL,
  labels = markers,
  vlines = markers,
  vline.colors = "red",
  sig.col = "red",
  models = c("MLM", "FarmCPU", "BLINK", "MLMM", "GLM", "CMLM", "SUPER"),
  model.colors = c("darkgreen", "darkorange3", "steelblue", "darkred", "darkorchid4",
    "burlywood4", "darkseagreen4"),
  facet = F,
  highlight.sig = F,
  highlight.marker.color = "red",
  legend.rows = 1,
  plotHBPvalues = F,
  skyline = NULL
)
}
\arguments{
\item{folder}{Folder containing GWAS results.}

\item{trait}{The trait to read.}

\item{chrom}{Chromosome to plot.}

\item{pos1}{Start position on chromosome.}

\item{pos2}{End position on chromosome.}

\item{title}{A title for the plot.}

\item{threshold}{Significant Threshold.}

\item{sug.threshold}{Suggested threshold.}

\item{markers}{Markers to be labelled.}

\item{labels}{Labels to be used for markers.}

\item{vlines}{Markers which will be used as a location for a vertical lines.}

\item{vline.colors}{colors for each vertical line.}

\item{sig.col}{Color for significant assoctiations.}

\item{models}{Models to read.}

\item{facet}{Logical, whether or not to produce a facetted or multi-model plot. Default is `facet = F`.}

\item{highlight.sig}{Logical, whether or not to highlight significant associations with a black circle. Used if `facet = F`.}

\item{legend.rows}{Number of rows for the legend.}

\item{plotHBPvalues}{Logical, should H.B.P.Values be uses.}

\item{skyline}{Which skyline type to use. Can be "NYC" or "Kansas". If left NULL, it will use the highest P.value.}
}
\value{
A manhattan plot.
}
\description{
Creates a manhattan plot zoomed in to a particular region.
}
