% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gg_GWAS_Hits.R
\name{gg_GWAS_Hits}
\alias{gg_GWAS_Hits}
\title{gg_GWAS_Hits}
\usage{
gg_GWAS_Hits(
  xx,
  myG,
  myTs,
  myR = 2e+06,
  myTitle = "",
  sigMin = 0,
  myCV = NULL,
  vlines = NULL,
  vline.colors = rep("red", length(vlines))
)
}
\arguments{
\item{xx}{Table of significant GWAS results. See ?table_GWAS_Results().}

\item{myG}{Genotype data.}

\item{myTs}{List of traits to use.}

\item{myR}{Range for binning GWAS hits.}

\item{myTitle}{Title for horizontal facet.}

\item{sigMin}{Minimum number of hits to plot.}

\item{myCV}{(optional) for filtering if you have a "CV" column.}

\item{vlines}{Markers to be labelled with a vertical red line.}

\item{vline.colors}{colors for each vertical line.}
}
\value{
A GWAS Hits plot.
}
\description{
Creates a summary GWAS plot of significant associations.
}
