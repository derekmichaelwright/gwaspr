# gg_QTL_Summary

Creates a summary QTL plot of significant associations. Note: this
function requires the GWAS results files to be ordered

## Usage

``` r
gg_QTL_Summary(
  xx = myQTLs,
  myG = myG,
  title = "Summary of QTL Results",
  caption = "derek was here"
)
```

## Arguments

- xx:

  Table of QTL results.

- myG:

  genotype file

- groups:

  Grouping for the traits. Should be equal length to \`traits\`.

## Value

A QTL summary plot.
