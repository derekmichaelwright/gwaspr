# gg_QTL_Summary

Creates a summary QTL plot of significant associations.

## Usage

``` r
gg_QTL_Summary(
  myG,
  myQ,
  title = "Summary of QTL Results",
  lodFill = F,
  fillColor = "darkgreen",
  fillColor_low = "grey50",
  xLab = "Pos",
  facetLab = "lg"
)
```

## Arguments

- myG:

  Genetic map. Must contain the following 3 columns
  c("Marker","Chr","Pos").

- myQ:

  Table of QTL results. Must contain the following 5 columns
  c("Trait","Chr","Pos","Pos_lo","Pos_hi","lod").

- title:

  Custom title for the plot.

- lodFill:

  Logical. If true, fill color will be a gradient based on lod values.

- fillColor:

  Color for filling points.

- fillColor_low:

  If `lodfill = T`, This will be the low end of the gradeint color scale
  for filling points.

- xLab:

  Custom x lab.

- facetLab:

  Custom label for facets.

## Value

A QTL summary plot.
