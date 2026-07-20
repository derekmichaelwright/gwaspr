# gg_QTL_Summary_Groups

Creates a summary QTL plot of significant associations.

## Usage

``` r
gg_QTL_Summary_Groups(
  myQ,
  myG,
  title = "Summary of QTL Results",
  yGroup,
  yLab = NULL,
  facetGroup,
  colorGroup,
  colorName = NULL,
  fillColors = gwaspr_Colors,
  xLab = "Pos",
  facetLab = "lg"
)
```

## Arguments

- myQ:

  Table of QTL results. Must contain the following 5 columns
  c("Trait","Chr","Pos","Pos_lo","Pos_hi","lod"). then it must also
  include grouping columns set by the user for `yGroup`,`facetGroup` &
  `colorGroup`.

- myG:

  Genetic map. Must contain the following 3 columns
  c("Marker","Chr","Pos").

- title:

  Custom title for the plot.

- yGroup:

  Name of column for y axis grouping.

- yLab:

  Custom y-axis label.

- facetGroup:

  Name of column for facet grouping.

- colorGroup:

  Name of column for color grouping.

- colorName:

  Title for the color legend.

- fillColors:

  Color for filling points.

- xLab:

  Custom x-axis label.

- facetLab:

  Custom label for facets.

## Value

A QTL summary plot.
