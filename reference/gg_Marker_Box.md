# gg_Marker_Box

Creates a marker plot with myG and myY objects.

## Usage

``` r
gg_Marker_Box(
  xG,
  xY,
  traits,
  markers,
  marker.colors = gwaspr_Colors,
  plot.violin = T,
  plot.points = T,
  plot.box = T,
  box.width = 0.1,
  point.size = 1,
  myncol = NULL,
  title = NULL,
  legend.rows = 1,
  subtitle = paste(markers, collapse = "\n"),
  yLab = traits,
  xCV = NULL,
  cv.colors = NULL,
  cv.label = NULL
)
```

## Arguments

- xG:

  GWAS genotype object. Note: needs to be in hapmap format.

- xY:

  GWAS phenotype object.

- traits:

  Traits to plot.

- markers:

  Markers to plot.

- marker.colors:

  Color palette.

- plot.violin:

  Logical, whether or not to plot violin.

- plot.points:

  Logical, whether or not to plot points.

- box.width:

  Width for the boxplot.

- point.size:

  Size for the points.

- myncol:

  Number of columns for facetting when plotting multiple traits.

- title:

  Title for the plot.

- legend.rows:

  Number of rows for the legend.

## Value

Marker plot.
