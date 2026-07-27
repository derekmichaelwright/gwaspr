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
  remove.hets = T,
  plot.violin = T,
  plot.box = T,
  plot.points = T,
  box.width = 0.1,
  point.size = 1,
  point.beeswarm = F,
  myncol = NULL,
  title = NULL,
  legend.rows = 1,
  subtitle = paste(markers, collapse = "\n"),
  yLab = traits,
  cv.source = "xG",
  cv.name = NULL,
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

  Colors to fill in the violin and boxplots.

- remove.hets:

  Logical, Whether to remove hets or not. advisded if plotting multiple
  markers.

- plot.violin:

  Logical, whether or not to plot violins.

- plot.box:

  Logical, whether or not to plot the boxplots.

- plot.points:

  Logical, whether or not to plot points.

- box.width:

  Width for the boxplot.

- point.size:

  Size for the points.

- point.beeswarm:

  Logical. If False (the default), will plot points with
  `geom_quasirandom`. If TRUE, will plot points with `geom_beeswarm`.

- myncol:

  Number of columns for facetting when plotting multiple traits.

- title:

  Title for the plot.

- legend.rows:

  Number of rows for the legend.

- subtitle:

  Subtitle for the plot. Defaults to the list of markers.

- yLab:

  Label for the y-axis.

- cv.source:

  Where to get your `cv.name` from. Default is "xG", while "xY" is the
  other option.

- cv.name:

  Covariable data for points.

- cv.colors:

  Covariable colors for filling points.

- cv.label:

  Label for the covariate.

## Value

Marker plot.
