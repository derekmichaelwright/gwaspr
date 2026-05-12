# gg_Marker_Bar

Creates a marker plot with myG and myY objects.

## Usage

``` r
gg_Marker_Bar(
  xG,
  xY,
  traits,
  markers,
  marker.colors = gwaspr_Colors,
  plot.histogram = T,
  plot.density = T,
  plot.counts = T,
  myncol = NULL,
  line.color = F,
  title = NULL
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

- plot.histogram:

  Logical, if true will plot histogram bars.

- plot.density:

  Logical, if true will plot density bands.

- plot.counts:

  Logical, if true will make a plot of counts, if false will make a
  density plot.

- myncol:

  Number of columns for facetting when plotting multiple traits.

- title:

  Title for the plot.

## Value

Marker plot.
