# gg_Marker_Scatter

[Creates a marker scatter plot with myG and myY
objects.](https://derekmichaelwright.github.io/gwaspr/articles/gg_Marker.html)

## Usage

``` r
gg_Marker_Scatter(
  xG,
  xY,
  trait.x,
  trait.y,
  markers,
  marker.colors = gwaspr_Colors,
  remove.hets = T,
  point.size = 1,
  title = NULL,
  legend.rows = 1,
  subtitle = paste(markers, collapse = "\n")
)
```

## Arguments

- xG:

  GWAS genotype object. Note: needs to be in hapmap format.

- xY:

  GWAS phenotype object.

- trait.x:

  Trait for the x-axis.

- trait.y:

  Trait for the y-axis.

- markers:

  Markers to plot.

- marker.colors:

  Colors to fill in the violin and boxplots.

- remove.hets:

  Logical, Whether to remove hets or not. advised if plotting multiple
  markers.

- point.size:

  Size for the points.

- title:

  Title for the plot.

- legend.rows:

  Number of rows for the legend.

- subtitle:

  Subtitle for the plot. Defaults to the list of markers.

- yLab:

  Label for the y-axis.

## Value

Marker plot.
