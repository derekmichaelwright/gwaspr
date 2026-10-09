# gg_Marker_Bubble

[Creates a marker scatter plot with myG and myY
objects.](https://derekmichaelwright.github.io/gwaspr/articles/gg_Marker_Bubble.html)

## Usage

``` r
gg_Marker_Bubble(
  xG,
  xY,
  trait.x = NULL,
  trait.y = "Alleles",
  facet.rows = NULL,
  facet.cols = NULL,
  markers,
  marker.colors = gwaspr_Colors,
  remove.hets = F,
  remove.N = T,
  legend.pos = "bottom",
  legend.rows = 1,
  title = NULL,
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

- facet.rows:

  Trait for the row facetting.

- facet.cols:

  Trait for the col facetting.

- markers:

  Markers to plot.

- marker.colors:

  Colors to fill the marker bubbles.

- remove.hets:

  Logical, Whether to remove hets or not. advised if plotting multiple
  markers.

- remove.N:

  Logical, Whether to remove N or NN.

- legend.pos:

  Position for the legend.

- legend.rows:

  Number of rows for the legend.

- title:

  Title for the plot.

- subtitle:

  Subtitle for the plot. Defaults to the list of markers.

## Value

Marker plot.
