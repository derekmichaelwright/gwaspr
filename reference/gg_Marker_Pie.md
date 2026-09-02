# gg_Marker_Pie

[Creates a marker plot with myG and myY
objects.](https://derekmichaelwright.github.io/gwaspr/articles/gg_Marker.html)

## Usage

``` r
gg_Marker_Pie(
  xG,
  xY,
  trait,
  trait.label = trait,
  markers,
  marker.colors = gwaspr_Colors,
  title = NULL,
  subtitle = paste(markers, collapse = "\n"),
  ncol = NULL
)
```

## Arguments

- xG:

  GWAS genotype object. Note: needs to be in hapmap format.

- xY:

  GWAS phenotype object.

- trait:

  Trait to plot.

- trait.label:

  Label for the Trait.

- markers:

  Markers to plot.

- marker.colors:

  Color palette.

- title:

  Title for the plot.

- subtitle:

  Subtitle for the plot. Defaults to the list of markers.

- ncol:

  number of columns for facetting.

## Value

Marker plot.
