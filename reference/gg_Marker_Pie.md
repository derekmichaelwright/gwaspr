# gg_Marker_Pie

Creates a marker plot with myG and myY objects.

## Usage

``` r
gg_Marker_Pie(
  xG,
  xY,
  trait,
  markers,
  marker.colors = gwaspr_Colors,
  title = NULL
)
```

## Arguments

- xG:

  GWAS genotype object. Note: needs to be in hapmap format.

- xY:

  GWAS phenotype object.

- trait:

  Trait to plot.

- markers:

  Markers to plot.

- marker.colors:

  Color palette.

- title:

  Title for the plot.

## Value

Marker plot.
