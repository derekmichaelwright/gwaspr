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
  trait.levels = NULL,
  markers,
  marker.colors = gwaspr_Colors,
  remove.hets = T,
  title = NULL,
  subtitle = paste(markers, collapse = "\n"),
  ncol = NULL,
  legend.rows = 1,
  removeHets = T,
  groupByTrait = F
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

  Label for the trait.

- trait.levels:

  Factor levels for the trait.

- markers:

  Markers to plot.

- marker.colors:

  Color palette.

- remove.hets:

  Logical, Whether to remove hets or not. advised if plotting multiple
  markers.

- title:

  Title for the plot.

- subtitle:

  Subtitle for the plot. Defaults to the list of markers.

- ncol:

  number of columns for facetting.

- legend.rows:

  number of rows in legend.

- groupByTrait:

  Logical, if TRUE, will make pies of each trait instead of each marker.

## Value

Marker plot.
