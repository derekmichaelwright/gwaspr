# gg_Marker_Bar

[Creates a marker plot with myG and myY
objects.](https://derekmichaelwright.github.io/gwaspr/articles/gg_Marker.html)

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
  ncol = NULL,
  line.color = F,
  title = NULL,
  subtitle = paste(markers, collapse = "\n"),
  asfactor = F
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

- ncol:

  Number of columns for facetting when plotting multiple traits.

- line.color:

  Color of the border lines.

- title:

  Title for the plot.

- subtitle:

  Subtile for the plot. Defaults to the list of markers.

- asfactor:

  Logical, whether or not to plot as a factor or numeric.

## Value

Marker plot.
