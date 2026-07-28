# gg_LD_Heatmap

[Creates a LD heatmap
plot.](https://derekmichaelwright.github.io/gwaspr/articles/gg_LD_Heatmap.html)

## Usage

``` r
gg_LD_Heatmap(
  xG = myG,
  chr,
  pos1,
  pos2,
  metric = "R^2",
  threshold = 0.9,
  myMs = NULL,
  myTitle = NULL,
  axisTextSize = NULL,
  nameTrim = NULL,
  color.low = "white",
  color.mid = "goldenrod1",
  color.high = "darkred"
)
```

## Arguments

- xG:

  GWAS genotype object. Note: needs to be in hapmap format.

- chr:

  Chromosome to plot.

- pos1:

  Start position within the selected chromosome.

- pos2:

  End position within the selected chromosome.

- metric:

  Which LD calculation to use. Default is "D'".

- threshold:

  Value for selecting linked markers. default is "0.9".

- myMs:

  Markers to highlight within the plot.

- myTitle:

  Title for the plot.

- axisTextSize:

  Text size for the axis labels (genotype names).

- nameTrim:

  String used to trim marker names.

- color.low:

  Color for gradient low.

- color.mid:

  Color for gradient mid.

- color.high:

  Color for gradient high.

## Value

A LD Heatmap plot.
