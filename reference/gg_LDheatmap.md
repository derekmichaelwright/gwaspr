# gg_LDHeatmap

Creates a manhattan plot.

## Usage

``` r
gg_LDheatmap(
  xg = myG,
  chr = 6,
  pos1 = 0,
  pos2 = 6e+06,
  myMs = NULL,
  myTitle = NULL,
  axisTextSize = NULL,
  nameTrim = NULL
)
```

## Arguments

- chr:

  Chromosome to plot.

- pos1:

  Start position within the selected chromosome.

- pos2:

  End position within the selected chromosome.

- myMs:

  Markers to highlight within the plot.

- myTitle:

  Title for the plot.

- axisTextSize:

  Text size for the axis labels (genotype names).

- nameTrim:

  String used to trim marker names.

- xG:

  GWAS genotype object. Note: needs to be in hapmap format.

## Value

A LD Heatmap plot.
