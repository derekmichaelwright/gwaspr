# gg_GWAS_Hits

[Creates a summary GWAS plot of significant
associations.](https://derekmichaelwright.github.io/gwaspr/articles/gg_GWAS_Hits.html)

## Usage

``` r
gg_GWAS_Hits(
  folder = "GWAS_Results/",
  traits = list_Traits(folder),
  models = c("BLINK", "FarmCPU", "MLMM", "MLM", "GLM", "CMLM", "SUPER"),
  model.colors = gwaspr_ModelColors,
  model.shapes = c(21:25, 21:22),
  groups = NULL,
  chr = NULL,
  pos1 = NULL,
  pos2 = NULL,
  range = 2e+06,
  title = "",
  sigMin = 0,
  vlines = NULL,
  vline.colors = rep("red", length(vlines)),
  vline.types = rep(1, length(vlines)),
  vline.legend = T,
  legend.rows = 1,
  chr.unit = "100 Mbp",
  facet = F,
  caption = T,
  skyline = "NYC"
)
```

## Arguments

- traits:

  List of traits to use.

- models:

  Models to read.

- model.colors:

  Colors for each model.

- range:

  Range for binning GWAS hits.

- title:

  Title for horizontal facet.

- sigMin:

  Minimum number of hits to plot.

- vlines:

  Markers to be labelled with a vertical red line.

- vline.colors:

  colors for each vertical line.

- vline.types:

  lty for each vertical line.

- legend.rows:

  Number of rows for the legend.

- facet:

  Logical. facet by GWAS model.

- caption:

  Logical. should a caption be added.

- xx:

  Table of significant GWAS results. See ?table_GWAS_Results().

- xG:

  Genotype data in hapmap format.

- xCV:

  Optional, for filtering if you have a "CV" column.

## Value

A GWAS Hits plot.
