# gg_GWAS_Hits

Creates a summary GWAS plot of significant associations.

## Usage

``` r
gg_GWAS_Hits(
  folder = "GWAS_Results/",
  xx,
  xG,
  xCV = NULL,
  traits = list_Traits(folder),
  groups = NULL,
  chr = NULL,
  pos1 = NULL,
  pos2 = NULL,
  range = 2e+06,
  title = "",
  sigMin = 0,
  models = c("BLINK", "FarmCPU", "MLMM", "MLM", "GLM", "CMLM", "SUPER"),
  model.colors = gwaspr_ModelColors,
  model.shapes = c(21:25, 21:22),
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

- xx:

  Table of significant GWAS results. See ?table_GWAS_Results().

- xG:

  Genotype data in hapmap format.

- xCV:

  Optional, for filtering if you have a "CV" column.

- traits:

  List of traits to use.

- range:

  Range for binning GWAS hits.

- title:

  Title for horizontal facet.

- sigMin:

  Minimum number of hits to plot.

- models:

  Models to read.

- model.colors:

  Colors for each model.

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

## Value

A GWAS Hits plot.
