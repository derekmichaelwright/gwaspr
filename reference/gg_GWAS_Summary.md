# gg_GWAS_Summary

Creates a summary GWAS plot of significant associations. Note: this
function requires the GWAS results files to be ordered

## Usage

``` r
gg_GWAS_Summary(
  folder = "GWAS_Results/",
  traits = list_Traits(folder),
  groups = NULL,
  threshold = round(-log10(5e-08), 1),
  sug.threshold = round(-log10(5e-06), 1),
  chr = NULL,
  pos1 = NULL,
  pos2 = NULL,
  models = c("MLM", "MLMM", "FarmCPU", "BLINK", "CMLM", "SUPER"),
  model.colors = gwaspr_Colors,
  shapes = 21:25,
  hlines = NULL,
  vlines = NULL,
  vline.colors = rep("red", length(vlines)),
  vline.types = rep(1, length(vlines)),
  vline.legend = T,
  threshold.legend = T,
  title = "Summary of Significant GWAS Results",
  caption = NULL,
  rowread = 2000,
  legend.position = "bottom",
  legend.rows = 1,
  plotHBPvalues = F,
  skyline = NULL
)
```

## Arguments

- folder:

  Folder containing GWAS results.

- traits:

  The traits to read.

- groups:

  Grouping for the traits. Should be equal length to `traits`.

- threshold:

  Significant threshold.

- sug.threshold:

  Suggestive threshold.

- chr:

  Chromosomes to plot.

- pos1:

  starting position to plot.

- pos2:

  ending position to plot.

- models:

  Models to read.

- model.colors:

  Colors for each model.

- shapes:

  The shape values to use for the different models. e.g., 21:25.

- hlines:

  Locations for horizontal lines. e.g., hlines = c(1.5,2.5).

- vlines:

  Markers to be labelled with a vertical red line.

- vline.colors:

  colors for each vertical line.

- vline.types:

  lty for each vertical line.

- vline.legend:

  Logical, display of vline color legend.

- threshold.legend:

  Logical, display of vline color legend.

- title:

  A title for the plot.

- caption:

  A caption for the plot.

- rowread:

  Number of rows to read for each GWAS results file.

- legend.rows:

  Number of rows for the legend.

- plotHBPvalues:

  Logical, should H.B.P.Values be uses.

- skyline:

  Which skyline type to use. Can be "NYC" or "Kansas". If left NULL, it
  will use the highest P.value.

## Value

A GWAS summary plot.
