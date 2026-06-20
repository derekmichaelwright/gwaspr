# gg_Manhattan_xModels

Creates a manhattan plot.

## Usage

``` r
gg_Manhattan_xModels(
  folder = "GWAS_Results/",
  traits = list_Traits(folder)[1:2],
  title = NULL,
  threshold = NULL,
  sug.threshold = NULL,
  chrom = NULL,
  markers = NULL,
  labels = markers,
  vlines = markers,
  vline.colors = rep("red", length(vlines)),
  vline.types = rep(1:6, length(vlines)),
  legend = F,
  legend.rows = 1,
  addQQ = T,
  pmax = NULL,
  pmin = 0,
  models = "MLM",
  trait.colors = gwaspr_Colors,
  chrom.unit = "100 Mbp",
  plotHBPvalues = F,
  skyline = "Kansas"
)
```

## Arguments

- folder:

  Folder containing GWAS results.

- traits:

  The traits to read.

- title:

  A title for the plot.

- threshold:

  Significant Threshold.

- sug.threshold:

  Suggested threshold.

- chrom:

  Chromosomes to plot. Use if you want to plot a single chromosome.

- markers:

  Markers to be labelled.

- labels:

  Labels to be used for markers.

- vlines:

  Markers which will be used as a location for a vertical lines.

- vline.colors:

  colors for each vertical line.

- vline.types:

  lty for each vertical line.

- legend:

  Logical, whether or not to add a legend.

- legend.rows:

  Number of rows for the legend.

- addQQ:

  Logical, whether or not to add a QQ plot.

- pmax:

  A max value for the y-axis.

- pmin:

  A min Value for plotting. Markers with lower values will be removed.

- models:

  Model to read.

- trait.colors:

  Colors for each trait.

- chrom.unit:

  Unit for the x-axis. Can be one of c("kbp","100 kbp","Mbp","100
  Mbp","Gbp").

- skyline:

  Which skyline type to use. Can be "NYC" or "Kansas". If left NULL, it
  will use the highest P.value.

## Value

A manhattan plot.
