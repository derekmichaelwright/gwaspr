# gg_Manhattan

Creates a manhattan plot.

## Usage

``` r
gg_Manhattan(
  folder = "GWAS_Results/",
  trait = list_Traits(folder)[1],
  title = trait,
  threshold = NULL,
  sug.threshold = NULL,
  chrom = NULL,
  markers = NULL,
  labels = markers,
  vlines = markers,
  vline.colors = rep("red", length(vlines)),
  vline.types = rep(1:6, length(vlines)),
  legend = T,
  legend.rows = 1,
  point.sizes = c(0.3, 1.25, 0.75),
  facet = F,
  addQQ = T,
  pmax = NULL,
  pmin = 0,
  models = c("MLM", "MLMM", "FarmCPU", "BLINK", "GLM", "CMLM", "SUPER"),
  model.colors = gwaspr_Colors,
  highlight.sig = F,
  sig.color = "darkred",
  chrom.colors = rep(c("darkgreen", "darkgoldenrod3"), 30),
  chrom.unit = "100 Mbp",
  plotHBPvalues = F,
  skyline = "Kansas"
)
```

## Arguments

- folder:

  Folder containing GWAS results.

- trait:

  The trait to read.

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

- point.sizes:

  Sizes for the points. c("Not Sig", "Sig", "Sug").

- facet:

  Logical, whether or not to produce a facetted or multi-model plot.
  Default is \`facet = F\`.

- addQQ:

  Logical, whether or not to add a QQ plot.

- pmax:

  A max value for the y-axis. Markers with higher values will be lowered
  to pmax..

- pmin:

  A min Value for plotting. Markers with lower values will be removed.

- models:

  Models to read.

- model.colors:

  Colors for each model. Used if \`facet = F\`.

- highlight.sig:

  Logical, whether or not to highlight significant associations with a
  black circle. Used if \`facet = F\`.

- sig.color:

  Color for significant assoctiations. Used if \`facet = T\`.

- chrom.colors:

  Colors for each chromosome. Used if \`facet = T\`.

- chrom.unit:

  Unit for the x-axis. Can be one of c("kbp","100 kbp","Mbp","100
  Mbp","Gbp").

- plotHBPvalues:

  Logical, if TRUE, H.B.P.values be uses.

- skyline:

  Which skyline type to use. Can be "NYC" or "Kansas". If NULL, it will
  use the highest P.value.

## Value

A manhattan plot.
