# gg_Manhattan_Zoom1

Creates a manhattan plot zoomed in to a particular region.

## Usage

``` r
gg_Manhattan_Zoom1(
  folder = "GWAS_Results/",
  traits = list_Traits(folder)[1],
  title = NULL,
  chrom,
  pos1,
  pos2,
  threshold = NULL,
  sug.threshold = NULL,
  markers = NULL,
  labels = markers,
  vlines = markers,
  vline.colors = rep("red", length(vlines)),
  vline.types = rep(1, length(vlines)),
  vline.legend = T,
  pmax = NULL,
  model = "MLM",
  model.colors = gwaspr_Colors,
  facet = F,
  highlight.sig = F,
  sig.color = "red",
  legend.rows = 1,
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

- chrom:

  Chromosome to plot.

- pos1:

  Start position on chromosome.

- pos2:

  End position on chromosome.

- threshold:

  Significant Threshold.

- sug.threshold:

  Suggested threshold.

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

- vline.legend:

  Logical, whether or not to add a legend for the vlines.

- pmax:

  A max value for the y-axis.

- model:

  Model to read.

- model.colors:

  Colors for each model. Used if \`facet = F\`.

- facet:

  Logical, whether or not to produce a facetted or multi-model plot.
  Default is \`facet = F\`.

- highlight.sig:

  Logical, whether or not to highlight significant associations with a
  black circle. Used if \`facet = F\`.

- sig.color:

  Color for significant assoctiations.

- legend.rows:

  Number of rows for the legend.

- plotHBPvalues:

  Logical, should H.B.P.Values be uses.

- skyline:

  Which skyline type to use. Can be "NYC" or "Kansas". If left NULL, it
  will use the highest P.value.

- addQQ:

  Logical, whether or not to add a QQ plot

- chrom.colors:

  Colors for each chromosome. Used if \`facet = T\`.

- chrom.unit:

  Unit for the x-axis. Can be one of c("kbp","100 kbp","Mbp","100
  Mbp","Gbp").

## Value

A manhattan plot.
