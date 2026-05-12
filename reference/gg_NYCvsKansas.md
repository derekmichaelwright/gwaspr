# gg_NYCvsKansas

Creates a manhattan plot comparing the NYC and Kansas results.

## Usage

``` r
gg_NYCvsKansas(
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
  vline.types = rep(1, length(vlines)),
  vline.legend = T,
  addQQ = T,
  pmax = NULL,
  models = c("FarmCPU", "BLINK"),
  sig.col = "darkred",
  chrom.colors = rep(c("darkgreen", "darkgoldenrod3"), 30),
  chrom.unit = "100 Mbp",
  legend.rows = 1,
  plotHBPvalues = F
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

- vline.legend:

  Logical, whether or not to add a legend for the vlines.

- addQQ:

  Logical, whether or not to add a QQ plot

- pmax:

  A max value for the y-axis.

- models:

  Models to read.

- sig.col:

  Color for significant assoctiations. Used if \`facet = T\`.

- chrom.colors:

  Colors for each chromosome. Used if \`facet = T\`.

- chrom.unit:

  Unit for the x-axis. Can be one of c("kbp","100 kbp","Mbp","100
  Mbp","Gbp").

- legend.rows:

  Number of rows for the legend.

- plotHBPvalues:

  Logical, if TRUE, H.B.P.Values be uses.

## Value

A manhattan plot.
