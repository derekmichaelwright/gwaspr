# gg_Manhattan_xTraits

Creates a manhattan plot.

## Usage

``` r
gg_Manhattan_xTraits(
  folder = "GWAS_Results/",
  traits = list_Traits(folder)[1:2],
  title = NULL,
  threshold = NULL,
  sug.threshold = NULL,
  chr = NULL,
  markers = NULL,
  labels = markers,
  vlines = markers,
  vline.colors = rep("red", length(vlines)),
  vline.types = rep(1, length(vlines)),
  vline.legend = F,
  addQQ = T,
  pmax = NULL,
  pmin = 0,
  models = c("MLM", "MLMM", "FarmCPU", "BLINK", "GLM", "CMLM", "SUPER"),
  trait.colors = c("darkgreen", "darkred", "darkorange3", "steelblue", "darkorchid4",
    "blue2", "magenta3"),
  chr.unit = "100 Mbp",
  legend.rows = 1,
  legend.box = "horizontal",
  point.sizes = c(0.3, 1, 0.75),
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

- chr:

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

  Logical, whether or not to add a QQ plot.

- pmax:

  A max value for the y-axis.

- pmin:

  A min Value for plotting. Markers with lower values will be removed.

- models:

  Model to read.

- trait.colors:

  Colors for each trait.

- chr.unit:

  Unit for the x-axis. Can be one of c("kbp","100 kbp","Mbp","100
  Mbp","Gbp").

- legend.rows:

  Number of rows for the legend.

- legend.box:

  Alignment of the legend. Default is "horizontal", but it can be
  changed to "vertical".

- point.sizes:

  Sizes for the points. c("Not Sig", "Sig", "Sug").

- skyline:

  Which skyline type to use. Can be "NYC" or "Kansas". If left NULL, it
  will use the highest P.value.

## Value

A manhattan plot.
