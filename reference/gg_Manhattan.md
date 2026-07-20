# gg_Manhattan

[Create manhattan plots from GAPIT GWAS
results.](https://derekmichaelwright.github.io/gwaspr/articles/04_gg_Manhattan.html)

## Usage

``` r
gg_Manhattan(
  folder = "GWAS_Results/",
  trait = list_Traits(folder)[1],
  title = trait,
  threshold = NULL,
  sug.threshold = NULL,
  chr = NULL,
  markers = NULL,
  labels = markers,
  vlines = markers,
  vline.colors = rep("red", length(vlines)),
  vline.types = rep(1, length(vlines)),
  legend = T,
  legend.rows = 1,
  legend.box = "horizontal",
  point.sizes = c(0.3, 1, 0.75),
  pmax = NULL,
  pmin = 0,
  models = c("MLM", "MLMM", "FarmCPU", "BLINK", "GLM", "CMLM", "SUPER"),
  model.colors = c("darkgreen", "darkred", "darkorange3", "steelblue", "darkorchid4",
    "blue2", "magenta3"),
  facet = F,
  highlight.sig = F,
  sig.color = "black",
  chr.colors = rep(c("darkgreen", "darkgoldenrod3"), 30),
  chr.unit = "100 Mbp",
  plotHBPvalues = F,
  skyline = "Kansas",
  addQQ = T
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

- legend:

  Logical, whether or not to add a legend.

- legend.rows:

  Number of rows for the legend.

- legend.box:

  Alignment of the legend. Default is "horizontal", but it can be
  changed to "vertical".

- point.sizes:

  Sizes for the points. c("Not Sig", "Sig", "Sug").

- pmax:

  A max value for the y-axis. Markers with higher values will be lowered
  to pmax..

- pmin:

  A min Value for plotting. Markers with lower values will be removed.

- models:

  Models to read.

- model.colors:

  Colors for each model. Used if `facet = F`.

- facet:

  Logical, whether or not to produce a facetted or multi-model plot.
  Default is `facet = F`.

- highlight.sig:

  Logical, whether or not to highlight significant associations with a
  black circle. Used if `facet = F`.

- sig.color:

  Color for significant assoctiations. Used if `facet = T`.

- chr.colors:

  Colors for each chromosome. Used if `facet = T`.

- chr.unit:

  Unit for the x-axis. Can be one of c("kbp","100 kbp","Mbp","100
  Mbp","Gbp").

- plotHBPvalues:

  Logical, if TRUE, H.B.P.values be uses.

- skyline:

  Which skyline type to use. Can be "NYC" or "Kansas". If NULL, it will
  use the highest P.value.

- addQQ:

  Logical, whether or not to add a QQ plot.

## Value

A manhattan plot.
