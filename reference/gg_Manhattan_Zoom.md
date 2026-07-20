# gg_Manhattan_Zoom

[Create manhattan plots from GAPIT GWAS results zoomed into a specific
region on a
chromosome.](https://derekmichaelwright.github.io/gwaspr/articles/05_gg_Manhattan_Zoom.html)

## Usage

``` r
gg_Manhattan_Zoom(
  folder = "GWAS_Results/",
  trait = list_Traits(folder)[1],
  title = trait,
  chr = 1,
  pos1 = NULL,
  pos2 = NULL,
  addGenome = T,
  threshold = NULL,
  sug.threshold = NULL,
  markers = NULL,
  labels = markers,
  vlines = markers,
  vline.colors = rep("red", length(vlines)),
  vline.types = rep(1, length(vlines)),
  vline.legend = T,
  facet = F,
  pmax = NULL,
  models = c("MLM", "MLMM", "FarmCPU", "BLINK", "GLM", "CMLM", "SUPER"),
  model.colors = c("darkgreen", "darkred", "darkorange3", "steelblue", "darkorchid4",
    "blue2", "magenta3"),
  sig.color = "black",
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

- trait:

  The trait to read.

- title:

  A title for the plot.

- chr:

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

- facet:

  Logical, whether or not to produce a facetted or multi-model plot.
  Default is `facet = F`.

- pmax:

  A max value for the y-axis.

- models:

  Models to read.

- model.colors:

  Colors for each model. Used if `facet = F`.

- sig.color:

  Color for significant assoctiations.

- legend.rows:

  Number of rows for the legend.

- legend.box:

  Alignment of the legend. Default is "horizontal", but it can be
  changed to "vertical".

- point.sizes:

  Sizes for the points. c("Not Sig", "Sig", "Sug").

- plotHBPvalues:

  Logical, should H.B.P.Values be uses.

- skyline:

  Which skyline type to use. Can be "NYC" or "Kansas". If left NULL, it
  will use the highest P.value.

## Value

A manhattan plot.
