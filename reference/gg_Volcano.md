# gg_Volcano

Creates a volcano plot.

## Usage

``` r
gg_Volcano(
  folder = "GWAS_Results/",
  trait = list_Traits(folder)[1],
  title = trait,
  markers = NULL,
  labels = markers,
  models = c("MLM", "FarmCPU", "BLINK", "MLMM", "GLM", "CMLM", "SUPER"),
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

- markers:

  Markers to be labelled.

- labels:

  Labels to be used for markers.

- models:

  Models to read.

- skyline:

  Which skyline type to use. Can be "NYC" or "Kansas". If left NULL, it
  will use the highest P.value.

## Value

A volcano plot.
