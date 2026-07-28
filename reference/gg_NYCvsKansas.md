# gg_NYCvsKansas

Creates a manhattan plot comparing the NYC and Kansas results.

## Usage

``` r
gg_NYCvsKansas(
  folder = "GWAS_Results/",
  trait = list_Traits(folder)[1],
  title = trait,
  threshold = NULL,
  chr.colors = rep(c("darkgreen", "darkgoldenrod3"), 30),
  chr.unit = "100 Mbp"
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

- chr.colors:

  Colors for each chromosome.

- chr.unit:

  Unit for the x-axis. Can be one of c("kbp","100 kbp","Mbp","100
  Mbp","Gbp").

## Value

A manhattan plot.
