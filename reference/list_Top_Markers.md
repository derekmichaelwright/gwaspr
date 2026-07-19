# list_Top_Markers

Finds the markers with the highest association on each chromosome.

## Usage

``` r
list_Top_Markers(
  folder = "GWAS_Results/",
  traits = list_Traits(folder),
  chr = 1:50,
  models = c("MLM", "FarmCPU", "BLINK", "MLMM", "GLM", "CMLM", "SUPER"),
  threshold = 5,
  n = 5
)
```

## Arguments

- folder:

  Folder containing GWAS results.

- traits:

  GWAS trait.

- chr:

  Chromosomes to include.

- models:

  GWAS models to include.

- threshold:

  filters results with -log10(p) below threshold.

- n:

  Number per trait and model.

## Value

Table of top results.
