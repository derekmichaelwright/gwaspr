# list_Top_Markers

Finds the markers with the highest association on each chromosome.

## Usage

``` r
list_Top_Markers(
  folder = "GWAS_Results/",
  traits = list_Traits(folder),
  models = c("MLM", "FarmCPU", "BLINK", "MLMM", "GLM", "CMLM", "SUPER"),
  threshold = 5,
  chroms = 1:50,
  n = 3
)
```

## Arguments

- folder:

  Folder containing GWAS results.

- traits:

  GWAS trait.

- models:

  GWAS models to include.

- threshold:

  filters results with -log10(p) below threshold.

- chroms:

  Chromosomes to include.

- n:

  Number per trait and model.

## Value

Table of top results.
