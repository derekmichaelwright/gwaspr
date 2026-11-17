# table_GWAS_Results

Create a table of significant GWAS results.

## Usage

``` r
table_GWAS_Results(
  folder = "GWAS_Results/",
  traits = list_Traits(folder),
  threshold = 6,
  sug.threshold = NULL,
  nrows = 10000,
  skyline = "NYC",
  models = c("MLM", "MLMM", "FarmCPU", "BLINK", "GLM", "CMLM", "SUPER")
)
```

## Arguments

- folder:

  Folder containing GWAS results.

- threshold:

  Significant threshold.

- sug.threshold:

  Suggestive threshold.

- nrows:

  Number of rows to read.

- skyline:

  Which skyline type to use. Can be "NYC" or "Kansas". If NULL, it will
  use the highest P.value.

- models:

  GWAS models to use.

## Value

A table of significant GWAS results.
