# table_GWAS_Results

Create a table of significant GWAS results.

## Usage

``` r
table_GWAS_Region(
  folder = "GWAS_Results/",
  chr = 1,
  pos1 = 3.5e+08,
  pos2 = 3.7e+08,
  fnames = list_Result_Files(folder),
  threshold = 6,
  sug.threshold = NULL,
  nrowstoread = 1000,
  skyline = NULL,
  models = c("MLM", "MLMM", "FarmCPU", "BLINK", "GLM", "CMLM", "SUPER")
)
```

## Arguments

- folder:

  Folder containing GWAS results.

- chr:

  Chromosome to plot.

- pos1:

  Start position on chromosome.

- pos2:

  End position on chromosome.

- fnames:

  The files to read.

- threshold:

  Significant threshold.

- sug.threshold:

  Suggestive threshold.

- nrowstoread:

  Number of rows to read.

- skyline:

  Which skyline type to use. Can be "NYC" or "Kansas". If left NULL, it
  will use the highest P.value.

- models:

  GWAS models to use.

- useHBPvalues:

  Logical, if TRUE, H.B.P.Values will be uses.

## Value

A table of significant GWAS results.
