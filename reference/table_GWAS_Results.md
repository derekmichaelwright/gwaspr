# table_GWAS_Results

Create a table of significant GWAS results.

## Usage

``` r
table_GWAS_Results(
  folder = "GWAS_Results/",
  fnames = list_Result_Files(folder),
  nrowstoread = 1000,
  threshold = 6,
  sug.threshold = NULL,
  skyline = NULL
)
```

## Arguments

- folder:

  Folder containing GWAS results.

- fnames:

  The files to read.

- nrowstoread:

  Number of rows to read.

- threshold:

  Significant threshold.

- sug.threshold:

  Suggestive threshold.

- skyline:

  Which skyline type to use. Can be "NYC" or "Kansas". If left NULL, it
  will use the highest P.value.

- useHBPvalues:

  Logical, if TRUE, H.B.P.Values will be uses.

## Value

A table of significant GWAS results.
