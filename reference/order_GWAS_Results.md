# order_GWAS_Results

Read in GWAS results, then arrange by P.value and write as csv..

## Usage

``` r
order_GWAS_Results(
  folder = "GWAS_Results/",
  files = list_Result_Files(folder),
  nrows = 10000
)
```

## Arguments

- folder:

  Folder containing GWAS results.

- files:

  The files to read. Default will read all files.

- nrows:

  Number of rows to read.

## Value

Rearranged files of GWAS results.
