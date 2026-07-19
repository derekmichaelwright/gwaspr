# clean_GWAS_Results

Delete extra GWAS files in the folder containing GWAS result files from
GAPIT.

## Usage

``` r
clean_GWAS_Results(
  folder = "GWAS_Results/",
  remove_nonResults = T,
  remove_Kansas = F,
  remove_NYC = F
)
```

## Arguments

- folder:

  Folder containing GWAS results.

- remove_nonResults:

  Delelets all non .csv files (all .pdf files)

- remove_Kansas:

  Delete Kansas files.

- remove_NYC:

  Delete NYC files.

## Value

A table of which models have been run on each trait.
