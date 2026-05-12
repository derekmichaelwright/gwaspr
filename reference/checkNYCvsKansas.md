# checkNYCvsKansas

Check if NYC and Kansas results are identical, if so, the Kansas file
will be deleted.

## Usage

``` r
checkNYCvsKansas(folder = "GWAS_Results/", deleteKansas = F)
```

## Arguments

- folder:

  Folder containing GWAS results.

- deteleKansas:

  Logical, if TRUE, will delete any \`Kansas\` files with no difference
  between the \`NYC\` files.

## Value

A table of which runs were identical and deleted.
