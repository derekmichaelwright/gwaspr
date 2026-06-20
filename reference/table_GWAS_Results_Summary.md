# table_GWAS_Results_Summary

Create a summary using the output from \`table_GWAS_Results()\`.

## Usage

``` r
table_GWAS_Results_Summary(xx, binMarkers = F, binSize = 5e+06, onlySig = F)
```

## Arguments

- xx:

  Object from \`table_GWAS_Results()\`.

- binMarkers:

  Logical, if TRUE, markers will be bined based on \`binSize\`.

- binSize:

  range on each side of marker to bin. default = 1,000,000.

- onlySig:

  Logical, if TRUE, any suggested associations will be removed.

## Value

A table of significant GWAS results.
