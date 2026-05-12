# table_GWAS_Results_Summary

Create a summary using the output from \`table_GWAS_Results()\`.

## Usage

``` r
table_GWAS_Results_Summary(xx, onlySig = F, binMarkers = F, binSize = 1e+06)
```

## Arguments

- xx:

  Object from \`table_GWAS_Results()\`.

- onlySig:

  Logical, if TRUE, any suggested associations will be removed.

- binMarkers:

  Logical, if TRUE, markers will be bined based on \`binSize\`.

- binSize:

  range on each side of marker to bin. default = 1,000,000.

## Value

A table of significant GWAS results.
