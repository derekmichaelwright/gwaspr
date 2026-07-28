# checkTaxaNames

Check to see if the names in your myY and myG file match.

## Usage

``` r
checkTaxaNames(xG = myG, xY = myY)
```

## Arguments

- xG:

  GWAS genotype object. Note: needs to be in hapmap format and read in
  with header = F.

- xY:

  GWAS phenotype object.

## Value

A list of names not present in both datasets.
