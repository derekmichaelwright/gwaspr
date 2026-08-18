# calc_LD_Decay

Calculates pairwise LD for markers within your myG genotype file.

## Usage

``` r
calc_LD_Decay(xG = myG, outputFolder, markerNum = 200)
```

## Arguments

- xG:

  GWAS genotype object. Note: needs to be in hapmap format.

- outputFolder:

  Folder to place RData files into.

- markerNum:

  Number of markers to randomly select for each chromosome.

## Value

RData files in the outputFolder location of LD data.
