# LD Decay

The function
[`gg_Volcano()`](https://derekmichaelwright.github.io/gwaspr/reference/gg_Volcano.md)
creates volcano plots from GAPIT GWAS results.

Specifying a `folder` and `trait` is all that is needed to create
manhattan plots.

``` r

# Calulate LD decay
calc_LD_Decay(
  # Load genotype data
  xG = myG,
  # Specify a folder with GWAS results
  outputFolder = "LD_Decay/",
  # Select the number of markers per chromosome to analyse
  markerNum = 200 )
```
