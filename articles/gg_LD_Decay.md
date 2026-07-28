# LD Decay

``` r

library("gwaspr")
```

The function
[`gg_Volcano()`](https://derekmichaelwright.github.io/gwaspr/reference/gg_Volcano.md)
creates volcano plots from GAPIT GWAS results.

Specifying a `folder` and `trait` is all that is needed to create
manhattan plots.

``` r

# Load genotype file (note: header = T)
myG <- read.csv("gwaspr_myG_hmp.csv", header = T)
```

``` r

# Calulate LD decay
calc_LD_Decay(
  # Load genotype data
  xG = myG,
  # Specify a folder with GWAS results
  outputFolder = "LD_Decay/",
  # Select the number of markers per chromosome to analyse
  markerNum = 1000 )
```

``` r

# Plot
mp <- gg_LD_Decay(
  # Load genotype data
  xG = myG,
  # Specify a folder with GWAS results
  outputFolder = "LD_Decay/",
  # Select the number of markers per chromosome to analyse
  markerNum = 1000 )
# Save
ggsave("figures/gg_LD_Decay.png", mp, width = 12, height = 4 )
```

![](figures/gg_LD_Decay.png)
