# gg_Marker\_\*()

``` r

library(gwaspr)
```

The functions
[`gg_Marker_Box()`](https://derekmichaelwright.github.io/gwaspr/reference/gg_Marker_Box.md),
[`gg_Marker_Bar()`](https://derekmichaelwright.github.io/gwaspr/reference/gg_Marker_Bar.md),
and
[`gg_Marker_Pie()`](https://derekmichaelwright.github.io/gwaspr/reference/gg_Marker_Pie.md)
create marker plots with your `myY` and `myG` GWAS input files.

## Load data

``` r

# Load genotype file (note: header = T)
myG <- read.csv("gwaspr_myG_hmp.csv", header = T)
```

``` r

# Load phenotype file
myY <- read.csv("gwaspr_myY.csv")
# Convert our nominal trait from numeric to factor.
#myY <- myY %>% 
#  mutate(Cotyledon_Color = mv(Cotyledon_RedvsYellow, c(1, 0, NA), c("Red", "Yellow", "Green")),
#         Cotyledon_Color = factor(Cotyledon_Color, levels = c("Red", "Yellow", "Green")))
myY[1:15,]
```

    ##                 Name DTF_Sask_2017 DTF_Nepal_2017 Cotyledon_RedvsYellow
    ## 1    CDC_Asterix_AGL          54.7          128.0                     0
    ## 2      CDC_Rosie_AGL          59.0          123.3                     1
    ## 3       X3156.11_AGL          60.7          125.3                     1
    ## 4  CDC_Greenstar_AGL          56.7          121.0                     0
    ## 5     CDC_Cherie_AGL          54.3          125.3                     1
    ## 6     CDC_Glamis_AGL          59.0          123.0                     0
    ## 7       CDC_Gold_AGL          54.0          125.3                     0
    ## 8       CDC_Imax_AGL          53.0          128.7                     1
    ## 9    CDC_Impower_AGL          57.0          123.3                     0
    ## 10      CDC_KR.1_AGL          57.0          121.7                     1
    ## 11     CDC_LeMay_AGL          55.7          126.7                     0
    ## 12     CDC_Maxim_AGL          54.3          125.7                     1
    ## 13      CDC_QG.1_AGL          57.7          126.7                    NA
    ## 14 CDC_Red_Rider_AGL          57.0          123.5                     1
    ## 15   CDC_Redcoat_AGL          57.7          126.7                     1

------------------------------------------------------------------------

## gg_Marker_Box

### Single marker, single trait

Specifying a `trait` and a `marker` along with your genotype and
phenotype data as `xG` and `xY` is needed to create marker plots.

``` r

# Plot
mp <- gg_Marker_Box(
  # Genotype data
  xG = myG, 
  # Phenotype data
  xY = myY,
  # Select traits to plot
  traits = "DTF_Sask_2017",
  # Select markers to plot
  markers = "Lcu.1GRN.Chr6p3269280" )
# Save
ggsave("figures/gg_Marker_Box_01.png", 
       mp, width = 6, height = 4 )
```

![](figures/gg_Marker_Box_01.png)
