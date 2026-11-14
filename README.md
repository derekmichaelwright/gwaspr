
# gwaspr <img src="man/figures/logo_gwaspr.png" align="right" width="90" height="100" />

`gwaspr`: an `R` package for plotting GWAS results from the `GAPIT`
package

------------------------------------------------------------------------

# Instal `gwaspr`

``` r
devtools::install_github("derekmichaelwright/gwaspr")
```

# Load `gwaspr`

``` r
library(gwaspr)
```

------------------------------------------------------------------------

# GWAS with `GAPIT`

`GAPIT`: an `R` package for performing Genome Wide Association Studies
(GWAS)

> - **GAPIT Website**: <a href="https://www.maizegenetics.net/gapit"
>   target="_blank">https://www.maizegenetics.net/gapit</a>
> - **GAPIT github**: <a href="https://github.com/jiabowang/GAPIT"
>   target="_blank">https://github.com/jiabowang/GAPIT</a>

``` r
devtools::install_github("jiabowang/GAPIT")
```

``` r
library(GAPIT)
```

------------------------------------------------------------------------

# `gwaspr` Github

- **github link**:
  <a href="https://github.com/derekmichaelwright/gwaspr/"
  target="_blank">https://github.com/derekmichaelwright/gwaspr/</a>

------------------------------------------------------------------------

# `gwapsr` Tutorial

- **dblogr GWAS Tutorial**: <a
  href="https://derekmichaelwright.github.io/dblogr/academic/gwas_tutorial/"
  target="_blank">https://derekmichaelwright.github.io/dblogr/academic/gwas_tutorial/</a>
- [Phenotype data (gwaspr_myY.csv)](vignettes/gwaspr_myY.csv)
- [Covariate data (gwaspr_myCV.csv)](vignettes/gwaspr_myCV.csv)
- [Genotype data (gwaspr_myG_hmp.zip)](vignettes/gwaspr_myG_hmp.zip)

------------------------------------------------------------------------

# Package Dependancies

`tidyverse`, `ggpubr`, `ggbeeswarm`, `ggrepel`, `ggtext`, `plotly`,
`htmlwidgets`

------------------------------------------------------------------------

# Example plots

## gg_Manhattan()

![](vignettes/figures/gg_Manhattan_02_2.png)

![](vignettes/figures/gg_Manhattan_03_4.png)

## gg_GWAS_Summary()

![](vignettes/figures/gg_GWAS_Summary_04.png)

## gg_Manhattan_Zoom

![](vignettes/figures/gg_Manhattan_Zoom_01_3.png)

## gg_Marker_Box()

![](vignettes/figures/gg_Marker_Box_03.png)

## gg_Marker_Pie()

![](vignettes/figures/gg_Marker_Pie_02.png)

## gg_Marker_Bubble()

![](vignettes/figures/gg_Marker_Bubble_03.png)

## gg_myG_Details()

![](vignettes/figures/gg_myG_Details_03_myG_Details.png)

## gg_LD_Heatmap()

![](vignettes/figures/gg_LD_decay_02.png)

![](vignettes/figures/gg_LD_Heatmap_01.png)

## gg_QTL_Summary()

![](vignettes/figures/gg_QTL_Summary_02.png)

------------------------------------------------------------------------

![](man/figures/logo_gwaspr.png)
