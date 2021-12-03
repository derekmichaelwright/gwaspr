gwaspr R Package
================

<img src="hex_gwaspr.png" align="right" width = "200px" />

`gwaspr`: an `R` package for plotting GWAS results from `GAPIT` package

`GAPIT`: and `R` package for performing Genome Wide Association Studies
(GWAS)

<https://github.com/jiabowang/GAPIT3>

# Installation

``` r
devtools::install_github("derekmichaelwright/gwaspr")
```

``` r
library(gwaspr)
```

# Dependancies

`tidyverse`, `ggpubr`, `ggbeeswarm`, `ggrepel`, `ggtext`, `plotly`,
`htmlwidgets`

# GWAS Tutorial

<https://dblogr.com/academic/gwas_tutorial/gwas_tutorial.html>

# Usage

``` r
myTraits <- list_Traits(folder = "Results/")
```

``` r
myResults <- list_Result_Files(folder = "Results/")
```

``` r
mp <- gg_GWAS_Summary(folder = "Results/", traits = myTraits)
ggsave("GWAS_Summary.png", mp, width = 12, height = 6)
```

``` r
for(i in myTraits) {
  mp <- gg_Manhattan(folder = "Results/", trait = myTraits[i], facet = F)
  ggsave(paste0("Multi_",i,".png"), mp, width = 12, height = 6)
}
```

``` r
for(i in myTraits) {
  mp <- gg_Manhattan(folder = "Results/", trait = myTraits[i], facet = T)
  ggsave(paste0("Facet_",i,".png"), mp, width = 12, height = 12)
}
```

------------------------------------------------------------------------

© Derek Michael Wright [www.dblogr.com/](https://dblogr.com/)
