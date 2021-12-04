gwaspr R Package
================

`gwaspr`: an `R` package for plotting GWAS results from the `GAPIT`
package

# Installation

``` r
devtools::install_github("derekmichaelwright/gwaspr")
```

``` r
library(gwaspr)
```

![](man/figures/hex_gwaspr.png)

# GWAS Tutorial

<a href="https://dblogr.com/academic/gwas_tutorial/gwas_tutorial.html">
<button class="btn btn-success"><i class="far fa-file-code"></i> https://dblogr.com/academic/gwas_tutorial/gwas_tutorial.html</button>
</a>

# Usage

For best practice, output from GAPIT should be in its own folder. In
this case, they are located in a folder called `Results/`.

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

# GAPIT

`GAPIT`: and `R` package for performing Genome Wide Association Studies
(GWAS)

<https://github.com/jiabowang/GAPIT3>

# Dependancies

`tidyverse`, `ggpubr`, `ggbeeswarm`, `ggrepel`, `ggtext`, `plotly`,
`htmlwidgets`

------------------------------------------------------------------------

© Derek Michael Wright [www.dblogr.com/](https://dblogr.com/)
