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
this case, they are located in a folder called `GWAS_Results/`. For this
example we will plot GWAS results from two traits in a lentil diversity
panel:

-   **Testa\_Pattern**: a *qualitative* trait describing the presence or
    absence of seed coat pigmentation.
-   **DTF\_Nepal\_2017**: a *quantitative* trait describing days from
    sowing to flowering in a 2017 Nepal field trial.

Note: for more info check out this [GWAS
tutorial](https://dblogr.com/academic/gwas_tutorial/gwas_tutorial.html%22).

``` r
myTraits <- list_Traits(folder = "GWAS_Results/")
myTraits
```

    ## [1] "DTF_Nepal_2017" "Testa_Pattern"

``` r
myResults <- list_Result_Files(folder = "GWAS_Results/")
myResults
```

    ## [1] "GAPIT.Blink.DTF_Nepal_2017.GWAS.Results.csv"  
    ## [2] "GAPIT.Blink.Testa_Pattern.GWAS.Results.csv"   
    ## [3] "GAPIT.FarmCPU.DTF_Nepal_2017.GWAS.Results.csv"
    ## [4] "GAPIT.FarmCPU.Testa_Pattern.GWAS.Results.csv" 
    ## [5] "GAPIT.MLM.DTF_Nepal_2017.GWAS.Results.csv"    
    ## [6] "GAPIT.MLM.Testa_Pattern.GWAS.Results.csv"     
    ## [7] "GAPIT.MLMM.DTF_Nepal_2017.GWAS.Results.csv"   
    ## [8] "GAPIT.MLMM.Testa_Pattern.GWAS.Results.csv"

``` r
list_Top_Markers(trait = "DTF_Nepal_2017", model = "MLMM", 
                 folder = "GWAS_Results/", 
                 threshold = 6.7, chroms = c(2,5), n = 1)
```

    ##                      SNP CHR      POS -log10(p)
    ## 1 Lcu.2RBY.Chr2p42543877   2 42543877     11.58
    ## 2  Lcu.2RBY.Chr5p1069654   5  1069654     16.71

## Summary Plot

``` r
mp <- gg_GWAS_Summary(folder = "GWAS_Results/", traits = myTraits,
                      threshold = 7.3, threshold2 = 6.7)
ggsave("man/figures/GWAS_Summary.png", mp, width = 10, height = 2)
```

![](man/figures/GWAS_Summary.png)

## Manhattan Plots

``` r
for(i in myTraits) {
  mp <- gg_Manhattan(folder = "GWAS_Results/", trait = i, facet = F,
                     threshold = 7.3, threshold2 = 6.7)
  ggsave(paste0("man/figures/Multi_",i,".png"), mp, width = 10, height = 4)
}
```

![](man/figures/Multi_DTF_Nepal_2017.png)

![](man/figures/Multi_Testa_Pattern.png)

``` r
for(i in myTraits) {
  mp <- gg_Manhattan(folder = "GWAS_Results/", trait = i, facet = T,
                     threshold = 7.3, threshold2 = 6.7)
  ggsave(paste0("man/figures/Facet_",i,".png"), mp, width = 10, height = 8)
}
```

![](man/figures/Facet_DTF_Nepal_2017.png)

![](man/figures/Facet_Testa_Pattern.png)

# GAPIT

`GAPIT`: and `R` package for performing Genome Wide Association Studies
(GWAS)

<https://github.com/jiabowang/GAPIT3>

# Dependancies

`tidyverse`, `ggpubr`, `ggbeeswarm`, `ggrepel`, `ggtext`, `plotly`,
`htmlwidgets`

------------------------------------------------------------------------

© Derek Michael Wright [www.dblogr.com/](https://dblogr.com/)
