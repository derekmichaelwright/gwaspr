# gwaspr R Package

`gwaspr`: an `R` package for plotting GWAS results from the `GAPIT`
package

![](reference/figures/logo_gwaspr.png)

------------------------------------------------------------------------

# GAPIT

`GAPIT`: and `R` package for performing Genome Wide Association Studies
(GWAS)

> - **GAPIT Website**:
>   [https://www.maizegenetics.net/gapit](https://www.maizegenetics.net/gapit)
> - **GAPIT github**:
>   [https://github.com/jiabowang/GAPIT](https://github.com/jiabowang/GAPIT)

------------------------------------------------------------------------

# Instal gwaspr

``` r

devtools::install_github("derekmichaelwright/gwaspr")
```

``` r

library(gwaspr)
```

------------------------------------------------------------------------

# Dependancies

`tidyverse`, `ggpubr`, `ggbeeswarm`, `ggrepel`, `ggtext`, `plotly`,
`htmlwidgets`

------------------------------------------------------------------------

# list_Traits()

List off the traits which have GWAS results files in the designated
folder.

``` r

data_path <- here::here()
setwd(data_path)
getwd()
list_Traits(folder = "vignettes/GWAS_Results/")
list_Traits(folder = data_path)
list_Traits(folder = "GWAS_Results/")
```

------------------------------------------------------------------------

# is_Ran()

Check to see which of your `myY` traits have GWAS results files int he
designated folder.

``` r

#myY <- read.csv("vignettes/myY.csv")
#is_Ran(folder = "GWAS_Results/", myY = myY)
```

------------------------------------------------------------------------

# run_Summary()

Checks which GWAS models have been run for each trait within the
designated folder.

``` r

run_Summary(folder = "vignettes/GWAS_Results/")
```

------------------------------------------------------------------------

# Order_GWAS_Results()

Reorders the result files if they are not already arranged by P.value.

``` r

order_GWAS_Results(folder = "vignettes/GWAS_Results/")
```

------------------------------------------------------------------------

# is_Ordered()

``` r

is_Ordered(folder = "vignettes/GWAS_Results/")
```

------------------------------------------------------------------------

# list_Top_Markers()

``` r

list_Top_Markers(folder = "GWAS_Results/", trait = "DTF_Nepal_2017", chroms = c(2,5))
```

``` r

list_Top_Markers(folder = "GWAS_Results/", trait = "DTF_Sask_2017", chroms = 6)
```

``` r

list_Top_Markers(folder = "GWAS_Results/", trait = "DTF_Sask_2017_b", chroms = 6)
```

``` r

list_Top_Markers(folder = "GWAS_Results/", trait = "Cotyledon_Color_RvsY", chroms = 1)
```

------------------------------------------------------------------------

If you want to remove some of the unnessesary files in the gwas results
folder, use the following function to delete all files except the ones
with GWAS result data.

``` r

clean_GWAS_Results(folder = "GWAS_Results/")
```

------------------------------------------------------------------------
