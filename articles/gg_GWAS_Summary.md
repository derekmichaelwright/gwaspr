# gg_GWAS_Summary()

``` r

library(gwaspr)
```

The function
[`gg_GWAS_Summary()`](https://derekmichaelwright.github.io/gwaspr/reference/gg_GWAS_Summary.md)
creates summary plots showing significant association from your GWAS
analyses.

Specifying a `folder` is all that is needed to create summary GWAS
plots.

``` r

# Plot
mp <- gg_GWAS_Summary(
  # Specify a folder with GWAS results
  folder = "GWAS_Results/" )
# Save
ggsave("figures/gg_GWAS_Summary_01.png",
       mp, width = 12, height = 4 )
```

![](figures/gg_GWAS_Summary_01.png)

------------------------------------------------------------------------

``` r

# Plot
mp <- gg_GWAS_Summary(
  # Specify a folder with GWAS results
  folder = "GWAS_Results/",
  # Specify a title
  title = "Summary of Significant GWAS Results",
  # Specify GWAS models to plot
  models = c("MLM","MLMM","FarmCPU","BLINK") )
# Save
ggsave("figures/gg_GWAS_Summary_02.png", 
       mp, width = 12, height = 4 )
```

![](figures/gg_GWAS_Summary_02.png)

------------------------------------------------------------------------

Make it an interactive plot with the following code

``` r

gg_plotly(mp, filename = "figures/gg_GWAS_Summary_02.html")
```

<https://derekmichaelwright.github.io/gwaspr/vignettes/figures/gg_GWAS_Summary_02.html>
