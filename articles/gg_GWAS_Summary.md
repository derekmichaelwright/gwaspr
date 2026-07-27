# gg_GWAS_Summary()

``` r

library("gwaspr")
```

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
