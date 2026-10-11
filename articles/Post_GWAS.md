# Post GWAS Analyses

``` r

library(gwaspr)
```

The function
[`run_Summary()`](https://derekmichaelwright.github.io/gwaspr/reference/run_Summary.md)
can be used to check which models have been ran on each trait.

``` r

run_Summary(folder = "GWAS_Results")
```

The function `is_ran()` can be used to check if all if your triats in
your myY object have been run successfully.

``` r

is_Ran(folder = "GWAS_Results", myY = read.csv("gwaspr_myY.csv"))
```

------------------------------------------------------------------------

## Vignettes

### Manhattan Plots

- [`gg_Manhattan()`](https://derekmichaelwright.github.io/gwaspr/reference/gg_Manhattan.md)

- <https://derekmichaelwright.github.io/gwaspr/articles/gg_Manhattan.html>

- [`gg_Manhattan_Zoom()`](https://derekmichaelwright.github.io/gwaspr/reference/gg_Manhattan_Zoom.md)

- <https://derekmichaelwright.github.io/gwaspr/articles/gg_Manhattan_Zoom.html>

### Summary Plots & Tables

- [`table_GWAS_Results()`](https://derekmichaelwright.github.io/gwaspr/reference/table_GWAS_Results.md)
  &
  [`table_GWAS_Results_Summary()`](https://derekmichaelwright.github.io/gwaspr/reference/table_GWAS_Results_Summary.md)

- <https://derekmichaelwright.github.io/gwaspr/articles/table_GWAS_Results.html>

- [`gg_GWAS_Summary()`](https://derekmichaelwright.github.io/gwaspr/reference/gg_GWAS_Summary.md)

- <https://derekmichaelwright.github.io/gwaspr/articles/gg_GWAS_Summary.html>

- [`gg_QTL_Summary()`](https://derekmichaelwright.github.io/gwaspr/reference/gg_QTL_Summary.md)

- <https://derekmichaelwright.github.io/gwaspr/articles/gg_QTL_Summary.html>

### Marker plots

- [`gg_Marker_Box()`](https://derekmichaelwright.github.io/gwaspr/reference/gg_Marker_Box.md)

- <https://derekmichaelwright.github.io/gwaspr/articles/gg_Marker_Box.html>

- [`gg_Marker_Box()`](https://derekmichaelwright.github.io/gwaspr/reference/gg_Marker_Box.md)

- <https://derekmichaelwright.github.io/gwaspr/articles/gg_Marker_Box.html>

- [`gg_Marker_Box()`](https://derekmichaelwright.github.io/gwaspr/reference/gg_Marker_Box.md)

- <https://derekmichaelwright.github.io/gwaspr/articles/gg_Marker_Box.html>

- [`gg_Marker_Box()`](https://derekmichaelwright.github.io/gwaspr/reference/gg_Marker_Box.md)

- <https://derekmichaelwright.github.io/gwaspr/articles/gg_Marker_Box.html>

### Other GWAS related plots & functions

- [`gg_myG_Details()`](https://derekmichaelwright.github.io/gwaspr/reference/gg_myG_Details.md)
- <https://derekmichaelwright.github.io/gwaspr/articles/gg_myG_Details.html>
