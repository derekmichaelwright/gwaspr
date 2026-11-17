# gwaspr_batch_Man

[Creates Manhattan plots for all results in your specified
folder.](https://derekmichaelwright.github.io/gwaspr/articles/gwaspr_batch_Man.html)

## Usage

``` r
gwaspr_batch_Man(
  folder = "vignettes/GWAS_Results/",
  output = "vignettes/figures/",
  traits = list_Traits(folder),
  models = c("MLM", "MLMM", "FarmCPU", "BLINK", "GLM", "CMLM", "SUPER")
)
```

## Arguments

- folder:

  Folder containing GWAS results.

- output:

  Output folder for Manhattan plots.

- traits:

  The traits to read.

- models:

  Models to read.

## Value

Manhattan plots.
