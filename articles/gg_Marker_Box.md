# gg_Marker_Box()

``` r

library(gwaspr)
```

The function
[`gg_Marker_Box()`](https://derekmichaelwright.github.io/gwaspr/reference/gg_Marker_Box.md)
creates marker plots with your `myY` and `myG` GWAS input files.

## Load data

> Note: `header = T`

``` r

# Load genotype file (note: header = T)
myG <- read.csv("gwaspr_myG_hmp.csv", header = T)
```

``` r

# Load phenotype file
myY <- read.csv("gwaspr_myY.csv")
# Convert our nominal trait from numeric to factor.
myY <- myY %>% 
  mutate(Cotyledon_Color = mv(Cotyledon_RedvsYellow, c(1, 0, NA), c("Red", "Yellow", "Green")),
         Cotyledon_Color = factor(Cotyledon_Color, levels = c("Red", "Yellow", "Green")))
```

``` r

myY[1:10,]
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
    ##    Cotyledon_Color
    ## 1           Yellow
    ## 2              Red
    ## 3              Red
    ## 4           Yellow
    ## 5              Red
    ## 6           Yellow
    ## 7           Yellow
    ## 8              Red
    ## 9           Yellow
    ## 10             Red

------------------------------------------------------------------------

## Single marker, single trait

Specifying `traits` and `markers` along with your genotype and phenotype
data as `xG` and `xY` is needed to create marker plots.

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

------------------------------------------------------------------------

## All Options

[`gg_Marker_Bar()`](https://derekmichaelwright.github.io/gwaspr/reference/gg_Marker_Bar.md)
contains a number of different options for customizing the marker plots.

``` r

# Prep data
myM <- c("Lcu.1GRN.Chr6p3269280")
# Plot
mp <- gg_Marker_Box(
  xG = myG, 
  xY = myY,
  traits = "DTF_Sask_2017",
  markers = myM,
  marker.colors = gwaspr_Colors,
  remove.hets = F,
  remove.N = T,
  plot.violin = T,
  plot.box = T,
  plot.points = T,
  box.width = 0.1,
  point.size = 1,
  point.beeswarm = F,
  ncol = NULL,
  legend.rows = 1,
  title = NULL,
  subtitle = paste(myM, collapse = "\n"),
  cv.source = "xG",
  cv.name = NULL,
  cv.colors = NULL,
  cv.label = cv.name,
  groupByCV = F )
# Save
ggsave("figures/gg_Marker_Box_02.png", 
       mp, width = 6, height = 4 )
```

![](figures/gg_Marker_Box_02.png)

------------------------------------------------------------------------

## Multiple markers, multiple traits

``` r

# Plot
mp <- gg_Marker_Box(
  # Genotype data
  xG = myG,
  # Phenotype data
  xY = myY,
  # Select traits to plot
  traits = c("DTF_Nepal_2017", "DTF_Sask_2017"),
  # Select markers to plot
  markers = c("Lcu.1GRN.Chr5p1658484", "Lcu.1GRN.Chr2p44545877"),
  # Remove heterozygotic markers  
  remove.hets = T )
# Save
ggsave("figures/gg_Marker_Box_03.png",
       mp, width = 8, height = 4 )
```

![](figures/gg_Marker_Box_03.png)

------------------------------------------------------------------------

## Boxplots

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
  markers = "Lcu.1GRN.Chr6p3269280",
  # Select marker colors
  marker.colors = c("darkorange3", "steelblue", "darkblue"),
  # Choose what should be plotted
  plot.violin = F,
  plot.box = T,
  plot.points = F,
  # Change the width of the boxplots
  box.width = 0.3 )
# Save
ggsave("figures/gg_Marker_Box_04.png", 
       mp, width = 6, height = 4 )
```

![](figures/gg_Marker_Box_04.png)

------------------------------------------------------------------------

## Violin + points

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
  markers = "Lcu.1GRN.Chr6p3269280",
  # Select marker colors
  marker.colors = c("darkorange3", "steelblue", "darkblue"),
  # Choose what should be plotted
  plot.violin = T,
  plot.box = F,
  plot.points = T,
  # Set the point size
  point.size = 2)
# Save
ggsave("figures/gg_Marker_Box_05.png", 
       mp, width = 6, height = 4 )
```

![](figures/gg_Marker_Box_05.png)

------------------------------------------------------------------------

## myG Covariable

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
  markers = "Lcu.1GRN.Chr6p3269280",
  # Select marker colors
  marker.colors = c("burlywood4", "maroon3", "purple3"),
  # Keep heterozygous markers
  remove.hets = F,
  # Choose what should be plotted
  plot.violin = T,
  plot.box = F,
  plot.points = T,
  # Set the point size
  point.size = 0.75,
  # Plot with geom_beeswarm instead of quasirandom
  point.beeswarm = T,
  # Select Covariable trait for points
  cv.name = "Lcu.1GRN.Chr2p44545877", 
  # Select colors for the covariable
  cv.colors = c("darkorange3", "steelblue", "darkblue") )
# Save
ggsave("figures/gg_Marker_Box_06.png", 
       mp, width = 6, height = 4 )
```

![](figures/gg_Marker_Box_06.png)

------------------------------------------------------------------------

## myY Covariable

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
  markers = "Lcu.1GRN.Chr6p3269280",
  # Select marker colors
  marker.colors = c("darkorange3", "steelblue", "darkblue"),
  # Choose what should be plotted
  plot.violin = T,
  plot.box = F,
  plot.points = T,
  # Set the point size
  point.size = 1.5,
  # Plot with geom_beeswarm instead of quasirandom
  point.beeswarm = T,
  # Select covariable source
  cv.source = "xY",
  # Select Covariable trait for points
  cv.name = "Cotyledon_Color", 
  # Select colors for the covariable
  cv.colors = c("darkred", "darkgoldenrod2", "darkgreen"),
  # Set a custom label for the covariable legend
  cv.label = "Cotyledon Color" )
# Save
ggsave("figures/gg_Marker_Box_07.png", 
       mp, width = 6, height = 4 )
```

![](figures/gg_Marker_Box_07.png)

------------------------------------------------------------------------

## Grouped by myY Covariable

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
  markers = "Lcu.1GRN.Chr6p3269280",
  # Select marker colors
  marker.colors = c("darkorange3", "steelblue", "darkblue"),
  # Choose what should be plotted
  plot.violin = T,
  plot.box = F,
  plot.points = T,
  # Set the point size
  point.size = 1.5,
  # Plot with geom_beeswarm instead of quasirandom
  point.beeswarm = T,
  # Select covariable source
  cv.source = "xY",
  # Select Covariable trait for points
  cv.name = "Cotyledon_Color", 
  # Select colors for the covariable
  cv.colors = c("darkred", "darkgoldenrod2", "darkgreen"),
  # Set a custom label for the covariable legend
  cv.label = "Cotyledon Color",
  # Plot CV on x-axis instead of marker
  groupByCV = T)
# Save
ggsave("figures/gg_Marker_Box_08.png", 
       mp, width = 6, height = 4 )
```

![](figures/gg_Marker_Box_08.png)
