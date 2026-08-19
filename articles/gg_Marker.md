# gg_Marker\_\*()

``` r

library(gwaspr)
```

The functions
[`gg_Marker_Box()`](https://derekmichaelwright.github.io/gwaspr/reference/gg_Marker_Box.md),
[`gg_Marker_Bar()`](https://derekmichaelwright.github.io/gwaspr/reference/gg_Marker_Bar.md),
and
[`gg_Marker_Pie()`](https://derekmichaelwright.github.io/gwaspr/reference/gg_Marker_Pie.md)
create marker plots with your `myY` and `myG` GWAS input files.

## Load data

``` r

# Load genotype file (note: header = T)
myG <- read.csv("gwaspr_myG_hmp.csv", header = T)
```

``` r

myG <- read.csv("gwaspr_myG_hmp_10.csv", header = F)
#write.csv(myG, "myG_hmp_10.csv", row.names = F)
```

``` r

# Load phenotype file
myY <- read.csv("gwaspr_myY.csv")
# Convert our nominal trait from numeric to factor.
myY <- myY %>% 
  mutate(Cotyledon_Color = mv(Cotyledon_RedvsYellow, c(1, 0, NA), c("Red", "Yellow", "Green")),
         Cotyledon_Color = factor(Cotyledon_Color, levels = c("Red", "Yellow", "Green")))
```

------------------------------------------------------------------------

## gg_Marker_Box

### Single marker, single trait

Specifying a `trait` and a `marker` along with your genotype and
phenotype data as `xG` and `xY` is needed to create marker plots.

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

### Customized Plots

#### Boxplots

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
  marker.colors = c("darkorange3", "steelblue"),
  # Choose what should be plotted
  plot.violin = F,
  plot.box = T,
  plot.points = F,
  # Change the width of the boxplots
  box.width = 0.3,
  # Create a custom label for the y-axis
  yLab = "Days from sowing to flower" )
# Save
ggsave("figures/gg_Marker_Box_02.png", 
       mp, width = 6, height = 4 )
```

![](figures/gg_Marker_Box_02.png)

------------------------------------------------------------------------

#### Violin + points

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
  marker.colors = c("darkorange3", "steelblue"),
  # Choose what should be plotted
  plot.violin = T,
  plot.box = F,
  plot.points = T,
  # Set the point size
  point.size = 2,
  # Create a custom label for the y-axis
  yLab = "Days from sowing to flower" )
# Save
ggsave("figures/gg_Marker_Box_03.png", 
       mp, width = 6, height = 4 )
```

![](figures/gg_Marker_Box_03.png)

------------------------------------------------------------------------

#### Covariable Point Color

##### myG covariable

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
  marker.colors = c("darkorange3", "darkseagreen4", "steelblue", "burlywood4"),
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
  cv.colors = c("darkslategray4", "maroon3", "purple3"),
  # Set a custom label for the covariable legend
  cv.label = "Chr2p44545877" )
# Save
ggsave("figures/gg_Marker_Box_04.png", 
       mp, width = 6, height = 4 )
```

![](figures/gg_Marker_Box_04.png)

------------------------------------------------------------------------

##### myY covariable

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
  marker.colors = c("darkorange3", "steelblue"),
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
ggsave("figures/gg_Marker_Box_05.png", 
       mp, width = 6, height = 4 )
```

![](figures/gg_Marker_Box_05.png)

------------------------------------------------------------------------

### Multiple markers, multiple traits

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
  markers = c("Lcu.1GRN.Chr5p1658484", "Lcu.1GRN.Chr2p44545877") )
# Save
ggsave("figures/gg_Marker_Box_06.png",
       mp, width = 8, height = 4 )
```

![](figures/gg_Marker_Box_06.png)

------------------------------------------------------------------------

## gg_Marker_Bar

### Single marker, single trait

``` r

# Plot
mp <- gg_Marker_Bar(
  # Genotype data
  xG = myG, 
  # Phenotype data
  xY = myY,
  # Select traits to plot
  traits = "DTF_Sask_2017",
  # Select markers to plot
  markers = "Lcu.1GRN.Chr6p3269280" )
# Save
ggsave("figures/gg_Marker_Bar_01.png", 
       mp, width = 6, height = 4 )
```

![](figures/gg_Marker_Bar_01.png)

------------------------------------------------------------------------

### Customized Plots

#### Histogram

``` r

# Plot
mp <- gg_Marker_Bar(
  # Genotype data
  xG = myG, 
  # Phenotype data
  xY = myY,
  # Select traits to plot
  traits = "DTF_Sask_2017",
  # Select markers to plot
  markers = "Lcu.1GRN.Chr6p3269280",
  # Select marker colors
  marker.colors = c("darkorange3", "steelblue"),
  # Choose what should be plotted
  plot.histogram = T,
  plot.density = F )
# Save
ggsave("figures/gg_Marker_Bar_02.png", 
       mp, width = 6, height = 4 )
```

![](figures/gg_Marker_Bar_02.png)

------------------------------------------------------------------------

#### Density

``` r

# Plot
mp <- gg_Marker_Bar(
  # Genotype data
  xG = myG, 
  # Phenotype data
  xY = myY,
  # Select traits to plot
  traits = "DTF_Sask_2017",
  # Select markers to plot
  markers = "Lcu.1GRN.Chr6p3269280",
  # Select marker colors
  marker.colors = c("darkorange3", "steelblue"),
  # Choose what should be plotted
  plot.histogram = F,
  plot.density = T )
# Save
ggsave("figures/gg_Marker_Bar_03.png", 
       mp, width = 6, height = 4)
```

![](figures/gg_Marker_Bar_03.png)

------------------------------------------------------------------------

### Multiple markers, multiple traits

``` r

# Plot
mp <- gg_Marker_Bar(
  # Genotype data
  xG = myG, 
  # Phenotype data
  xY = myY,
  # Select traits to plot
  traits = c("DTF_Nepal_2017", "DTF_Sask_2017"),
  # Select markers to plot
  markers = c("Lcu.1GRN.Chr2p44545877", "Lcu.1GRN.Chr5p1658484") )
# Save
ggsave("figures/gg_Marker_Bar_04.png", 
       mp, width = 8, height = 4 )
```

![](figures/gg_Marker_Bar_04.png)

------------------------------------------------------------------------

### Factor data

#### Numeric format

``` r

# Plot 
mp <- gg_Marker_Bar(
  # Genotype data
  xG = myG, 
  # Phenotype data
  xY = myY,
  # Select traits to plot
  traits = "Cotyledon_RedvsYellow",
  # Select markers to plot
  markers = "Lcu.1GRN.Chr1p365986872",
  # Select marker colors
  marker.colors = c("darkred", "darkgoldenrod2"),
  # Choose what should be plotted
  plot.histogram = T,
  plot.density = F )
# Save
ggsave("figures/gg_Marker_Bar_05.png", 
       mp, width = 6, height = 4 )
```

![](figures/gg_Marker_Bar_05.png)

------------------------------------------------------------------------

#### Factor format

``` r

# Plot 
mp <- gg_Marker_Bar(
  # Genotype data
  xG = myG, 
  # Phenotype data
  xY = myY,
  # Select traits to plot
  traits = "Cotyledon_Color",
  # Select markers to plot
  markers = "Lcu.1GRN.Chr1p365986872",
  # Select marker colors
  marker.colors = c("darkred", "darkgoldenrod2", "darkgreen"),
  # Choose what should be plotted
  plot.histogram = T,
  plot.density = F )
# Save
ggsave("figures/gg_Marker_Bar_06.png", 
       mp, width = 6, height = 4 )
```

![](figures/gg_Marker_Bar_06.png)

------------------------------------------------------------------------

## gg_Marker_Pie

``` r

# Plot 
mp <- gg_Marker_Pie(
  # Genotype data
  xG = myG, 
  # Phenotype data
  xY = myY,
  # Select traits to plot
  trait = "Cotyledon_Color",
  # Select markers to plot
  markers = "Lcu.1GRN.Chr1p365986872",
  # Select marker colors
  marker.colors = c("darkred", "darkgoldenrod2", "darkgreen") )
# Save
ggsave("figures/gg_Marker_Pie_01.png", 
       mp, width = 6, height = 4 )
```

![](figures/gg_Marker_Pie_01.png)
