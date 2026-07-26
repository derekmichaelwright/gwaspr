# gg_LDheatmap()

``` r

library("gwaspr")
```

``` r

# Load genotype file (note: header = T)
myG <- read.csv("myG_hmp.csv", header = T)
# Plot
mp <- gg_LDheatmap(
  xG = myG, 
  chr = 6, 
  pos1 = 2500000, 
  pos2 = 4000000,
  myMs = "Lcu.1GRN.Chr6p3269280", 
  axisTextSize = 3, 
  nameTrim = "Lcu.1GRN.Chr6",
  myTitle = "Lentil Diversity Panel"
  )
# Save
ggsave("figures/gg_LDheatmap_01.png", 
       mp, width = 12, height = 12, dpi = 600)
```

![](figures/gg_LDheatmap_01.png)
