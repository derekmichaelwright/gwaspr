# gg_LDheatmap()

``` r

library("gwaspr")
```

``` r

# Load genotype file (note: header = T)
myG <- read.csv("myG_hmp.csv", header = T)
# Plot
mp <- gg_LDheatmap(
  # Load genotype data
  xG = myG, 
  # Select chromosome and region to analyse
  chr = 6, 
  pos1 = 2500000, 
  pos2 = 4000000,
  # Select a LD metric
  metric = "R^2",
  # Set the threshold for identifying linked markers
  threshold = 0.9,
  # Select markers to highlight
  myMs = "Lcu.1GRN.Chr6p3269280",
  # Set the size of the marker names on the x and y axes
  axisTextSize = 3, 
  # prefix to trim from the marker name
  nameTrim = "Lcu.1GRN.Chr6",
  # Set custom title for plot
  myTitle = "Lentil Diversity Panel",
  # Select colors for the gradient scale
  color.low = "white",
  color.mid = "goldenrod1",
  color.high = "darkred" )
# Save Plot
ggsave("figures/gg_LDheatmap_01.png", 
       mp[[1]], width = 12, height = 12, dpi = 600 )
# Save data
write.csv(mp[[2]], "figures/LDheatmap_data.csv", row.names = F)
```

![](figures/gg_LDheatmap_01.png)

``` r

read.csv("figures/LDheatmap_data.csv")[1:30,]
```

    ##        SNP1     SNP2          R.2 Chr  SNP1_d  SNP2_d
    ## 1  p2525330 p2525330 0.0000000000   6 2525330 2525330
    ## 2  p2525330 p2556449 0.0009017462   6 2525330 2556449
    ## 3  p2525330 p2556737 0.4896666611   6 2525330 2556737
    ## 4  p2525330 p2557053 0.0050843582   6 2525330 2557053
    ## 5  p2525330 p2557113 0.1313557757   6 2525330 2557113
    ## 6  p2525330 p2558106 0.4561908009   6 2525330 2558106
    ## 7  p2525330 p2624430 0.3253504820   6 2525330 2624430
    ## 8  p2525330 p2624453 0.3463576884   6 2525330 2624453
    ## 9  p2525330 p2747121 0.2089890347   6 2525330 2747121
    ## 10 p2525330 p2747122 0.2089330650   6 2525330 2747122
    ## 11 p2525330 p2747161 0.0226736004   6 2525330 2747161
    ## 12 p2525330 p2747170 0.2914167831   6 2525330 2747170
    ## 13 p2525330 p2747202 0.0218454389   6 2525330 2747202
    ## 14 p2525330 p2747230 0.2890620222   6 2525330 2747230
    ## 15 p2525330 p2747250 0.0240439548   6 2525330 2747250
    ## 16 p2525330 p2747291 0.2394054943   6 2525330 2747291
    ## 17 p2525330 p2747305 0.3787211156   6 2525330 2747305
    ## 18 p2525330 p2747308 0.0203400873   6 2525330 2747308
    ## 19 p2525330 p2747315 0.2690796078   6 2525330 2747315
    ## 20 p2525330 p2747329 0.0026214735   6 2525330 2747329
    ## 21 p2525330 p2747332 0.2658085703   6 2525330 2747332
    ## 22 p2525330 p2747377 0.0043181052   6 2525330 2747377
    ## 23 p2525330 p2747395 0.2726070560   6 2525330 2747395
    ## 24 p2525330 p2747397 0.2198628404   6 2525330 2747397
    ## 25 p2525330 p2747428 0.0144815126   6 2525330 2747428
    ## 26 p2525330 p2747444 0.2404931749   6 2525330 2747444
    ## 27 p2525330 p2747446 0.2393698489   6 2525330 2747446
    ## 28 p2525330 p2747487 0.2296461371   6 2525330 2747487
    ## 29 p2525330 p2747494 0.3048660208   6 2525330 2747494
    ## 30 p2525330 p2747548 0.0227952146   6 2525330 2747548
