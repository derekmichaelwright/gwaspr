# Summary Tables

``` r

library("gwaspr")
```

## table_GWAS_Results()

The function
[`table_GWAS_Results()`](https://derekmichaelwright.github.io/gwaspr/reference/table_GWAS_Results.md)
creates a table of GWAS results.

Specifying a `folder` is all that is required to produce a GWAS results
table.

``` r

myResults <- table_GWAS_Results(
  # Specify a folder with GWAS results
  folder = "GWAS_Results/" )
```

Customization

``` r

myResults <- table_GWAS_Results(
  # Specify a folder with GWAS results
  folder = "GWAS_Results/",
  # Set threshold level
  threshold = 6.7, 
  # Set suggestive threshold level
  sug.threshold = 5 )
# Save
write.csv(myResults, "GWAS_Results_Table.csv", row.names = F)
```

``` r

read.csv("GWAS_Results_Table.csv")[1:20,]
```

    ##                        SNP Chr       Pos       P.value        MAF effect
    ## 1  Lcu.1GRN.Chr1p365986872   1 365986872 1.754402e-175 0.47096774     NA
    ## 2  Lcu.1GRN.Chr1p365986872   1 365986872 1.164052e-154 0.47096774     NA
    ## 3  Lcu.1GRN.Chr1p365986872   1 365986872 3.281941e-128 0.47096774     NA
    ## 4  Lcu.1GRN.Chr1p361840399   1 361840399  2.178943e-84 0.04193548     NA
    ## 5  Lcu.1GRN.Chr1p361840399   1 361840399  6.963600e-47 0.04193548     NA
    ## 6  Lcu.1GRN.Chr1p361856257   1 361856257  3.377658e-40 0.04516129     NA
    ## 7  Lcu.1GRN.Chr1p365986872   1 365986872  2.845683e-37 0.47096774     NA
    ## 8  Lcu.1GRN.Chr1p365318023   1 365318023  1.769760e-35 0.47096774     NA
    ## 9  Lcu.1GRN.Chr1p365318027   1 365318027  2.803987e-35 0.47258065     NA
    ## 10   Lcu.1GRN.Chr5p1658484   5   1658484  1.050791e-34 0.12345679     NA
    ## 11 Lcu.1GRN.Chr1p368962767   1 368962767  2.535833e-34 0.43225806     NA
    ## 12  Lcu.1GRN.Chr2p44546658   2  44546658  9.063107e-34 0.06790123     NA
    ## 13   Lcu.1GRN.Chr6p3269280   6   3269280  1.213338e-33 0.28105590     NA
    ## 14 Lcu.1GRN.Chr1p437374598   1 437374598  2.783736e-31 0.47670807     NA
    ## 15 Lcu.1GRN.Chr1p365786633   1 365786633  6.742839e-31 0.43225806     NA
    ## 16 Lcu.1GRN.Chr1p365987296   1 365987296  2.173474e-30 0.42741935     NA
    ## 17 Lcu.1GRN.Chr1p365785855   1 365785855  2.180688e-30 0.43064516     NA
    ## 18 Lcu.1GRN.Chr1p361407757   1 361407757  3.259584e-30 0.15000000     NA
    ## 19 Lcu.1GRN.Chr1p367671439   1 367671439  3.283698e-30 0.42741935     NA
    ## 20 Lcu.1GRN.Chr1p365986872   1 365986872  6.530328e-30 0.47096774     NA
    ##      H.B.P.Value   Model Type                 Trait negLog10_P negLog10_HBP
    ## 1  5.901229e-170   BLINK  NYC Cotyledon_RedvsYellow  174.75587    169.22906
    ## 2  3.915488e-149    MLMM  NYC Cotyledon_RedvsYellow  153.93403    148.40721
    ## 3  1.103937e-122 FarmCPU  NYC Cotyledon_RedvsYellow  127.48387    121.95706
    ## 4   3.664623e-79   BLINK  NYC Cotyledon_RedvsYellow   83.66175     78.43597
    ## 5   1.171163e-41 FarmCPU  NYC Cotyledon_RedvsYellow   46.15717     40.93138
    ## 6   5.680663e-35    MLMM  NYC Cotyledon_RedvsYellow   39.47138     34.24560
    ## 7   9.571939e-32     GLM  NYC Cotyledon_RedvsYellow   36.54581     31.01900
    ## 8   2.976445e-30     GLM  NYC Cotyledon_RedvsYellow   34.75209     29.52630
    ## 9   3.143895e-30     GLM  NYC Cotyledon_RedvsYellow   34.55222     29.50253
    ## 10  3.534516e-29   BLINK  NYC        DTF_Nepal_2017   33.97848     28.45167
    ## 11  2.132426e-29     GLM  NYC Cotyledon_RedvsYellow   33.59588     28.67113
    ## 12  3.048530e-28 FarmCPU  NYC        DTF_Nepal_2017   33.04272     27.51591
    ## 13  4.081267e-28 FarmCPU  NYC       DTF_Sask_2017_b   32.91602     27.38920
    ## 14  4.681785e-26 FarmCPU  NYC       DTF_Sask_2017_b   30.55537     25.32959
    ## 15  4.536137e-26     GLM  NYC Cotyledon_RedvsYellow   30.17116     25.34331
    ## 16  1.047873e-25     GLM  NYC Cotyledon_RedvsYellow   29.66285     24.97969
    ## 17  1.047873e-25     GLM  NYC Cotyledon_RedvsYellow   29.66141     24.97969
    ## 18  3.654722e-25 FarmCPU  NYC Cotyledon_RedvsYellow   29.48684     24.43715
    ## 19  1.380660e-25     GLM  NYC Cotyledon_RedvsYellow   29.48364     24.85991
    ## 20  2.196587e-24     MLM  NYC Cotyledon_RedvsYellow   29.18507     23.65825
    ##      Threshold      Effect
    ## 1  Significant -0.46837845
    ## 2  Significant -0.49056078
    ## 3  Significant -0.49543836
    ## 4  Significant -0.46035890
    ## 5  Significant -0.47155436
    ## 6  Significant -1.49473301
    ## 7  Significant -0.45583432
    ## 8  Significant -0.44709570
    ## 9  Significant -0.44512589
    ## 10 Significant 12.56191501
    ## 11 Significant -0.43780205
    ## 12 Significant -1.31597036
    ## 13 Significant -2.25531798
    ## 14 Significant -0.19434407
    ## 15 Significant -0.41587018
    ## 16 Significant -0.41397885
    ## 17 Significant -0.41358916
    ## 18 Significant -0.02021127
    ## 19 Significant -0.40080355
    ## 20 Significant -0.43389870

------------------------------------------------------------------------

## table_GWAS_Results_Summary()

The function
[`table_GWAS_Results_Summary()`](https://derekmichaelwright.github.io/gwaspr/reference/table_GWAS_Results_Summary.md)
creates a summary table of GWAS results.

Specifying a `folder` is all that is required to produce a GWAS summary
results table.

``` r

mySummary <- table_GWAS_Results_Summary(
  # Load table_GWAS_Results() object 
  xx = myResults )
```

Customization

``` r

mySummary <- table_GWAS_Results_Summary(
  # Load table_GWAS_Results() object 
  xx = myResults, 
  # Should markers be binned
  binMarkers = T, 
  # Set the size of the bin
  binSize = 5000000 )
# Save
write.csv(mySummary, "GWAS_Results_Summary_Table.csv", row.names = F)
```

``` r

read.csv("GWAS_Results_Summary_Table.csv")[1:20,]
```

    ##                        SNP Chr       Pos Hits        MAF max_negLog10_P
    ## 1  Lcu.1GRN.Chr1p365986872   1 365986872  697 0.47096774      174.75587
    ## 2  Lcu.1GRN.Chr1p361840399   1 361840399  179 0.04193548       83.66175
    ## 3    Lcu.1GRN.Chr5p1658484   5   1658484   88 0.12345679       33.97848
    ## 4  Lcu.1GRN.Chr1p368962767   1 368962767  109 0.43225806       33.59588
    ## 5   Lcu.1GRN.Chr2p44546658   2  44546658  229 0.06790123       33.04272
    ## 6    Lcu.1GRN.Chr6p3269280   6   3269280   54 0.28105590       32.91602
    ## 7  Lcu.1GRN.Chr1p437374598   1 437374598  107 0.47670807       30.55537
    ## 8  Lcu.1GRN.Chr1p446411579   1 446411579   45 0.27950311       28.18039
    ## 9  Lcu.1GRN.Chr1p357509147   1 357509147   86 0.36129032       27.09940
    ## 10 Lcu.1GRN.Chr1p351845237   1 351845237   41 0.34677419       27.01084
    ## 11 Lcu.1GRN.Chr1p371692181   1 371692181   87 0.33870968       24.49782
    ## 12 Lcu.1GRN.Chr1p381165558   1 381165558   42 0.37258065       23.48246
    ## 13 Lcu.1GRN.Chr4p416272481   4 416272481   28 0.27795031       23.48129
    ## 14  Lcu.1GRN.Chr1p27007485   1  27007485   15 0.07741935       23.11123
    ## 15 Lcu.1GRN.Chr4p423043571   4 423043571   22 0.28881988       22.58392
    ## 16 Lcu.1GRN.Chr3p239208186   3 239208186   56 0.46118012       22.55951
    ## 17 Lcu.1GRN.Chr4p441072457   4 441072457  107 0.39751553       22.40887
    ## 18 Lcu.1GRN.Chr4p431243285   4 431243285    3 0.26863354       22.34247
    ## 19 Lcu.1GRN.Chr2p414512712   2 414512712   23 0.22670807       22.24036
    ## 20 Lcu.1GRN.Chr4p418804416   4 418804416   15 0.28726708       22.22975
    ##    min_negLog10_P                     Models
    ## 1        5.002856 BLINK;MLMM;FarmCPU;GLM;MLM
    ## 2        5.021118 BLINK;FarmCPU;MLMM;GLM;MLM
    ## 3        5.068647 BLINK;FarmCPU;GLM;MLMM;MLM
    ## 4        5.053575                    GLM;MLM
    ## 5        5.056161 FarmCPU;GLM;BLINK;MLMM;MLM
    ## 6        5.028903 FarmCPU;GLM;BLINK;MLMM;MLM
    ## 7        5.028706 FarmCPU;GLM;BLINK;MLMM;MLM
    ## 8        6.528411                FarmCPU;GLM
    ## 9        5.081874                    GLM;MLM
    ## 10       5.150784                    GLM;MLM
    ## 11       5.001179      GLM;MLM;FarmCPU;BLINK
    ## 12       5.218290                    GLM;MLM
    ## 13       6.283576                        GLM
    ## 14       6.132473     MLMM;GLM;FarmCPU;BLINK
    ## 15       5.147086                GLM;FarmCPU
    ## 16      13.972620                  GLM;BLINK
    ## 17       6.544255                        GLM
    ## 18      18.687908                        GLM
    ## 19      18.627747                        GLM
    ## 20       6.211708                        GLM
    ##                                          Traits max_negLog10_HBP
    ## 1                         Cotyledon_RedvsYellow        169.22906
    ## 2                         Cotyledon_RedvsYellow         78.43597
    ## 3  DTF_Nepal_2017;DTF_Sask_2017_b;DTF_Sask_2017         28.45167
    ## 4                         Cotyledon_RedvsYellow         28.67113
    ## 5                  DTF_Nepal_2017;DTF_Sask_2017         27.51591
    ## 6  DTF_Sask_2017_b;DTF_Nepal_2017;DTF_Sask_2017         27.38920
    ## 7  DTF_Sask_2017_b;DTF_Sask_2017;DTF_Nepal_2017         25.32959
    ## 8  DTF_Sask_2017_b;DTF_Nepal_2017;DTF_Sask_2017         23.13070
    ## 9                         Cotyledon_RedvsYellow         22.95280
    ## 10                        Cotyledon_RedvsYellow         22.88197
    ## 11         Cotyledon_RedvsYellow;DTF_Nepal_2017         20.66120
    ## 12                        Cotyledon_RedvsYellow         19.80690
    ## 13                DTF_Sask_2017_b;DTF_Sask_2017         18.54795
    ## 14                        Cotyledon_RedvsYellow         18.06154
    ## 15                DTF_Sask_2017_b;DTF_Sask_2017         18.12774
    ## 16                              DTF_Sask_2017_b         18.12774
    ## 17                DTF_Sask_2017_b;DTF_Sask_2017         18.12774
    ## 18                              DTF_Sask_2017_b         18.12774
    ## 19                              DTF_Sask_2017_b         18.12774
    ## 20                DTF_Sask_2017_b;DTF_Sask_2017         18.12774
    ##    min_negLog10_HBP         min_P        max_P       min_HBP      max_HBP
    ## 1         2.2949280 1.754402e-175 9.934454e-06 5.901229e-170 5.070747e-03
    ## 2         0.9188869  2.178943e-84 9.525370e-06  3.664623e-79 1.205350e-01
    ## 3         1.2595187  1.050791e-34 8.537931e-06  3.534516e-29 5.501503e-02
    ## 4         2.3403431  2.535833e-34 8.839436e-06  2.132426e-29 4.567273e-03
    ## 5         1.2595187  9.063107e-34 8.786970e-06  3.048530e-28 5.501503e-02
    ## 6         0.6316665  1.213338e-33 9.356142e-06  4.081267e-28 2.335250e-01
    ## 7         0.6316665  2.783736e-31 9.360395e-06  4.681785e-26 2.335250e-01
    ## 8         3.8754994  6.600970e-29 2.962025e-07  7.401162e-24 1.331989e-04
    ## 9         2.3679743  7.954236e-28 8.281817e-06  1.114809e-23 4.285738e-03
    ## 10        2.4274280  9.753502e-28 7.066686e-06  1.312302e-23 3.737421e-03
    ## 11        1.7421589  3.178199e-25 9.972885e-06  2.181717e-21 1.810677e-02
    ## 12        2.4894361  3.292641e-24 6.049371e-06  1.559910e-20 3.240141e-03
    ## 13        3.7152165  3.301480e-24 5.205045e-07  2.831713e-19 1.926564e-04
    ## 14        1.7196030  7.740496e-24 7.371009e-07  8.678825e-19 1.907203e-02
    ## 15        0.7813179  2.606656e-23 7.127122e-06  7.451844e-19 1.654559e-01
    ## 16        8.9229276  2.757344e-23 1.065075e-14  7.451844e-19 1.194187e-09
    ## 17        3.8884305  3.900628e-23 2.855913e-07  7.451844e-19 1.292914e-04
    ## 18       16.1385144  4.544939e-23 2.051598e-19  7.451844e-19 7.269183e-17
    ## 19       16.0965685  5.749635e-23 2.356423e-19  7.451844e-19 8.006293e-17
    ## 20        3.6778994  5.891782e-23 6.141754e-07  7.451844e-19 2.099426e-04
