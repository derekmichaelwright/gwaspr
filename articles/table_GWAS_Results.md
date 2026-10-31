# Summary Tables

``` r

library(gwaspr)
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

    ##                        SNP Chr       Pos negLog10_P                 Trait
    ## 1  Lcu.1GRN.Chr1p365986872   1 365986872  174.75587 Cotyledon_RedvsYellow
    ## 2  Lcu.1GRN.Chr1p365986872   1 365986872  153.93403 Cotyledon_RedvsYellow
    ## 3  Lcu.1GRN.Chr1p365986872   1 365986872  127.48387 Cotyledon_RedvsYellow
    ## 4  Lcu.1GRN.Chr1p361840399   1 361840399   83.66175 Cotyledon_RedvsYellow
    ## 5  Lcu.1GRN.Chr1p361840399   1 361840399   46.15717 Cotyledon_RedvsYellow
    ## 6  Lcu.1GRN.Chr1p361856257   1 361856257   39.47138 Cotyledon_RedvsYellow
    ## 7  Lcu.1GRN.Chr1p365986872   1 365986872   36.54581 Cotyledon_RedvsYellow
    ## 8  Lcu.1GRN.Chr1p365318023   1 365318023   34.75209 Cotyledon_RedvsYellow
    ## 9  Lcu.1GRN.Chr1p365318027   1 365318027   34.55222 Cotyledon_RedvsYellow
    ## 10   Lcu.1GRN.Chr5p1658484   5   1658484   33.97848        DTF_Nepal_2017
    ## 11 Lcu.1GRN.Chr1p368962767   1 368962767   33.59588 Cotyledon_RedvsYellow
    ## 12  Lcu.1GRN.Chr2p44546658   2  44546658   33.04272        DTF_Nepal_2017
    ## 13   Lcu.1GRN.Chr6p3269280   6   3269280   32.91602       DTF_Sask_2017_b
    ## 14 Lcu.1GRN.Chr1p437374598   1 437374598   30.55537       DTF_Sask_2017_b
    ## 15 Lcu.1GRN.Chr1p365786633   1 365786633   30.17116 Cotyledon_RedvsYellow
    ## 16 Lcu.1GRN.Chr1p365987296   1 365987296   29.66285 Cotyledon_RedvsYellow
    ## 17 Lcu.1GRN.Chr1p365785855   1 365785855   29.66141 Cotyledon_RedvsYellow
    ## 18 Lcu.1GRN.Chr1p361407757   1 361407757   29.48684 Cotyledon_RedvsYellow
    ## 19 Lcu.1GRN.Chr1p367671439   1 367671439   29.48364 Cotyledon_RedvsYellow
    ## 20 Lcu.1GRN.Chr1p365986872   1 365986872   29.18507 Cotyledon_RedvsYellow
    ##      Model   Threshold       P.value        MAF effect   H.B.P.Value Type
    ## 1    BLINK Significant 1.754402e-175 0.47096774     NA 5.901229e-170  NYC
    ## 2     MLMM Significant 1.164052e-154 0.47096774     NA 3.915488e-149  NYC
    ## 3  FarmCPU Significant 3.281941e-128 0.47096774     NA 1.103937e-122  NYC
    ## 4    BLINK Significant  2.178943e-84 0.04193548     NA  3.664623e-79  NYC
    ## 5  FarmCPU Significant  6.963600e-47 0.04193548     NA  1.171163e-41  NYC
    ## 6     MLMM Significant  3.377658e-40 0.04516129     NA  5.680663e-35  NYC
    ## 7      GLM Significant  2.845683e-37 0.47096774     NA  9.571939e-32  NYC
    ## 8      GLM Significant  1.769760e-35 0.47096774     NA  2.976445e-30  NYC
    ## 9      GLM Significant  2.803987e-35 0.47258065     NA  3.143895e-30  NYC
    ## 10   BLINK Significant  1.050791e-34 0.12345679     NA  3.534516e-29  NYC
    ## 11     GLM Significant  2.535833e-34 0.43225806     NA  2.132426e-29  NYC
    ## 12 FarmCPU Significant  9.063107e-34 0.06790123     NA  3.048530e-28  NYC
    ## 13 FarmCPU Significant  1.213338e-33 0.28105590     NA  4.081267e-28  NYC
    ## 14 FarmCPU Significant  2.783736e-31 0.47670807     NA  4.681785e-26  NYC
    ## 15     GLM Significant  6.742839e-31 0.43225806     NA  4.536137e-26  NYC
    ## 16     GLM Significant  2.173474e-30 0.42741935     NA  1.047873e-25  NYC
    ## 17     GLM Significant  2.180688e-30 0.43064516     NA  1.047873e-25  NYC
    ## 18 FarmCPU Significant  3.259584e-30 0.15000000     NA  3.654722e-25  NYC
    ## 19     GLM Significant  3.283698e-30 0.42741935     NA  1.380660e-25  NYC
    ## 20     MLM Significant  6.530328e-30 0.47096774     NA  2.196587e-24  NYC
    ##    negLog10_HBP      Effect
    ## 1     169.22906 -0.46837845
    ## 2     148.40721 -0.49056078
    ## 3     121.95706 -0.49543836
    ## 4      78.43597 -0.46035890
    ## 5      40.93138 -0.47155436
    ## 6      34.24560 -1.49473301
    ## 7      31.01900 -0.45583432
    ## 8      29.52630 -0.44709570
    ## 9      29.50253 -0.44512589
    ## 10     28.45167 12.56191501
    ## 11     28.67113 -0.43780205
    ## 12     27.51591 -1.31597036
    ## 13     27.38920 -2.25531798
    ## 14     25.32959 -0.19434407
    ## 15     25.34331 -0.41587018
    ## 16     24.97969 -0.41397885
    ## 17     24.97969 -0.41358916
    ## 18     24.43715 -0.02021127
    ## 19     24.85991 -0.40080355
    ## 20     23.65825 -0.43389870

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

    ##                        SNP Chr       Pos Hits max_negLog10_P
    ## 1  Lcu.1GRN.Chr1p365986872   1 365986872  774      174.75587
    ## 2  Lcu.1GRN.Chr1p361840399   1 361840399  252       83.66175
    ## 3    Lcu.1GRN.Chr5p1658484   5   1658484  300       33.97848
    ## 4  Lcu.1GRN.Chr1p368962767   1 368962767  142       33.59588
    ## 5   Lcu.1GRN.Chr2p44546658   2  44546658  334       33.04272
    ## 6    Lcu.1GRN.Chr6p3269280   6   3269280  206       32.91602
    ## 7  Lcu.1GRN.Chr1p437374598   1 437374598  235       30.55537
    ## 8  Lcu.1GRN.Chr1p446411579   1 446411579  396       28.18039
    ## 9  Lcu.1GRN.Chr1p357509147   1 357509147  128       27.09940
    ## 10 Lcu.1GRN.Chr1p351845237   1 351845237  227       27.01084
    ## 11 Lcu.1GRN.Chr1p371692181   1 371692181  161       24.49782
    ## 12 Lcu.1GRN.Chr1p381165558   1 381165558  164       23.48246
    ## 13 Lcu.1GRN.Chr4p416272481   4 416272481   64       23.48129
    ## 14  Lcu.1GRN.Chr1p27007485   1  27007485   62       23.11123
    ## 15 Lcu.1GRN.Chr4p423043571   4 423043571   53       22.58392
    ## 16 Lcu.1GRN.Chr3p239208186   3 239208186  160       22.55951
    ## 17 Lcu.1GRN.Chr4p441072457   4 441072457  185       22.40887
    ## 18 Lcu.1GRN.Chr4p431243285   4 431243285   14       22.34247
    ## 19 Lcu.1GRN.Chr2p414512712   2 414512712  136       22.24036
    ## 20 Lcu.1GRN.Chr4p418804416   4 418804416   60       22.22975
    ##                                                                Traits
    ## 1                                               Cotyledon_RedvsYellow
    ## 2                                               Cotyledon_RedvsYellow
    ## 3  DTF_Nepal_2017;DTF_Sask_2017_b;Cotyledon_RedvsYellow;DTF_Sask_2017
    ## 4                                               Cotyledon_RedvsYellow
    ## 5                                        DTF_Nepal_2017;DTF_Sask_2017
    ## 6  DTF_Sask_2017_b;Cotyledon_RedvsYellow;DTF_Nepal_2017;DTF_Sask_2017
    ## 7  DTF_Sask_2017_b;DTF_Sask_2017;Cotyledon_RedvsYellow;DTF_Nepal_2017
    ## 8  DTF_Sask_2017_b;Cotyledon_RedvsYellow;DTF_Nepal_2017;DTF_Sask_2017
    ## 9                                               Cotyledon_RedvsYellow
    ## 10               Cotyledon_RedvsYellow;DTF_Sask_2017_b;DTF_Nepal_2017
    ## 11                               Cotyledon_RedvsYellow;DTF_Nepal_2017
    ## 12               Cotyledon_RedvsYellow;DTF_Sask_2017_b;DTF_Nepal_2017
    ## 13                DTF_Sask_2017_b;DTF_Sask_2017;Cotyledon_RedvsYellow
    ## 14               Cotyledon_RedvsYellow;DTF_Sask_2017_b;DTF_Nepal_2017
    ## 15 DTF_Sask_2017_b;Cotyledon_RedvsYellow;DTF_Nepal_2017;DTF_Sask_2017
    ## 16               DTF_Sask_2017_b;Cotyledon_RedvsYellow;DTF_Nepal_2017
    ## 17                                      DTF_Sask_2017_b;DTF_Sask_2017
    ## 18                                      DTF_Sask_2017_b;DTF_Sask_2017
    ## 19 DTF_Sask_2017_b;Cotyledon_RedvsYellow;DTF_Sask_2017;DTF_Nepal_2017
    ## 20                                      DTF_Sask_2017_b;DTF_Sask_2017
    ##                        Models
    ## 1  BLINK;MLMM;FarmCPU;GLM;MLM
    ## 2  BLINK;FarmCPU;MLMM;GLM;MLM
    ## 3  BLINK;FarmCPU;GLM;MLMM;MLM
    ## 4                     GLM;MLM
    ## 5  FarmCPU;GLM;BLINK;MLMM;MLM
    ## 6  FarmCPU;GLM;BLINK;MLMM;MLM
    ## 7  FarmCPU;GLM;BLINK;MLMM;MLM
    ## 8                 FarmCPU;GLM
    ## 9                     GLM;MLM
    ## 10                    GLM;MLM
    ## 11      GLM;MLM;FarmCPU;BLINK
    ## 12                    GLM;MLM
    ## 13                        GLM
    ## 14     MLMM;GLM;FarmCPU;BLINK
    ## 15                GLM;FarmCPU
    ## 16                  GLM;BLINK
    ## 17                        GLM
    ## 18                        GLM
    ## 19                        GLM
    ## 20                        GLM
