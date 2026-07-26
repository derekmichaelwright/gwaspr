# Run GWAS with GAPIT

``` r

library(gwaspr)
```

Before the `gwaspr` package can be used, we first need to run GWAS with
`GAPIT`.

## Genotype data (myG)

Our population, the LDP, was genotyped using an *exome caputure array*
and filtered with the following criteria to yield a data set of 336,367
SNPs, stored in our object `myG`.

- Only Bi-allelic = TRUE
- Minimum Read Depth = 3
- Minor Allele Frequency = \>5%
- Maximum Missing Frequency = 25%
- Heterozygous = \<25%

``` r

# Load our genotype file in hapmap format (note: header = F)
myG <- read.csv("myG_hmp.csv", header = F)
```

``` r

# Map + marker Info
myG[1:10,1:11]
```

    ##                      V1      V2    V3     V4     V5       V6     V7       V8
    ## 1                    V1      V2    V3     V4     V5       V6     V7       V8
    ## 2                    rs alleles chrom    pos strand assembly center protLSID
    ## 3  Lcu.1GRN.Chr1p853882     A/G     1 853882   <NA>     <NA>   <NA>     <NA>
    ## 4  Lcu.1GRN.Chr1p854117     G/A     1 854117   <NA>     <NA>   <NA>     <NA>
    ## 5  Lcu.1GRN.Chr1p854159     A/C     1 854159   <NA>     <NA>   <NA>     <NA>
    ## 6  Lcu.1GRN.Chr1p854174     C/G     1 854174   <NA>     <NA>   <NA>     <NA>
    ## 7  Lcu.1GRN.Chr1p870836     A/G     1 870836   <NA>     <NA>   <NA>     <NA>
    ## 8  Lcu.1GRN.Chr1p870898     C/T     1 870898   <NA>     <NA>   <NA>     <NA>
    ## 9  Lcu.1GRN.Chr1p870903     G/A     1 870903   <NA>     <NA>   <NA>     <NA>
    ## 10 Lcu.1GRN.Chr1p871108     A/G     1 871108   <NA>     <NA>   <NA>     <NA>
    ##           V9   V10    V11
    ## 1         V9   V10    V11
    ## 2  assayLSID panel QCcode
    ## 3       <NA>  <NA>   <NA>
    ## 4       <NA>  <NA>   <NA>
    ## 5       <NA>  <NA>   <NA>
    ## 6       <NA>  <NA>   <NA>
    ## 7       <NA>  <NA>   <NA>
    ## 8       <NA>  <NA>   <NA>
    ## 9       <NA>  <NA>   <NA>
    ## 10      <NA>  <NA>   <NA>

``` r

# genotrype calls
myG[1:10,12:17]
```

    ##             V12             V13            V14            V15          V16
    ## 1           V12             V13            V14            V15          V16
    ## 2  X3156.11_AGL CDC_Asterix_AGL CDC_Cherie_AGL CDC_Glamis_AGL CDC_Gold_AGL
    ## 3             G               A              G              A            G
    ## 4             G               G              G              A            G
    ## 5             C               C              C              C            C
    ## 6             G               G              G              G            G
    ## 7             G               G              N              G            G
    ## 8             T               C              N              C            T
    ## 9             G               G              N              A            G
    ## 10            A               A              A              A            A
    ##                  V17
    ## 1                V17
    ## 2  CDC_Greenstar_AGL
    ## 3                  A
    ## 4                  G
    ## 5                  A
    ## 6                  C
    ## 7                  A
    ## 8                  C
    ## 9                  G
    ## 10                 A

## Phenotype data (myY)

``` r

# Load our phenotype file
myY <- read.csv("myY.csv")
```

``` r

myY[1:20,]
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
    ## 11     CDC_LeMay_AGL          55.7          126.7                     0
    ## 12     CDC_Maxim_AGL          54.3          125.7                     1
    ## 13      CDC_QG.1_AGL          57.7          126.7                    NA
    ## 14 CDC_Red_Rider_AGL          57.0          123.5                     1
    ## 15   CDC_Redcoat_AGL          57.7          126.7                     1
    ## 16   CDC_Redwing_AGL          57.3          128.7                     1
    ## 17     CDC_Robin_AGL          52.0          124.7                     1
    ## 18   CDC_Rosebud_AGL          54.7          121.3                     1
    ## 19  CDC_Rosetown_AGL          57.7          123.7                     1
    ## 20   CDC_Rouleau_AGL          55.3          123.7                     1

## Covariate data (myCV)

``` r

# Load our covariate file
myCV <- read.csv("myCV.csv")
```

``` r

myCV[1:20,]
```

    ##                 Name        b
    ## 1    CDC_Asterix_AGL 0.000470
    ## 2      CDC_Rosie_AGL 0.000335
    ## 3       X3156.11_AGL 0.000369
    ## 4  CDC_Greenstar_AGL 0.000521
    ## 5     CDC_Cherie_AGL 0.000416
    ## 6     CDC_Glamis_AGL 0.000498
    ## 7       CDC_Gold_AGL 0.000543
    ## 8       CDC_Imax_AGL 0.000259
    ## 9    CDC_Impower_AGL 0.000497
    ## 10      CDC_KR.1_AGL 0.000746
    ## 11     CDC_LeMay_AGL 0.000411
    ## 12     CDC_Maxim_AGL 0.000405
    ## 13      CDC_QG.1_AGL 0.000468
    ## 14 CDC_Red_Rider_AGL 0.000251
    ## 15   CDC_Redcoat_AGL 0.000185
    ## 16   CDC_Redwing_AGL 0.000225
    ## 17     CDC_Robin_AGL 0.000461
    ## 18   CDC_Rosebud_AGL 0.000730
    ## 19  CDC_Rosetown_AGL 0.000230
    ## 20   CDC_Rouleau_AGL 0.000348

------------------------------------------------------------------------

## Run GWAS

Run GWAS on the 3 traits `myY`

``` r

library(GAPIT)
```

``` r

# Run GWAS
myGAPIT <- GAPIT( 
  # Phenotype data
  Y = myY,
  # Genotype data
  G = myG, 
  # Set PCA number
  PCA.total = 4,
  # Select Models
  model = c("GLM","MLM","MLMM","FarmCPU","BLINK"),
  # Extra functions to prevent errors
  Random.model = F, Phenotype.View = F )
```

------------------------------------------------------------------------

## GWAS with Covariate

Run GWAS on the `DTF_Sask_2017` trait with `b` as a covariate.

``` r

# Prep data
myY2 <- myY[,c("Name","DTF_Sask_2017")]
# Rename the trait to be rerun with a covariate
colnames(myY2)[2] <- "DTF_Sask_2017_b"
# Run GWAS with a covariate
myGAPIT <- GAPIT(
  # Phenotype data
  Y = myY2,
  # Genotype data
  G = myG,
  # Covariate data
  CV = myCV,
  # Set PCA number
  PCA.total = 0,
  # Select Models
  model = c("GLM","MLM","MLMM","FarmCPU","BLINK"),
  # Extra functions to prevent errors
  Random.model = F, Phenotype.View = F )
```
