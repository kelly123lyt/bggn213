---
title: "class13"
output: 
  html_document: 
    keep_md: yes
---




#Load genome project data

```r
genotype <- read.csv("373531-SampleGenotypes-Homo_sapiens_Variation_Sample_rs8067378.csv")
```


```r
table(genotype[,2])/nrow(genotype)
```

```
## 
##      A|A      A|G      G|A      G|G 
## 0.343750 0.328125 0.187500 0.140625
```

```r
#G/G is 14% so 14% are homozygous
```

Base quality scores from fastqsanger
Q8: Does the first sequence have good quality? 

```r
library(seqinr)
library(gtools)
phred <- asc( s2c("DDDDCDEDCDDDDBBDDDCC@") ) - 33
phred 
```

```
##  D  D  D  D  C  D  E  D  C  D  D  D  D  B  B  D  D  D  C  C  @ 
## 35 35 35 35 34 35 36 35 34 35 35 35 35 33 33 35 35 35 34 34 31
```

#Populational scale analysis



```r
populationgeno <- read.table("rs8067378_ENSG00000172057.6.txt")
```



```r
summary(populationgeno$exp[populationgeno$geno =="A/A"])
```

```
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##   11.40   27.02   31.25   31.82   35.92   51.52
```

```r
summary(populationgeno$exp[populationgeno$geno =="A/G"])
```

```
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##   7.075  20.626  25.065  25.397  30.552  48.034
```

```r
summary(populationgeno$exp[populationgeno$geno =="G/G"])
```

```
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##   6.675  16.903  20.074  20.594  24.457  33.956
```


```r
boxplot(exp~geno, data=populationgeno)
```

![](class13_files/figure-html/unnamed-chunk-6-1.png)<!-- -->

