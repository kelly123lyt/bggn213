---
title: "lecture12"
output: 
  html_document: 
    keep_md: yes
---



#Lecture 12


```r
library(bio3d)
```

Set up HIV for docking study

Get the protein first

```r
file.name <- get.pdb("1hsg")
```

```
## Warning in get.pdb("1hsg"): ./1hsg.pdb exists. Skipping download
```

```r
hiv <- read.pdb(file.name)
protein <- trim.pdb(hiv, "protein")
ligand <- trim.pdb(hiv, "ligand")
write.pdb(protein, file="1hsg_protein.pdb")
write.pdb(ligand, "1hsg_ligand.pdb")
```

process docking results for viewing in VMD


```r
library(bio3d)
res <- read.pdb("all.pdbqt", multi=TRUE)
write.pdb(res, "results.pdb")
```



```r
res <- read.pdb("all.pdbqt", multi=TRUE)
ori <- read.pdb("ligand.pdbqt")
rmsd(ori, res)
```

```
##  [1]  0.697  4.195 11.146 10.606 10.852 10.945 10.945  3.844  5.473  4.092
## [11] 10.404  5.574  3.448 11.396  6.126  3.848  8.237 11.196 10.981 11.950
```

