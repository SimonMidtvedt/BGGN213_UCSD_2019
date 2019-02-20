Lecture 12. Structural Bioinformatics (Part 2)
================

### Get PDB file

``` r
library(bio3d)
file.name <- get.pdb("1hsg")
```

    ## Warning in get.pdb("1hsg"): ./1hsg.pdb exists. Skipping download

``` r
hiv <- read.pdb(file.name)
```

### Extract protein only and ligand only segments and write a new PDB file

``` r
prot <- trim.pdb(hiv, "protein")
#prot.filename <- paste(file.name, "_protein.pdb", sep="")
lig <- trim.pdb(hiv, "ligand")

write.pdb(prot, "1hsg_protein.pdb")
write.pdb(lig, "1hsg_ligand.pdb")
```

### Convert docking results for viewing in VMD

``` r
res <- read.pdb("all.pdbqt", multi = TRUE)
write.pdb(res, "results.pdb")
```

### Calculate RMSD for qualitative analysis of results

``` r
ori <- read.pdb("1hsg_ligand.pdbqt")
rmsd(ori, res)
```

    ##  [1]  0.590 11.163 10.531  4.364 11.040  3.682  5.741  3.864  5.442 10.920
    ## [11]  4.318  6.249 11.084  8.929
