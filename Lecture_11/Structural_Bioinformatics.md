Lecture 11. Structural Bioinformatics
================

Read CSV file and determine values
----------------------------------

``` r
df <- read.csv("Data Export Summary.csv", row.names = 1)

# Percentage of X-ray determined structures
xray <- df$Total[1]/sum(df$Total)*100
em <- df$Total[2]/sum(df$Total)*100

# All percentages
percent <- df$Total/sum(df$Total)*100
names(percent) <- row.names(df)
percent
```

    ##               X-Ray                 NMR Electron Microscopy 
    ##         89.43068692          8.40846082          1.88696977 
    ##               Other        Multi Method 
    ##          0.18325960          0.09062288

``` r
# Percentage of proteins
proteins <- sum(df$Proteins)/sum(df$Total)*100

# Print
round(xray, 1)
```

    ## [1] 89.4

``` r
round(em, 1)
```

    ## [1] 8.4

``` r
round(proteins, 1)
```

    ## [1] 92.8

Get data from anywhere using datapasta

``` r
tmp <- data.frame(stringsAsFactors=FALSE,
   Experimental.Method = c("X-Ray", "NMR", "Electron Microscopy", "Other",
                           "Multi Method", "Total"),
              Proteins = c(124770, 10988, 2057, 250, 127, 138192),
         Nucleic.Acids = c(1993, 1273, 31, 4, 5, 3306),
    ProteinNA.Complex = c(6451, 257, 723, 6, 2, 7439),
                 Other = c(10, 8, 0, 13, 1, 32),
                 Total = c(133224, 12526, 2811, 273, 135, 148969)
)
```

Bio3D
-----

``` r
library(bio3d)
```

``` r
pdb <- read.pdb("1hsg")
```

    ##   Note: Accessing on-line PDB file

``` r
library(bio3d.view)
library(rgl)
view(pdb, "overview", col = "sse")
```

    ## Computing connectivity from coordinates...

Extract the protein only portion of this PDB structure and wrie it out to a new PDB file.

``` r
prot <- atom.select(pdb, "protein")
head(pdb$atom[prot$atom, ])
```

    ##   type eleno elety  alt resid chain resno insert      x      y     z o
    ## 1 ATOM     1     N <NA>   PRO     A     1   <NA> 29.361 39.686 5.862 1
    ## 2 ATOM     2    CA <NA>   PRO     A     1   <NA> 30.307 38.663 5.319 1
    ## 3 ATOM     3     C <NA>   PRO     A     1   <NA> 29.760 38.071 4.022 1
    ## 4 ATOM     4     O <NA>   PRO     A     1   <NA> 28.600 38.302 3.676 1
    ## 5 ATOM     5    CB <NA>   PRO     A     1   <NA> 30.508 37.541 6.342 1
    ## 6 ATOM     6    CG <NA>   PRO     A     1   <NA> 29.296 37.591 7.162 1
    ##       b segid elesy charge
    ## 1 38.10  <NA>     N   <NA>
    ## 2 40.62  <NA>     C   <NA>
    ## 3 42.64  <NA>     C   <NA>
    ## 4 43.40  <NA>     O   <NA>
    ## 5 37.87  <NA>     C   <NA>
    ## 6 38.40  <NA>     C   <NA>

``` r
head(pdb$xyz[ ,prot$xyz])
```

    ## [1] 29.361 39.686  5.862 30.307 38.663  5.319

``` r
prot.pdb <- trim.pdb(pdb, prot)
write.pdb(prot.pdb, file = "prot.pdb")
```

Extract the ligand (i.e. drug) and write out to a seperate file.

``` r
lig <- atom.select(pdb, "ligand")
head(pdb$atom[lig$atom, ])
```

    ##        type eleno elety  alt resid chain resno insert      x      y     z
    ## 1515 HETATM  1517    N1 <NA>   MK1     B   902   <NA>  9.280 23.763 3.004
    ## 1516 HETATM  1518    C1 <NA>   MK1     B   902   <NA>  9.498 23.983 4.459
    ## 1517 HETATM  1519    C2 <NA>   MK1     B   902   <NA> 10.591 24.905 4.962
    ## 1518 HETATM  1520    C3 <NA>   MK1     B   902   <NA> 10.591 24.864 6.466
    ## 1519 HETATM  1521    O1 <NA>   MK1     B   902   <NA> 10.937 23.849 7.057
    ## 1520 HETATM  1522    N2 <NA>   MK1     B   902   <NA> 10.193 25.953 7.094
    ##      o     b segid elesy charge
    ## 1515 1 28.25  <NA>     N   <NA>
    ## 1516 1 30.30  <NA>     C   <NA>
    ## 1517 1 27.27  <NA>     C   <NA>
    ## 1518 1 28.85  <NA>     C   <NA>
    ## 1519 1 29.59  <NA>     O   <NA>
    ## 1520 1 22.29  <NA>     N   <NA>

``` r
head(pdb$xyz[ ,lig$xyz])
```

    ## [1]  9.280 23.763  3.004  9.498 23.983  4.459

``` r
lig.pdb <- trim.pdb(pdb, lig)
write.pdb(lig.pdb, file = "lig.pdb")
```

Access residues

``` r
pdb$seqres[50:60]
```

    ##     A     A     A     A     A     A     A     A     A     A     A 
    ## "ILE" "GLY" "GLY" "PHE" "ILE" "LYS" "VAL" "ARG" "GLN" "TYR" "ASP"
