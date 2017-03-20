# PCA and MDS
*see https://www.cog-genomics.org/plink/1.9/strat*


Create a new directory called `pca` in `grants_package/processed` to store output files

```
mkdir pca
```


## Localities

```
plink --bfile locality_grants --pca --within ../locality_clusters.txt -out pca/pca_locality
plink --bfile locality_grants --cluster  --within ../locality_clusters.txt -mds-plot 2  -out pca/mds_locality
```



```
cd pca
R
```

in `R`

```
pca <- read.table("mds_locality.mds", header=T)
plot(pca$C1 , pca$C2 , pch = 20 , cex = 1.5 , col = pca$IID)
```




## Species

```
plink --bfile species_grants --pca --within ../species_clusters.txt -out pca/pca_species
plink --bfile species_grants --cluster  --within ../species_clusters.txt -mds-plot 2  -out pca/mds_species
```



```
cd pca
R
```

in `R`

```
pca_species <- read.table("mds_species.mds", header=T)
plot(pca_species$C1 , pca_species$C2 , pch = 20 , cex = 1.5 , col = pca_species$IID)
```
