# Pre-processing PLINK files

Create a new directory called `processed` to store processed files

```
cd grants_package
mkdir processed
cd processed
```

## Set minimum minor allele frequency (maf) to 0.05

```
plink --bfile ../grantsThomsons2017_maxMissing0.2 --maf 0.05 --no-sex --make-bed --out grants_maf0.05
```

## Change within family ID from 1 to locality/species

Use given information from `gazelle_popinfo_withPLINKid.txt` to set the within family IDs (2nd column in .fam) to either locality or species.

```
plink --bfile grants_maf0.05 --update-ids ../locality_as_within_fam_ID.txt --make-bed -out locality_grant
plink --bfile grants_maf0.05 --update-ids ../species_as_within_fam_ID.txt --make-bed -out species_grant
```
