---
title: "Untitled"
output: github_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

from https://github.com/na399/popgen2017-pink_salmon/blob/master/scripts/5_estimate_Ne.r

```{r}

source('R_functions.r')

# Make an empty data frame to store the results
DF <- data.frame(Pop=rep(NA, 6), Site=rep(NA, 6), Lineage=rep(NA, 6), Ne_est=rep(NA, 6))

# Calculate 
Pops = c('Nome_ODD','Nome_EVEN', 'Koppen_ODD', 'Koppen_EVEN', 'Puget_ODD', 'Puget_EVEN')

for (index in 1:6){
    POP = Pops[index]
    site = strsplit(POP, split = '_')[[1]][1]
    lin = strsplit(POP, split = '_')[[1]][2]
    res = get_Ne(base_path = paste("./work/", POP, sep = ''))
    DF[index, ] = c(POP, site, lin, res$Ne_est)
}

# take a look at the results 
DF

write.table(DF, "./work/Ne_estimates.txt", sep = '\t', quote = FALSE, row.names = FALSE)

```