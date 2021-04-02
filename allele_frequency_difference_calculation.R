### Calculating absolute allele frequency difference (AFD; Berner, 2019)
### K. Andres 
### Last updated April 2021

library(dplyr)
afd <- function(x,y,z){ # x is a data frame with 4 columns: locus, allele, Pop1, Pop2
  x %>%
    mutate(diff = abs(x[,y]-x[,z])) %>%
    group_by(locus) %>% 
    summarise(afd = sum(diff, na.rm=TRUE)/2)
}

# Example provided in Berner, 2019
dat <- data.frame(locus=rep("Nmel1103", 4), allele=c(1,2,3,4),
                  Pop1=c(0.175,0.575,0.25,0), Pop2=c(0.094,0.125,0.5,0.281))
afd(dat,"Pop1","Pop2") # it works!

# Version to compute pairwise matrix across multiple populations and loci
dat <- data.frame(locus=rep(c("Locus_1","Locus_2"), each=4), allele=rep(c(1,2,3,4),2),
                  Pop1=c(0.175,0.575,0.25,0,0.25,0.25,0.25,0.25), Pop2=c(0.094,0.125,0.5,0.281,0,0,0.5,0.5),
                  Pop3=c(0.25,0.25,0.25,0.25,1,0,0,0), Pop4=c(0,0,0.5,0.5,0.5,0.25,0.12,0.13))

pairwise <- combn(colnames(dat[,-c(1:2)]), 2, # all combinations of populations
                  FUN = function(x) mean(afd(dat,x[1],x[2])$afd)) # calculate mean AFD across all loci
mat <- matrix(nrow = ncol(dat[,-c(1:2)]), ncol = ncol(dat[,-c(1:2)])) # turn into pairwise matrix
mat[lower.tri(mat)] <- pairwise
mat <- as.dist(mat)
