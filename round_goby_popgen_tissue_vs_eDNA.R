### Round goby population genetics: tissue vs. eDNA analyses in field samples 
# This script accomplishes the following:
# PART 1: Tissue-based population genetic analysis
# PART 2: eDNA-based population genetic analysis
# PART 3: Comparison of approaches

#####################################################################################
############# PART 1: Tissue-based population genetic analysis #####################
#####################################################################################

# import tissue data as genind object 
library(adegenet)
library(pegas)
library(poppr)
library(devtools)
library(hierfstat)

# Import read count and genotype matrix 
# read in MultAmp output: col = samples, row = loci, cells = allele1/allele2:reads1,reads2
hap_genotype <- read.csv("/Users/kbja10/Documents/Cornell/Research/Round goby/Populations genetics: eDNA/Field_experiment_2019/data_analysis/ref_out_9nov2020/hap_genotype.csv", header=TRUE)
hap_genotype <- data.frame(lapply(hap_genotype[,-2], gsub, pattern=".*:", replacement=""))   # Remove haplotype column, remove everything before ":" in each cell
hap_genotype <- hap_genotype[,-grep("e_", colnames(hap_genotype))] # remove eDNA samples
hap_genotype_mat <- t(as.matrix(hap_genotype[,-1]))   # Transpose, turn into matrix, remove col and row headers

# Re-code genotypes for missing data 
# Create a new matrix w/ the total number of reads per locus per sample
counts_total <- matrix(NA, nrow = nrow(hap_genotype_mat), ncol = ncol(hap_genotype_mat)) # nrow = samples, ncol = loci
for (i in 1:(nrow(hap_genotype_mat)*ncol(hap_genotype_mat))){
  counts_total[i] <- sum(as.numeric(strsplit(hap_genotype_mat[i], ",")[[1]]))
}
# Create a new matrix that treats total # reads/locus <=10 as missing data
counts_less10 <- matrix(NA, nrow = nrow(hap_genotype_mat), ncol = ncol(hap_genotype_mat))
for (i in 1:(nrow(hap_genotype_mat)*ncol(hap_genotype_mat))){
  if (counts_total[i]<=10) {
    counts_less10[i]=NA} else {
      counts_less10[i]=counts_total[i]
    }
}
# Read in genotypes
genotypes <- read.csv("/Users/kbja10/Documents/Cornell/Research/Round goby/Populations genetics: eDNA/Field_experiment_2019/data_analysis/ref_out_9nov2020/hap_genotype_matrix.csv", header=TRUE)
genotypes <- genotypes[-grep("e_", genotypes$X),] # remove eDNA samples
rownames(genotypes) <- genotypes$X # rownames are sample names
genotypes <- as.matrix(genotypes[,!names(genotypes) %in% "X"]) # remove column of sample names
# Loop through to re-code alleles appropriately
#If total count <10, locus gets "0/0"
genotypes_edited <-  matrix(NA, nrow = nrow(hap_genotype_mat), ncol = ncol(hap_genotype_mat))
for (i in 1:(nrow(hap_genotype_mat)*ncol(hap_genotype_mat))){
  if (is.na(counts_less10[i]==TRUE)) {
    genotypes_edited[i] <- paste0(0, "|", 0)} else {
      genotypes_edited[i] = genotypes[i]
    }
}

rownames(genotypes_edited) <- rownames(genotypes)
colnames(genotypes_edited) <- colnames(genotypes)
genotypes_edited <- as.data.frame(genotypes_edited)
genotypes_edited[genotypes_edited == "0|0"] <- "NA|NA" # replace missing data with NAs
patterns <- c(3,8,19) # loci to remove
genotypes_edited <- genotypes_edited[,-patterns]
pops <- gsub(rownames(genotypes), pattern="t_", replacement="")
pops <- gsub(pops, pattern="_.*", replacement="")
geno.obj <- df2genind(genotypes_edited, sep="\\|", ploidy=2, type="codom", pop=pops,
                      loc.names=names(genotypes_edited), NA.char="NA|NA") # "|" must be preceded by double backslashes
geno.obj$tab
geno.cent <- scaleGen(geno.obj, center= T, scale= F, NA.method="mean")
geno.clusters <- find.clusters(geno.obj) #Retained 200 Axes (>90% of variance)
# looks to be 3-4 clusters

# Testing for Hardy-Weinberg equilibrium
hw.output <- data.frame(hw.test(geno.obj, 10000))
hw.output$observed_heterozygosity <- summary(geno.obj)[[6]]
hw.output$expected_heterozygosity <- summary(geno.obj)[[7]]
hw.output[hw.output$Pr.exact< 0.05/35 & hw.output$observed>hw.output$expected,] # bonferroni corrected pvalue

# look at PCA plot
x.geno.obj <- tab(geno.obj, freq=TRUE, NA.method="mean")
pca.geno.obj <- dudi.pca(x.geno.obj, center=TRUE, scale=FALSE, scannf = FALSE, nf = 200)
s.class(pca.geno.obj$li, fac=pop(geno.obj),col=transp(funky(15),.6),axesel=FALSE, cstar=0, cpoint=3)
add.scatter.eig(pca.geno.obj$eig[1:50],3,1,2, ratio=.2)
loadingplot(pca.geno.obj$c1^2)

# pcoa 
pca.geno.obj <- dudi.pco(dist(x.geno.obj))
s.class(pca.geno.obj$li, fac=pop(geno.obj),col=transp(funky(15),.6),axesel=FALSE, cstar=0, cpoint=3)
add.scatter.eig(pca.geno.obj$eig[1:50],3,1,2, ratio=.2)

# dapc
dapc1 <- dapc(geno.obj, pop=pops, n.pca=200, n.da=12)
scatter(dapc1)

# pairwise Fst
# geno.obj2 <- genind2hierfstat(geno.obj) 
# basic.stats(geno.obj2)
pairwise_fst_mat <- as.matrix(genet.dist(geno.obj, method="Ds"))
rownames(pairwise_fst_mat) <- levels(geno.obj$pop)
colnames(pairwise_fst_mat) <- levels(geno.obj$pop)
pairwise_fst_mat <- as.dist(pairwise_fst_mat)

# pairwise AFD
genpop.obj <- genind2genpop(geno.obj) # turn into genpop object
allele_freqs <- as.data.frame(t(makefreq(genpop.obj))) # calculate population allele frequencies
allele_freqs$locus <- data.frame(do.call('rbind', strsplit(as.character(rownames(allele_freqs)),'.',fixed=TRUE)))[,1] # create column of loci names
allele_freqs$allele <- data.frame(do.call('rbind', strsplit(as.character(rownames(allele_freqs)),'.',fixed=TRUE)))[,2] # create column of alleles
allele_freqs <- allele_freqs %>% select(locus, allele, everything()) # move locus and allele columns to the front

# AFD function
library(dplyr)
afd <- function(x,y,z){ # x is a data frame with 4 columns: locus, allele, Pop1, Pop2
  x %>%
    mutate(diff = abs(x[,y]-x[,z])) %>%
    group_by(locus) %>% 
    summarise(afd = sum(diff, na.rm=TRUE)/2)
}

# calculate pairwise AFD for all loci and populations 
pairwise_AFD <- combn(colnames(allele_freqs[,-c(1:2)]), 2, # all combinations of populations
                  FUN = function(x) mean(afd(allele_freqs,x[1],x[2])$afd)) # calculate mean AFD across all loci
mat <- matrix(nrow = ncol(allele_freqs[,-c(1:2)]), ncol = ncol(allele_freqs[,-c(1:2)])) # turn into pairwise matrix
mat[lower.tri(mat)] <- pairwise_AFD
rownames(mat) <- levels(geno.obj$pop)
colnames(mat) <- levels(geno.obj$pop)
pairwise_AFD_mat <- as.dist(mat)

# correlation between Fst and AFD pairwise distances 
library(ggplot2)
library(spaa)
library(ade4)
pairwise_df <- data.frame(dist2list(pairwise_fst_mat),afd_dist=dist2list(pairwise_AFD_mat)$value)
names(pairwise_df)[3] <- "fst_dist"
ggplot(pairwise_df,aes(x=fst_dist,y=afd_dist)) + 
  geom_point(size=3, color="darkgray") +
  xlab("Fst pairwise distance") + ylab("AFD pairwise distance") +
  theme_bw() +
  theme(text = element_text(size=20))

r1 <- mantel.rtest(pairwise_fst_mat,pairwise_AFD_mat,nrepet = 999)
plot(r1)
r1
# highly correlated! 


#####################################################################################
############# PART 2: eDNA-based population genetic analysis #####################
#####################################################################################
# Total normalized read frequencies for all counts > 10 and frequencies > 0.01
# Import read counts for all alleles at all loci
dataFiles <- lapply(Sys.glob("/Users/kbja10/Documents/Cornell/Research/Round goby/Populations genetics: eDNA/Field_experiment_2019/data_analysis/ref_out_9nov2020/haplotype2sample_raw/Nmel*.haplotype2sample.txt"), read.delim)
filenames <- Sys.glob("/Users/kbja10/Documents/Cornell/Research/Round goby/Populations genetics: eDNA/Field_experiment_2019/data_analysis/ref_out_9nov2020/haplotype2sample_raw/Nmel*.haplotype2sample.txt") # make list of file names
filenames <- gsub(".haplotype2sample\\.txt", "", filenames) # take out extra characters in file names
filenames <- gsub("/Users/kbja10/Documents/Cornell/Research/Round goby/Populations genetics: eDNA/Field_experiment_2019/data_analysis/ref_out_9nov2020/haplotype2sample_raw/", "", filenames)
names(dataFiles) <- filenames # rename files for each locus
site_codes <- c("LSH","LMM","LMK","ONO","ROC","ONE","OSW","CAY","ERIW","ECAN","ERIE","LSH","CRO")

# Begin loop
# create empty dataframe to populate with allele frequencies for each mesocosm
edna_allele_freqs <- data.frame(NULL)
j = 1 # individual locus indexed from "filenames"
# Create matrix with total number of reads/haplotype <=1 as missing data
for (file in dataFiles) { # for the reads at each locus for all samples 
  file_mat <- as.matrix(file[,-1]) # remove "haplotype" column, turn into matrix
  file_mat <- file_mat[,grep("e_", colnames(file_mat))] # subset to eDNA samples
  file_mat <- file_mat[,-grep("_BL", colnames(file_mat))] # remove blank samples
  file_mat[file_mat <= 1] <- 0 # replace counts <=1 with 0
  file_mat_norm <- as.data.frame(apply(file_mat, 2, function(x) x/sum(x, na.rm = TRUE))) # read count normalized per sample
  file_mat_norm[file_mat_norm<0.01] <- 0 # remove alleles below threshold of 0.01 per sample
  file_mat_norm <- as.data.frame(apply(file_mat_norm, 2, function(x) x/sum(x, na.rm = TRUE))) # overwrite with alleles < threshold removed
  edna_allele_freqs_temp <- data.frame(locus=rep(filenames[[j]]), 
                                       allele=file[,1])
  for (i in site_codes){ # for each site (15)
    subset <- file_mat_norm[,grep(i, colnames(file_mat_norm))]
    subset[is.na(subset)] <- 0
    freq <- rowSums(subset, na.rm = TRUE)/sum(subset, na.rm = TRUE) # read count normalized per sample
    edna_allele_freqs_temp <- cbind(edna_allele_freqs_temp, freq)
  }
  colnames(edna_allele_freqs_temp) <- c("locus","allele",site_codes)
  edna_allele_freqs <- rbind(edna_allele_freqs, edna_allele_freqs_temp)
  j <- j+1
}

edna_allele_freqs[is.na(edna_allele_freqs)] <- 0
edna_allele_freqs <- edna_allele_freqs[rowSums(edna_allele_freqs[,-c(1:2)])>0,]
nrow(edna_allele_freqs)
