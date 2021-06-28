### Round goby population genetics: tissue vs. eDNA analyses in field samples 
# This script accomplishes the following:
# PART 1: Genotype round goby tissues, remove bad loci, 
#         calculate tissue-based population allele frequencies
# PART 2: Calculate eDNA-based population allele frequencies

rm(list = ls())
setwd("/Users/kbja10/Github/eDNA_goby_popgen/")
library(splitstackshape)
library(dplyr)
library(adegenet)
library(RColorBrewer)
library(ggplot2)
library(hierfstat)
library(reshape2)
library(ape)
library(tidyr)
library(pegas)

#########################################################################
################## PART 1: Genotype round goby tissues ##################
#########################################################################
# Input: DADA2 ASV-by-sample table for each locus
# Output: : data.frame containing allele data (hap_genotype_matrix.csv)

# Import read counts for all error-corrected ASVs (alleles)
dataFiles <- lapply(Sys.glob("/Users/kbja10/Documents/Cornell/Research/Round goby/Populations genetics: eDNA/Field_experiment_2019/data_analysis/out_dada2/Nmel*/seqtab.nochim.csv"), read.csv)
filenames <- Sys.glob("/Users/kbja10/Documents/Cornell/Research/Round goby/Populations genetics: eDNA/Field_experiment_2019/data_analysis/out_dada2/Nmel*") # make list of file names
filenames <- gsub(".*(Nmel)", "\\1", filenames) # take out extra characters in file names
names(dataFiles) <- filenames # rename files for each locus

# Create empty dataframe to populate with genotypes at each locus
hap_genotype_matrix <- data.frame(matrix(nrow=285, ncol=35)) # 285 individuals, 35 loci
sample_names <- colnames(dataFiles[[1]])
sample_names <- sample_names[-grep("_e_", sample_names)] # subset to tissue samples
sample_names <- gsub(".*JBL7W.", "\\1", sample_names) # remove sequence lane info from sample names
sample_names <- gsub("(\\d)[^0-9]+$", "\\1", sample_names) # remove barcode info from sample names
rownames(hap_genotype_matrix) <- sample_names[-1]
colnames(hap_genotype_matrix) <- filenames

# Create empty dataframe to calculate read depth per locus 
locus_read_depth <- data.frame(matrix(nrow=285, ncol=35))
rownames(locus_read_depth) <- sample_names[-1]
colnames(locus_read_depth) <- filenames

# Loop thru locus files to determine genotypes and read depth  
j = 1 # individual locus indexed from "filenames"
allele_diff <- 0.2 # proportion between 1st and 2nd allele
# Create matrix with total number of reads/haplotype <=1 as missing data
for (file in dataFiles) { # for the reads at each locus for all samples 
  file_mat <- file[,-1] # remove sequence column
  file_mat <- file_mat[,-grep(".e.", colnames(file_mat))] # subset to tissue samples
  colnames(file_mat) <- gsub(".*JBL7W.", "\\1", colnames(file_mat)) # remove sequence lane info from sample names
  colnames(file_mat) <- gsub("(\\d)[^0-9]+$", "\\1", colnames(file_mat)) # remove barcode info from sample names
  colnames(file_mat) <- gsub("\\.","_",colnames(file_mat)) # change . to _ 
  file_mat <- data.frame(allele=1:nrow(file_mat), file_mat) # add haplotype column
  for (i in rownames(hap_genotype_matrix)) { # for each sample
    if(i %in% colnames(file_mat)){ # if sample exists in file_mat
      locus_read_depth[i,j] <- sum(file_mat[,i])
      genotype_ind <- file_mat[,c("allele",i)] # subset to just alleles and reads
      genotype_order <- genotype_ind[order(genotype_ind[,2], decreasing = TRUE),] # top 2 alleles based on read counts
      if (genotype_order[1,2]<10) { # if top allele has < 10 reads
        hap_genotype_matrix[i,j] <- paste("NA|NA") # code as missing data 
      } else if (genotype_order[2,2]>(allele_diff*genotype_order[1,2])) { # if 2nd allele is <20% of reads of top allele
        hap_genotype_matrix[i,j] <- paste(genotype_order[1,1], genotype_order[2,1], sep="|") # code as heterozygote
      } else hap_genotype_matrix[i,j] <- paste(genotype_order[1,1], genotype_order[1,1], sep="|") # otherwise homozygote for top allele
    } else {
      hap_genotype_matrix[i,j] <- paste("NA|NA") # if sample not in file_mat assign missing data
      locus_read_depth[i,j] <- 0
    }
  }
  j=j+1
}

# Check out locus read depth
apply(locus_read_depth, 2, mean)
ggplot(melt(locus_read_depth,value.name="reads"), aes(x=variable, y=reads)) + 
  geom_boxplot() + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# Remove loci with >50% missing data
sapply(hap_genotype_matrix, function(x) sum(grepl("NA|NA",x))) # number of individuals out of 285 w/ missing data 
patterns <- "Nmel185" 
hap_genotype_matrix <- hap_genotype_matrix[ , -which(names(hap_genotype_matrix) %in% patterns)]

# Convert genotype data.frame to a genind object
# Add site codes to genotype matrix and reorder W to E
hap_genotype_matrix$site_codes <- sub(rownames(hap_genotype_matrix), 
                                      pattern="_.*", replacement="")
hap_genotype_matrix$site_codes <- factor(hap_genotype_matrix$site_codes, 
                                         levels=c("LMK","LMM","LHS","LSC","ERIW","ERIE",
                                                  "ROC","OSW","ECAN","CAY","CRO","ONO","ONE"))
hap_genotype_matrix <- hap_genotype_matrix %>%
  arrange(hap_genotype_matrix$site_codes) 

# Import tissue data as genind object 
geno.obj <- df2genind(hap_genotype_matrix[,-ncol(hap_genotype_matrix)], sep="\\|", # "|" must be preceded by double backslashes
                      ploidy=2, type="codom", pop=hap_genotype_matrix$site_codes,
                      loc.names=names(hap_genotype_matrix[,-ncol(hap_genotype_matrix)]), 
                      NA.char="NA|NA")
geno.obj # 285 individuals; 34 loci; 597 alleles

# Exclude loci not in HWE and with heterozygote excess 
hw.output <- data.frame(hw.test(geno.obj, 1000))
hw.output$observed_heterozygosity <- summary(geno.obj)[[6]]
hw.output$expected_heterozygosity <- summary(geno.obj)[[7]]
hw.output[hw.output$Pr.exact< 0.05/34 & hw.output$observed>hw.output$expected,] # bonferroni corrected pvalue
# remove loci 1103, 1531, 155, 262, 361, 660 -- deviation from HWE and heterzygote excess
patterns <- c("Nmel1103","Nmel1531","Nmel155","Nmel262","Nmel361","Nmel660") # loci to remove
hap_genotype_matrix <- hap_genotype_matrix[ , -which(names(hap_genotype_matrix) %in% patterns)]
ncol(hap_genotype_matrix) # 28 loci remaining

# Turn df back into genind object
geno.obj <- df2genind(hap_genotype_matrix[,-ncol(hap_genotype_matrix)], sep="\\|", # "|" must be preceded by double backslashes
                      ploidy=2, type="codom", pop=hap_genotype_matrix$site_codes,
                      loc.names=names(hap_genotype_matrix[,-ncol(hap_genotype_matrix)]), 
                      NA.char="NA|NA")
geno.obj # 285 individuals; 28 loci; 438 alleles

# Calculate allelic richness per locus
hap_genotype_matrix_split <- as.data.frame(cSplit(hap_genotype_matrix, 1:ncol(hap_genotype_matrix), '|'))  # split each locus by "|"
rownames(hap_genotype_matrix_split) <- rownames(hap_genotype_matrix)
allelic_richness <- data.frame(allele = NA, allelic_rich = NA)
for (i in seq(1,ncol(hap_genotype_matrix_split)-1,2)){
  uniq_alleles <- c(hap_genotype_matrix_split[,i], hap_genotype_matrix_split[,i+1])
  allele <- colnames(hap_genotype_matrix_split)[i]
  allelic_rich <- length(unique(uniq_alleles[!is.na(uniq_alleles)]))
  allelic_richness_temp <- cbind(allele, allelic_rich)
  allelic_richness <- rbind(allelic_richness, allelic_richness_temp) 
}
allelic_richness <- allelic_richness[-1,]
allelic_richness$allelic_rich <- as.numeric(allelic_richness$allelic_rich)
mean(allelic_richness$allelic_rich, na.rm = TRUE) # mean allelic richness = 15.6
min(allelic_richness$allelic_rich, na.rm = TRUE) # min allelic richness = 5
max(allelic_richness$allelic_rich, na.rm = TRUE) # max allelic richness = 28

# Calculate tissue-based population allele frequencies
genpop.obj <- genind2genpop(geno.obj) # turn into genpop object
allele_freqs <- as.data.frame(t(makefreq(genpop.obj))) # calculate population allele frequencies

# Write file for output
# write.csv(hap_genotype_matrix, "datasets/hap_genotype_matrix.csv",row.names=FALSE)
# write.csv(allele_freqs,  "datasets/tissue_allele_freqs.csv",row.names=TRUE)

#########################################################################
##### PART 2: Calculate eDNA-based population allele frequencies ########
#########################################################################
# Total normalized read frequencies for all counts >= 1 and frequencies > 0.01
# Import read counts for all error-corrected ASVs (alleles)
dataFiles <- lapply(Sys.glob("/Users/kbja10/Documents/Cornell/Research/Round goby/Populations genetics: eDNA/Field_experiment_2019/data_analysis/out_dada2/Nmel*/seqtab.nochim.csv"), read.csv)
filenames <- Sys.glob("/Users/kbja10/Documents/Cornell/Research/Round goby/Populations genetics: eDNA/Field_experiment_2019/data_analysis/out_dada2/Nmel*") # make list of file names
filenames <- gsub(".*(Nmel)", "\\1", filenames) # take out extra characters in file names
names(dataFiles) <- filenames # rename files for each locus

# Create empty dataframes to populate with read counts and frequencies for each eDNA sample
sample_names <- colnames(dataFiles[[1]])
sample_names <- sample_names[grep("_e_", sample_names)] # subset to tissue samples
sample_names <- gsub(".*_e_", "\\1", sample_names) # remove sequence lane info from sample names
sample_names <- gsub("(.*)_\\w+_\\w+", "\\1", sample_names) # remove barcode info from sample names
sample_names <- gsub("\\w+_(\\w+_.*)", "\\1", sample_names) # remove sequence run number

edna_read_counts <- NULL
edna_read_freqs <- NULL

# Begin loop
j = 1 # individual locus indexed from "filenames"
for (file in dataFiles) { # for the reads at each locus for all samples 
  file_mat <- file[,-1] # remove sequence column
  file_mat <- file_mat[,grep(".e.", colnames(file_mat))] # subset to tissue samples
  colnames(file_mat) <- gsub(".*.e.", "\\1", colnames(file_mat)) # remove sequence lane info from sample names
  colnames(file_mat) <- gsub("\\.","_",colnames(file_mat)) # change . to _ 
  colnames(file_mat) <- gsub("(.*)_\\w+_\\w+", "\\1", colnames(file_mat)) # remove barcode info from sample names
  colnames(file_mat) <- gsub("\\w+_(\\w+_.*)", "\\1", colnames(file_mat)) # remove sequence run number
  file_mat <- data.frame(allele=1:nrow(file_mat), file_mat) # add allele column
  edna_read_counts_temp <- data.frame(locus=rep(filenames[[j]]), allele=file_mat$allele)
  for (i in unique(sample_names)) { # for each sample
    if(length(grep(i,colnames(file_mat)))>0){ # if sample exists in file_mat
      sample_reads <- file_mat[,c(grep(i,colnames(file_mat)))] # subset to just alleles and reads per sample 
      if (class(sample_reads)=="data.frame"){
        edna_read_counts_temp <- cbind(edna_read_counts_temp, rowSums(sample_reads))
      } else edna_read_counts_temp <- cbind(edna_read_counts_temp, sample_reads)
    } else edna_read_counts_temp <- cbind(edna_read_counts_temp, rep(0,nrow(edna_read_counts_temp)))
  }
  colnames(edna_read_counts_temp) <- c("locus","allele",unique(sample_names))
  #edna_read_counts_temp[,-c(1:2)][edna_read_counts_temp[,-c(1:2)] <= 1] <- 0 # replace counts <=1 with 0
  edna_read_freqs_temp <- data.frame(locus=rep(filenames[[j]]), allele=file_mat$allele, apply(edna_read_counts_temp[,-c(1:2)], 2, function(x) x/sum(x)))
  edna_read_freqs_temp[is.na(edna_read_freqs_temp)] <- 0
  edna_read_counts <- rbind(edna_read_counts, edna_read_counts_temp)
  edna_read_freqs <- rbind(edna_read_freqs, edna_read_freqs_temp)
  j=j+1
}

# Change "LSH" sites to "LHS"
colnames(edna_read_counts) <- gsub("LSH", "LHS", colnames(edna_read_counts))
colnames(edna_read_freqs) <- gsub("LSH", "LHS", colnames(edna_read_freqs))

# Move locus, allele to rownames
rownames(edna_read_counts) <- paste(edna_read_counts$locus, 
                                     edna_read_counts$allele, sep=".")
rownames(edna_read_freqs) <- paste(edna_read_freqs$locus, 
                                   edna_read_freqs$allele, sep=".")
# edna_read_counts <- edna_read_counts[,-c(1:2)]
# edna_read_freqs <- edna_read_freqs[,-c(1:2)]

### Cleaning the eDNA data
# Decisions to be made:
# 1) Remove entire loci: same loci removed from tissues
# 2) Remove samples: min. number of reads, max. number of missing loci
# 3) Remove sites: Poor amplification in all samples

# Check out reads in eDNA blanks (field, lab, PCR)
edna_blanks <- edna_read_counts[,grep("BL", colnames(edna_read_counts))]
edna_blanks <- edna_blanks[rowSums(edna_blanks)>0,]
edna_blanks
colSums(edna_blanks) # low read counts in all blanks

# Remove blanks and loci that we removed from tissue samples
edna_read_counts <- edna_read_counts[,-grep("BL", colnames(edna_read_counts))]
edna_read_freqs <- edna_read_freqs[,-grep("BL", colnames(edna_read_freqs))]
edna_read_counts <- edna_read_counts[!grepl("Nmel185|Nmel1103|Nmel1531|Nmel155|Nmel262|Nmel361|Nmel660", edna_read_counts$locus),]
edna_read_freqs <- edna_read_freqs[!grepl("Nmel185|Nmel1103|Nmel1531|Nmel155|Nmel262|Nmel361|Nmel660", edna_read_freqs$locus),]

# Output file with replicate samples separate
# write.csv(edna_read_counts, "datasets/edna_read_counts_separate.csv", row.names=TRUE)

# Total number of reads per locus
options(scipen=999) # no scientific notation
par(mar=c(6,6,4,2)+0.1, mgp=c(5,1,0)) # set margins and space btw axis label
locus_reads <- rowsum(edna_read_counts[,-(1:2)], edna_read_counts$locus) # reads per locus per sample
barplot(rowSums(locus_reads), log="y", ylab="Number of reads",xlab="Locus",
        names.arg=rownames(locus_reads), las=2, col="#69b3a2",
        main="Total number of reads per locus")

# Histogram of number of reads per sample (63 samples: 15 sites, 3 triplicates, 18 blanks)
dev.off()
hist(colSums(edna_read_counts[,-(1:2)]), breaks=50, 
     xlab="Number of reads", ylab="Number of eDNA samples", col="tomato",
     main="Total number of reads per sample \n (singletons removed)")
abline(v=mean(colSums(edna_read_counts[,-(1:2)])), col="green", lty=3, lwd=2)
abline(v=median(colSums(edna_read_counts[,-(1:2)])), col="blue", lty=3, lwd=2)
legend(125000, 15, legend=c("Mean", "Median"),
       col=c("green", "blue"), lty=3, lwd=2)

# Histogram of number of reads at low read counts -- most samples here have < 50 reads
hist(colSums(edna_read_counts[,-(1:2)][,colSums(edna_read_counts[,-(1:2)])<3000]), breaks=50, 
     xlab="Number of reads", ylab="Number of eDNA samples", col="tomato",
     main="Total number of reads per sample \n (subset to samples < 3,000 reads)")

# Histogram of total number of missing loci per sample
hist(colSums(locus_reads<=1), breaks=20, 
     xlab="Number of loci with missing data", ylab="Number of eDNA samples",
     main="Total number of missing loci (out of 35) per sample", col="darkred")
abline(v=mean(colSums(locus_reads<=1)), col="green", lty=3, lwd=2)
abline(v=median(colSums(locus_reads<=1)), col="blue", lty=3, lwd=2)
legend(20, 15, legend=c("Mean", "Median"),
       col=c("green", "blue"), lty=3, lwd=2)

# Contributor estimation function
MixtureLikelihood_Parsed <- function(x,p.v){
  require(RcppAlgos)
  if(length(p.v)>(2*x)) {sum.p <- 0} 
  else {
    counter <- 0
    if (length(p.v)>0){
      sum.p <- sum(p.v)^(2*x)
      if(length(p.v)>1){ # only conduct if > 1 allele is observed
        for(i in 1:(length(p.v)-1)){
          counter <- counter+1
          temp.combo <- comboGeneral(v=length(p.v), m = length(p.v)-i, repetition = FALSE, Parallel=TRUE,nThreads=10) #; i; dim(temp.combo);i<-i+1
          if(nrow(temp.combo)>15e6){
            ix <- 999999		
            row.sums.pow2x <- 1:(floor(nrow(temp.combo)/ix)+1)
            for(m in 1:length(row.sums.pow2x)){
              temp.m <- temp.combo[ ((m-1)*ix+1) : min(m*ix,nrow(temp.combo)), ]
              temp.combo.v <- c(temp.m)
              p.v.m <- matrix(p.v[temp.combo.v],nr=dim(temp.m)[1],nc=dim(temp.m)[2],byrow=F)
              row.sums.pow2x[m] <- sum(rowSums(p.v.m)^(2*x)) # reference as folows
            } # end m loop
            sum.p <- sum.p+((-1)^counter)*sum(row.sums.pow2x) } else 
            {
              # non-parsed
              temp.combo.v <- c(temp.combo)
              p.v.m <- matrix(p.v[temp.combo.v],nr=dim(temp.combo)[1],nc=dim(temp.combo)[2],byrow=F)
              row.sums.pow2x <- rowSums(p.v.m)^(2*x) # reference as folows
              sum.p <- sum.p+((-1)^counter)*sum(row.sums.pow2x) 
            } # end parsing option
        }# end i loop
      } # end if
    } else sum.p <- NA # if no alleles observed, assign NA
  }
  return(sum.p)
} # end function
MixtureLikelihood2 <- function(y, loci){
  Likelihood_df <- data.frame()
  for (i in 1:length(loci)){ # for each locus
    for (j in 1:y) { # for each of 1:y putative contributors
      Likelihood_df[j,i] <- MixtureLikelihood_Parsed(x=j,p.v=unlist(loci[i])) # function
    }  # end j loop
  } # end i loop
  Likelihood_df$Product <- apply(Likelihood_df, 1, function(x) prod(x, na.rm = TRUE)) # product across all loci, exclude NAs
  Likelihood_df$Contributors <- c(1:y) # putative numbers of contributors
  return(Likelihood_df[Likelihood_df$Product==max(Likelihood_df$Product),]) # find row with max likelihood across all loci
} # end function

# Contributor estimation when allele freqs specified per population
pop_allele_freqs_tissue <- cbind(rownames(allele_freqs),allele_freqs)
names(pop_allele_freqs_tissue)[1]<- "locus_allele"
edna_read_freqs$locus_allele <- rownames(edna_read_freqs)
combined_allele_freqs <- merge(edna_read_freqs, pop_allele_freqs_tissue, by = "locus_allele")
contrib_estim_mat <- matrix(nrow=39,ncol=2)
z <- 1
for (i in colnames(pop_allele_freqs_tissue[,-1])){
  site_i <- cbind(combined_allele_freqs[,1:3], combined_allele_freqs[,grep(i, colnames(combined_allele_freqs))])
  for (j in colnames(site_i[,4:(ncol(site_i)-1)])){
    sample_i <- site_i[,c("locus_allele","locus","allele",j,i)]
    sample_i <- sample_i[sample_i[,j]>0 & sample_i[,i]>0,] #  
    loci.list <- split(sample_i, sample_i$locus)
    loci.list <- lapply(loci.list, function(x) {x[,i]})
    contrib_estim_mat[z,1] <- j
    if (length(loci.list)==0) {contrib_estim_mat[z,2] <- 0
    } else {contrib_estim_mat[z,2] <- MixtureLikelihood2(y = 50, loci = loci.list)$Contributors}
    z <- z+1
  }
} 

contrib_estim <- as.data.frame(contrib_estim_mat)
colnames(contrib_estim) = c("sample", "estimation")
contrib_estim$estimation <- as.numeric(contrib_estim$estimation)

barplot(contrib_estim$estimation, ylab="Contributors",xlab="Sample",
        col="#AA4371", main="Estimated number of contributors per sample")
abline(h=mean(contrib_estim$estimation), col="green", lty=3, lwd=2)
abline(h=median(contrib_estim$estimation), col="blue", lty=3, lwd=2)
legend(40, 20, legend=c("Mean", "Median"),
       col=c("green", "blue"), lty=3, lwd=2)

# Remove samples that don't exceed thresholds for: 
# total # of reads, max. # of missing loci, number of estimated contributors
# min_reads <- 50
# max_missing_loci <- 15
to_remove <- names(edna_read_counts[,-(1:2)][colSums(edna_read_counts[,-(1:2)])<50])
to_remove <- c(to_remove, names(locus_reads[colSums(locus_reads<=1)>15]))
edna_read_counts <- edna_read_counts[, -which(colnames(edna_read_counts) %in% to_remove)]

# Make final files for output
site_codes <- c("LMK","LMM","LHS","LSC","ERIW","ERIE","ROC","OSW","ECAN","CAY","CRO","ONO","ONE")
edna_read_counts_final <- data.frame(matrix(nrow=nrow(edna_read_counts),ncol=length(site_codes)+2))
rownames(edna_read_counts_final) <- rownames(edna_read_counts)
colnames(edna_read_counts_final) <- c("locus","allele",site_codes)
edna_read_counts_final[,1:2] <- edna_read_counts[,1:2]

for (i in site_codes){ # for each site (13)
  if (length(grep(i,colnames(edna_read_counts)))>0){
  site_i <- edna_read_counts[,grep(i, colnames(edna_read_counts))]
  edna_read_counts_final[,i] <- rowSums(site_i)
  }
}

edna_read_counts_final <- edna_read_counts_final[rowSums(edna_read_counts_final[,-c(1:2)],na.rm=TRUE)>0,] # remove alleles with no reads at any site
nrow(edna_read_counts_final) # 2648
# write.csv(edna_read_counts_final, "datasets/edna_read_counts.csv", row.names=TRUE)

edna_read_freqs_final <- NULL
for (i in unique(edna_read_counts_final$locus)){
  temp <- edna_read_counts_final[edna_read_counts_final$locus==i,]
  temp <- apply(temp[,-(1:2)], 2, function(x) x/sum(x, na.rm = TRUE))
  edna_read_freqs_final <- rbind(edna_read_freqs_final, temp)
}

edna_read_freqs_final <- cbind(edna_read_counts_final[,1:2], edna_read_freqs_final)
nrow(edna_read_freqs_final) # 2648
# write.csv(edna_read_freqs_final, "datasets/edna_allele_freqs.csv", row.names=TRUE)
