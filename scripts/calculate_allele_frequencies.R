### Round goby population genetics: tissue vs. eDNA analyses in field samples 
# This script accomplishes the following:
# PART 1: Genotype round goby tissues, remove bad loci, 
#         calculate tissue-based population allele frequencies
# PART 2: Calculate eDNA-based population allele frequencies
# Last updated: 1 November 2021

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
library(vegan)
library(ggpubr)

#########################################################################
################## PART 1: Genotype round goby tissues ##################
#########################################################################
# Input: DADA2 ASV-by-sample table for each locus
# Output: : data.frame containing allele data (hap_genotype_matrix.csv)

# Import read counts for all error-corrected ASVs (alleles)
dataFiles <- lapply(Sys.glob("/Users/kbja10/Documents/Cornell/Research/Round goby/Populations genetics: eDNA/Field_experiment_2019/data_analysis/out_dada2_28oct2021/Nmel*/*mod.csv"), read.csv)
filenames <- Sys.glob("/Users/kbja10/Documents/Cornell/Research/Round goby/Populations genetics: eDNA/Field_experiment_2019/data_analysis/out_dada2_28oct2021/Nmel*") # make list of file names
filenames <- gsub(".*(Nmel)", "\\1", filenames) # take out extra characters in file names
names(dataFiles) <- filenames # rename files for each locus

# Create empty dataframe to populate with genotypes at each locus
hap_genotype_matrix <- data.frame(matrix(nrow=285, ncol=35)) # 285 individuals, 35 loci
sample_names <- colnames(dataFiles[[1]])
sample_names <- sample_names[grep("_t_", sample_names)] # subset to tissue samples
sample_names <- gsub(".*__", "", sample_names) # remove duplicate sample names
colnames(hap_genotype_matrix) <- filenames
rownames(hap_genotype_matrix) <- sample_names

# Create empty dataframe to calculate read depth per locus 
locus_read_depth <- data.frame(matrix(nrow=285, ncol=35))
rownames(locus_read_depth) <- sample_names
colnames(locus_read_depth) <- filenames

# Loop thru locus files to determine genotypes and read depth  
j <-  1 # individual locus indexed from "filenames"
allele_diff <- 0.2 # proportion between 1st and 2nd allele
# Create matrix with total number of reads/haplotype <=1 as missing data
for (file in dataFiles) { # for the reads at each locus for all samples 
  file_mat <- file[,grep("_t_", colnames(file))] # subset to tissue samples
  # colnames(file_mat) <- gsub(".*JBL7W.", "\\1", colnames(file_mat)) # remove sequence lane info from sample names
  # colnames(file_mat) <- gsub("(\\d)[^0-9]+$", "\\1", colnames(file_mat)) # remove barcode info from sample names
  # colnames(file_mat) <- gsub("\\.","_",colnames(file_mat)) # change . to _ 
  colnames(file_mat) <- gsub(".*__", "", colnames(file_mat)) # remove duplicate sample names
  file_mat <- data.frame(allele=1:nrow(file_mat), file_mat) # add haplotype column
  # file_mat <- data.frame(allele=file$Id, file_mat) # add haplotype column
  for (i in rownames(hap_genotype_matrix)) { # for each sample
    if(i %in% colnames(file_mat)){ # if sample exists in file_mat
      locus_read_depth[i,j] <- sum(file_mat[,i]) # sample read depth at this locus
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

# Convert genotype data.frame to a genind object
# Add site codes to genotype matrix and reorder W to E
hap_genotype_matrix$site_codes <- gsub("t_(.+)_.*", "\\1", rownames(hap_genotype_matrix))
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
geno.obj # 285 individuals; 34 loci; 520 alleles

# Remove same alleles as mesocosm paper
patterns <- c("Nmel1103","Nmel1531","Nmel155","Nmel248","Nmel361","Nmel726") # loci to remove

# Check out locus read depth
apply(locus_read_depth, 2, mean)
ggplot(melt(locus_read_depth,value.name="reads"), aes(x=variable, y=reads)) + 
  geom_boxplot() + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# Remove loci with >50% missing data (>142 individuals missing data)
missing <- sapply(hap_genotype_matrix, function(x) sum(grepl("NA|NA",x))) # number of individuals out of 285 w/ missing data 
missing[missing>142]
patterns <- c("Nmel185","Nmel660",patterns)
hap_genotype_matrix <- hap_genotype_matrix[ , -which(names(hap_genotype_matrix) %in% patterns)]
ncol(hap_genotype_matrix) # 27 loci remaining

# Turn df back into genind object
geno.obj <- df2genind(hap_genotype_matrix[,-ncol(hap_genotype_matrix)], sep="\\|", # "|" must be preceded by double backslashes
                      ploidy=2, type="codom", pop=hap_genotype_matrix$site_codes,
                      loc.names=names(hap_genotype_matrix[,-ncol(hap_genotype_matrix)]), 
                      NA.char="NA|NA")
geno.obj # 285 individuals; 27 loci; 389 alleles

# Calculate allelic richness per locus
hap_genotype_matrix_split <- as.data.frame(cSplit(hap_genotype_matrix[,-ncol(hap_genotype_matrix)], 1:(ncol(hap_genotype_matrix)-1), '|'))  # split each locus by "|"
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
mean(allelic_richness$allelic_rich, na.rm = TRUE) # mean allelic richness = 14.41
min(allelic_richness$allelic_rich, na.rm = TRUE) # min allelic richness = 3
max(allelic_richness$allelic_rich, na.rm = TRUE) # max allelic richness = 29

# Calculate tissue-based population allele frequencies
genpop.obj <- genind2genpop(geno.obj) # turn into genpop object
allele_freqs <- as.data.frame(t(makefreq(genpop.obj))) # calculate population allele frequencies

# Write file for output
# write.csv(hap_genotype_matrix, "datasets/hap_genotype_matrix.csv",row.names=FALSE)
# write.csv(allele_freqs,  "datasets/tissue_allele_freqs.csv",row.names=TRUE)

# Convert file to something compatible with STRUCTURE
genind2structure <- function(obj, file="", pops=TRUE){
  pl <- max(obj@ploidy)   # get the max ploidy of the dataset
  S <- adegenet::nInd(obj)   # get the number of individuals
  tab <- data.frame(ind=rep(indNames(obj), each=pl))  # column of individual names to write; set up data.frame
  if(pops){   # column of pop ids to write
    popnums <- 1:adegenet::nPop(obj)
    names(popnums) <- as.character(unique(adegenet::pop(obj)))
    popcol <- rep(popnums[as.character(adegenet::pop(obj))], each=pl)
    tab <- cbind(tab, data.frame(pop=popcol))
  }
  loci <- adegenet::locNames(obj) # add columns for genotypes
  tab <- cbind(tab, matrix(-9, nrow=dim(tab)[1], ncol=adegenet::nLoc(obj),
                           dimnames=list(NULL,loci)))
  # begin going through loci
  for(L in loci){
    thesegen <- obj@tab[,grep(paste("^", L, "\\.", sep=""), 
                              dimnames(obj@tab)[[2]]), 
                        drop = FALSE] # genotypes by locus
    al <- 1:dim(thesegen)[2] # numbered alleles
    for(s in 1:S){
      if(all(!is.na(thesegen[s,]))){
        tabrows <- (1:dim(tab)[1])[tab[[1]] == indNames(obj)[s]] # index of rows in output to write to
        tabrows <- tabrows[1:sum(thesegen[s,])] # subset if this is lower ploidy than max ploidy
        tab[tabrows,L] <- rep(al, times = thesegen[s,])
      }
    }
  }
  write.table(tab, file=file, sep="\t", quote=FALSE, row.names=FALSE)
}
# genind2structure(geno.obj, "structure/genotypes_edited.txt", pops=TRUE)

#########################################################################
##### PART 2: Calculate eDNA-based population allele frequencies ########
#########################################################################
# Import read counts for all error-corrected ASVs (alleles)
dataFiles <- lapply(Sys.glob("/Users/kbja10/Documents/Cornell/Research/Round goby/Populations genetics: eDNA/Field_experiment_2019/data_analysis/out_dada2_28oct2021/Nmel*/*mod.csv"), read.csv)
filenames <- Sys.glob("/Users/kbja10/Documents/Cornell/Research/Round goby/Populations genetics: eDNA/Field_experiment_2019/data_analysis/out_dada2_28oct2021/Nmel*") # make list of file names
filenames <- gsub(".*(Nmel)", "\\1", filenames) # take out extra characters in file names
names(dataFiles) <- filenames # rename files for each locus

# Create empty dataframes to populate with read counts and frequencies for each eDNA sample
sample_names <- colnames(dataFiles[[1]])
sample_names <- sample_names[grep("e_", sample_names)] # subset to tissue samples
sample_names <- gsub(".*__", "", sample_names) # remove duplicate sample names
sample_names <- gsub("e_([^_]+_)", "", sample_names) # remove sequence run number
length(unique(sample_names)) # 62 samples
edna_read_counts <- NULL
edna_read_freqs <- NULL

# Read counts and normalized read frequencies 
# Combine reads from separate runs, keep triplicate samples per site separate
# Begin loop
j = 1 # individual locus indexed from "filenames"
LminLmax <- data.frame()
for (file in dataFiles) { # for the reads at each locus for all samples   file_mat <- file[,grep("e_", colnames(file))] # subset to tissue samples
  LminLmax <- rbind(LminLmax,data.frame(locus=filenames[[j]],
                   Lmin=min(nchar(file$alleleSequence)),
                   Lmax=max(nchar(file$alleleSequence))))
  file_mat <- file[,grep("e_", colnames(file))] # subset to eDNA samples
  file_mat[file_mat<=2] <- 0
  colnames(file_mat) <- gsub(".*__", "", colnames(file_mat)) # remove duplicate sample names
  file_mat <- data.frame(allele=1:nrow(file_mat), file_mat) # add allele column
  edna_read_counts_temp <- data.frame(locus=rep(filenames[[j]]), allele=file_mat$allele)
  edna_read_freqs_temp <- data.frame(locus=rep(filenames[[j]]), allele=file_mat$allele)
  for (i in unique(sample_names)) { # for each sample
    if(length(grep(i,colnames(file_mat)))>0){ # if sample exists in file_mat
      sample_reads <- file_mat[,c(grep(i,colnames(file_mat)))] # subset to just alleles and reads per sample 
      if (class(sample_reads)=="data.frame"){
        edna_read_counts_temp <- cbind(edna_read_counts_temp, rowSums(sample_reads))
        sample_reads_norm <- as.data.frame(apply(sample_reads, 2, function(x) x/sum(x, na.rm = TRUE))) # read count normalized per sample
        sample_reads_norm[is.na(sample_reads_norm)] <- 0 # replace NAs with 0
        if (sum(rowSums(sample_reads_norm))>0){
          edna_read_freqs_temp <-  cbind(edna_read_freqs_temp, rowSums(sample_reads_norm)/sum(rowSums(sample_reads_norm)))
        } else edna_read_freqs_temp <-  cbind(edna_read_freqs_temp, rowSums(sample_reads_norm))
      } else edna_read_counts_temp <- cbind(edna_read_counts_temp, sample_reads)
    } else edna_read_counts_temp <- cbind(edna_read_counts_temp, rep(0, nrow(edna_read_counts_temp)))
  }
  colnames(edna_read_counts_temp) <- c("locus","allele",unique(sample_names))
  edna_read_counts <- rbind(edna_read_counts, edna_read_counts_temp)
  colnames(edna_read_freqs_temp) <- c("locus","allele",unique(sample_names))
  edna_read_freqs <- rbind(edna_read_freqs, edna_read_freqs_temp)
  j=j+1
}

# min and max allele sizes for each locus
# write.csv(LminLmax, "datasets/allele_sizes.csv")

# Change "LSH" sites to "LHS"
colnames(edna_read_counts) <- gsub("LSH", "LHS", colnames(edna_read_counts))
colnames(edna_read_freqs) <- gsub("LSH", "LHS", colnames(edna_read_freqs))

# Move locus, allele to rownames
rownames(edna_read_counts) <- paste(edna_read_counts$locus, 
                                    edna_read_counts$allele, sep=".")
rownames(edna_read_freqs) <- paste(edna_read_freqs$locus, 
                                   edna_read_freqs$allele, sep=".")

# Remove alleles (rows) with no reads in any sample 
edna_read_counts <- edna_read_counts[rowSums(edna_read_counts[,-c(1:2)])>0,]
edna_read_freqs <- edna_read_freqs[rowSums(edna_read_freqs[,-c(1:2)])>0,]

# Check out reads in eDNA blanks (field, lab, PCR)
edna_blanks <- edna_read_counts[,grep("BL", colnames(edna_read_counts))]
edna_blanks <- edna_blanks[rowSums(edna_blanks)>0,]
edna_blanks
max(edna_blanks) # max 9 reads in blank
mean(colSums(edna_blanks != 0)) # avg # alleles w/ reads per blank
mean(colSums(edna_blanks)) # low read counts in all blanks

# Remove blanks and loci that we removed from tissue samples
edna_read_counts <- edna_read_counts[,-grep("BL|SEN|BUF", colnames(edna_read_counts))]
edna_read_freqs <- edna_read_freqs[,-grep("BL|SEN|BUF", colnames(edna_read_freqs))]
edna_read_counts <- edna_read_counts[!grepl("Nmel185|Nmel660|Nmel1103|Nmel1531|Nmel155|Nmel248|Nmel361|Nmel746", edna_read_counts$locus),]
edna_read_freqs <- edna_read_freqs[!grepl("Nmel185|Nmel660|Nmel1103|Nmel1531|Nmel155|Nmel248|Nmel361|Nmel746", edna_read_freqs$locus),]
nrow(edna_read_counts) # 627
nrow(edna_read_freqs) # 627

# Output file with replicate samples separate
# write.csv(edna_read_counts, "datasets/edna_read_counts_separate.csv", row.names=TRUE)
# write.csv(edna_read_freqs, "datasets/edna_read_freqs_separate.csv", row.names=TRUE)

### Cleaning the eDNA data
# Decisions to be made:
# 1) Remove entire loci: same loci removed from tissues
# 2) Remove samples: min. number of reads, max. number of missing loci
# 3) Remove sites: Poor amplification in all samples

# Total number of reads per locus
options(scipen=999) # no scientific notation
par(mar=c(6,6,4,2)+0.1, mgp=c(5,1,0)) # set margins and space btw axis label
locus_reads <- rowsum(edna_read_counts[,-(1:2)], edna_read_counts$locus) # reads per locus per sample
barplot(rowSums(locus_reads), log="y", ylab="Number of reads",xlab="Locus",
        names.arg=rownames(locus_reads), las=2, col="#69b3a2",
        main="Total number of reads per locus")

### Number of reads per sample 
edna_read_counts_long <- as.data.frame(pivot_longer(edna_read_counts, cols=ONE_01:LMK_03, names_to="sample", values_to="reads"))
edna_read_counts_sum <- as.data.frame(edna_read_counts_long %>% 
  group_by(sample) %>% 
  summarize(Reads=sum(reads), Alleles=sum(reads>0)))
n_loci <- edna_read_counts_long[edna_read_counts_long$reads>0,]
n_loci <- as.data.frame(n_loci %>% 
                          group_by(locus, sample) %>% 
                          count() %>%
                          group_by(sample) %>% 
                          count())
edna_read_counts_sum$Loci <- n_loci$n
edna_read_counts_sum$pop <- gsub("_.*", "", edna_read_counts_sum$sample)
edna_read_counts_sum$pop <- factor(edna_read_counts_sum$pop, levels=c("LMK","LMM","LHS","LSC","ERIW","ERIE","ROC","OSW","ECAN","CAY","CRO","ONO","ONE"))
edna_read_counts_sum <- edna_read_counts_sum[order(edna_read_counts_sum$pop),]
edna_read_counts_sum$sample <- factor(edna_read_counts_sum$sample, levels=unique(edna_read_counts_sum$sample))
edna_read_counts_sum_long <- edna_read_counts_sum %>%
  pivot_longer(Reads:Loci, names_to="group", values_to="reads")
edna_read_counts_sum_long$group <- factor(edna_read_counts_sum_long$group, levels=c("Reads","Alleles","Loci"))

supp.labs <- c("Lake Michigan W","Lake Michigan E","Lake Huron","Lake St. Clair",
               "Lake Erie W","Lake Erie E","Lake Ontario W","Lake Ontario E",
               "Erie Canal","Cayuga Lake","Cross Lake","Onondaga Lake", "Oneida Lake")
p1 <- ggplot(data=edna_read_counts_sum_long, aes(x=sample, y=reads, fill=pop)) +
  geom_col() + xlab("Sample") + ylab("Number") +
  scale_fill_manual(name="Population", values=cols, labels=supp.labs) +
  theme_bw() +
  theme(axis.text.y = element_text(size = 12), axis.text.x=element_blank()) +
  facet_grid(rows = vars(group), scales="free_y")
p1
# ggsave(filename = paste("figures/FigS5_sample_reads.pdf"), plot=p1, dpi=330, width=8, height=8, units="in")   

# Histogram of total number of missing loci per sample
hist(colSums(locus_reads<=10), breaks=20, 
     xlab="Number of loci with missing data", ylab="Number of eDNA samples",
     main="Total number of missing loci (out of 27) per sample", col="darkred")
abline(v=mean(colSums(locus_reads<=1)), col="green", lty=3, lwd=2)
abline(v=median(colSums(locus_reads<=1)), col="blue", lty=3, lwd=2)
legend(18, 10, legend=c("Mean", "Median"),
       col=c("green", "blue"), lty=3, lwd=2)

# Remove samples that don't exceed thresholds for max. # of missing loci
# max_missing_loci=5 (need 80% of loci) 27-(27*0.8)
# avg. allelic richness>=3 
to_remove <- edna_read_counts_sum[edna_read_counts_sum$Loci<=6,]$sample #  samples
to_remove <- c(to_remove, edna_read_counts_sum[edna_read_counts_sum$Alleles<=54,]$sample) #  samples
edna_read_counts <- edna_read_counts[, -which(colnames(edna_read_counts) %in% to_remove)]
length(unique(to_remove)) # 24 removed
ncol(edna_read_counts) # 15 remaining

# Make final files for output
site_codes <- c("LMK","LMM","LHS","LSC","ERIW","ERIE","ROC","OSW","ECAN","CAY","CRO","ONO","ONE")
edna_read_counts_final <- data.frame(matrix(nrow=nrow(edna_read_counts),ncol=length(site_codes)+2))
rownames(edna_read_counts_final) <- rownames(edna_read_counts)
colnames(edna_read_counts_final) <- c("locus","allele",site_codes)
edna_read_counts_final[,1:2] <- edna_read_counts[,1:2]

for (i in site_codes){ # for each site (13)
  if (length(grep(i,colnames(edna_read_counts)))>0){
    site_i <- edna_read_counts[,grep(i, colnames(edna_read_counts))]
    site_i <- data.frame(locus=edna_read_counts$locus, site_i)
    if (length(grep(i,colnames(edna_read_counts)))==1){
      edna_read_counts_final[,i] <- site_i[,-1]
    } else if (length(grep(i,colnames(edna_read_counts)))>1){
      edna_read_counts_final[,i] <- rowSums(site_i[,-1], na.rm=TRUE)
    }
  }
}

edna_read_counts_final <- edna_read_counts_final[rowSums(edna_read_counts_final[,-c(1:2)],na.rm=TRUE)>0,] # remove alleles with no reads at any site
nrow(edna_read_counts_final) # 543
# write.csv(edna_read_counts_final, "datasets/edna_read_counts.csv", row.names=TRUE)

edna_read_freqs_final <- NULL
for (i in unique(edna_read_counts_final$locus)){
  temp <- edna_read_counts_final[edna_read_counts_final$locus==i,]
  temp <- apply(temp[,-(1:2)], 2, function(x) x/sum(x, na.rm = TRUE))
  edna_read_freqs_final <- rbind(edna_read_freqs_final, temp)
}

edna_read_freqs_final <- cbind(edna_read_counts_final[,1:2], edna_read_freqs_final)
nrow(edna_read_freqs_final) # 543
# write.csv(edna_read_freqs_final, "datasets/edna_allele_freqs.csv", row.names=TRUE)

### Accumulation curves of eDNA vs. tissues
accum_curve_edna <- t(edna_read_counts[,-c(1:2)])
accum_curve_tissue <- geno.obj$tab
accum_curve_tissue[is.na(accum_curve_tissue)] <- 0
sites <- c("CAY","LHS","LMM","LSC","ONE","ONO","OSW")
accum_curve_tissue <- accum_curve_tissue[grepl(paste(sites,collapse="|"),rownames(accum_curve_tissue)),]
cols <- rev(colorRampPalette(brewer.pal(name="Spectral", n = 11)[c(1:5,8:11)])(13))
# pdf("figures/FigS4_alllele_accum_curve.pdf", width=12, height=6)
plot(specaccum(accum_curve_edna), ci.type="poly", 
     col=cols[2], lwd=3, ci.lty=0, ci.col=alpha(cols[2], 0.5),
     xlab="Number of samples", ylab="Alleles", xlim=c(0,nrow(accum_curve_tissue)))
plot(specaccum(accum_curve_tissue), ci.type="poly", 
     col=cols[1], lwd=3, ci.lty=0, ci.col=alpha(cols[1], 0.5),
     xlab="Number of samples", ylab="Alleles", add=TRUE)
legend("topright", legend=c("eDNA","Tissues"), col=cols[2:1], lwd=3,bty="n")
# dev.off()
# pdf("figures/FigS4_allele_accum_curve_locus.pdf", width=12, height=6)
for(i in unique(edna_read_counts$locus)){
  temp <- accum_curve_edna[,grep(i, colnames(accum_curve_edna))]
  plot(specaccum(temp), col=cols[2], lwd=3, ci = 0, ylim=c(0,60),
       xlab="Number of samples", ylab="Alleles", xlim=c(0,nrow(accum_curve_tissue)), add=TRUE)
  temp2 <- accum_curve_tissue[,grep(i, colnames(accum_curve_tissue))]
  plot(specaccum(temp2), col=cols[1], lwd=3, ci = 0,
       xlab="Number of samples", ylab="Alleles", add=TRUE)
  legend("topright", legend=c("eDNA","Tissues"), col=cols[2:1], lwd=3,bty="n")
}
# dev.off()
