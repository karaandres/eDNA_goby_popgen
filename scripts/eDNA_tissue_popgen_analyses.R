### Round goby population genetics: tissue vs. eDNA analyses in field samples 
# This script accomplishes the following:
# PART 1: Tissue-based population genetic analysis
# PART 2: eDNA-based population genetic analysis
# PART 3: Comparison of approaches

rm(list = ls())
setwd("/Users/kbja10/Github/eDNA_goby_popgen/")

library(adegenet)
library(ggplot2) 
library(poppr)
library(RColorBrewer)
library(hierfstat) 
library(reshape2)
library(ape)
library(dplyr)
library(tidyr)
library(ade4)
library(viridis)

##########################################################################
########### PART 1: Tissue-based population genetic analysis #############
##########################################################################
hap_genotype_matrix <- read.csv("datasets/hap_genotype_matrix.csv")
hap_genotype_matrix$site_codes <- factor(hap_genotype_matrix$site_codes, 
                                         levels=c("LMK","LMM","LHS","LSC","ERIW","ERIE",
                                                  "ROC","OSW","ECAN","CAY","CRO","ONO","ONE"))

# Turn tissue data as genind object 
geno.obj <- df2genind(hap_genotype_matrix[,-ncol(hap_genotype_matrix)], sep="\\|", # "|" must be preceded by double backslashes
                      ploidy=2, type="codom", pop=hap_genotype_matrix$site_codes,
                      loc.names=names(hap_genotype_matrix[,-ncol(hap_genotype_matrix)]), 
                      NA.char="NA|NA")
geno.obj # 285 individuals; 28 loci; 438 alleles

# K-means clustering: looks to be 2-3 clusters
# geno.clusters <- find.clusters(geno.obj, n.pca = 120) # retained 120 Axes (>90% of variance)
# pdf(file = "markdown_images/goby_find_clusters.pdf", width = 8, height = 8)
# plot(geno.clusters$Kstat,type="o", col="#408DBF", lwd=2, 
#       xlab="Number of clusters", ylab="BIC",cex.axis=1.5,cex.lab=1.5)
# dev.off()

# Look at PCA plot
x.geno.obj <- tab(geno.obj, freq=TRUE, NA.method="mean")
pca.geno.obj <- dudi.pca(x.geno.obj, center=TRUE, scale=FALSE, scannf = FALSE, nf = 120)
s.class(pca.geno.obj$li, fac=pop(geno.obj),col=transp(funky(15),.6),axesel=FALSE, cstar=0, cpoint=3)
add.scatter.eig(pca.geno.obj$eig[1:50],3,1,2, ratio=.2)
# loadingplot(pca.geno.obj$c1^2)

# Look at PCOA plot 
# pca.geno.obj <- dudi.pco(d = dist(x.geno.obj), nf = 120, scannf = FALSE)
# s.class(pca.geno.obj$li, fac=pop(geno.obj),col=transp(funky(15),.6),axesel=FALSE, cstar=0, cpoint=3)
# add.scatter.eig(pca.geno.obj$eig[1:50],3,1,2, ratio=.2)

# Look at DAPC plot
set.seed(867)
dapc1 <- dapc(geno.obj, pop=hap_genotype_matrix$site_codes, n.pca=120, n.da=12)
supp.labs <- c("Lake Michigan W","Lake Michigan E","Lake Huron","Lake St. Clair",
               "Lake Erie W","Lake Erie E","Lake Ontario W","Lake Ontario E",
               "Erie Canal","Cayuga Lake","Cross Lake","Onondaga Lake", "Oneida Lake")
names(supp.labs) <- unique(geno.obj$pop)
cols <- rev(colorRampPalette(brewer.pal(name="Spectral", n = 11)[c(1:5,8:11)])(13))
# pdf(file = "markdown_images/goby_popgen_dapc.pdf", width=8, height=6)
scatter(dapc1, bg="white", pch=20,  col=cols, 
        cex=3, clab=0, leg=TRUE, txt.leg=supp.labs, posi.leg="bottomleft")
# dev.off()

# Calculate and plot expected heterozygosity (Hs) within populations
hs_df <- data.frame(pop=factor(unique(hap_genotype_matrix$site_codes)), Hs=Hs(geno.obj, pop=hap_genotype_matrix$site_codes))
hs_df$pop <- factor(hs_df$pop, levels=c("LMK","LMM","LHS","LSC","ERIW","ERIE","ROC","OSW","ECAN","CAY","CRO","ONO","ONE"))
p <- ggplot(data=hs_df, aes(x=pop, y=Hs)) +
  geom_bar(stat="identity", fill="#408DBF") + ylim(0,1) +
  scale_x_discrete(labels=supp.labs) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 11, hjust = 1))
p

# Calculate pairwise Fst and plot heatmap
pairwise_fst_mat <- as.matrix(genet.dist(geno.obj, method="Ds"))
rownames(pairwise_fst_mat) <- levels(geno.obj$pop)
colnames(pairwise_fst_mat) <- levels(geno.obj$pop)
pairwise_fst_mat.tri <- pairwise_fst_mat
pairwise_fst_mat.tri[upper.tri(pairwise_fst_mat.tri, diag=TRUE)] <- NA # only keep upper triangle
pairwise_fst_long <- melt(pairwise_fst_mat.tri, na.rm =TRUE) # turn into long df

fst_heatmap <- ggplot(data = pairwise_fst_long, aes(Var2, Var1, fill = value))+ 
  geom_tile(color = "white") + 
  scale_fill_gradientn(colours = brewer.pal(9,"YlGnBu"), name="FST")  + 
  scale_x_discrete(labels= supp.labs) +
  scale_y_discrete(labels= supp.labs) +
  theme_bw() +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, size = 11, hjust = 1),
        axis.text.y = element_text(size = 12)) + coord_fixed()
fst_heatmap
# ggsave(filename = paste("markdown_images/fst_heatmap_tissue.pdf"), plot=fst_heatmap, dpi = 300,  width = 6, height = 6, units = "in")   

# Plot neighbor joining tree from distance matrix
theTree <- pairwise_fst_mat %>%
  nj() %>%    # calculate neighbor-joining tree
  ladderize() # organize branches by clade
supp.labs_df <- as.data.frame(supp.labs)
theTree$tip.label <- supp.labs_df[,1][match(theTree$tip.label, rownames(supp.labs_df))]
# pdf(file = "markdown_images/goby_popgen_tree.pdf", width = 5, height = 5)
plot(theTree)
add.scale.bar()
# dev.off() # save the tree

# pairwise AFD
genpop.obj <- genind2genpop(geno.obj) # turn into genpop object
allele_freqs <- as.data.frame(t(makefreq(genpop.obj))) # calculate population allele frequencies
allele_freqs$locus <- data.frame(do.call('rbind', strsplit(as.character(rownames(allele_freqs)),'.',fixed=TRUE)))[,1] # create column of loci names
allele_freqs$allele <- data.frame(do.call('rbind', strsplit(as.character(rownames(allele_freqs)),'.',fixed=TRUE)))[,2] # create column of alleles
allele_freqs <- allele_freqs %>% dplyr::select(locus, allele, everything()) # move locus and allele columns to the front

# AFD function
afd <- function(x,y,z){ # x is a data frame with 4 columns: locus, allele, Pop1, Pop2
  x %>%
    mutate(diff = abs(x[,y]-x[,z])) %>%
    group_by(locus) %>% 
    summarise(afd = sum(diff, na.rm=TRUE)/2)
}

# calculate pairwise AFD for all loci and populations 
pairwise_AFD <- data.frame(Var1=combn(colnames(allele_freqs[,-c(1:2)]), 2)[1,],
                           Var2=combn(colnames(allele_freqs[,-c(1:2)]), 2)[2,],
                           AFD=combn(colnames(allele_freqs[,-c(1:2)]), 2, # all combinations of populations
                                     FUN = function(x) mean(afd(allele_freqs,x[1],x[2])$afd))) # calculate mean AFD across all loci
pairwise_AFD_mat <- data.matrix(pivot_wider(pairwise_AFD, names_from=Var1, values_from=AFD))
pairwise_AFD_mat <- pairwise_AFD_mat[,-1] # remove column with rownames
pairwise_AFD_mat <- cbind(pairwise_AFD_mat, rep(NA,12))
pairwise_AFD_mat <- rbind(rep(NA,13), pairwise_AFD_mat)

# turn into symmetrical matrix so we can reorder populations E-W
makeSymm <- function(m) {
  m[upper.tri(m)] <- t(m)[upper.tri(m)]
  return(m)
}
pairwise_AFD_mat <- makeSymm(pairwise_AFD_mat)
rownames(pairwise_AFD_mat) <- levels(geno.obj$pop)
colnames(pairwise_AFD_mat) <- levels(geno.obj$pop)
pairwise_AFD_mat.tri <- pairwise_AFD_mat
pairwise_AFD_mat.tri[upper.tri(pairwise_AFD_mat.tri, diag=TRUE)] <- NA # only keep lower triangle
pairwise_AFD_long <- melt(pairwise_AFD_mat.tri, na.rm =TRUE) # turn into long df

# correlation between Fst and AFD pairwise distances
pairwise_df <- data.frame(pairwise_fst_long,pairwise_AFD_long)
names(pairwise_df)[c(3,6)] <- c("fst_dist", "afd_dist")
p <- ggplot(pairwise_df,aes(x=fst_dist,y=afd_dist)) + 
  geom_point(size=3, col="#408DBF") +
  xlab("Nei's standard genetic distance") + ylab("AFD pairwise distance") +
  xlim(0,0.5) + ylim(0,0.5) +
  theme_bw() +
  theme(text = element_text(size=20))
p
# ggsave(filename = paste("markdown_images/AFD_Fst_correlation.eps"), plot=p, dpi = 300,  width = 6, height = 6, units = "in")   

r1 <- mantel.rtest(as.dist(pairwise_fst_mat),as.dist(pairwise_AFD_mat),nrepet = 999)
plot(r1)
r1
# highly correlated! 


##########################################################################
############# PART 2: eDNA-based population genetic analysis #############
##########################################################################

edna_allele_freqs <- read.csv("datasets/edna_allele_freqs.csv", header=TRUE)
colnames(edna_allele_freqs)[1] <- "locus_allele"
rownames(edna_allele_freqs) <- edna_allele_freqs[,1]

# Subset to only alleles found in round goby tissues and recalculate frequencies
edna_allele_freqs <- edna_allele_freqs[edna_allele_freqs$locus_allele %in% colnames(geno.obj$tab), ]
edna_freqs <- edna_allele_freqs[,-1]
nrow(edna_freqs) # 303
edna_freqs_adjusted <- data.frame()
for (locus in unique(edna_freqs$locus)){
  sub <- edna_freqs[edna_freqs$locus==locus,] # subset to each locus
  sub_new <- cbind(sub[,(1:2)], apply(sub[,-(1:2)], 2, function(x) x/sum(x)))
  edna_freqs_adjusted <- rbind(edna_freqs_adjusted, sub_new)
}
edna_freqs_adjusted[is.na(edna_freqs_adjusted)] <- 0

# Reorder sites east to west
edna_freqs_adjusted <- edna_freqs_adjusted[c("locus","allele","LMK","LMM","LHS","LSC","ERIW","ERIE",
                                             "ROC","OSW","ECAN","CAY","CRO","ONO","ONE")]

# calculate pairwise AFD for all loci and populations 
pairwise_AFD_eDNA <- data.frame(Var1=combn(colnames(edna_freqs_adjusted[,-c(1:2)]), 2)[1,],
                                Var2=combn(colnames(edna_freqs_adjusted[,-c(1:2)]), 2)[2,],
                                AFD=combn(colnames(edna_freqs_adjusted[,-c(1:2)]), 2, # all combinations of populations
                                          FUN = function(x) mean(afd(edna_freqs_adjusted,x[1],x[2])$afd))) # calculate mean AFD across all loci
pairwise_AFD_eDNA_mat <- data.matrix(pivot_wider(pairwise_AFD_eDNA, names_from=Var1, values_from=AFD))
pairwise_AFD_eDNA_mat <- pairwise_AFD_eDNA_mat[,-1] # remove column with rownames
pairwise_AFD_eDNA_mat <- cbind(pairwise_AFD_eDNA_mat, rep(NA,12))
pairwise_AFD_eDNA_mat <- rbind(rep(NA,13), pairwise_AFD_eDNA_mat)

# turn into symmetrical matrix so we can reorder populations E-W
pairwise_AFD_eDNA_mat <- makeSymm(pairwise_AFD_eDNA_mat)
rownames(pairwise_AFD_eDNA_mat) <- colnames(edna_freqs_adjusted[,-(1:2)])
pairwise_AFD_eDNA_mat.tri <- pairwise_AFD_eDNA_mat
pairwise_AFD_eDNA_mat.tri[upper.tri(pairwise_AFD_eDNA_mat.tri, diag=TRUE)] <- NA # only keep lower triangle
pairwise_AFD_eDNA_long <- melt(pairwise_AFD_eDNA_mat.tri, na.rm =TRUE) # turn into long df

# correlation between tissue-based and AFD-based pairwise distances
pairwise_df <- data.frame(pairwise_AFD_long, pairwise_AFD_eDNA_long)
sites <- c("CAY","ERIE","LHS","LMM","LSC","ONE","ONO","OSW")
pairwise_df <- filter(pairwise_df, grepl(paste(sites,collapse="|"),Var1))
pairwise_df <- filter(pairwise_df, grepl(paste(sites,collapse="|"),Var2))
names(pairwise_df)[c(3,6)] <- c("afd_dist_tissue", "afd_dist_eDNA")
cols <- rev(brewer.pal(name="Spectral", n=8))
p <- ggplot(pairwise_df,aes(x=afd_dist_tissue,y=afd_dist_eDNA,colour=Var2,shape=Var1)) + 
  geom_point(size=5, stroke=1.5) +
  xlab("Tissue-based AFD") + ylab("eDNA-based AFD") +
  #xlim(0,1) + ylim(0,1) +
  scale_color_manual(name="Pop 2", values=cols, labels=supp.labs) +
  scale_shape_manual(name="Pop 1", values=c(15:18,3,4,8), labels=supp.labs) +
  theme_bw() +
  theme(text = element_text(size=20))
p
# ggsave(filename = paste("markdown_images/eDNA_tissue_AFD_correlation.eps"), plot=p, dpi = 300,  width = 9, height = 6, units = "in")   
cor.test(pairwise_df$afd_dist_tissue, pairwise_df$afd_dist_eDNA)
r1 <- mantel.rtest(as.dist(pairwise_AFD_mat), as.dist(pairwise_AFD_eDNA_mat),nrepet = 999)
plot(r1)
r1

#############################################################################################
############## Correlation between population allele/read frequencies ###############
#############################################################################################
# Combine eDNA and tissue allele freqs
tissue_freqs_long <- pivot_longer(allele_freqs,cols=LMK:ONE,names_to="population")
edna_freqs_long <- pivot_longer(edna_freqs_adjusted,cols=LMK:ONE,names_to="population")
total_allele_freqs <- merge(edna_freqs_long, tissue_freqs_long, by = c("locus","allele","population"))
names(total_allele_freqs) <- c("locus","allele","population","edna_freq","tissue_freq")
total_allele_freqs$population <- factor(total_allele_freqs$population, 
                                      levels=c("LMK","LMM","LHS","LSC","ERIW","ERIE","ROC","OSW","ECAN","CAY","CRO","ONO","ONE"))
total_allele_freqs <- total_allele_freqs %>%
  arrange(total_allele_freqs$population) 

# Subset to only the alleles present in tissues from each population and recalculate frequencies
total_allele_freqs <-  total_allele_freqs[total_allele_freqs$tissue_freq>0,]
nrow(total_allele_freqs) # 3132
total_allele_freqs_adjusted <- data.frame()
for (locus in unique(total_allele_freqs$locus)){
  sub <- total_allele_freqs[total_allele_freqs$locus==locus,] # subset to each locus
  for (pop in unique(sub$population)){
    sub_pop <- sub[sub$population==pop,]
    sub_pop$edna_freq <- sub_pop$edna_freq/sum(sub_pop$edna_freq)
    total_allele_freqs_adjusted <- rbind(total_allele_freqs_adjusted, sub_pop)
  }
}
total_allele_freqs_adjusted[is.na(total_allele_freqs_adjusted)] <- 0
sites <- c("CAY","ERIE","LHS","LMM","LSC","ONE","ONO","OSW")
total_allele_freqs_adjusted <- filter(total_allele_freqs_adjusted, 
                                      grepl(paste(sites,collapse="|"),population))
unique(total_allele_freqs_adjusted$population)

# Correlation coefficient between tissue and eDNA allele frequencies per site
cors <- total_allele_freqs_adjusted %>%
  group_by(population) %>%
  dplyr::summarize(cor=round(cor(tissue_freq, edna_freq, use="complete.obs"), 2))
# write.csv(cors, "/Users/kbja10/Github/eDNA_goby_popgen/allele_freq_cors.csv")

p <- ggplot(total_allele_freqs_adjusted, aes(x=tissue_freq, y=edna_freq, color=locus)) +
  geom_point() + ylab("Field eDNA allele frequency\n") + xlab("\nTissue allele frequency") + 
  scale_x_continuous(breaks = seq(0.00, 1.00, 0.25)) + scale_y_continuous(breaks = seq(0.00, 1.00, 0.25)) +
  scale_color_manual(values=viridis(32)) + theme_bw() + 
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold"))
p <- p + facet_wrap(vars(population), labeller = labeller(population = supp.labs)) +
  geom_text(data=cors, aes(label=paste("r=", cor, sep="")), x=0.2, y=0.9, inherit.aes = FALSE) +
  theme(strip.text.x = element_text(size = 10, face="bold")) +
  theme(legend.position="none")
p
# ggsave(filename=paste("markdown_images/allele_freq_correlation.eps"), plot=p, dpi=300,  width=7, height=6, units="in")   

p <- ggplot(total_allele_freqs, aes(x=tissue_freq, y=edna_freq, color=locus)) +
  geom_point() + ylab("Field eDNA allele frequency\n") + xlab("\nTissue allele frequency") + 
  scale_x_continuous(breaks = seq(0.00, 1.00, 0.25)) + scale_y_continuous(breaks = seq(0.00, 1.00, 0.25)) +
  scale_color_manual(values=viridis(32)) + theme_bw() + 
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold"))
p <- p + facet_wrap(vars(population), labeller = labeller(population = supp.labs)) +
  geom_text(data=cors, aes(label=paste("r=", cor, sep="")), x=0.2, y=0.9, inherit.aes = FALSE) +
  theme(strip.text.x = element_text(size = 10, face="bold")) +
  theme(legend.position="none")
p

#############################################################################################
######################## PCA of population read frequencies ###########################
#############################################################################################
# PCA
data_long <- data.frame(locus_allele=paste(total_allele_freqs$locus,"_", total_allele_freqs$allele,sep = ""),
                        total_allele_freqs[,(3:4)])
data_long <- filter(data_long, grepl(paste(sites,collapse="|"),population))
data_wide <- data_long %>% 
  pivot_wider(names_from = c(population), values_from=edna_freq)
PCA_data <- as.data.frame(t(subset(data_wide, select = -locus_allele)))
PCA_data[is.na(PCA_data)] <- 0
colnames(PCA_data) <- data_wide$locus_allele # name the columns: each allele is a column 
my.pca <- prcomp(PCA_data, center = TRUE)
pca_summary <- summary(my.pca)
PC1_and_PC2 <- data.frame(PC1=my.pca$x[,1], PC2= my.pca$x[,2], site=rownames(my.pca$x))
PC1_and_PC2$site <- factor(PC1_and_PC2$site,levels=PC1_and_PC2$site)
my_plot <- ggplot(PC1_and_PC2, aes(x=PC1, y=PC2, col=site)) +
  scale_color_manual(name = "Population", values = cols, labels=supp.labs) +
  theme_bw() +
  labs(x = paste("PC1 (",round(pca_summary$importance[,1][2]*100,1),"%)", sep = ""), 
       y = paste("PC2 (",round(pca_summary$importance[,2][2]*100,1),"%)", sep = "")) +
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=20)) +
  geom_point(size = 6, stroke = 2)
my_plot
# ggsave(filename = paste("markdown_images/field_samples_PCA1PC2.eps"), plot=my_plot, dpi = 300,  width = 7, height = 5, units = "in")   

PC1_and_PC2 <- data.frame(PC1=my.pca$x[,3], PC2= my.pca$x[,4], site=rownames(my.pca$x))
PC1_and_PC2$site <- factor(PC1_and_PC2$site,levels=PC1_and_PC2$site)
my_plot <- ggplot(PC1_and_PC2, aes(x=PC1, y=PC2, col=site)) +
  scale_color_manual(name = "Population", values = cols, labels=supp.labs) +
  theme_bw() +
  labs(x = paste("PC3 (",round(pca_summary$importance[,3][2]*100,1),"%)", sep = ""), 
       y = paste("PC4 (",round(pca_summary$importance[,4][2]*100,1),"%)", sep = "")) +
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=20)) +
  geom_point(size = 6, stroke = 2)
my_plot
# ggsave(filename = paste("markdown_images/field_samples_PC3PC4.eps"), plot=my_plot, dpi = 300,  width = 7, height = 5, units = "in")   

#############################################################################################
############## Correlation between population allele/read frequencies ###############
#############################################################################################
cors <- ddply(total_allele_freqs_HWE, "population", summarise, cor = round(cor(tissue_freq, edna_freq, use="complete.obs"), 2))
# write.csv(cors, "/Users/kbja10/Github/eDNA_goby_popgen/allele_freq_cors.csv")
supp.labs <- c("Cayuga Lake","Cross Lake","Erie Canal","Lake Erie (east)","Lake Erie (west)",
               "Lake Huron","Lake Michigan (west)","Lake Michigan (east)","Lake St. Clair",
               "Oneida Lake","Onondaga Lake","Lake Ontario (east)","Lake Ontario (west)")
names(supp.labs) <- c("CAY","CRO","ECAN","ERIE","ERIW","LHS","LMK","LMM",
                      "LSC","ONE","ONO","OSW","ROC")
p <- ggplot(total_allele_freqs_HWE, aes(x=tissue_freq, y=edna_freq, color=locus)) +
  geom_point() + ylab("Field eDNA allele frequency") + xlab("Genotyped tissue allele frequency") + 
  scale_x_continuous(breaks = seq(0.00, 1.00, 0.25)) + scale_y_continuous(breaks = seq(0.00, 1.00, 0.25)) +
  scale_color_manual(values=viridis(29)) + theme_bw() 
p <- p + facet_wrap(vars(population), labeller = labeller(population = supp.labs)) +
  geom_text(data=cors, aes(label=paste("r=", cor, sep="")), x=0.2, y=0.9, inherit.aes = FALSE)
p
# ggsave(filename = paste("/Users/kbja10/Downloads/field_samples_correlation.eps"), plot=p, dpi = 300,  width = 8, height = 6, units = "in")   

#############################################################################################
############################# Subset to just sites with good read recovery #######################################
#############################################################################################
# correlation
sites <- c("CAY","ERIE","LHS","LMM","LSC","ONE","ONO","OSW")
subset_total_allele_freqs_HWE <- filter(total_allele_freqs_HWE, grepl(paste(sites,collapse="|"),population))
cors <- ddply(subset_total_allele_freqs_HWE, "population", summarise, cor = round(cor(tissue_freq, edna_freq, use="complete.obs"), 2))
p <- ggplot(subset_total_allele_freqs_HWE, aes(x=tissue_freq, y=edna_freq, color=locus)) +
  geom_point() + ylab("Field eDNA allele frequency") + xlab("Genotyped tissue allele frequency") + 
  scale_x_continuous(breaks = seq(0.00, 1.00, 0.25)) + scale_y_continuous(breaks = seq(0.00, 1.00, 0.25)) +
  scale_color_manual(values=viridis(29)) + theme_bw() 
p <- p + facet_wrap(vars(population), labeller = labeller(population = supp.labs)) +
  geom_text(data=cors, aes(label=paste("r=", cor, sep="")), x=0.15, y=0.9, inherit.aes = FALSE)
# ggsave(filename = paste("/Users/kbja10/Downloads/field_samples_correlation.eps"), plot=p, dpi = 300,  width = 8, height = 6, units = "in")   

# PCA
data_long <- data.frame(locus_allele=paste(subset_total_allele_freqs_HWE$locus,"_", subset_total_allele_freqs_HWE$allele,sep = ""),
                        population=subset_total_allele_freqs_HWE$population,
                        edna_freq=subset_total_allele_freqs_HWE$edna_freq,
                        tissue_freq=subset_total_allele_freqs_HWE$tissue_freq)
data_wide <- data_long %>% 
  pivot_wider(names_from = c(population), values_from = c(edna_freq,tissue_freq))
PCA_data <- as.data.frame(t(subset(data_wide, select = -locus_allele)))
PCA_data[is.na(PCA_data)] <- 0
colnames(PCA_data) <- data_wide$locus_allele # name the columns: each allele is a column 
my.pca <- prcomp(PCA_data, center = TRUE)
pca_summary <- summary(my.pca)
PC1_and_PC2 <- data.frame(PC1=my.pca$x[,1], PC2= my.pca$x[,2], type = rownames(my.pca$x), 
                          site = rep(sites, 2), sample = rep(c("eDNA","tissue"), each=8))
my_palette <- viridis(8)
my_plot <- ggplot(PC1_and_PC2, aes(x=PC1, y=PC2, col=factor(site), shape = factor(sample))) +
  scale_color_manual(name = "site", values = my_palette) +
  scale_shape_manual(name = "sample", values=c(24,21)) +
  theme_bw() +
  labs(x = paste("PC1 (",round(pca_summary$importance[,1][2]*100,1),"%)", sep = ""), 
       y = paste("PC2 (",round(pca_summary$importance[,2][2]*100,1),"%)", sep = "")) +
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=20)) +
  geom_point(size = 5, stroke = 2)
my_plot
# ggsave(filename = paste("/Users/kbja10/Downloads/field_samples_PCA.eps"), plot=my_plot, dpi = 300,  width = 7, height = 5, units = "in")   



# PCA
data_long <- data.frame(locus_allele=paste(total_allele_freqs_HWE$locus,"_", total_allele_freqs_HWE$allele,sep = ""),
                        population=total_allele_freqs_HWE$population,
                        edna_freq=total_allele_freqs_HWE$edna_freq)
data_wide <- data_long %>% 
  pivot_wider(names_from = c(population), values_from = c(edna_freq))
PCA_data <- as.data.frame(t(subset(data_wide, select = -locus_allele)))
PCA_data[is.na(PCA_data)] <- 0
colnames(PCA_data) <- data_wide$locus_allele # name the columns: each allele is a column 
my.pca <- prcomp(PCA_data)
pca_summary <- summary(my.pca)
PC1_and_PC2 <- data.frame(PC1=my.pca$x[,1], PC2= my.pca$x[,2], site = rownames(my.pca$x))
my_palette <- viridis(13)
my_plot <- ggplot(PC1_and_PC2, aes(x=PC1, y=PC2, col=factor(site))) +
  scale_color_manual(name = "site", values = my_palette) +
  theme_bw() +
  labs(x = paste("PC1 (",round(pca_summary$importance[,1][2]*100,1),"%)", sep = ""), 
       y = paste("PC2 (",round(pca_summary$importance[,2][2]*100,1),"%)", sep = "")) +
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=20)) +
  geom_point(size = 5, stroke = 2)
my_plot

