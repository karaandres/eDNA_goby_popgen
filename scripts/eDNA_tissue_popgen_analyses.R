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
library(plyr)
library(dplyr)
library(tidyr)
library(ade4)
library(viridis)
library(pegas)
library(VennDiagram)
library(devtools)
library(pophelper)
library(gridExtra)
library(geodist)
library(ggpubr)

##########################################################################
########### PART 1: Tissue-based population genetic analysis #############
##########################################################################
hap_genotype_matrix <- read.csv("datasets/hap_genotype_matrix.csv")
hap_genotype_matrix$site_codes <- factor(hap_genotype_matrix$site_codes, levels=c("LMK","LMM","LHS","LSC","ERIW","ERIE",
                                                  "ROC","OSW","ECAN","CAY","CRO","ONO","ONE"))
goby_strata <- read.csv("datasets/goby_strata.csv")
goby_strata$site <- factor(goby_strata$site, levels=c("LMK","LMM","LHS","LSC","ERIW","ERIE",
                                    "ROC","OSW","ECAN","CAY","CRO","ONO","ONE"))
goby_strata <- goby_strata[order(goby_strata$site),]

# Turn tissue data as genind object 
geno.obj <- df2genind(hap_genotype_matrix[,-ncol(hap_genotype_matrix)], sep="\\|", # "|" must be preceded by double backslashes
                      ploidy=2, type="codom", pop=hap_genotype_matrix$site_codes,
                      loc.names=names(hap_genotype_matrix[,-ncol(hap_genotype_matrix)]), 
                      NA.char="NA|NA", strata=goby_strata)
geno.obj # 285 individuals; 27 loci; 389 alleles

# K-means clustering: looks to be 3 clusters
# geno.clusters <- find.clusters(geno.obj, n.pca = 120) # retained 120 Axes (>90% of variance)
# pdf(file = "figures/goby_find_clusters.pdf", width = 8, height = 8)
# plot(geno.clusters$Kstat,type="o", col="#408DBF", lwd=2, 
#       xlab="Number of clusters", ylab="BIC",cex.axis=1.5,cex.lab=1.5)
# dev.off()

# Look at PCA plot
cols <- rev(colorRampPalette(brewer.pal(name="Spectral", n = 11)[c(1:5,8:11)])(13))
x.geno.obj <- tab(geno.obj, freq=TRUE, NA.method="mean")
pca.geno.obj <- dudi.pca(x.geno.obj, center=TRUE, scale=FALSE, scannf=FALSE, nf=3)
s.class(pca.geno.obj$li, fac=pop(geno.obj),col=cols,axesel=FALSE, cstar=0, cpoint=3)
add.scatter.eig(pca.geno.obj$eig[1:50],3,1,2, ratio=.2)
# loadingplot(pca.geno.obj$c1^2)

# Look at PCOA plot 
# pca.geno.obj <- dudi.pco(d = dist(x.geno.obj), nf = 120, scannf = FALSE)
# s.class(pca.geno.obj$li, fac=pop(geno.obj),col=transp(funky(15),.6),axesel=FALSE, cstar=0, cpoint=3)
# add.scatter.eig(pca.geno.obj$eig[1:50],3,1,2, ratio=.2)

# Look at DAPC plot
set.seed(867)
dapc1 <- dapc(geno.obj, pop=hap_genotype_matrix$site_codes, 
              n.pca=120, n.da=12)
supp.labs <- c("Lake Michigan west","Lake Michigan east","Lake Huron","Lake St. Clair",
               "Lake Erie west","Lake Erie east","Lake Ontario west","Lake Ontario east",
               "Erie Canal","Cayuga Lake","Cross Lake","Onondaga Lake", "Oneida Lake")
names(supp.labs) <- unique(geno.obj$pop)
cols <- rev(colorRampPalette(brewer.pal(name="Spectral", n = 11)[c(1:5,8:11)])(13))
# pdf(file = "figures/goby_popgen_dapc.pdf", width=8, height=6)
scatter(dapc1, bg="white", pch=20,  col=cols, 
        cex=3, clab=0, leg=TRUE, txt.leg=supp.labs, posi.leg="bottomleft")
# dev.off()

# STRUCTURE plot: Delta K
sfiles <- list.files(path="structure/output_4Jul2022/", full.names=T)
slist <- readQ(files=sfiles,filetype="structure")
tr1 <- tabulateQ(qlist=slist)
sr1 <- summariseQ(tr1)
p <- evannoMethodStructure(data=sr1)
p1 <- ggplot(p, aes(x=k, y=deltaK)) +
  geom_line(col="#408DBF") +
  geom_point(col="#408DBF", size=3) +
  scale_x_continuous(labels=2:12, breaks=2:12) +
  theme_bw() +
  theme(axis.text=element_text(size=15), axis.title=element_text(size=20))
p1
# ggsave(filename = paste("figures/deltaK.pdf"), plot=p1, dpi = 300,  width = 6, height = 6, units = "in")   

### AMOVA 
amva <- poppr.amova(geno.obj, hier=~structure/site, method ="ade4", within = FALSE)
amva_test <- randtest(amva, nrepet=1000)
plot(amva_test)
amva

### ALLELIC RICHNESS PER SITE
genpop.obj <- genind2genpop(geno.obj) # turn into genpop object
allele_freqs <- as.data.frame(t(makefreq(genpop.obj))) # calculate population allele frequencies
allele_freqs$locus <- gsub("\\..*", "", rownames(allele_freqs))
Nal_df <- data.frame()
Nal_anova <- data.frame()
for (site in 1:(ncol(allele_freqs)-1)){
  Nal_per_site=c()
  for (locus in unique(allele_freqs$locus)){
    temp <- allele_freqs[allele_freqs$locus==locus, site]
    Nal_per_site=c(Nal_per_site, length(temp[temp>0]))
  }
  Nal_df_temp <- data.frame(pop=names(allele_freqs)[site],Nal=mean(Nal_per_site), 
                            sd=sd(Nal_per_site))
  Nal_df <- rbind(Nal_df, Nal_df_temp)
  Nal_anova_temp <- data.frame(pop=names(allele_freqs)[site], Nal=Nal_per_site)
  Nal_anova <- rbind(Nal_anova, Nal_anova_temp) 
}
mean(Nal_df$Nal) # avg = 6.63 avg # alleles per locus per pop
sd(Nal_df$Nal) # sd = 0.81 avg # alleles per locus per pop
one.way <- aov(Nal ~ pop, data = Nal_anova) # anova of allelic richness per site
summary(one.way) 

### PRIVATE ALLELES PER SITE
PA <- rowSums(private_alleles(geno.obj, count.alleles=FALSE)) # poppr function
allele_freqs2 <- allele_freqs[,1:13]
colSums(allele_freqs2[apply(allele_freqs2, 1, function(x) length(x[x>0]))==1,] !=0) # my function

### HETEROZYGOSITY & Fis PER SITE
geno.obj_sep <- seppop(geno.obj)
mean.hobs <- do.call("c", lapply(geno.obj_sep, function(x) mean(summary(x)$Hobs))) 
mean.hexp <- do.call("c", lapply(geno.obj_sep, function(x) mean(summary(x)$Hexp))) 
basic_goby <- basic.stats(geno.obj)
apply(basic_goby$Fis, MARGIN = 2, FUN = mean, na.rm = TRUE) %>%
  round(digits = 3)

### EXPECTED HETEROZYGOSITY
# Calculate and plot expected heterozygosity within populations
HE <- function(x) 1-(sum(x^2))
He_df <- data.frame()
for (column in 1:(ncol(allele_freqs)-1)){
  He_per_locus=c()
  for (locus in unique(allele_freqs$locus)){
    temp <- allele_freqs[allele_freqs$locus==locus,column]
    He_per_locus=c(He_per_locus, HE(temp))
  }
  He_df_temp <- data.frame(pop=names(allele_freqs)[column],He=mean(He_per_locus))
  He_df <- rbind(He_df, He_df_temp)
}

He_df$pop <- factor(He_df$pop, levels=c("LMK","LMM","LHS","LSC","ERIW","ERIE","ROC","OSW","ECAN","CAY","CRO","ONO","ONE"))

### Gene diversity -- slightly different than heterozygosity
hs_df <- data.frame(pop=factor(unique(hap_genotype_matrix$site_codes)), Hs=Hs(geno.obj))

### Plot allelic richness per locus and per population
temp <- summary(geno.obj)
barplot(temp$loc.n.all, las=3,
        xlab="Population", ylab="Number of alleles")
mean(temp$loc.n.all) # mean = 14.4 alleles/locus
barplot(temp$pop.n.all, col=cols, las=3,
        xlab="Population", ylab="Number of alleles")
mean(temp$pop.n.all) # mean = 178.9 alleles/pop
sd(temp$pop.n.all) # sd = 21.8 alleles/pop
plot(temp$Hexp, temp$Hobs, pch=20, cex=3, xlim=c(0.3,1), ylim=c(0.3,1),
     xlab="Expected heterozygosity", ylab="Observed heterozygosity")
abline(0,1,lty=2)

# Test effect of residence time on genetic diversity metrics 
yrs_since <- data.frame()
for (i in unique(goby_strata$site)){
  df <- goby_strata[goby_strata$site==i,]
  yrs_since <- rbind(yrs_since, df[1,])
}
yrs_since$Nal <- Nal_df$Nal
yrs_since$PA <- PA 
yrs_since$He <- He_df$He
yrs_since$Na <- temp$pop.n.all
lmNa <- lm(Na~since_invasion, data=yrs_since) # Create the linear regression
summary(lmNa)
lmNal <- lm(Nal~since_invasion, data=yrs_since) # Create the linear regression
summary(lmNal)
lmPA <- lm(PA~since_invasion, data=yrs_since) #Create the linear regression
summary(lmPA)
lmHe <- lm(He~since_invasion, data=yrs_since) #Create the linear regression
summary(lmHe)
# plot normality
plot(predict(lmPA), residuals(lmPA))
qqnorm(residuals(lmPA)); qqline(residuals(lmPA))
hist(residuals(lmPA))

# Calculate pairwise Fst and plot heatmap
pairwise_fst_mat <- as.matrix(genet.dist(geno.obj, method="WC84"))
rownames(pairwise_fst_mat) <- levels(geno.obj$pop)
colnames(pairwise_fst_mat) <- levels(geno.obj$pop)
pairwise_fst_mat.tri <- pairwise_fst_mat
pairwise_fst_mat.tri[lower.tri(pairwise_fst_mat.tri, diag=TRUE)] <- NA # only keep lower triangle
pairwise_fst_long <- melt(pairwise_fst_mat.tri, na.rm =TRUE) # turn into long df

fst_heatmap <- ggplot(data = pairwise_fst_long, aes(Var2, Var1, fill = value))+ 
  geom_tile(color = "white") + 
  geom_text(aes(label = round(value, 2)), size=4, color="bisque4", fontface = "bold") +
  scale_fill_gradientn(colours = brewer.pal(9,"YlGnBu"), name="FST")  + 
  scale_x_discrete(labels= supp.labs) +
  scale_y_discrete(labels= supp.labs) +
  theme_bw() +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, size = 11, hjust = 1),
        axis.text.y = element_text(size = 12)) + coord_fixed()
fst_heatmap
# ggsave(filename = paste("figures/fst_heatmap_tissue.pdf"), plot=fst_heatmap, dpi=300,  width=8, height=8, units="in")   

# CI of pairwise Fst estimates
boot.ppfst(geno.obj,nboot=100,quant=c(0.025,0.975),diploid=TRUE)
min(pairwise_fst_mat)
max(pairwise_fst_mat)

# Get upper tri
pairwise_fst_mat.tri <- pairwise_fst_mat
pairwise_fst_mat.tri[upper.tri(pairwise_fst_mat.tri, diag=TRUE)] <- NA # only keep upper triangle
pairwise_fst_long <- melt(pairwise_fst_mat.tri, na.rm =TRUE) # turn into long df

# Mantel test (IBD)
Dgen <- genet.dist(geno.obj, method="WC84")
Dgeo <- as.dist(geodist(cbind(unique(goby_strata$longitude), unique(goby_strata$latitude)), measure="geodesic"))
ibd <- mantel.randtest(Dgen, Dgeo)
ibd
plot(ibd)
plot(Dgeo, Dgen)

# Plot neighbor joining tree from distance matrix
theTree <- pairwise_fst_mat %>%
  nj() %>%    # calculate neighbor-joining tree
  ladderize() # organize branches by clade
supp.labs_df <- as.data.frame(supp.labs)
theTree$tip.label <- supp.labs_df[,1][match(theTree$tip.label, rownames(supp.labs_df))]
plot(theTree, font=2, tip.col=cols)
add.scale.bar()

# Subset to just the sites with good recovery 
sites <- c("CAY","LHS","LMM","LSC","ONE","ONO","OSW")
col_names = grepl(paste(sites,collapse="|"),colnames(pairwise_fst_mat))
row_names = grepl(paste(sites,collapse="|"),rownames(pairwise_fst_mat))
pairwise_fst_mat_sub <- pairwise_fst_mat[row_names, col_names]
cols <- rev(colorRampPalette(brewer.pal(name="Spectral", n = 11)[c(1:5,8:11)])(13))
cols <- cols[c(2,3,4,8,10,12,13)]
theTree <- pairwise_fst_mat_sub %>%
  nj() %>%    # calculate neighbor-joining tree
  ladderize() # organize branches by clade
supp.labs_df <- as.data.frame(supp.labs)
theTree$tip.label <- supp.labs_df[,1][match(theTree$tip.label, rownames(supp.labs_df))]
plot(theTree, type="unr", tip.col=cols, font=2)
annot <- round(theTree$edge.length,2)
add.scale.bar()

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

pairwise_AFD_mat_sub <- pairwise_AFD_mat[row_names, col_names]
theTree <- pairwise_AFD_mat_sub %>%
  nj() %>%    # calculate neighbor-joining tree
  ladderize() # organize branches by clade
supp.labs_df <- as.data.frame(supp.labs)
theTree$tip.label <- supp.labs_df[,1][match(theTree$tip.label, rownames(supp.labs_df))]
# pdf(file = "figures/goby_tissues_tree.pdf", width = 5, height = 5)
plot(theTree, type="unr", tip.col=cols, font=2)
annot <- round(theTree$edge.length,2)
add.scale.bar()
# dev.off() # save the tree

# correlation between Fst and AFD pairwise distances
pairwise_df <- data.frame(pairwise_fst_long,pairwise_AFD_long)
names(pairwise_df)[c(3,6)] <- c("fst_dist", "afd_dist")
p <- ggplot(pairwise_df,aes(x=fst_dist,y=afd_dist)) + 
  geom_point(size=3, col="#408DBF") +
  xlab("Fst") + ylab("AFD") +
  xlim(0,0.2) + ylim(0,0.5) +
  theme_bw() +
  theme(text = element_text(size=20))
p
# ggsave(filename = paste("figures/AFD_Fst_correlation.eps"), plot=p, dpi = 300,  width = 6, height = 6, units = "in")   

r1 <- mantel.rtest(as.dist(pairwise_fst_mat),as.dist(pairwise_AFD_mat),nrepet = 999)
plot(r1)
r1
# highly correlated

##########################################################################
############# PART 2: eDNA-based population genetic analysis #############
##########################################################################

edna_allele_freqs <- read.csv("datasets/edna_allele_freqs.csv", header=TRUE)
colnames(edna_allele_freqs)[1] <- "locus_allele"
rownames(edna_allele_freqs) <- edna_allele_freqs[,1]
edna_allele_freqs[is.na(edna_allele_freqs)] <- 0
edna_allele_freqs <- edna_allele_freqs[rowSums(edna_allele_freqs[,4:16])>0,]
nrow(edna_allele_freqs)
ncol(edna_allele_freqs)

### TOTAL ALLELES PER SITE 
Na_eDNA <- apply(edna_allele_freqs[4:16], 2, function(x) length(x[x>0]))
Na_eDNA

### AVERAGE ALLELIC RICHNESS PER SITE
Nal_df_eDNA <- data.frame()
for (site in 4:(ncol(edna_allele_freqs))){
  Nal_per_site=c()
  for (locus in unique(edna_allele_freqs$locus)){
    temp <- edna_allele_freqs[edna_allele_freqs$locus==locus, site]
    Nal_per_site=c(Nal_per_site, length(temp[temp>0]))
  }
  Nal_df_temp <- data.frame(pop=names(edna_allele_freqs)[site],Nal=mean(Nal_per_site), 
                            sd=sd(Nal_per_site))
  Nal_df_eDNA <- rbind(Nal_df_eDNA, Nal_df_temp)
}
Nal_df_eDNA

# DF of summary stats
sites <- c("LMM","LHS","LSC","OSW","CAY","ONO","ONE")
summary_comparison <- data.frame(yrs_since, Na_eDNA, Nal_df_eDNA)
summary_comparison <- summary_comparison[summary_comparison$pop %in% sites,]

### Remove samples that don't meet threshold
sites <- c("locus_allele","locus","allele","LMM","LHS","LSC","OSW","CAY","ONO","ONE")
edna_allele_freqs <- edna_allele_freqs[, names(edna_allele_freqs) %in% sites]

### PRIVATE ALLELES PER SITE
allele_freqs2 <- edna_allele_freqs[,4:ncol(edna_allele_freqs)]
PA_eDNA <- colSums(allele_freqs2[apply(allele_freqs2, 1, function(x) length(x[x>0]))==1,] !=0) # my function
PA_eDNA

### EXPECTED HETEROZYGOSITY
# Calculate and plot expected heterozygosity within populations
HE <- function(x) 1-(sum(x^2))
He_df_eDNA <- data.frame()
for (column in 4:(ncol(edna_allele_freqs))){
  He_per_locus=c()
  for (locus in unique(edna_allele_freqs$locus)){
    temp <- edna_allele_freqs[edna_allele_freqs$locus==locus,column]
    He_per_locus=c(He_per_locus, HE(temp))
  }
  He_df_temp <- data.frame(pop=names(edna_allele_freqs)[column],He=mean(He_per_locus, na.rm=TRUE))
  He_df_eDNA <- rbind(He_df_eDNA, He_df_temp)
}
He_df_eDNA

### Tests/plots between groups
summary_comparison <- data.frame(summary_comparison, PA_eDNA, He_df_eDNA$He)
# Na
cor.test(summary_comparison$Na, summary_comparison$Na_eDNA)
t.test(summary_comparison$Na, summary_comparison$Na_eDNA, paired=TRUE)

# Nal
cor.test(summary_comparison$Nal, summary_comparison$Nal.1)
t.test(summary_comparison$Nal, summary_comparison$Nal.1, paired=TRUE)
p1 <- ggpaired(summary_comparison %>% pivot_longer(cols=c(Nal,Nal.1), names_to="sample", values_to="Nal"),
         x="sample", y="Nal", palette=c("#5E4FA2","#4075B4"), xlab="",
         color="sample", line.color="gray", line.size=0.4, ylab="Allelic richness") +
         theme(legend.position="none") +
         stat_compare_means(method="t.test", paired=TRUE)
# PA
cor.test(summary_comparison$PA, summary_comparison$PA_eDNA)
t.test(summary_comparison$PA, summary_comparison$PA_eDNA, paired=TRUE)
p2 <-ggpaired(summary_comparison %>% pivot_longer(cols=c(PA,PA_eDNA), names_to="sample", values_to="PA"),
         x="sample", y="PA", palette=c("#5E4FA2","#4075B4"), xlab="",
         color="sample", line.color="gray", line.size=0.4, ylab="Private alleles") +
         theme(legend.position="none") +
         stat_compare_means(method="t.test", paired=TRUE)
# He
cor.test(summary_comparison$He, summary_comparison$He_df_eDNA.He)
t.test(summary_comparison$He, summary_comparison$He_df_eDNA.He, paired=TRUE)
p3 <- ggpaired(summary_comparison %>% pivot_longer(cols=c(He,He_df_eDNA.He), names_to="sample", values_to="He"),
         x="sample", y="He", palette=c("#5E4FA2","#4075B4"), xlab="",
         color="sample", line.color="gray", line.size=0.4, ylab="Expected heterozygosity") +
         theme(legend.position="none") +
         stat_compare_means(method="t.test", paired=TRUE)
arrange_plots <- ggarrange(p1,p2,p3, labels=c("(A)","(B)","(C)"), nrow=1)
arrange_plots
# ggsave(filename=paste("figures/FigS6_diversity_stats_comparison.pdf"), plot=arrange_plots, dpi=300, width=10, height=4, units="in")

# eDNA private alleles as a function of time since invasion 
lmPA_eDNA <- lm(PA_eDNA~since_invasion, data=summary_comparison) #Create the linear regression
summary(lmPA_eDNA)

### Venn diagram of alleles in eDNA vs. tissues 
tissue_allele_freqs <- read.csv("datasets/tissue_allele_freqs.csv", header=TRUE)
cols <- rev(colorRampPalette(brewer.pal(name="Spectral", n = 11)[c(1:5,8:11)])(13))

venn.diagram(x = list(tissue_allele_freqs$X, edna_allele_freqs$locus_allele),
  category.names = c("Tissue alleles" , "eDNA alleles"),
  filename = "figures/Fig2_eDNA_alleles_venn_diagram.png",cex = 1.5, cat.cex = 1.5,
  lty=0,fill=cols[1:2], alpha=0.7, fontfamily="sans", 
  cat.pos = c(-20, 20), cat.dist=c(0.05, 0.05))

# Subset to only alleles found in round goby tissues
edna_allele_freqs <- edna_allele_freqs[edna_allele_freqs$locus_allele %in% colnames(geno.obj$tab), ]
nrow(edna_allele_freqs) # 260
nrow(edna_allele_freqs)/nrow(tissue_allele_freqs) # 66.8% of all known alleles
tissue_allele_freqs_0.1 <- tissue_allele_freqs[,-1]
rownames(tissue_allele_freqs_0.1) <- tissue_allele_freqs[,1]
tissue_allele_freqs_0.1[tissue_allele_freqs_0.1<0.1] <- 0
tissue_allele_freqs_0.1 <- tissue_allele_freqs_0.1[rowSums(tissue_allele_freqs_0.1)>0,]
edna_allele_freqs_0.1 <- edna_allele_freqs[edna_allele_freqs$locus_allele %in% rownames(tissue_allele_freqs_0.1),]
nrow(edna_allele_freqs_0.1)/nrow(tissue_allele_freqs_0.1) # and 83.7% of all alleles with frequencies greater than 10% at a locus

edna_freqs <- edna_allele_freqs[,-1]

# calculate pairwise AFD for all loci and populations 
pairwise_AFD_eDNA <- data.frame(Var1=combn(colnames(edna_freqs[,-c(1:2)]), 2)[1,],
                                Var2=combn(colnames(edna_freqs[,-c(1:2)]), 2)[2,],
                                AFD=combn(colnames(edna_freqs[,-c(1:2)]), 2, # all combinations of populations
                                          FUN = function(x) mean(afd(edna_freqs,x[1],x[2])$afd))) # calculate mean AFD across all loci
pairwise_AFD_eDNA_mat <- data.matrix(pivot_wider(pairwise_AFD_eDNA, names_from=Var1, values_from=AFD))
pairwise_AFD_eDNA_mat <- pairwise_AFD_eDNA_mat[,-1] # remove column with rownames
pairwise_AFD_eDNA_mat <- cbind(pairwise_AFD_eDNA_mat, rep(NA,nrow(pairwise_AFD_eDNA_mat)))
pairwise_AFD_eDNA_mat <- rbind(rep(NA,ncol(pairwise_AFD_eDNA_mat)), pairwise_AFD_eDNA_mat)

# turn into symmetrical matrix so we can reorder populations E-W
pairwise_AFD_eDNA_mat <- makeSymm(pairwise_AFD_eDNA_mat)
rownames(pairwise_AFD_eDNA_mat) <- colnames(edna_freqs[,-(1:2)])
colnames(pairwise_AFD_eDNA_mat) <- colnames(edna_freqs[,-(1:2)])
pairwise_AFD_eDNA_mat.tri <- pairwise_AFD_eDNA_mat
pairwise_AFD_eDNA_mat.tri[upper.tri(pairwise_AFD_eDNA_mat.tri, diag=TRUE)] <- NA # only keep lower triangle
pairwise_AFD_eDNA_long <- melt(pairwise_AFD_eDNA_mat.tri, na.rm =TRUE) # turn into long df

# correlation between tissue-based and AFD-based pairwise distances
pairwise_AFD_long <- filter(pairwise_AFD_long, grepl(paste(colnames(pairwise_AFD_eDNA_mat),collapse="|"),Var1))
pairwise_AFD_long <- filter(pairwise_AFD_long, grepl(paste(colnames(pairwise_AFD_eDNA_mat),collapse="|"),Var2))
pairwise_df <- data.frame(pairwise_AFD_long, pairwise_AFD_eDNA_long)
names(pairwise_df)[c(3,6)] <- c("afd_dist_tissue", "afd_dist_eDNA")
cols <- rev(colorRampPalette(brewer.pal(name="Spectral", n = 11)[c(1:5,8:11)])(13))
cols <- cols[c(2,3,4,8,10,12,13)]
p <- ggplot(pairwise_df,aes(x=afd_dist_tissue,y=afd_dist_eDNA,colour=Var2,shape=Var1)) + 
  geom_point(size=5, stroke=1.5) +
  xlab("Tissue-based AFD") + ylab("eDNA-based AFD") +
  #xlim(0,1) + ylim(0,1) +
  scale_color_manual(name="Pop 2", values=cols, labels=supp.labs) +
  scale_shape_manual(name="Pop 1", values=c(15:18,3,4,8), labels=supp.labs) +
  theme_bw() +
  theme(text = element_text(size=20))
p
# ggsave(filename = paste("figures/Fig4_eDNA_tissue_AFD_correlation.pdf"), plot=p, dpi = 300,  width = 9, height = 6, units = "in")   
cor.test(pairwise_df$afd_dist_tissue, pairwise_df$afd_dist_eDNA)
# r1 <- mantel.rtest(as.dist(pairwise_AFD_mat), as.dist(pairwise_AFD_eDNA_mat),nrepet = 999)
# plot(r1)
# r1

afd_heatmap <- ggplot(data = pairwise_AFD_eDNA_long, aes(Var2, Var1, fill = value))+ 
  geom_tile(color = "white") + 
  scale_fill_gradientn(colours = brewer.pal(9,"YlGnBu"), name="FST")  + 
  scale_x_discrete(labels= supp.labs) +
  scale_y_discrete(labels= supp.labs) +
  theme_bw() +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, size = 11, hjust = 1),
        axis.text.y = element_text(size = 12)) + coord_fixed()
afd_heatmap
# ggsave(filename = paste("figures/afd_heatmap_eDNA.pdf"), plot=afd_heatmap, dpi = 300,  width = 6, height = 6, units = "in")   

# Plot neighbor joining tree from distance matrix
# turn into symmetrical matrix so we can reorder populations E-W
col_names = grepl(paste(sites,collapse="|"),colnames(pairwise_AFD_eDNA_mat))
row_names = grepl(paste(sites,collapse="|"),rownames(pairwise_AFD_eDNA_mat))
pairwise_AFD_eDNA_mat_sub <- pairwise_AFD_eDNA_mat[row_names, col_names]
cols <- rev(colorRampPalette(brewer.pal(name="Spectral", n = 11)[c(1:5,8:11)])(13))
cols <- cols[c(2,3,4,8,10,12,13)]
theTree <- as.dist(pairwise_AFD_eDNA_mat_sub) %>%
  nj() %>%    # calculate neighbor-joining tree
  ladderize() # organize branches by clade
supp.labs_df <- as.data.frame(supp.labs)
theTree$tip.label <- supp.labs_df[,1][match(theTree$tip.label, rownames(supp.labs_df))]
# pdf(file = "figures/goby_edna_tree.pdf", width = 5, height = 5)
plot(theTree, type="unr", font=2, tip.col=cols)
add.scale.bar()
# dev.off() # save the tree

### AMOVA on AFD -- can't figure out bc no genetic distance among alleles
# amova_eDNA <- ade4::amova(edna_freqs[,-c(1:2)], as.dist(pairwise_AFD_eDNA_mat_sub), data.frame(structure=c("west","west","west","east","east","east","east")))
# edna_dist <- as.dist(pairwise_AFD_eDNA_mat_sub)
# edna_strata <- data.frame(structure=as.factor(c("west","west","west","east","east","east","east")),
#                           pop=as.factor(colnames(edna_freqs[,-c(1:2)])))
# pegas::amova(edna_dist ~ structure, data=edna_strata, nperm = 10)

# Mantel test (IBD)
pairwise_AFD_eDNA_mat <- data.matrix(pivot_wider(pairwise_AFD_eDNA, names_from=Var1, values_from=AFD))
pairwise_AFD_eDNA_mat <- pairwise_AFD_eDNA_mat[,-1] # remove column with rownames
pairwise_AFD_eDNA_mat <- cbind(pairwise_AFD_eDNA_mat, rep(NA,nrow(pairwise_AFD_eDNA_mat)))
pairwise_AFD_eDNA_mat <- rbind(rep(NA,ncol(pairwise_AFD_eDNA_mat)), pairwise_AFD_eDNA_mat)
rownames(pairwise_AFD_eDNA_mat) <- colnames(edna_freqs[,-(1:2)])
colnames(pairwise_AFD_eDNA_mat) <- colnames(edna_freqs[,-(1:2)])
Dgen_eDNA <-  as.dist(pairwise_AFD_eDNA_mat)
goby_strata_eDNA <- goby_strata[goby_strata$site %in% sites,]
Dgeo_eDNA <- as.dist(geodist(cbind(unique(goby_strata_eDNA$longitude), unique(goby_strata_eDNA$latitude)), measure="geodesic"))
ibd <- mantel.randtest(Dgen_eDNA, Dgeo_eDNA, nrepet = 999)
ibd
plot(ibd)
plot(Dgeo_eDNA, Dgen_eDNA)

#############################################################################################
############## Correlation between population allele/read frequencies ###############
#############################################################################################
# Overall allele frequencies
geno.obj$tab

# Combine eDNA and tissue allele freqs
tissue_freqs_long <- pivot_longer(allele_freqs,cols=LMK:ONE,names_to="population")
edna_freqs_long <- pivot_longer(edna_freqs,cols=LMM:ONE,names_to="population")
total_allele_freqs <- merge(edna_freqs_long, tissue_freqs_long, by = c("locus","allele","population"))
names(total_allele_freqs) <- c("locus","allele","population","edna_freq","tissue_freq")
total_allele_freqs$population <- factor(total_allele_freqs$population, 
                                      levels=c("LMK","LMM","LHS","LSC","ERIW","ERIE","ROC","OSW","ECAN","CAY","CRO","ONO","ONE"))
total_allele_freqs <- total_allele_freqs %>%
  arrange(total_allele_freqs$population) 

# Correlation coefficient between tissue and eDNA allele frequencies per site
cor.test(total_allele_freqs$edna_freq, total_allele_freqs$tissue_freq) # r = 0.69
cors <- total_allele_freqs %>%
  group_by(population) %>%
  dplyr::summarize(cor=round(cor(tissue_freq, edna_freq, use="complete.obs"), 2))
# write.csv(cors, "/Users/kbja10/Github/eDNA_goby_popgen/allele_freq_cors.csv")

p <- ggplot(total_allele_freqs, aes(x=tissue_freq, y=edna_freq, color=locus)) +
  geom_point() + ylab("eDNA allele frequency\n") + xlab("\nTissue allele frequency") + 
  scale_x_continuous(breaks = seq(0.00, 1.00, 0.25)) + scale_y_continuous(breaks = seq(0.00, 1.00, 0.25)) +
  scale_color_manual(values=viridis(27)) + theme_bw() + 
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
  theme(legend.position="none")
p <- p + facet_wrap(vars(population), labeller = labeller(population = supp.labs)) +
  geom_text(data=cors, aes(label=paste("r=", cor, sep="")), x=0.2, y=0.9, inherit.aes = FALSE) +
  theme(strip.text.x = element_text(size = 10, face="bold")) +
  theme(legend.position="none")
p
# ggsave(filename=paste("figures/FigS7_allele_freq_correlation.pdf"), plot=p, dpi=300,  width=7, height=6, units="in")   

#############################################################################################
######################## PCA of population read frequencies ###########################
#############################################################################################
# PCA of eDNA allele frequencies
data_long <- data.frame(locus_allele=paste(total_allele_freqs$locus,"_", total_allele_freqs$allele,sep = ""),
                        total_allele_freqs[,(3:4)])
data_long <- filter(data_long, grepl(paste(sites,collapse="|"),population))
data_wide <- data_long %>% 
  pivot_wider(names_from = c(population), values_from=edna_freq)
PCA_data <- as.data.frame(t(subset(data_wide, select = -locus_allele)))
PCA_data[is.na(PCA_data)] <- 0
colnames(PCA_data) <- data_wide$locus_allele # name the columns: each allele is a column 
PCA_data <- PCA_data[,colSums(PCA_data, na.rm = TRUE)>0]
my.pca <- prcomp(PCA_data, center=TRUE)
pca_summary <- summary(my.pca)
PC1_and_PC2 <- data.frame(PC1=my.pca$x[,1], PC2= my.pca$x[,2], site=rownames(my.pca$x))
PC1_and_PC2$site <- factor(PC1_and_PC2$site,levels=PC1_and_PC2$site)
cols <- rev(colorRampPalette(brewer.pal(name="Spectral", n = 11)[c(1:5,8:11)])(13))
cols <- cols[c(2,3,4,8,10,12,13)]
my_plot <- ggplot(PC1_and_PC2, aes(x=PC1, y=PC2, fill=site)) +
  theme_bw() +
  labs(x = paste("PC1 (",round(pca_summary$importance[,1][2]*100,1),"%)", sep = ""), 
       y = paste("PC2 (",round(pca_summary$importance[,2][2]*100,1),"%)", sep = "")) +
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=20)) +
  geom_point(shape=21, size=6, color="grey23") +
  scale_fill_manual(name = "Population", values=cols, labels=supp.labs)
my_plot

# PCA of tissue allele frequencies
data_long <- data.frame(locus_allele=paste(total_allele_freqs$locus,"_", total_allele_freqs$allele,sep = ""),
                        total_allele_freqs[,c(3,5)])
data_long <- filter(data_long, grepl(paste(sites,collapse="|"),population))
data_wide <- data_long %>% 
  pivot_wider(names_from = c(population), values_from=tissue_freq)
PCA_data <- as.data.frame(t(subset(data_wide, select = -locus_allele)))
PCA_data[is.na(PCA_data)] <- 0
colnames(PCA_data) <- data_wide$locus_allele # name the columns: each allele is a column 
PCA_data <- PCA_data[,colSums(PCA_data, na.rm = TRUE)>0]
my.pca <- prcomp(PCA_data, center=TRUE)
pca_summary <- summary(my.pca)
PC1_and_PC2 <- data.frame(PC1=my.pca$x[,1], PC2= my.pca$x[,2], site=rownames(my.pca$x))
PC1_and_PC2$site <- factor(PC1_and_PC2$site,levels=PC1_and_PC2$site)
cols <- rev(colorRampPalette(brewer.pal(name="Spectral", n = 11)[c(1:5,8:11)])(13))
cols <- cols[c(2,3,4,8,10,12,13)]
my_plot2 <- ggplot(PC1_and_PC2, aes(x=PC1, y=PC2, fill=site)) +
  theme_bw() +
  theme(legend.position = "none") +
  labs(x = paste("PC1 (",round(pca_summary$importance[,1][2]*100,1),"%)", sep = ""), 
       y = paste("PC2 (",round(pca_summary$importance[,2][2]*100,1),"%)", sep = "")) +
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=20)) +
  geom_point(shape=21, size=6, color="grey23") +
  scale_fill_manual(name = "Population", values=cols, labels=supp.labs)
my_plot2

arrange_plots <- ggarrange(my_plot2, my_plot, ncol=2, widths=c(3,4))
# ggsave(filename=paste("figures/PCA_tissues_eDNA.pdf"), plot=arrange_plots, dpi=300, width=9, height=4, units="in")   
