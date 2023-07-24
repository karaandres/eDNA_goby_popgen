### Round goby population genetics: tissue vs. eDNA analyses in field samples 
# This script accomplishes the following:
# PART 1: mtDNA COI qPCR results,
#         nuDNA microsatellite qPCR results
#         mt:nu DNA ratio

rm(list = ls())
setwd("/Users/kbja10/Github/eDNA_goby_popgen/")

library(ggplot2)
library(ggpubr)
library(dplyr)
library(tidyr)
library(RColorBrewer)
library(ggpmisc)

# qPCR efficiency and R2
mean(c(86.691, 92.918, 99.429, 99.406)) # efficiency
sd(c(86.691, 92.918, 99.429, 99.406))
mean(c(0.996, 0.996, 0.996, 0.992)) # R2
sd(c(0.996, 0.996, 0.996, 0.992))

# Plot mt:nu DNA ratio
qpcr_dat <- read.csv("datasets/2L_COI_msat_qPCR_results.csv", header = TRUE)
qpcr_dat <- qpcr_dat[-grep("STD_", qpcr_dat$Sample),]
qpcr_dat <- qpcr_dat[-grep("BUF", qpcr_dat$Sample),]
qpcr_dat <- qpcr_dat[-grep("SEN", qpcr_dat$Sample),]
qpcr_dat <- qpcr_dat[-grep("NEG", qpcr_dat$Sample),]
qpcr_dat[is.na(qpcr_dat)] <- 0

# Transform from CN/µL (elution) to CN/L (sample) -- 
qpcr_dat$COI_Quantity <- qpcr_dat$COI_Quantity*50/2
qpcr_dat$msat_Quantity <- qpcr_dat$msat_Quantity*50/2

# Calculate mean of site replicates (3)
qpcr_dat_summary <- as.data.frame(qpcr_dat %>% 
                          group_by(Sample) %>% 
                          summarise(COI_mean=mean(COI_Quantity, na.rm=TRUE),
                                    msat_mean=mean(msat_Quantity, na.rm=TRUE)))
qpcr_dat_summary$ratio <- qpcr_dat_summary$COI_mean/qpcr_dat_summary$msat_mean
qpcr_dat_summary$site <- gsub(qpcr_dat_summary$Sample, pattern="_.*", replacement="")
qpcr_dat_summary$ratio[is.infinite(qpcr_dat_summary$ratio)] <- NA
site_levels <- factor(qpcr_dat_summary$site, levels=c("LMK","LMM","LHS","LSC","ERIW","ERIE", "ROC","OSW","ECAN","CAY","CRO","ONO","ONE"))
cols <- rev(colorRampPalette(brewer.pal(name="Spectral", n=11)[c(1:5,8:11)])(13))
supp.labs <- c("Cayuga Lake","Cross Lake","Erie Canal","Lake Erie E","Lake Erie W","Lake Huron","Lake Michigan W","Lake Michigan E","Lake St. Clair", "Oneida Lake","Onondaga Lake","Lake Ontario E","Lake Ontario W","Seneca Lake", "Buffalo")
names(supp.labs) <- c("CAY","CRO","ECAN","ERIE","ERIW","LHS","LMK","LMM","LSC","ONE","ONO","OSW","ROC","SEN","BUF")

# Proportion of positive detections 
length(qpcr_dat_summary$COI_mean[qpcr_dat_summary$COI_mean>0])/nrow(qpcr_dat_summary) # 100% mtDNA
length(qpcr_dat_summary$msat_mean[qpcr_dat_summary$msat_mean>0])/nrow(qpcr_dat_summary) # 92% msat

# Avg reads and ratio of mtDNA and nuDNA per sample 
mean(qpcr_dat_summary$COI_mean)
sd(qpcr_dat_summary$COI_mean)
mean(qpcr_dat_summary$msat_mean)
sd(qpcr_dat_summary$msat_mean)
mean(qpcr_dat_summary$ratio, na.rm=TRUE)
sd(qpcr_dat_summary$ratio, na.rm=TRUE)
min(qpcr_dat_summary$ratio, na.rm=TRUE)
max(qpcr_dat_summary$ratio, na.rm=TRUE)

# Plot mean copy number of mtDNA & nuDNA (log scale)
p1 <- qpcr_dat_summary %>% 
  pivot_longer(c(COI_mean,msat_mean)) %>%
  ggplot(aes(x=name, y=value, fill=name)) +
  geom_boxplot() +
  geom_point() +
  scale_fill_manual(values=cols, labels=c("COI", "Microsatellite"), name="") +
  scale_x_discrete(labels=c("COI_mean"="mtDNA", "msat_mean"="nuDNA")) + 
  xlab("DNA marker type") + ylab("eDNA concentrationn (CN/L)") + scale_y_continuous(trans='log10') +
  theme_bw() +
  theme(legend.position = "none") +
  theme(text = element_text(size=16))
p1

# ANOVA of nuDNA vs mtDNA
library(lme4)
library(lmerTest)
library(emmeans)
qpcr_dat_summary_long <- as.data.frame(pivot_longer(qpcr_dat_summary,
                                                    cols=COI_mean:msat_mean,
                                                    names_to="marker",
                                                    values_to="quantity"))
lm2 <- lmer(log(quantity+1) ~ marker + (1|Sample)+ (1|site), data=qpcr_dat_summary_long)
summary(lm2)
anova(lm2)
# normality tests
plot(predict(lm2), residuals(lm2))
qqnorm(residuals(lm2)); qqline(residuals(lm2))
hist(residuals(lm2))
summary(emmeans(lm2,  ~ marker), type= "response") # estimated means

# Correlation analysis and linear regression of mtDNA & nuDNA copy numbers
cor.test(qpcr_dat_summary$COI_mean, qpcr_dat_summary$msat_mean)
lm3 <- lm(log(qpcr_dat_summary$COI_mean+1)~log(qpcr_dat_summary$msat_mean+1))
summary(lm3)

p2 <- ggplot(data=qpcr_dat_summary, aes(x=COI_mean+1, y=msat_mean+1)) +
  stat_poly_line() +
  stat_poly_eq(aes(label = after_stat(eq.label))) +
  stat_poly_eq(label.y = 0.9) +
  geom_point() + scale_x_log10() + scale_y_log10() +
  xlab("mtDNA concentration (CN/L)") + ylab("nuDNA concentration (CN/L)") +
  theme_bw() +
  theme(text = element_text(size=16))
p2

p3 <- ggplot(qpcr_dat_summary, aes(x=site_levels, y=ratio, color=site_levels)) + 
  geom_boxplot() +
  geom_hline(yintercept=mean(qpcr_dat_summary$ratio, na.rm=TRUE), linetype="dashed", color = "gray") +
  scale_color_manual(labels=supp.labs, values=cols, name="") +
  xlab("Population") + ylab("mtDNA:nuDNA ratio \n(qPCR copy number)") +
  theme_bw() +
  theme(axis.text.x=element_blank()) +
  theme(text = element_text(size=16))
p3

arrange_plots <- ggarrange(ggarrange(p1, p2, ncol = 2, labels = c("b", "c")))
# ggsave(filename=paste("figures/Fig1_COI_msat_qPCR_results.pdf"), plot=arrange_plots, dpi=300, width=8, height=4, units="in")
# ggsave(filename=paste("figures/FigS1_COI_msat_qPCR_ratio.pdf"), plot=p3, dpi=300, width=8, height=4, units="in")

### Plot nuDNA copy number vs. read frequency & allelic richness
edna_read_counts <- read.csv("datasets/edna_read_counts_separate.csv")
total_read_counts <- data.frame(Sample=colnames(edna_read_counts[,-c(1:3)]),
                                Counts=colSums(edna_read_counts[,-c(1:3)]))
total_read_counts$Sample <-  gsub(total_read_counts$Sample, pattern="0", replacement="")
dat <- merge(total_read_counts,qpcr_dat_summary,by="Sample")
dat$site <- gsub(dat$Sample, pattern="_.*", replacement="")
site_levels <- factor(dat$site, levels=c("LMK","LMM","LHS","LSC","ERIW","ERIE", "ROC","OSW","ECAN","CAY","CRO","ONO","ONE"))

p1 <- ggplot(dat, aes(x=msat_mean, y=Counts, color=site_levels)) +
  geom_point(size=6) +
  scale_color_manual(labels=supp.labs, values=cols, name="") +
  xlab("qPCR nuDNA copies (CN/L)") + ylab("NGS read count") +
  theme_bw() +
  theme(text = element_text(size=20))
p1
# ggsave(filename=paste("figures/qpcr_vs_read_count.pdf"), plot=p1, dpi=300, width=8, height=5, unit="in")

# Allelic richness per locus
allelic_richness <- data.frame(Sample=names(edna_read_counts[,-c(1:3)]),
                               richness=apply(edna_read_counts[,-c(1:3)], 2, function(c) sum(c!=0, na.rm=TRUE)))
allelic_richness$Sample <-  gsub(allelic_richness$Sample, pattern="0", replacement="")
dat <- merge(dat,allelic_richness,by="Sample")
p2 <- ggplot(dat, aes(x=Counts, y=richness, color=site_levels)) +
  geom_point(size=6) +
  scale_color_manual(labels=supp.labs, values=cols, name="") +
  xlab("NGS read count") + ylab("Allelic richness") +
  theme_bw() +
  theme(text = element_text(size=20))
p2
# ggsave(filename=paste("figures/read_count_vs_allelic_richness.pdf"), plot=p2, dpi=300, width=8, height=5, unit="in")
