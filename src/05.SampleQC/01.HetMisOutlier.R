setwd("/scratch/richards/tomoko.nakanishi/09.COVID19/05.BQC/01.genotypeQC/v5.0/")
library(data.table)
f <- fread("het-vs-imiss.txt")
s <- read.table("all.sample.ancestry", sep = " ", header=T)
s$UMAP_POP <- as.character(s$UMAP_POP)
#s <- s %>% mutate(final_POP = ifelse(confident == 1, ADMIXTURE_POP, "oth"))
fs <- merge(f, s, by.x=c("V1"), by.y=c("FID"))
png(file="/project/richards/tomoko.nakanishi/repo/BQC19_genotype_pipeline/results/05.SampleQC/01_misshet.png", width=800, height=600, res=150)
#fs1 <- fs[fs$clusterUMAP == 0,]
plot(log10(fs$V4), fs$V3, col=as.factor(fs$final_POP),
     xaxt = "n", cex=0.6, pch=16,
     xlab="Proportion of missing genotype", ylab="Heterozygosity Rate")
xvalues <- c(-3,-2.5,-2,-1.5)
xstrings <- c(10^(-3),0.003,10^(-2),0.03)
axis(side=1,　# X軸に挿入。Y軸はside=2 
     at=xvalues, # 上記で作成したxvaluesを利用。
     labels=xstrings)

mean(fs$V3[fs$final_POP == "nfe"])
mean(fs$V3[fs$final_POP == "nfe"]) + 5*sd(fs$V3[fs$final_POP == "nfe"])
abline(v=log10(0.04), lty="dashed")
#abline(h=mean(fs$V3.x[fs$final_POP  == "AFR"]) + 3*sd(fs$V3.x[fs$final_POP  == "AFR"]), lty="dashed")
#abline(h=mean(fs$V3.x[fs$final_POP  == "EAS"]) - 3*sd(fs$V3.x[fs$final_POP  == "EAS"]), lty="dashed")
legend("bottomright", legend = unique(fs$final_POP), pch=16, cex=0.7,
       col = as.factor(unique(fs$final_POP)))

dev.off()

##correct herterozygousity by ancestry
LM <- glm(V3 ~ PC1 + PC2 + PC3 + PC4, data=fs, family="gaussian")
fs$heterozygousity_corrected <- residuals(LM) 
png(file="/project/richards/tomoko.nakanishi/repo/BQC19_genotype_pipeline/results/05.SampleQC/02_misshet_corrected.png", width=800, height=600, res=150)
plot(log10(fs$V4), fs$heterozygousity_corrected, col=as.factor(fs$final_POP),
     xaxt = "n", cex=0.6, pch=16,
     xlab="Proportion of missing genotype", ylab="Corrected Heterozygosity Rate")
xvalues <- c(-3,-2.5,-2,-1.5)
xstrings <- c(10^(-3),0.003,10^(-2),0.03)
axis(side=1,　# X軸に挿入。Y軸はside=2 
     at=xvalues, # 上記で作成したxvaluesを利用。
     labels=xstrings)
mean(fs$heterozygousity_corrected)
abline(v=log10(0.04), lty="dashed")
abline(h=mean(fs$heterozygousity_corrected) + 5*sd(fs$heterozygousity_corrected), lty="dashed")
abline(h=mean(fs$heterozygousity_corrected) - 5*sd(fs$heterozygousity_corrected), lty="dashed")
legend("bottomright", legend = unique(fs$final_POP), pch=16, cex=0.7,
       col = as.factor(unique(fs$final_POP)))

dev.off()

roh <- fread("bqc19-v5.0-qc4.hom.indiv") %>% select(c("FID","KB"))
fs <- fs %>% merge(roh, by.x="V1", by.y="FID")
png(file="/project/richards/tomoko.nakanishi/repo/BQC19_genotype_pipeline/results/05.SampleQC/03_het_roh.png", width=800, height=600, res=150)
plot(log10(fs$KB), fs$heterozygousity_corrected, col=as.factor(fs$final_POP),
     xaxt = "n", cex=0.6, pch=16,
     xlab="total length of ROH segments (kb)", ylab="Corrected Heterozygosity Rate")
xvalues <- c(3.5,4,4.5,5,5.5)
xstrings <- c("10^(3.5)", "10^4","10^(4.5)", "10^5", "10^(5.5)")
axis(side=1,　# X軸に挿入。Y軸はside=2 
     at=xvalues, # 上記で作成したxvaluesを利用。
     labels=xstrings)
dev.off()

fs %>% filter(heterozygousity_corrected > mean(fs$heterozygousity_corrected) + 3*sd(fs$heterozygousity_corrected))

fs$remove <- 0
fs$remove[fs$V4 > 0.04] <- 1
fs$remove[fs$heterozygousity_corrected > mean(fs$heterozygousity_corrected) + 3*sd(fs$heterozygousity_corrected) & fs$final_POP == "eur"] <- 1
fs1 <- fs[fs$remove == 1,]
write.table(fs1[,c(1,2)], file="hetmis.outlier", sep = " ", row.names = F, col.names = F, quote = F)

fs1 <- fs[fs$remove == 0,]
write.table(fs, file="Ethnicity", sep = " ", row.names = F, col.names = F, quote = F)

