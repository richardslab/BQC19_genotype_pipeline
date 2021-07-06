setwd("/home/richards/tomoko.nakanishi/scratch/09.COVID19/05.BQC/01.genotypeQC/v4.0/")

library(data.table)
#library(ggplot2)
pop <- fread("/project/richards/public/hgdp_tgp/gnomad.genomes.v3.1.hgdp_1kg_subset.sample_meta.tsv.gz")
pop <- pop %>% select(population_inference.pop,labeled_subpop)
pcs <- fread("hgdp_1kg.projections", header=TRUE)

pcs <- bind_cols(pcs, pop)
setkey(pcs, IID)
proj <- fread('bqc19-v4.0-qc4.projections', header=TRUE)
setkey(proj, IID)
proj$population_inference.pop <- "0BQC"
proj$labeled_subpop <- "0BQC"

dat <- rbind(proj, pcs)
# K-means clustering
#install.packages("umap")
set.seed(1234)
library(umap)
dat.umap <- umap(dat[,3:12])
dat$UMAP1 <- dat.umap$layout[,1]
dat$UMAP2 <- dat.umap$layout[,2]
library("dbscan")
x <- rep(0,96)
for(i in seq(5,200)){
  cl <- hdbscan(dat[,.(UMAP1, UMAP2)], minPts = i)
  dat[, clusterUMAP := cl$cluster]
  x[i-4] <- sum(dat$clusterUMAP == 0 & dat$population_inference.pop == "0BQC") 
}
which(x == min(x))
cl <- hdbscan(dat[,.(UMAP1, UMAP2)], minPts = 90)
colors <- mapply(function(col, i) adjustcolor(col, alpha.f = cl$membership_prob[i]), 
                 palette()[cl$cluster+1], seq_along(cl$cluster))
dat[, clusterUMAP := cl$cluster]
dat[, table(clusterUMAP, population_inference.pop)]
dat[, table(clusterUMAP, labeled_subpop)]

tmp <- dat %>% filter(clusterUMAP == 15)
table(tmp[,c("population_inference.pop", "labeled_subpop")])

write.table(dat, file="sample.all", quote = F, sep=" ", row.names = F, col.names = T)


batch <- fread("/project/richards/restricted/bqc19/data/array-genotype-data/Export_BQC-AX001toAX028_JGH-CRCHUM/ALL DATA AX001 to AX031/SampleReport_AX001toAX031.txt")
batch <- batch[,c(1,7)]
colnames(batch) <- c("FID", "batch")
dat2 <- merge(dat, batch, by="FID")
dim(dat2[dat2$clusterUMAP %in% c(16:17),])/dim(dat2[dat2$population_inference.pop == "0BQC",])
#0.699579

dat2 <- dat2 %>% mutate(UMAP_POP = case_when(clusterUMAP %in% c(16:17) ~ "nfe",
                                             clusterUMAP %in% c(9:10,13,14) ~ "afr",
                                             clusterUMAP %in% c(2,3,4) ~ "eas",
                                             clusterUMAP %in% c(5) ~ "fin",
                                             clusterUMAP %in% c(15) ~ "mid",
                                             clusterUMAP %in% c(11,12) ~ "amr",
                                             clusterUMAP %in% c(6,7,8) ~ "sas",
                                             TRUE ~ "oth"
                                             ))

png(file="popstrat.fig1.png", width=1200, height=1200, res=150)
par(mfrow=c(2,2))

dat1 <- dat %>% filter(population_inference.pop != "0BQC")
p1 <- ggplot(dat1, aes(PC1, PC2, col=population_inference.pop)) + geom_point() + xlim(min(dat1$PC1)-0.01, max(dat1$PC1)+0.01) +
  ylim(min(dat1$PC2)-0.01, max(dat1$PC2)+0.01) + theme_minimal()
p1

p2 <- ggplot(dat2, aes(PC1, PC2, col=UMAP_POP)) + geom_point() + xlim(min(dat1$PC1)-0.01, max(dat1$PC1)+0.01) +
  ylim(min(dat1$PC2)-0.01, max(dat1$PC2)+0.01) + theme_minimal()
p2

p1dat1[, plot(PC1, PC2, pch=16, cex=0.8,col="#00000033", bg="#00000006", main="HGDP + 1000 Genomes (PCA plot)", xlim=c(-0.2,0.45), ylim=c(-0.3,0.25))]
dat1[population_inference.pop=="nfe", points(PC1, PC2, pch=16, cex=0.8,col="#ff0000",  xlim=c(-0.2,0.45), ylim=c(-0.3,0.25))]
legend("bottomright", legend = "nfe", pch=16, 
       col = c("#ff0000"))

dat1 <- dat[dat$population_inference.pop == "0BQC",]
dat1[, plot(PC1, PC2, pch=16, cex=0.8,col="#00000033", bg="#00000006", main="BQC19 - JGH (PCA plot)", xlim=c(-0.2,0.45), ylim=c(-0.3,0.25))]
dat1[clusterUMAP %in% c(8:10), points(PC1, PC2, pch=16, cex=0.8,col="#ff0000", xlim=c(-0.2,0.45), ylim=c(-0.3,0.25))]
legend("bottomright", legend = paste0(round(dim(dat1[dat1$clusterUMAP %in% c(8:10),])[1]/dim(dat1)[1], 3)*100, " %"), pch=16, 
       col = c("#ff0000"))

dat1 <- dat[dat$population_inference.pop != "BQC",]
dat1[, plot(UMAP1, UMAP2, pch=16, cex=0.8,col="#00000033", bg="#00000006", main="1000 Genomes (UMAP plot)", xlim=c(-20,25), ylim=c(-10,25))]
dat1[SPOP=="EUR", points(UMAP1, UMAP2, pch=16, cex=0.8,col="#ff0000", xlim=c(-20,25), ylim=c(-10,25))]
legend("topright", legend = "EUR", pch=16, 
       col = c("#ff0000"))

dat1 <- dat[dat$SPOP == "BQC",]
dat1[, plot(UMAP1, UMAP2, pch=16, cex=0.8,col="#00000033", bg="#00000006",main="BQC19 - JGH (UMAP plot)", xlim=c(-20,25), ylim=c(-10,25))]
dat1[clusterUMAP %in% c(2,18,19), points(UMAP1, UMAP2, pch=16, cex=0.8,col="#ff0000", xlim=c(-20,25), ylim=c(-10,25))]
legend("topright", legend = paste0(round(dim(dat1[dat1$clusterUMAP %in% c(2,18,19),])[1]/dim(dat1)[1], 3)*100, " %"), pch=16, 
       col = c("#ff0000"))

dev.off()

write.table(dat2[dat2$clusterUMAP %in% c(2,18,19), c(1:2,18)], file="../../scratch/v3.1/bqc19-v3.1-qc4.EUR.IID_IID.batch", quote=F, col.names = F, row.names = F)
