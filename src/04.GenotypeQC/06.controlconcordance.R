setwd("/scratch/richards/tomoko.nakanishi/09.COVID19/05.BQC/01.genotypeQC/v5.0/")

library(data.table)
png(file="/project/richards/tomoko.nakanishi/repo/BQC19_genotype_pipeline/results/04.genotypeQC/control.fig1.png", width=1200, height=600, res=150)
par(mfrow=c(1,2))
c1 <- fread("CEPH1463.frqx")
c1$d <- apply(c1[,5:7], 1, function(x){1- max(x, na.rm=T)/sum(x, na.rm=T)})
write.table(c1$SNP[c1$d > 0.05 & !is.na(c1$d)], file="CEPH1463.fail.SNP", quote = F, row.names = F, col.names = F)
c1 <- c1[!is.na(c1$d),]
mydata_hist <- hist(c1$d, breaks=50, plot=FALSE)
plot(mydata_hist$mids, mydata_hist$count, log="y", type='h', yaxt = "n", xlab="Genotype discordance", ylab="Frequency", main="CEPH1463-02")
yvalues <- c(0,10,10^2,10^3,10^4,10^5,10^6)
ystrings <- c(expression(0), expression(10^1), expression(10^2), expression(10^3), expression(10^4), expression(10^5),
              expression(10^6))
axis(side=2,　# X軸に挿入。Y軸はside=2 
     at=yvalues, # 上記で作成したxvaluesを利用。
     labels=ystrings)
abline(v=0.05, col="red")

c1 <- fread("NA24385.frqx")
c1$d <- apply(c1[,5:7], 1, function(x){1- max(x, na.rm=T)/sum(x, na.rm=T)})
write.table(c1$SNP[c1$d > 0.05 & !is.na(c1$d)], file="NA24385.fail.SNP", quote = F, row.names = F, col.names = F)
c1 <- c1[!is.na(c1$d),]

mydata_hist <- hist(c1$d, breaks=50, plot=FALSE)
plot(mydata_hist$mids, mydata_hist$count, log="y", type='h', yaxt = "n", xlab="Genotype discordance", ylab="Frequency", main="NA24385")
yvalues <- c(0,10,10^2,10^3,10^4,10^5,10^6)
ystrings <- c(expression(0), expression(10^1), expression(10^2), expression(10^3), expression(10^4), expression(10^5),
              expression(10^6))
axis(side=2,　# X軸に挿入。Y軸はside=2 
     at=yvalues, # 上記で作成したxvaluesを利用。
     labels=ystrings)
abline(v=0.05, col="red")
dev.off()
