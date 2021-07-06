setwd("/home/richards/tomoko.nakanishi/scratch/09.COVID19/05.BQC/01.genotypeQC/v4.0/")
library(data.table)
w <- fread("~/scratch/09.COVID19/05.BQC/BQC_phenotype/JGH/20210615/withdrawal.list", header=F)
w$id <- paste0("JGH_",w$V1)
f <- fread("bqc19-v4.0-qc1.fam")
map <- fread("/home/richards/tomoko.nakanishi/09.COVID19/src/05.BQC/01.genotypeQC/data/v4.0/01.batch/sampleid.studyid.sex.map", header=F)
with_list <- map$V1[map$V2 %in% w$id]
f$withdrawal <- 0
f$withdrawal[f$V1 %in% with_list] <- 1

write.table(f[f$withdrawal == 1,1:2], file="/home/richards/tomoko.nakanishi/09.COVID19/src/05.BQC/01.genotypeQC/data/v4.0/02.withdrawal/withdrawal.remove", row.names = F, col.names = F, quote = F, sep=" ", append=T)
