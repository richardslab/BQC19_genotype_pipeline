setwd("/home/richards/tomoko.nakanishi/09.COVID19/data/05.BQC/01.genotype/v4.0/09.CHUM")

library(tidyr)
library(data.table)
i <- 1
tmp <- fread(paste0("all_",i,".sscore")) %>%
    dplyr::select(`IID`, SCORE1_SUM)
colnames(tmp) <- c("IID",paste0("SCORE_",i))

for(i in c(2:22)){
    TMP <- fread(paste0("all_",i,".sscore")) %>%
      dplyr::select(`IID`, SCORE1_SUM)
    colnames(TMP) <- c("IID",paste0("SCORE_",i))
    tmp <- merge(tmp, TMP, by="IID")
  }
  tmp$prs <- apply(tmp[,c(-1)], 1, sum) 
  tmp <- tmp %>% dplyr::select(IID, prs)
  saveRDS(tmp, file=paste0("all_prs.rds"))
