setwd("/scratch/richards/tomoko.nakanishi/09.COVID19/05.BQC/01.genotypeQC/v5.0/")

library(data.table)
f <- fread("bqc19-v5.0-qc4.sexcheck")
png("/project/richards/tomoko.nakanishi/repo/BQC19_genotype_pipeline/results/05.SampleQC/XChromosomeHomogozyosityEstimate.png", width=600, height = 400)
hist(f$F, breaks=100, xlab="X Chromosome Homogozyosity Estimate", main="All sample")
dev.off()
f1 <- fread("bqc19-v5.0-qc5.sexcheck")
f1 <- f1[,c(1,7)]
f <- f[,c(1:6)]
final <- merge(f,f1,by="FID")
final$mLOY <- 0
final$mLOY[final$YCOUNT == 0 & final$F > 0.8] <- 1
final$mLOY[final$YCOUNT > 0 & final$F > 0.8] <- 2

final[final$YCOUNT == 0 & final$F > 0.8,]

map <- fread("/project/richards/tomoko.nakanishi/09.COVID19/data/05.BQC/01.genotype/v5.0/01.batch/sampleid.studyid.sex.map", header=F)
final <- final %>% merge(map, by.x="FID", by.y="V1", all.x=T)

#discordance check with self-reported
final <- final %>% mutate(sex_discordant = case_when(SNPSEX == 2 & V3 == "F" ~ 0,
                                                     SNPSEX == 1 & V3 == "M" ~ 0,
                                                     TRUE ~ 1
))
##list from Zaman, 2022/01/07
final <- final %>% mutate(sex_discordant = ifelse(V2 %in% c("BQC_JGH0099", "BQC_JGH0701", "BQC_JGH1189", "BQC_JGH1273"), 0, sex_discordant))

tmp <- final %>% filter(STATUS == "OK" & SNPSEX == 2) 

final %>% filter(sex_discordant == 1) %>% #filter(grepl("CRCHUM", V2)) %>%
  rename(ID = V2) %>%
  mutate(genotypeSEX = ifelse(SNPSEX == 2, "F", "M"),
         selfreportedSEX = ifelse(V3 == "M", "M", "F")) %>%
  select("ID", "genotypeSEX", "selfreportedSEX", "F", "YCOUNT") %>% write.xlsx("discordantsex.xlsx")

final %>% filter(sex_discordant == 1)
final <- final %>% mutate(remove = ifelse(sex_discordant == 1, 1, 0))
#final <- final %>% filter(remove == 0)

write.table(final[final$remove == 1,1:2], file="sex.discordance.remove", row.names = F, col.names = F, quote = F, sep=" ")



