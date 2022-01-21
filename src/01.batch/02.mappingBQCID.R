setwd("/project/richards/tomoko.nakanishi/09.COVID19/data/05.BQC/01.genotype/v5.0/01.batch/")

master <- read.xlsx("Master_list_WHOLE_BLOOD_CORRECTED_10112021.xlsx", sheet = "concilie")
master <- master %>% select(Individual.id, Alias.BQCid) %>% unique()
map <- fread("sampleid.studyid.sex.map", header=F)

library(stringr)
colnames(map) <- c("GenotypeID", "Individual.id", "Sex")

#############NOT DONE BELOW

file <- merge(map, master, by="Individual.id", all.x=T)

map1 <- map %>% merge(master, by.x="VACcode", by.y="Sample.name.VAC.code", all.x=T)
map1 <- map1 %>% drop_na(VACcode) %>% drop_na(Alias.BQCid)
map1 <- unique(map1)

for(i in seq(1, dim(map1)[1])){
  map$Individual.id[map$V1 == map1$V1[i]] <- map1$Individual.id[i]
}

map <- map %>% mutate(Individual.id = case_when(!is.na(Individual.id) ~ Individual.id,
                                                grepl("CHUQ_", V2) ~ gsub("CHUQ_", "CHUQ", paste0("BQC_",V2)),
                                                grepl("CHUS_", V2) ~ gsub("CHUS_", "CHUS", paste0("BQC_",V2)),
                                                grepl("CHUSJ_", V2) ~ gsub("CHUSJ_", "CHUSJ", paste0("BQC_",V2)),
                                                grepl("CHUM_", V2) ~ gsub("CHUM_", "CRCHUM", paste0("BQC_",V2)),
                                                grepl("HSCM_", V2) ~ gsub("HSCM_", "HSCM", paste0("BQC_",V2)),
                                                grepl("JGH_", V2) ~ gsub("JGH_", "JGH", paste0("BQC_JGH",str_pad(tmp, 4, pad = "0"))),
                                                grepl("MUHC_MED", V2) ~ gsub("MUHC_", "MUHC", paste0(V2)),
                                                grepl("MUHC_DVA", V2) ~ gsub("MUHC_", "MUHC", paste0(V2)),
                                                grepl("MUHC_", V2) ~ gsub("MUHC_", "MUHC", paste0("BQC_MUHC",tmp)),
                                                grepl("SLSJ_", V2) ~ gsub("SLSJ_", "SLSJ", paste0("BQC_",V2))))


map1 <- map %>% merge(master, by="Individual.id", all.x=T)
map1 <- map1 %>% select(Individual.id, V1, V3, Sex.submit, Sex.BTRSR, Sex.REDCap, Alias.BQCid)
map1 <- map1 %>% rename(genotypeID = V1, 
                        genotypeSex = V3)
map1 <- unique(map1)
map1 <- map1 %>% filter(!(Individual.id == "BQC_HSCM49" & Alias.BQCid == "BQC10016"))

write.table(map1, file="mapping_BQCID_20211105.tsv", quote=F, sep="\t", row.names = F)
