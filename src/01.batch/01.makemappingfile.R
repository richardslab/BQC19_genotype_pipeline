setwd("/home/richards/tomoko.nakanishi/scratch/09.COVID19/05.BQC/01.genotypeQC/v4.0/")
library(data.table)
library(dplyr)
s <- fread("/project/richards/restricted/bqc19/data/array-genotype-data/Export_BQC-AX001toAX028_JGH-CRCHUM/ALL DATA AX001 to AX031/SampleReport_AX001toAX031.txt")
#b <- fread("../../data/v2.1/01.batch/previous.batchid", header = F)
#s$batch <- ifelse(s$`affymetrix-plate-barcode` %in% b$V1, 1, 0)
#s$`affymetrix-plate-barcode` <- as.factor(s$`affymetrix-plate-barcode`)

ids2 <- fread("/project/richards/restricted/bqc19/data/array-genotype-data/Export_BQC-AX001toAX028_JGH-CRCHUM/ALL DATA AX001 to AX031/REFERENCE List ALL Samples 2021-06-28.txt")
ids2 <- ids2 %>% select(`Axiom Sample ID`, `Individual id`, Sex)

fam <- fread("/project/richards/restricted/bqc19/data/array-genotype-data/Export_BQC-AX001toAX028_JGH-CRCHUM/ALL DATA AX001 to AX031/SampleReport_AX001toAX031.txt")
fam <- fam %>% select(`Sample Filename`)
fam <- fam %>% merge(ids2, by.x="Sample Filename", by.y="Axiom Sample ID", all.x=T)

ids1 <- read.csv("/scratch/richards/tomoko.nakanishi/09.COVID19/05.BQC/01.genotypeQC/v3.2/BQC_AX001toAX028_Individuals.csv")
ids1 <- ids1 %>% mutate(Individual.id = gsub("'","", X.Individual.id),
                        sex = gsub("'","", X.Sex),
                        ID = gsub("'","", X.Sample.name))
ids1 <- ids1 %>% select(c("Individual.id","ID", "sex"))
ids1 <- unique(ids1)
fam <- fam %>% mutate(ID = case_when(grepl("CRCHUM",`Sample Filename`) ~ paste0(str_split(`Sample Filename`, pattern = "_", simplify = TRUE)[,1],"_",str_split(`Sample Filename`, pattern = "_", simplify = TRUE)[,2]),
                                       TRUE ~ str_split(`Sample Filename`, pattern = "_", simplify = TRUE)[,1]))

fam <- fam %>% merge(ids1, by="ID", all.x=T)
fam <- fam %>% filter(!grepl("NA24385", ID))
fam <- fam %>% filter(!grepl("Negative", ID))
fam <- fam %>% filter(!grepl("CEPH1463-02", ID))
fam <- fam %>% filter(!grepl("BLANK", ID))

dim(fam)
map <- fread("/home/richards/tomoko.nakanishi/09.COVID19/src/05.BQC/01.genotypeQC/data/v3.1/01.batch/sampleid.studyid.sex.map", header=F)

fam <- fam %>% merge(map, by.x="Sample Filename",by.y="V1", all.x=T)
dim(fam)

fam %>% filter(is.na(sex) & is.na(Sex))

fam <- fam %>% mutate(SampleID = `Sample Filename`,
                      anonymized_patient_id = case_when(!is.na(V2) ~ V2,
                                                        grepl("BQC_JGH", Individual.id) ~ paste0("JGH_",as.numeric(gsub("BQC_JGH","", Individual.id))),
                                                        grepl("BQC_CRCHUM", Individual.id) ~ paste0("CHUM_",gsub("BQC_CRCHUM","", Individual.id)),
                                                        grepl("BQC_MUHC", Individual.id) ~ paste0("MUHC_",gsub("BQC_MUHC","", Individual.id)),
                                                        grepl("BQC_CHUQ", Individual.id) ~ paste0("CHUQ_",gsub("BQC_CHUQ","", Individual.id)),
                                                        grepl("BQC_CHUS", Individual.id) ~ paste0("CHUS_",gsub("BQC_CHUS","", Individual.id)),
                                                        grepl("BQC_HSCM", Individual.id) ~ paste0("HSCM_",gsub("BQC_HSCM","", Individual.id)),
                                                        grepl("BQC_SLSJ", Individual.id) ~ paste0("SLSJ_",gsub("BQC_SLSJ","", Individual.id)),
                                                        grepl("BQC_CHUSJ", Individual.id) ~ paste0("CHUSJ_",gsub("BQC_CHUSJ","", Individual.id)),
                                                        grepl("BQC_JGH", `Individual id`) ~ paste0("JGH_",as.numeric(gsub("BQC_JGH","", `Individual id`))),
                                                        grepl("BQC_CRCHUM", `Individual id`) ~ paste0("CHUM_",gsub("BQC_CRCHUM","", `Individual id`)),
                                                        grepl("BQC_MUHC", `Individual id`) ~ paste0("MUHC_",gsub("BQC_MUHC","", `Individual id`)),
                                                        grepl("BQC_CHUQ", `Individual id`) ~ paste0("CHUQ_",gsub("BQC_CHUQ","", `Individual id`)),
                                                        grepl("BQC_CHUS", `Individual id`) ~ paste0("CHUS_",gsub("BQC_CHUS","", `Individual id`)),
                                                        grepl("BQC_HSCM", `Individual id`) ~ paste0("HSCM_",gsub("BQC_HSCM","", `Individual id`)),
                                                        grepl("BQC_SLSJ", `Individual id`) ~ paste0("SLSJ_",gsub("BQC_SLSJ","", `Individual id`)),
                                                        grepl("BQC_CHUSJ", `Individual id`) ~ paste0("CHUSJ_",gsub("BQC_CHUSJ","", `Individual id`))
                                                        ))

fam <- fam %>% mutate(sex = case_when(!is.na(sex) ~ sex,
                                      !is.na(Sex) ~ Sex,
                                      !is.na(V3) ~ V3
                                      ))

fam <- fam %>% select(SampleID, anonymized_patient_id, sex)
write.table(fam, "/home/richards/tomoko.nakanishi/09.COVID19/src/05.BQC/01.genotypeQC/data/v4.0/01.batch/sampleid.studyid.sex.map", quote=F, col.names = F, row.names = F, sep="\t")

s$remove <- 0
s <- s %>% mutate(remove = case_when(`Pass/Fail` == "Fail" ~ 1,
                                     TRUE ~ remove))
s <- s %>% merge(fam, by.x="Sample Filename", by.y="SampleID", all.x=T)

s <- s[order(s$`QC call_rate`, decreasing = T),]
s$remove[duplicated(s$anonymized_patient_id) & !is.na(s$anonymized_patient_id)] <- 1
s %>% filter(remove == 1)
s %>% filter(`affymetrix-plate-barcode` == "5507784388583022821343" & remove == 0 & duplicated(anonymized_patient_id))
s <- s %>% mutate(remove = case_when(grepl("CEPH", `Sample Filename`) ~ 1,
                                     grepl("NA24385", `Sample Filename`) ~ 1,
                                     grepl("BLANK", `Sample Filename`) ~ 1,
                                     grepl("Negative", `Sample Filename`) ~ 1,
                                     TRUE ~ remove
))
s1 <- s %>% filter(remove == 0)

write.table(s1[,c("Sample Filename","anonymized_patient_id","sex")], file="/home/richards/tomoko.nakanishi/09.COVID19/src/05.BQC/01.genotypeQC/data/v4.0/01.batch/sampleid.studyid.sex.map", quote = F, col.names = F, row.names = F)

write.table(cbind(s1$`Sample Filename`,s1$`Sample Filename`), file="/home/richards/tomoko.nakanishi/09.COVID19/src/05.BQC/01.genotypeQC/data/v4.0/01.batch/all.sample.removed.duplicated", quote = F, col.names = F, row.names = F)



