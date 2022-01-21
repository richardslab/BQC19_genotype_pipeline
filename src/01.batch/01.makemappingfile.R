setwd("/scratch/richards/tomoko.nakanishi/09.COVID19/05.BQC/01.genotypeQC/v5.0/")
library(data.table)
library(dplyr)
#~/projects/richards/restricted/bqc19/data/array-genotype-data/BQC_AX001_AX032R_BPW/

s <- fread("/project/richards/restricted/bqc19/data/array-genotype-data/BQC_AX001_AX032R_BPW/EXPORT/SampleReport.txt")
#b <- fread("../../data/v2.1/01.batch/previous.batchid", header = F)
#s$batch <- ifelse(s$`affymetrix-plate-barcode` %in% b$V1, 1, 0)
#s$`affymetrix-plate-barcode` <- as.factor(s$`affymetrix-plate-barcode`)

ids2 <- read.xlsx("/project/richards/restricted/bqc19/data/array-genotype-data/BQC_AX001_AX032R_BPW/REFERENCE List ALL Samples 2021-11-18.xlsx", sheet = "ALL Sites - Genotyping Infos")
ids2 <- ids2 %>% select(Axiom.Sample.ID, Individual.id, Sex)
ids2 <- unique(ids2)
colnames(ids2) <- c("SampleID",  "anonymized_patient_id", "sex")

master <- read.xlsx("/project/richards/tomoko.nakanishi/09.COVID19/data/05.BQC/01.genotype/v5.0/01.batch/Master_list_WHOLE_BLOOD_CORRECTED_10112021.xlsx")
master <- master %>% select(Sample.name, Individual.id, Collection.site, Sex)
master1 <- unique(master)
master <- read.xlsx("/project/richards/tomoko.nakanishi/09.COVID19/data/05.BQC/01.genotype/v5.0/01.batch/Master_list_WHOLE_BLOOD_CORRECTED_10112021.xlsx", sheet = "concilie")
master <- master %>% select(Sample.name.VAC.code, Individual.id, Alias.BQCid, Collection.site, Sex.submit, Sex.BTRSR, Sex.REDCap)
master2 <- unique(master)

fam <- s
fam <- fam %>% select(`Sample Filename`)
fam <- fam %>% merge(ids2, by.x="Sample Filename", by.y="SampleID", all.x=T)
fam1 <- fam %>% filter(!is.na(sex))


# ids1 <- read.csv("/scratch/richards/tomoko.nakanishi/09.COVID19/05.BQC/01.genotypeQC/v3.2/BQC_AX001toAX028_Individuals.csv")
# ids1 <- ids1 %>% mutate(Individual.id = gsub("'","", X.Individual.id),
#                         sex = gsub("'","", X.Sex),
#                         ID = gsub("'","", X.Sample.name))
# ids1 <- ids1 %>% select(c("Individual.id","ID", "sex"))
# ids1 <- unique(ids1)

fam2 <- fam %>% filter(is.na(sex)) %>% filter(!grepl("NA24385", `Sample Filename`)) %>%
  filter(!grepl("Negative", `Sample Filename`)) %>%
  filter(!grepl("CEPH1463-02", `Sample Filename`)) %>%
  filter(!grepl("BLANK", `Sample Filename`)) 
  
fam2 <- fam2 %>% mutate(anonymized_patient_id = case_when(grepl("VAC14472", `Sample Filename`) ~ "BQC_HSCM7",
                                                  grepl("CRCHUM",`Sample Filename`) ~ paste0("BQC_",str_split(`Sample Filename`, pattern = "_", simplify = TRUE)[,1],"",str_split(`Sample Filename`, pattern = "_", simplify = TRUE)[,2]),
                                                  grepl("918",`Sample Filename`) ~ paste0("BQC_JGH0",str_split(`Sample Filename`, pattern = "_", simplify = TRUE)[,1]),
                                                  TRUE ~ paste0("BQC_JGH",str_split(`Sample Filename`, pattern = "_", simplify = TRUE)[,1])))
# fam2 <- fam2 %>% filter(!(Individual.id %in% fam1$Individual.id))
# fam2 <- fam2 %>% merge(ids1, by="Individual.id", all.x=T)

fam2 <- fam2 %>% merge(master, by.x="anonymized_patient_id", by.y="Individual.id", all.x=T)
fam2 <- fam2 %>% mutate(sex = Sex.submit)
fam2 <- fam2 %>% select(`Sample Filename`, anonymized_patient_id, sex)

fam <- bind_rows(fam1, fam2)
#colnames(FAM) <- c("SampleID", "anonymized_patient_id", "sex")

fam <- unique(fam)
write.table(fam, "/project/richards/tomoko.nakanishi/09.COVID19/data/05.BQC/01.genotype/v5.0/01.batch/sampleid.studyid.sex.map", quote=F, col.names = F, row.names = F, sep="\t")

s <- fread("/project/richards/restricted/bqc19/data/array-genotype-data/BQC_AX001_AX032R_BPW/EXPORT/SampleReport.txt")

#CRCHUM_160_2-924672.CEL
s1 <- s %>% merge(fam1, by="Sample Filename")
s2 <- s %>% merge(fam2, by="Sample Filename")
s <- bind_rows(s1, s2)
s$remove <- 0
s <- s %>% mutate(remove = case_when(`Pass/Fail` == "Fail" ~ 1,
                                     TRUE ~ remove))
s %>% filter(remove == 1)

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
s1 <- s %>% filter(remove == 0) %>% filter(!is.na(sex))

write.table(cbind(s1$`Sample Filename`,s1$`Sample Filename`), file="/project/richards/tomoko.nakanishi/09.COVID19/data/05.BQC/01.genotype/v5.0/01.batch/all.sample.removed.duplicated", quote = F, col.names = F, row.names = F)



