setwd("/home/richards/tomoko.nakanishi/09.COVID19/data/05.BQC/01.genotype/v4.0/")
library(data.table)
library(openxlsx)
ans <- fread("08.PC/EUR.pc")
#ans <- fread("03.Ancestry/all.sample.ancestry")
map <- fread("01.batch/sampleid.studyid.sex.map", header=F)

data <- merge(ans, map, by.x="FID", by.y="V1")

list <- read.xlsx("/home/richards/tomoko.nakanishi/scratch/09.COVID19/05.BQC/BQC_phenotype/20211002_Clusters_for_Tomoko.xlsx")
list <- list %>% mutate(anonymized_patient_id = paste0(Center,"_",Local_Code))

list <- list %>% merge(data, by.x="anonymized_patient_id", by.y="V2", all.x=T)

prs <- readRDS("09.CHUM/all_prs.rds")
list <- list %>% merge(prs, by="IID", all.x=T)

list <- list %>% mutate(prs = ifelse(UMAP_POP == "eur", prs, NA))
list <- list %>% mutate(prs = scale(prs))
list$PHATE_Cluster <- as.factor(list$PHATE_Cluster)
list$`DECEASED.<DSO60` <- as.factor(list$`DECEASED.<DSO60`)

write.xlsx(list, "09.CHUM/20211002_Clusters_for_Tomoko_appendPRS.xlsx")
ggplot(list, aes(x=prs)) +
  geom_density(aes(y = ..density.., fill = PHATE_Cluster), adjust = 1, alpha = 0.5) 

ggplot(list, aes(x=prs)) +
  geom_density(aes(y = ..density.., fill = `DECEASED.<DSO60`), adjust = 1, alpha = 0.5) 

LM <- glm(`DECEASED.<DSO60` ~ PHATE_Cluster + Age + Sex + prs, data=list, family="binomial")
summary(LM)

colnames(data)[4] <- "anonymized_patient_id"
data %>% write.xlsx("../../data/v3.1/12.CHUM_data/CHUM_ancestry_rs10774671_20210426.xlsx")