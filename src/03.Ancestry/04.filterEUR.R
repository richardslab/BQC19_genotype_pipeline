setwd("/home/richards/tomoko.nakanishi/scratch/09.COVID19/05.BQC/01.genotypeQC/v4.0/")
library(data.table)

tbl=read.table("bqc19_hgdp_1kg.5.Q")
colnames(tbl) <- paste0("ADMIXTURE",1:5)
all <- fread("bqc19_hgdp_1kg.fam", header = F)
all <- all %>% select(V1, V2)
tbl1 <- cbind(all, tbl)
tbl1$id <- paste0(tbl1$V1, tbl1$V2)
pop <- fread("sample.all", header=TRUE)
pop$id <- paste0(pop$FID,pop$IID)

res <- merge(tbl1, pop, by = "id")
res <- res %>% mutate(UMAP_POP = case_when(clusterUMAP %in% c(16:17) ~ "eur",
                                           clusterUMAP %in% c(9:10,13,14) ~ "afr",
                                           clusterUMAP %in% c(2,3,4) ~ "eas",
                                           clusterUMAP %in% c(5) ~ "eur",#fin
                                           clusterUMAP %in% c(15) ~ "mid",
                                           clusterUMAP %in% c(11,12) ~ "amr",
                                           clusterUMAP %in% c(6,7,8) ~ "sas",
                                           TRUE ~ "oth"
))

library(dplyr)
library(tidyr)

res_long <- gather(res, ADMIXTURE, value, ADMIXTURE1:ADMIXTURE5, factor_key=TRUE)
tmp <- res_long %>% filter(UMAP_POP %in% c("eur"))
tmp <- tmp %>% arrange(ADMIXTURE, value,labeled_subpop)

png("ADMIXTURE.nfe.fin.png", width=1200)
ggplot(tmp, aes(x=id, y=value, fill=ADMIXTURE)) + geom_bar(stat="identity") +
  facet_grid(~fct_inorder(labeled_subpop),scales = "free_x",space = "free_x") +
  theme_minimal() + 
  theme(
    axis.text.x = element_blank(),
    panel.grid = element_blank(), panel.spacing = unit(0, "lines"),
    strip.text.x = element_text(angle = 90, size=10)
  )  + scale_fill_brewer(palette="Dark2")
dev.off()

tmp_nfe <- tmp %>% filter(ADMIXTURE == "ADMIXTURE2")
min_nfe <- min(tmp_nfe$value[tmp_nfe$labeled_subpop %in% c("ceu", "french", "gbr","basque", "ibs","italian","tsi","tuscan", "orcadian")])

tmp_afr <- res_long %>% filter(UMAP_POP %in% c("afr"))
tmp_afr <- tmp_afr %>% arrange(ADMIXTURE, value,labeled_subpop)

png("ADMIXTURE.afr.png", width=1200)
ggplot(tmp_afr, aes(x=id, y=value, fill=ADMIXTURE)) + geom_bar(stat="identity") +
  facet_grid(~fct_inorder(labeled_subpop),scales = "free_x",space = "free_x") +
  theme_minimal() + 
  theme(
    axis.text.x = element_blank(),
    panel.grid = element_blank(), panel.spacing = unit(0, "lines"),
    strip.text.x = element_text(angle = 90, size=10)
  )  + scale_fill_brewer(palette="Dark2")
dev.off()

tmp_afr <- tmp_afr %>% filter(ADMIXTURE == "ADMIXTURE1")
min_afr <- min(tmp_afr$value[tmp_afr$labeled_subpop %in% c("esn","gwd","mandenka","msl","yoruba","yri")])

tmp_eas <- res_long %>% filter(UMAP_POP %in% c("eas"))
tmp_eas <- tmp_eas %>% arrange(ADMIXTURE, value,labeled_subpop)

png("ADMIXTURE.eas.png", width=1200)
ggplot(tmp_eas, aes(x=id, y=value, fill=ADMIXTURE)) + geom_bar(stat="identity") +
  facet_grid(~fct_inorder(labeled_subpop),scales = "free_x",space = "free_x") +
  theme_minimal() + 
  theme(
    axis.text.x = element_blank(),
    panel.grid = element_blank(), panel.spacing = unit(0, "lines"),
    strip.text.x = element_text(angle = 90, size=10)
  )  + scale_fill_brewer(palette="Dark2")
dev.off()

tmp_eas <- tmp_eas %>% filter(ADMIXTURE == "ADMIXTURE4")
min_eas <- min(tmp_eas$value)

tmp_sas <- res_long %>% filter(UMAP_POP %in% c("sas"))
tmp_sas <- tmp_sas %>% arrange(ADMIXTURE, value,labeled_subpop)

png("ADMIXTURE.sas.png", width=1200)
ggplot(tmp_sas, aes(x=id, y=value, fill=ADMIXTURE)) + geom_bar(stat="identity") +
  facet_grid(~fct_inorder(labeled_subpop),scales = "free_x",space = "free_x") +
  theme_minimal() + 
  theme(
    axis.text.x = element_blank(),
    panel.grid = element_blank(), panel.spacing = unit(0, "lines"),
    strip.text.x = element_text(angle = 90, size=10)
  )  + scale_fill_brewer(palette="Dark2")
dev.off()

tmp_sas <- tmp_sas %>% filter(ADMIXTURE == "ADMIXTURE3")
min_sas <- min(tmp_sas$value)

tmp_amr <- res_long %>% filter(UMAP_POP %in% c("amr"))
tmp_amr <- tmp_amr %>% arrange(ADMIXTURE, value,labeled_subpop)

png("ADMXITURE.amr.png", width=1200)
ggplot(tmp_amr, aes(x=id, y=value, fill=ADMIXTURE)) + geom_bar(stat="identity") +
  facet_grid(~fct_inorder(labeled_subpop),scales = "free_x",space = "free_x") +
  theme_minimal() + 
  theme(
    axis.text.x = element_blank(),
    panel.grid = element_blank(), panel.spacing = unit(0, "lines"),
    strip.text.x = element_text(angle = 90, size=10)
  )  + scale_fill_brewer(palette="Dark2")
dev.off()


res <- res %>% mutate(ADMIXTURE_POP = case_when(ADMIXTURE2 >= min_nfe ~ "eur",
                                                ADMIXTURE1 >= min_afr ~ "afr",
                                                ADMIXTURE4 >= min_eas ~ "eas",
                                                ADMIXTURE3 >= min_sas ~ "sas",
                                                UMAP_POP == "0BQC" ~ "amr",
                                                TRUE ~ "oth"))

res <- res %>% mutate(confident = ifelse(UMAP_POP == ADMIXTURE_POP, 1, 0))
res <- res %>% select(FID, IID, ADMIXTURE1, ADMIXTURE2, ADMIXTURE3, ADMIXTURE4, ADMIXTURE5,
                      PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10, population_inference.pop, labeled_subpop,
                      UMAP1, UMAP2, clusterUMAP, UMAP_POP, ADMIXTURE_POP, confident)

batch <- fread("~/projects/richards/restricted/bqc19/data/array-genotype-data/Export_BQC-AX001toAX028_JGH-CRCHUM/ALL DATA AX001 to AX031/SampleReport_AX001toAX031.txt")
batch <- batch[,c(1,7)]
colnames(batch) <- c("ID", "batch")
dat2 <- merge(res, batch, by.x="FID", by.y="ID")
dat2 <- dat2 %>% filter(UMAP_POP == "eur" & confident == 1)

write.table(dat2, file="bqc19-v4.0-qc4.EUR.IID_IID.batch", quote=F, col.names = F, row.names = F)

write.table(res, file="all.sample.ancestry", col.names = T, row.names = F, quote = F, sep=" ")

