setwd("/scratch/richards/tomoko.nakanishi/09.COVID19/05.BQC/01.genotypeQC/v5.0/")

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
res <- res %>% mutate(UMAP_POP = case_when(clusterUMAP %in% c(15:16) ~ "nfe",
                                             clusterUMAP %in% c(2, 3, 12, 13) ~ "afr",
                                             clusterUMAP %in% c(4,5,6,9) ~ "eas",
                                             clusterUMAP %in% c(14) ~ "mid",
                                             clusterUMAP %in% c(11) ~ "amr",
                                             clusterUMAP %in% c(7,9,10) ~ "sas",
                                             TRUE ~ "oth"
))


library(dplyr)
library(tidyr)

res_long <- gather(res, ADMIXTURE, value, ADMIXTURE1:ADMIXTURE5, factor_key=TRUE)
tmp <- res_long %>% filter(UMAP_POP %in% c("nfe"))
tmp <- tmp %>% arrange(ADMIXTURE, value,labeled_subpop)

png("/project/richards/tomoko.nakanishi/repo/BQC19_genotype_pipeline/results/03.Ancestry/ADMIXTURE.nfe.png", width=1200)
tmp %>% mutate(name = fct_reorder(id, value)) %>%
  ggplot(aes(x=name, y=value, fill=ADMIXTURE)) + geom_bar(stat="identity") +
  facet_grid(~fct_inorder(labeled_subpop),scales = "free_x",space = "free_x") +
  theme_minimal() + 
  theme(
    axis.text.x = element_blank(),
    panel.grid = element_blank(), panel.spacing = unit(0, "lines"),
    strip.text.x = element_text(angle = 90, size=10)
  )  + scale_fill_brewer(palette="Dark2")
dev.off()

tmp_nfe <- tmp %>% filter(ADMIXTURE == "ADMIXTURE2")
min_nfe <- min(tmp_nfe$value[tmp_nfe$labeled_subpop %in% c("basque",  "ceu", "french",  "gbr",  "ibs", "italian", "orcadian", "tsi", "tuscan")])
#0.922999

tmp_afr <- res_long %>% filter(UMAP_POP %in% c("afr"))
tmp_afr <- tmp_afr %>% arrange(ADMIXTURE, value,labeled_subpop)

png("/project/richards/tomoko.nakanishi/repo/BQC19_genotype_pipeline/results/03.Ancestry/ADMIXTURE.afr.png", width=1200)
tmp_afr %>% mutate(name = fct_reorder(id, value)) %>%
  ggplot(aes(x=name, y=value, fill=ADMIXTURE)) + geom_bar(stat="identity") +
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
#0.914902

tmp_eas <- res_long %>% filter(UMAP_POP %in% c("eas"))
tmp_eas <- tmp_eas %>% arrange(ADMIXTURE, value,labeled_subpop)

png("/project/richards/tomoko.nakanishi/repo/BQC19_genotype_pipeline/results/03.Ancestry/ADMIXTURE.eas.png", width=1200)
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
dat <- tmp_eas %>% group_by(labeled_subpop) %>% summarise(mean = mean(value, na.rm=T))

eas_list <- c("chs", "cdx", "chs", "dai", "han", "hezhen", "jpt", "khv", "lahu", "miaozu",
              "naxi", "she", "tujia", "yizu")
min_eas <- min(tmp_eas$value[tmp_eas$labeled_subpop %in% eas_list])
#0.891172

tmp_sas <- res_long %>% filter(UMAP_POP %in% c("sas"))
tmp_sas <- tmp_sas %>% arrange(ADMIXTURE, value,labeled_subpop)

png("/project/richards/tomoko.nakanishi/repo/BQC19_genotype_pipeline/results/03.Ancestry/ADMIXTURE.sas.png", width=1200)
ggplot(tmp_sas, aes(x=id, y=value, fill=ADMIXTURE)) + geom_bar(stat="identity") +
  facet_grid(~fct_inorder(labeled_subpop),scales = "free_x",space = "free_x") +
  theme_minimal() + 
  theme(
    axis.text.x = element_blank(),
    panel.grid = element_blank(), panel.spacing = unit(0, "lines"),
    strip.text.x = element_text(angle = 90, size=10)
  )  + scale_fill_brewer(palette="Dark2")
dev.off()

tmp_sas <- tmp_sas %>% filter(ADMIXTURE == "ADMIXTURE5")
dat <- tmp_sas %>% group_by(labeled_subpop) %>% summarise(mean = mean(value, na.rm=T))

sas_list <- c("itu", "stu", "beb", "pjl", "gih")
min_sas <- min(tmp_sas$value[tmp_sas$labeled_subpop %in% sas_list])
#0.715746

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

tmp_amr <- tmp_amr %>% filter(ADMIXTURE == "ADMIXTURE3")
dat <- tmp_amr %>% group_by(labeled_subpop) %>% summarise(mean = mean(value, na.rm=T))

amr_list <- c("mxl", "pur", "clm", "pel")
#min_amr <- min(tmp_amr$value[tmp_amr$labeled_subpop %in% amr_list])

res <- res %>% mutate(ADMIXTURE_POP = case_when(ADMIXTURE2 >= min_nfe ~ "nfe",
                                                ADMIXTURE1 >= min_afr ~ "afr",
                                                ADMIXTURE4 >= min_eas ~ "eas",
                                                ADMIXTURE5 >= min_sas ~ "sas",
                                                UMAP_POP == "amr" ~ "amr",
                                                TRUE ~ "oth")) %>%
  mutate(final_POP = ifelse(UMAP_POP == ADMIXTURE_POP, UMAP_POP, "oth"))

res <- res %>% mutate(nfe_UMAP = ifelse(UMAP_POP == "nfe", TRUE, FALSE),
                      nfe_ADMIXTURE = ifelse(ADMIXTURE_POP == "nfe", TRUE, FALSE),
                      nfe_final = ifelse(final_POP == "nfe", TRUE, FALSE))

library(gridExtra)
png(file="/project/richards/tomoko.nakanishi/repo/BQC19_genotype_pipeline/results/03.Ancestry/popstrat.PC12.png", width=1500, height=1200, res=150)
p1 <- ggplot(res, aes(PC1, PC2, col=population_inference.pop)) + geom_point() + xlim(min(res$PC1)-0.01, max(res$PC1)+0.01) +
  ylim(min(res$PC2)-0.01, max(res$PC2)+0.01) + theme_minimal() + theme(legend.title = element_text(size=8))
bqc <- res %>% filter(population_inference.pop == "0BQC")
p2 <- ggplot(bqc, aes(PC1, PC2, col=nfe_UMAP)) + geom_point() + xlim(min(res$PC1)-0.01, max(res$PC1)+0.01) +
  ylim(min(res$PC2)-0.01, max(res$PC2)+0.01) + theme_minimal() + scale_color_manual(values=c("#999999", "#E69F00"))
p3 <- ggplot(bqc, aes(PC1, PC2, col=nfe_ADMIXTURE)) + geom_point() + xlim(min(res$PC1)-0.01, max(res$PC1)+0.01) +
  ylim(min(res$PC2)-0.01, max(res$PC2)+0.01) + theme_minimal() + scale_color_manual(values=c("#999999", "#E69F00"))
p4 <- ggplot(bqc, aes(PC1, PC2, col=nfe_final)) + geom_point() + xlim(min(res$PC1)-0.01, max(res$PC1)+0.01) +
  ylim(min(res$PC2)-0.01, max(res$PC2)+0.01) + theme_minimal() + scale_color_manual(values=c("#999999", "#E69F00"))
grid.arrange(
  p1,
  p2,
  p3,
  p4,
  nrow = 2)
dev.off() 

png(file="/project/richards/tomoko.nakanishi/repo/BQC19_genotype_pipeline/results/03.Ancestry/popstrat.PC34.png", width=1500, height=1200, res=150)
p1 <- ggplot(res, aes(PC3, PC4, col=population_inference.pop)) + geom_point() + xlim(min(res$PC3)-0.01, max(res$PC3)+0.01) +
  ylim(min(res$PC4)-0.01, max(res$PC4)+0.01) + theme_minimal() + theme(legend.title = element_text(size=8))
bqc <- res %>% filter(population_inference.pop == "0BQC")
p2 <- ggplot(bqc, aes(PC3, PC4, col=nfe_UMAP)) + geom_point() + xlim(min(res$PC3)-0.01, max(res$PC3)+0.01) +
  ylim(min(res$PC4)-0.01, max(res$PC4)+0.01) + theme_minimal() + scale_color_manual(values=c("#999999", "#E69F00"))
p3 <- ggplot(bqc, aes(PC3, PC4, col=nfe_ADMIXTURE)) + geom_point() + xlim(min(res$PC3)-0.01, max(res$PC3)+0.01) +
  ylim(min(res$PC4)-0.01, max(res$PC4)+0.01) + theme_minimal() + scale_color_manual(values=c("#999999", "#E69F00"))
p4 <- ggplot(bqc, aes(PC3, PC4, col=nfe_final)) + geom_point() + xlim(min(res$PC3)-0.01, max(res$PC3)+0.01) +
  ylim(min(res$PC4)-0.01, max(res$PC4)+0.01) + theme_minimal() + scale_color_manual(values=c("#999999", "#E69F00"))
grid.arrange(
  p1,
  p2,
  p3,
  p4,
  nrow = 2)
dev.off()  

png(file="/project/richards/tomoko.nakanishi/repo/BQC19_genotype_pipeline/results/03.Ancestry/popstrat.PC56.png", width=1500, height=1200, res=150)
p1 <- ggplot(res, aes(PC5, PC6, col=population_inference.pop)) + geom_point() + xlim(min(res$PC5)-0.01, max(res$PC5)+0.01) +
  ylim(min(res$PC6)-0.01, max(res$PC6)+0.01) + theme_minimal() + theme(legend.title = element_text(size=8))
bqc <- res %>% filter(population_inference.pop == "0BQC")
p2 <- ggplot(bqc, aes(PC5, PC6, col=nfe_UMAP)) + geom_point() + xlim(min(res$PC5)-0.01, max(res$PC5)+0.01) +
  ylim(min(res$PC6)-0.01, max(res$PC6)+0.01) + theme_minimal() + scale_color_manual(values=c("#999999", "#E69F00"))
p3 <- ggplot(bqc, aes(PC5, PC6, col=nfe_ADMIXTURE)) + geom_point() + xlim(min(res$PC5)-0.01, max(res$PC5)+0.01) +
  ylim(min(res$PC6)-0.01, max(res$PC6)+0.01) + theme_minimal() + scale_color_manual(values=c("#999999", "#E69F00"))
p4 <- ggplot(bqc, aes(PC5, PC6, col=nfe_final)) + geom_point() + xlim(min(res$PC5)-0.01, max(res$PC5)+0.01) +
  ylim(min(res$PC6)-0.01, max(res$PC6)+0.01) + theme_minimal() + scale_color_manual(values=c("#999999", "#E69F00"))
grid.arrange(
  p1,
  p2,
  p3,
  p4,
  nrow = 2)
dev.off()  

batch <- fread("/project/richards/restricted/bqc19/data/array-genotype-data/BQC_AX001_AX032R_BPW/EXPORT/SampleReport.txt")
batch <- batch[,c(1,7)]
colnames(batch) <- c("ID", "batch")
res1 <- res %>% 
  select(FID, IID, ADMIXTURE1, ADMIXTURE2, ADMIXTURE3, ADMIXTURE4, ADMIXTURE5, 
         PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10, population_inference.pop, labeled_subpop, UMAP1, UMAP2, clusterUMAP, UMAP_POP,
         ADMIXTURE_POP, final_POP)

dat2 <- merge(res1, batch, by.x="FID", by.y="ID")
dat2 <- dat2 %>% filter(final_POP == "nfe")

write.table(dat2, file="bqc19-v5.0-qc4.EUR.IID_IID.batch", quote=F, col.names = F, row.names = F)

write.table(res1, file="/project/richards/tomoko.nakanishi/09.COVID19/data/05.BQC/01.genotype/v5.0/all.sample.ancestry", col.names = T, row.names = F, quote = F, sep=" ")
