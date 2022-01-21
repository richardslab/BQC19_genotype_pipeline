setwd("/project/richards/tomoko.nakanishi/09.COVID19/data/05.BQC/01.genotype/v5.0/08.PC")

res <- fread("EUR.pc")

png(file="/project/richards/tomoko.nakanishi/repo/BQC19_genotype_pipeline/results/08.PC/EUR.PCs.png", width=1500, height=1200, res=150)
p1 <- ggplot(res, aes(PC1, PC2)) + geom_point() + theme_minimal() + theme(legend.title = element_text(size=8))
p2 <- ggplot(res, aes(PC3, PC4)) + geom_point() + theme_minimal() + scale_color_manual(values=c("#999999", "#E69F00"))
p3 <- ggplot(res, aes(PC5, PC6)) + geom_point() + theme_minimal() + scale_color_manual(values=c("#999999", "#E69F00"))
p4 <- ggplot(res, aes(PC7, PC8)) + geom_point() + theme_minimal() + scale_color_manual(values=c("#999999", "#E69F00"))
grid.arrange(
  p1,
  p2,
  p3,
  p4,
  nrow = 2)
dev.off()
