setwd("/scratch/richards/tomoko.nakanishi/09.COVID19/05.BQC/01.genotypeQC/v5.0/")
library(foreach)
library(doMC)
library(data.table)

run.fisher <- function(gt) {
  f <- foreach(i=1:nrow(gt), .combine=rbind) %dopar% tryCatch({ 
    c(gt[i][,1],gt[i][,2], 
      fisher.p=fisher.test(matrix(as.numeric(gt[i][,3:ncol(gt)]), nrow=2), 
                           workspace = 2000000)$p.value) }, 
    warning = function(w) { NA },
    error = function(e) { NA },
    finally = { NA })
  as.data.table(f)
}

get.prefix <- function(batch, suffix){
  paste0("bqc-b",batch,"-qc3-EUR.",suffix)
}

args <- commandArgs(trailingOnly=TRUE)
q <- args[1] # 1
t <- args[2] # 2345

message("SETTING OUT TO SEA")
frqx <- merge(fread(get.prefix(q, "frqx"), header=TRUE), 
              fread(get.prefix(t, "frqx"), header=TRUE), 
              by="SNP", suffixes = c(".1", ".2"))
setkey(frqx, SNP)

frqx.males <- merge(fread(get.prefix(q,"males.frqx"), header=TRUE), 
                    fread(get.prefix(t,"males.frqx"), header=TRUE), 
                    by="SNP", suffixes = c(".1", ".2"))
setkey(frqx.males, SNP)

frqx.females <- merge(fread(get.prefix(q,"females.frqx"), header=TRUE), 
                      fread(get.prefix(t,"females.frqx"), header=TRUE), 
                      by="SNP", suffixes = c(".1", ".2"))
setkey(frqx.females, SNP)

message("CASTING NET")
# diploid chromosomes, but include X (diploid in females), and PAR XY.
gt.diploid.autoXXY <- frqx[CHR.1 %in% c(1:22,23,25), 
                           .(CHR.1, SNP, `C(HOM A1).1`, `C(HOM A1).2`, `C(HET).1`, 
                             `C(HET).2`, `C(HOM A2).1`, `C(HOM A2).2`)]

# diploid X in females only
gt.diploid.females.X <- frqx.females[CHR.1 %in% c(23), 
                                     .(CHR.1, SNP, `C(HOM A1).1`, `C(HOM A1).2`, `C(HET).1`, 
                                       `C(HET).2`, `C(HOM A2).1`, `C(HOM A2).2`)]

# Haploid mitochondria
#gt.haploid.mito <- frqx[CHR.1 %in% c(26), 
#                        .(CHR.1, SNP, `C(HOM A1).1`, `C(HOM A1).2`, `C(HOM A2).1`, 
#                          `C(HOM A2).2`)]

# Haploid male sex chromosomes X and Y
gt.haploid.males.sex <- frqx.males[CHR.1 %in% c(23,24), 
                                   .(CHR.1, SNP, `C(HAP A1).1`, `C(HAP A1).2`, `C(HAP A2).1`, 
                                     `C(HAP A2).2`)]

message("GONE FISHER'ing")
set.seed(1234)
registerDoMC(8)
f <- list(run.fisher(gt.diploid.autoXXY),
          run.fisher(gt.diploid.females.X),
          run.fisher(gt.haploid.males.sex))

message("REELING IN CATCH")
r <- rbindlist(f, use.names = TRUE)

message("COOKIN' DINNER")
fwrite(r, file=paste0("bqc-b",q,"vs",t,"-qc3-EUR.het.fisher.p.txt"), quote=FALSE, sep="\t")
