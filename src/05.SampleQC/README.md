
```{bash}
path=/scratch/richards/tomoko.nakanishi/09.COVID19/05.BQC/01.genotypeQC/v5.0/

awk 'BEGEIN{FS="\t"}{print $2}' $path/bqc-v5.0-qc5-EUR.bim > $path/qc.snplist

plink --bfile $path/bqc19-v5.0-qc3 \
--extract $path/qc.snplist \
--exclude ../../data/LdRegion-AbecasisHg38.txt \
--snps-only 'just-acgt' \
--missing \
--het \
--out $path/bqc19-v5.0-qc4

plink --bfile $path/bqc19-v5.0-qc3 \
--extract $path/qc.snplist \
--exclude ../../data/LdRegion-AbecasisHg38.txt \
--snps-only 'just-acgt' \
--homozyg \
--homozyg-kb 1000 \
--out $path/bqc19-v5.0-qc4

paste $path/bqc19-v5.0-qc4.het $path/bqc19-v5.0-qc4.imiss | tail -n +2 | awk '{ print $1, $2, ($5-$3)/$5, $12 }' > $path/het-vs-imiss.txt

#hetmis.outlier remove sample 01.HetMisOutlier.R 

cat $path/hetmis.outlier > $path/all.remove

plink --bfile $path/bqc19-v5.0-qc3 \
--extract $path/qc.snplist \
--chr 23,24 \
--geno 0.05 \
--maf 0.05 \
--indep-pairphase 20000 2000 0.5 \
--out $path/bqc19-v5.0-qc3

plink --bfile $path/bqc19-v5.0-qc3 \
--check-sex 0.4 0.7 \
--extract $path/bqc19-v5.0-qc3.prune.in \
--out $path/bqc19-v5.0-qc4

plink --bfile $path/bqc19-v5.0-qc3 \
--chr 23,24 \
--check-sex ycount 0.4 0.7 \
--out $path/bqc19-v5.0-qc5

#sexdiscordance check 02.Sexdiscordance.R

cat $path/sex.discordance.remove >> $path/all.remove

awk 'FNR==NR { m[$1]=$4; next } $1 in m { print $1, $2, $3, $4, m[$1], $6 }' $path/bqc19-v5.0-qc5.sexcheck \
$path/bqc19-v5.0-qc3.fam > tmp
mv tmp $path/bqc19-v5.0-qc3.fam

plink --bfile $path/bqc19-v5.0-qc3 \
--remove $path/all.remove \
--extract <(awk 'BEGEIN{FS="\t"}{print $2}' $path/bqc-v5.0-qc5-EUR.bim) \
--make-bed \
--out $path/bqc19-v5.0-qc6

plink --bfile $path/bqc19-v5.0-qc6 \
--biallelic-only strict \
--chr 1-22 \
--exclude range ../../data/LdRegion-AbecasisHg38.txt \
--geno 0.05 \
--hwe 1E-6 midp \
--indep-pairwise 50 5 0.05 \
--keep-allele-order \
--mac 5 \
--maf 0.01 \
--out $path/bqc19-v5.0-qc6

plink --bfile $path/bqc19-v5.0-qc6 \
--extract $path/bqc19-v5.0-qc6.prune.in \
--make-bed \
--out $path/bqc19-v5.0-qc7

king -b $path/bqc19-v5.0-qc7.bed --duplicate 
cut -f3 king.con | sort | uniq | grep -v "FID" | awk '$1 != " "{print $1, $1}' >> $path/all.remove
mv king* $path

plink --bfile $path/bqc19-v5.0-qc3 \
--remove $path/all.remove \
--extract <(awk 'BEGEIN{FS="\t"}{print $2}' $path/bqc-v5.0-qc5-EUR.bim) \
--make-bed \
--out $path/bqc19-v5.0-qc6

#final QCed sample and genotype!!!
```
