```{bash}
path=/scratch/richards/tomoko.nakanishi/09.COVID19/05.BQC/01.genotypeQC/v4.0
data=/project/richards/tomoko.nakanishi/09.COVID19/data/05.BQC/01.genotype/v4.0/08.PC/
plink --bfile $path/bqc19-v4.0-qc6 \
--biallelic-only strict \
--chr 1-22 \
--exclude range ../data/LdRegion-AbecasisHg38.txt \
--geno 0.05 \
--hwe 1E-6 midp \
--indep-pairwise 50 5 0.05 \
--keep-allele-order \
--mac 5 \
--maf 0.01 \
--out $path/bqc19-v4.0-qc6

plink --bfile $path/bqc19-v4.0-qc6 \
--extract $path/bqc19-v4.0-qc6.prune.in \
--make-bed --out $data/grm_final

awk -F" " '($23 == "eur" && $25 ==1){print $1,$1}' $data/../03.Ancestry/all.sample.ancestry > $data/EUR.sample
awk -F" " '($23 == "afr" && $25 ==1){print $1,$1}' $data/../03.Ancestry/all.sample.ancestry > $data/AFR.sample
awk -F" " '($23 == "sas" && $25 ==1){print $1,$1}' $data/../03.Ancestry/all.sample.ancestry > $data/SAS.sample
awk -F" " '($23 == "eas" && $25 ==1){print $1,$1}' $data/../03.Ancestry/all.sample.ancestry > $data/EAS.sample
awk -F" " '($23 == "amr" && $25 ==1){print $1,$1}' $data/../03.Ancestry/all.sample.ancestry > $data/AMR.sample

for POP in EUR SAS AFR AMR EAS
do
plink --bfile  $data/grm_final \
--keep $data/${POP}.sample \
--remove $path/kingunrelated_toberemoved.txt \
--make-bed --out $path/unrelated_${POP}

plink --bfile  $data/grm_final \
--keep <(cat $data/${POP}.sample) \
--remove $path/kingunrelated_toberemoved.txt \
--make-bed --out $path/all_${POP}

flashpca --bfile $path/unrelated_${POP} --outload $path/${POP}-loadings.pop \
--outmeansd $path/${POP}-meansd.pop --outpc $path/${POP}.pc --suffix .pop --ndim 10
flashpca --bfile $path/all_${POP} --project --inmeansd $path/${POP}-meansd.pop --outproj $path/all_${POP}.projections --inload $path/${POP}-loadings.pop -v
 
##for regenie
plink2 \
  --bfile $path/bqc19-v4.0-qc6 \
  --keep $data/${POP}.sample \
  --maf 0.01 --mac 100 --geno 0.1 --hwe 1e-15 \
  --chr 1-22 \
  --mind 0.1 \
  --write-snplist --write-samples --no-id-header \
  --out $data/${POP}_qc_pass

done

```
