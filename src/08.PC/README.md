```{bash}
path=/scratch/richards/tomoko.nakanishi/09.COVID19/05.BQC/01.genotypeQC/v5.0/
data=/project/richards/tomoko.nakanishi/09.COVID19/data/05.BQC/01.genotype/v5.0

plink --bfile $path/bqc19-v5.0-qc6 \
--biallelic-only strict \
--chr 1-22 \
--exclude range ../data/LdRegion-AbecasisHg38.txt \
--geno 0.05 \
--hwe 1E-6 midp \
--indep-pairwise 50 5 0.05 \
--keep-allele-order \
--mac 5 \
--maf 0.01 \
--out $path/bqc19-v5.0-qc6

plink --bfile $path/bqc19-v5.0-qc6 \
--extract $path/bqc19-v5.0-qc6.prune.in \
--make-bed --out $data/08.PC/grm_final

king -b $data/08.PC/grm_final.bed --unrelated --degree 2

mv *.txt $path

awk -F" " '($25 == "nfe"){print $1,$1}' $data/03.Ancestry/all.sample.ancestry > $data/08.PC/EUR.sample
awk -F" " '($25 == "afr"){print $1,$1}' $data/03.Ancestry/all.sample.ancestry > $data/08.PC/AFR.sample
awk -F" " '($25 == "sas"){print $1,$1}' $data/03.Ancestry/all.sample.ancestry > $data/08.PC/SAS.sample
awk -F" " '($25 == "eas"){print $1,$1}' $data/03.Ancestry/all.sample.ancestry > $data/08.PC/EAS.sample
awk -F" " '($25 == "amr"){print $1,$1}' $data/03.Ancestry/all.sample.ancestry > $data/08.PC/AMR.sample

for POP in EUR SAS AFR AMR EAS
do
plink --bfile  $data/08.PC/grm_final \
--keep $data/08.PC/${POP}.sample \
--remove $path/kingunrelated_toberemoved.txt \
--make-bed --out $path/unrelated_${POP}

plink --bfile  $data/08.PC/grm_final \
--keep $data/08.PC/${POP}.sample \
--make-bed --out $path/all_${POP}

flashpca --bfile $path/unrelated_${POP} --outload $path/${POP}-loadings.pop \
--outmeansd $path/${POP}-meansd.pop --outpc $data/08.PC/${POP}.pc --suffix .pop --ndim 10
flashpca --bfile $path/all_${POP} --project --inmeansd $path/${POP}-meansd.pop --outproj $data/08.PC/all_${POP}.projections --inload $path/${POP}-loadings.pop -v
 
##for regenie
plink2 \
  --bfile $path/bqc19-v5.0-qc6 \
  --keep $data/08.PC/${POP}.sample \
  --maf 0.01 --mac 100 --geno 0.1 --hwe 1e-15 \
  --chr 1-22 \
  --mind 0.1 \
  --write-snplist --write-samples --no-id-header \
  --out $data/08.PC/${POP}_qc_pass

done

mv *pop $path
```
