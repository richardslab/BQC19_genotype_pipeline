path=/scratch/richards/tomoko.nakanishi/09.COVID19/05.BQC/01.genotypeQC/v5.0
data=/project/richards/tomoko.nakanishi/09.COVID19/data/05.BQC/01.genotype/v5.0

#remove 5507784388583022821343
cut -f7 /project/richards/restricted/bqc19/data/array-genotype-data/BQC_AX001_AX032R_BPW/EXPORT/SampleReport.txt | tail -n+2 | sort | uniq | grep -v "5507784388583022821343" > $path/batch

counter=0
for b in {A..Z} {a..e}; {
    PREFIX=$path/bqc-b${b}-qc3
    counter=$(($counter+1))
    B=$(sed -n ${counter}p $path/batch)
    awk 'BEGIN{FS=" "}($26 == "'${B}'"){print $1,$2}' $path/bqc19-v5.0-qc4.EUR.IID_IID.batch > ${PREFIX}.EUR.IID_IID
    # All
    plink --bfile $path/bqc19-v5.0-qc3 \
    --a1-allele $path/bqc19-v5.0-qc3.bim 5 2 \
    --keep ${PREFIX}.EUR.IID_IID \
    --make-bed \
    --out ${PREFIX}-EUR
}

for i in {A..Z} {a..e}; { 
BGRP=$(for j in {A..Z} {a..e}; do printf $j; done | sed -e "s/${i}//g")
rm -f $path/b${BGRP}.txt
touch $path/b${BGRP}.txt
for b in $(for j in {A..Z} {a..e}; do echo $j; done | grep -v "$i"); { 
PREFIX=$path/bqc-b${b}-qc3-EUR ; echo $PREFIX >> $path/b${BGRP}.txt 
}
plink --merge-list $path/b${BGRP}.txt --make-bed --out $path/bqc-b${BGRP}-qc3-EUR
}

cat $path/bqc-b*-qc3.EUR.IID_IID > $path/bqc-bAtoe-qc3.EUR.IID_IID

for j in {A..Z} {a..e}; {
for b in $j $(for i in {A..Z} {a..e}; do printf $i; done | sed -e "s/${j}//g") ; {
PREFIX=$path/bqc-b${b}-qc3-EUR
ls ${PREFIX}.bed ${PREFIX}.bim ${PREFIX}.fam
# All
plink \
--bfile ${PREFIX} \
--a1-allele ${PREFIX}.bim 5 2 \
--freqx \
--hardy midp \
--missing \
--out ${PREFIX}
# Males
plink \
--bfile ${PREFIX} \
--a1-allele ${PREFIX}.bim 5 2 \
--freqx \
--hardy midp \
--missing \
--filter-males \
--out ${PREFIX}.males
plink \
--bfile ${PREFIX} \
--a1-allele ${PREFIX}.bim 5 2 \
--freqx \
--hardy midp \
--missing \
--filter-females \
--out ${PREFIX}.females
# All
plink \
--bfile ${PREFIX} \
--freq \
--out ${PREFIX}
# Males
plink \
--bfile ${PREFIX} \
--freq \
--filter-males \
--out ${PREFIX}.males
# Females
plink \
--bfile ${PREFIX} \
--freq \
--filter-females \
--out ${PREFIX}.females
}
}


for j in {A..Z} {a..e}; {
for b in $(for i in {A..Z} {a..e}; do printf $i; done | sed -e "s/${j}//g") ; {
Rscript 02_generate-gt-freq-fisher-sex.R $j $b &
} 
}
