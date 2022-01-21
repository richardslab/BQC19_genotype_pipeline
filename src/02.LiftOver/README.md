
```{bash}
path=/scratch/richards/tomoko.nakanishi/09.COVID19/05.BQC/01.genotypeQC/v5.0
data=/project/richards/tomoko.nakanishi/09.COVID19/data/05.BQC/01.genotype/v5.0
```

##BQC19
```{bash}
plink --file /project/richards/restricted/bqc19/data/array-genotype-data/BQC_AX001_AX032R_BPW/EXPORT/BQC-AX-001toAX032R_PLINK \
--make-bed --out $path/bqc19-v5.0-raw --no-fid --no-parents --no-sex --no-pheno
cat $path/bqc19-v5.0-raw.bim | awk -F "\t" '($1 != 25){print $0}($1 == 25){OFS="\t"; print 23,$2,$3,$4,$5,$6}' > $path/bqc19-v5.0-raw1.bim
plink --bed $path/bqc19-v5.0-raw.bed --bim $path/bqc19-v5.0-raw1.bim --fam $path/bqc19-v5.0-raw.fam --split-x 'hg19' --make-bed --out $path/bqc19-v5.0-raw2

awk 'FNR==NR { m[$1]=$6; next } ($1 in m) { print $1, $2, $3, $4, m[$1], $6 }' /project/richards/restricted/bqc19/data/array-genotype-data/BQC_AX001_AX032R_BPW/EXPORT/SampleReport.txt \
$path/bqc19-v5.0-raw2.fam | sed -e "s/female/2/g" | sed -e "s/male/1/g" | sed -e "s/unknown/-9/g" | grep -v "Negative" | grep -v "NA24385" | grep -v "CEPH1463" | grep -v "BLANK" \
> $path/bqc19-v5.0-qc1.fam

awk 'BEGIN{FS="\t"} FNR==NR {m[$1]=$2; next} ($2 in m){OFS="\t"; print $1,m[$2],$3,$4,$5,$6} !($2 in m){OFS="\t"; print $0}' \
<(cat ../../data/Axiom_PMRA.na35.annot.tbl | sed -e "s/ /\t/g") \
$path/bqc19-v5.0-raw2.bim > $path/bqc19-v5.0-qc1.bim
##867939 SNPs

#liftover
cat $path/bqc19-v5.0-qc1.bim | awk -F"\t" '($1 != 23 && $1 != 24 && $1 != 25 && $1 != 26){OFS="\t"; print "chr"$1,$4,$4+1,$2}($1 == 23 || $1 == 25){OFS="\t"; print "chrX",$4,$4+1,$2}($1 == 24){OFS="\t"; print "chrY",$4,$4+1,$2}($1 == 26){OFS="\t"; print "chrMT",$4,$4+1,$2}' > $path/b37.bed
liftOver $path/b37.bed ~/projects/richards/share/scratch/ucsc/hg19/hg19ToHg38.over.chain $path/b38.bed $path/unMapped

plink --bed $path/bqc19-v5.0-raw2.bed --fam $path/bqc19-v5.0-raw2.fam --bim $path/bqc19-v5.0-qc1.bim \
--keep $data/01.batch/all.sample.removed.duplicated \
--extract <(awk '{print $4}' $path/b38.bed) \
--make-bed --out $path/bqc19-v5.0-qc2

awk -F "\t" '(FNR=NR){m[$4] = $1 && n[$4] = $2}($2 in m){OFS="\t"; print $1,$2,$3,n[$2],$5,$6}' $path/b38.bed $path/bqc19-v5.0-qc2.bim > $path/tmp
awk -F "\t" '($1 != 25){OFS="\t"; print $0} ($1 == 25){OFS="\t"; print 23,$2,$3,$4,$5,$6}' $path/tmp > $path/bqc19-v5.0-qc2.bim

plink --bfile $path/bqc19-v5.0-qc2 --split-x 'hg38' --make-bed --out $path/bqc-v5.0-qc2

```

##SweCovid
attempted to combine, but there were very few SNPs overlapped.
BQC19 used Axiom Precision Medicine array (N=867698 SNPs, b37) and Swecovid used GSA (N=730059 SNPs, b38).
However, there were only 117481 SNPs that were overlapped btw two chips.
