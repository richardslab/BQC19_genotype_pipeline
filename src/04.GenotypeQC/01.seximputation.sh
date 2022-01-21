path=/scratch/richards/tomoko.nakanishi/09.COVID19/05.BQC/01.genotypeQC/v5.0
data=/project/richards/tomoko.nakanishi/09.COVID19/data/05.BQC/01.genotype/v5.0

plink --bfile $path/bqc19-v5.0-qc3 \
--chr 23,24 \
--geno 0.05 \
--maf 0.05 \
--indep-pairphase 20000 2000 0.5 \
--out $path/bqc19-v5.0-qc3

plink --bfile $path/bqc19-v5.0-qc3 \
--check-sex 0.4 0.8 \
--extract $path/bqc19-v5.0-qc3.prune.in \
--out $path/bqc19-v5.0-qc4

awk 'FNR==NR { m[$1]=$4; next } $1 in m { print $1, $2, $3, $4, m[$1], $6 }' $path/bqc19-v5.0-qc4.sexcheck \
$path/bqc19-v5.0-qc3.fam > tmp
mv tmp $path/bqc19-v5.0-qc3.fam
