path=/scratch/richards/tomoko.nakanishi/09.COVID19/05.BQC/01.genotypeQC/v5.0/

plink --bfile $path/bqc-v5.0-qc4-EUR \
--chr 23 \
--maf 0.001 \
--geno 0.02 \
--filter-females \
--exclude <(cat $path/CEPH1463.fail.SNP $path/NA24385.fail.SNP) \
--hwe 1e-6 midp \
--make-bed \
--out $path/bqc-v5.0-qc5-EUR.X

plink --bfile $path/bqc-v5.0-qc4-EUR \
--extract <(awk '{print $2}' $path/bqc-v5.0-qc5-EUR.X.bim) \
--make-bed \
--out $path/bqc-v5.0-qc5-EUR.X.rev

plink --bfile $path/bqc-v5.0-qc4-EUR \
--chr 1-22,25 \
--maf 0.001 \
--geno 0.02 \
--exclude <(cat $path/CEPH1463.fail.SNP $path/NA24385.fail.SNP) \
--hwe 1e-6 midp \
--make-bed \
--out $path/bqc-v5.0-qc5-EUR.1-22

echo "$path/bqc-v5.0-qc5-EUR.X.rev" > $path/mergelist1

plink --bfile $path/bqc-v5.0-qc5-EUR.1-22 \
--merge-list $path/mergelist1 --make-bed --out $path/bqc-v5.0-qc5-EUR

#503693 (v4 507114) > final QCed snp
