#LiftOver from hg19 to Hg38
cat $path/bqc19-v3.1-qc6.bim | awk -F"\t" '($1 != 23 && $1 != 25){OFS="\t"; print "chr"$1,$4,$4+1,$2}($1 == 23 || $1 == 25){OFS="\t"; print "chrX",$4,$4+1,$2}' > $path/b37.bed
liftOver $path/b37.bed ~/projects/richards/share/scratch/ucsc/hg19/hg19ToHg38.over.chain $path/b38.bed $path/unMapped

plink --bfile $path/bqc19-v3.1-qc6 \
--extract <(awk '{print $4}' $path/b38.bed) \
--make-bed --out $path/bqc19-v3.1-qc7

awk -F "\t" '(FNR=NR){m[$4] = $1 && n[$4] = $2}($2 in m){OFS="\t"; print $1,$2,$3,n[$2],$5,$6}' $path/b38.bed $path/bqc19-v3.1-qc7.bim > $path/tmp
awk -F "\t" '($1 != 25){OFS="\t"; print $0} ($1 == 25){OFS="\t"; print 23,$2,$3,$4,$5,$6}' $path/tmp > $path/bqc19-v3.1-qc7.rev.bim
cd $path/
ln -s bqc19-v3.1-qc7.bed bqc19-v3.1-qc7.rev.bed
ln -s bqc19-v3.1-qc7.fam bqc19-v3.1-qc7.rev.fam

plink --bfile $path/bqc19-v3.1-qc7.rev \
--freq \
--out $path/bqc19-v3.1-qc8

# Run HRC check tool to prepare PLINK input:
perl HRC-1000G-check-bim.pl \
-b $path/bqc19-v3.1-qc7.rev.bim \
-f $path/bqc19-v3.1-qc8.frq \
-r ../../data/v2.1/07.prepare_imputation/ALL.TOPMed_freeze5_hg38_dbSNP.PASS.sites.tab -h

# Convert PLINK to VCF running some checks:
mv $path/Run-plink.sh .
bash Run-plink.sh

for i in {1..22}
do
grep "#" $path/bqc19-v3.1-qc7.rev-updated-chr${i}.vcf > $path/tmp
grep -v "#" $path/bqc19-v3.1-qc7.rev-updated-chr${i}.vcf | awk '{OFS=""; print "chr"$0}' >> $path/tmp
bgzip -c $path/tmp > $path/bqc19-v3.1-qc7.rev-updated-chr${i}.vcf.gz
done

for i in 23
do
grep "#" $path/bqc19-v3.1-qc7.rev-updated-chr${i}.vcf > $path/tmp
grep -v "#" $path/bqc19-v3.1-qc7.rev-updated-chr${i}.vcf | awk '{OFS=""; print "chr"$0}' | sed -e "s/chr23/chrX/g" >> $path/tmp
bgzip -c $path/tmp > $path/bqc19-v3.1-qc7.rev-updated-chr${i}.vcf.gz
done

##upload bqc19-v3.1-qc7.rev-updated-chr*.vcf.gz
