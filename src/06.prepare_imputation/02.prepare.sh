path=/home/richards/tomoko.nakanishi/scratch/09.COVID19/05.BQC/01.genotypeQC/v4.0

plink --bfile $path/bqc19-v4.0-qc6 \
--freq \
--out $path/bqc19-v4.0-qc7

# Run HRC check tool to prepare PLINK input:
perl HRC-1000G-check-bim.pl \
-b $path/bqc19-v4.0-qc6.bim \
-f $path/bqc19-v4.0-qc7.frq \
-r ../../data/ALL.TOPMed_freeze5_hg38_dbSNP.PASS.sites.tab -h

# Convert PLINK to VCF running some checks:
mv $path/Run-plink.sh .
bash Run-plink.sh

for i in {1..22}
do
grep "#" $path/bqc19-v4.0-qc6-updated-chr${i}.vcf > $path/tmp
grep -v "#" $path/bqc19-v4.0-qc6-updated-chr${i}.vcf | awk '{OFS=""; print "chr"$0}' >> $path/tmp
bgzip -c $path/tmp > $path/bqc19-v4.0-qc6-updated-chr${i}.vcf.gz
done

for i in 23
do
grep "#" $path/bqc19-v4.0-qc6-updated-chr${i}.vcf > $path/tmp
grep -v "#" $path/bqc19-v4.0-qc6-updated-chr${i}.vcf | awk '{OFS=""; print "chr"$0}' | sed -e "s/chr23/chrX/g" >> $path/tmp
bgzip -c $path/tmp > $path/bqc19-v4.0-qc6-updated-chr${i}.vcf.gz
done

##upload bqc19-v4.0-qc6-updated-chr*.vcf.gz
