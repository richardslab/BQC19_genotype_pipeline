plink --bfile /scratch/richards/tomoko.nakanishi/09.COVID19/05.BQC/01.genotypeQC/v5.0/bqc19-v5.0-qc7 --exclude /scratch/richards/tomoko.nakanishi/09.COVID19/05.BQC/01.genotypeQC/v5.0/Exclude-bqc19-v5.0-qc7-HRC.txt --make-bed --out /scratch/richards/tomoko.nakanishi/09.COVID19/05.BQC/01.genotypeQC/v5.0/TEMP1
plink --bfile /scratch/richards/tomoko.nakanishi/09.COVID19/05.BQC/01.genotypeQC/v5.0/TEMP1 --update-map /scratch/richards/tomoko.nakanishi/09.COVID19/05.BQC/01.genotypeQC/v5.0/Chromosome-bqc19-v5.0-qc7-HRC.txt --update-chr --make-bed --out /scratch/richards/tomoko.nakanishi/09.COVID19/05.BQC/01.genotypeQC/v5.0/TEMP2
plink --bfile /scratch/richards/tomoko.nakanishi/09.COVID19/05.BQC/01.genotypeQC/v5.0/TEMP2 --update-map /scratch/richards/tomoko.nakanishi/09.COVID19/05.BQC/01.genotypeQC/v5.0/Position-bqc19-v5.0-qc7-HRC.txt --make-bed --out /scratch/richards/tomoko.nakanishi/09.COVID19/05.BQC/01.genotypeQC/v5.0/TEMP3
plink --bfile /scratch/richards/tomoko.nakanishi/09.COVID19/05.BQC/01.genotypeQC/v5.0/TEMP3 --flip /scratch/richards/tomoko.nakanishi/09.COVID19/05.BQC/01.genotypeQC/v5.0/Strand-Flip-bqc19-v5.0-qc7-HRC.txt --make-bed --out /scratch/richards/tomoko.nakanishi/09.COVID19/05.BQC/01.genotypeQC/v5.0/TEMP4
plink --bfile /scratch/richards/tomoko.nakanishi/09.COVID19/05.BQC/01.genotypeQC/v5.0/TEMP4 --a2-allele /scratch/richards/tomoko.nakanishi/09.COVID19/05.BQC/01.genotypeQC/v5.0/Force-Allele1-bqc19-v5.0-qc7-HRC.txt --make-bed --out /scratch/richards/tomoko.nakanishi/09.COVID19/05.BQC/01.genotypeQC/v5.0/bqc19-v5.0-qc7-updated
plink --bfile /scratch/richards/tomoko.nakanishi/09.COVID19/05.BQC/01.genotypeQC/v5.0/bqc19-v5.0-qc7-updated --real-ref-alleles --make-bed --chr 1 --out /scratch/richards/tomoko.nakanishi/09.COVID19/05.BQC/01.genotypeQC/v5.0/bqc19-v5.0-qc7-updated-chr1
plink --bfile /scratch/richards/tomoko.nakanishi/09.COVID19/05.BQC/01.genotypeQC/v5.0/bqc19-v5.0-qc7-updated --real-ref-alleles --recode vcf --chr 1 --out /scratch/richards/tomoko.nakanishi/09.COVID19/05.BQC/01.genotypeQC/v5.0/bqc19-v5.0-qc7-updated-chr1
plink --bfile /scratch/richards/tomoko.nakanishi/09.COVID19/05.BQC/01.genotypeQC/v5.0/bqc19-v5.0-qc7-updated --real-ref-alleles --make-bed --chr 2 --out /scratch/richards/tomoko.nakanishi/09.COVID19/05.BQC/01.genotypeQC/v5.0/bqc19-v5.0-qc7-updated-chr2
plink --bfile /scratch/richards/tomoko.nakanishi/09.COVID19/05.BQC/01.genotypeQC/v5.0/bqc19-v5.0-qc7-updated --real-ref-alleles --recode vcf --chr 2 --out /scratch/richards/tomoko.nakanishi/09.COVID19/05.BQC/01.genotypeQC/v5.0/bqc19-v5.0-qc7-updated-chr2
plink --bfile /scratch/richards/tomoko.nakanishi/09.COVID19/05.BQC/01.genotypeQC/v5.0/bqc19-v5.0-qc7-updated --real-ref-alleles --make-bed --chr 3 --out /scratch/richards/tomoko.nakanishi/09.COVID19/05.BQC/01.genotypeQC/v5.0/bqc19-v5.0-qc7-updated-chr3
plink --bfile /scratch/richards/tomoko.nakanishi/09.COVID19/05.BQC/01.genotypeQC/v5.0/bqc19-v5.0-qc7-updated --real-ref-alleles --recode vcf --chr 3 --out /scratch/richards/tomoko.nakanishi/09.COVID19/05.BQC/01.genotypeQC/v5.0/bqc19-v5.0-qc7-updated-chr3
plink --bfile /scratch/richards/tomoko.nakanishi/09.COVID19/05.BQC/01.genotypeQC/v5.0/bqc19-v5.0-qc7-updated --real-ref-alleles --make-bed --chr 4 --out /scratch/richards/tomoko.nakanishi/09.COVID19/05.BQC/01.genotypeQC/v5.0/bqc19-v5.0-qc7-updated-chr4
plink --bfile /scratch/richards/tomoko.nakanishi/09.COVID19/05.BQC/01.genotypeQC/v5.0/bqc19-v5.0-qc7-updated --real-ref-alleles --recode vcf --chr 4 --out /scratch/richards/tomoko.nakanishi/09.COVID19/05.BQC/01.genotypeQC/v5.0/bqc19-v5.0-qc7-updated-chr4
plink --bfile /scratch/richards/tomoko.nakanishi/09.COVID19/05.BQC/01.genotypeQC/v5.0/bqc19-v5.0-qc7-updated --real-ref-alleles --make-bed --chr 5 --out /scratch/richards/tomoko.nakanishi/09.COVID19/05.BQC/01.genotypeQC/v5.0/bqc19-v5.0-qc7-updated-chr5
plink --bfile /scratch/richards/tomoko.nakanishi/09.COVID19/05.BQC/01.genotypeQC/v5.0/bqc19-v5.0-qc7-updated --real-ref-alleles --recode vcf --chr 5 --out /scratch/richards/tomoko.nakanishi/09.COVID19/05.BQC/01.genotypeQC/v5.0/bqc19-v5.0-qc7-updated-chr5
plink --bfile /scratch/richards/tomoko.nakanishi/09.COVID19/05.BQC/01.genotypeQC/v5.0/bqc19-v5.0-qc7-updated --real-ref-alleles --make-bed --chr 6 --out /scratch/richards/tomoko.nakanishi/09.COVID19/05.BQC/01.genotypeQC/v5.0/bqc19-v5.0-qc7-updated-chr6
plink --bfile /scratch/richards/tomoko.nakanishi/09.COVID19/05.BQC/01.genotypeQC/v5.0/bqc19-v5.0-qc7-updated --real-ref-alleles --recode vcf --chr 6 --out /scratch/richards/tomoko.nakanishi/09.COVID19/05.BQC/01.genotypeQC/v5.0/bqc19-v5.0-qc7-updated-chr6
plink --bfile /scratch/richards/tomoko.nakanishi/09.COVID19/05.BQC/01.genotypeQC/v5.0/bqc19-v5.0-qc7-updated --real-ref-alleles --make-bed --chr 7 --out /scratch/richards/tomoko.nakanishi/09.COVID19/05.BQC/01.genotypeQC/v5.0/bqc19-v5.0-qc7-updated-chr7
plink --bfile /scratch/richards/tomoko.nakanishi/09.COVID19/05.BQC/01.genotypeQC/v5.0/bqc19-v5.0-qc7-updated --real-ref-alleles --recode vcf --chr 7 --out /scratch/richards/tomoko.nakanishi/09.COVID19/05.BQC/01.genotypeQC/v5.0/bqc19-v5.0-qc7-updated-chr7
plink --bfile /scratch/richards/tomoko.nakanishi/09.COVID19/05.BQC/01.genotypeQC/v5.0/bqc19-v5.0-qc7-updated --real-ref-alleles --make-bed --chr 8 --out /scratch/richards/tomoko.nakanishi/09.COVID19/05.BQC/01.genotypeQC/v5.0/bqc19-v5.0-qc7-updated-chr8
plink --bfile /scratch/richards/tomoko.nakanishi/09.COVID19/05.BQC/01.genotypeQC/v5.0/bqc19-v5.0-qc7-updated --real-ref-alleles --recode vcf --chr 8 --out /scratch/richards/tomoko.nakanishi/09.COVID19/05.BQC/01.genotypeQC/v5.0/bqc19-v5.0-qc7-updated-chr8
plink --bfile /scratch/richards/tomoko.nakanishi/09.COVID19/05.BQC/01.genotypeQC/v5.0/bqc19-v5.0-qc7-updated --real-ref-alleles --make-bed --chr 9 --out /scratch/richards/tomoko.nakanishi/09.COVID19/05.BQC/01.genotypeQC/v5.0/bqc19-v5.0-qc7-updated-chr9
plink --bfile /scratch/richards/tomoko.nakanishi/09.COVID19/05.BQC/01.genotypeQC/v5.0/bqc19-v5.0-qc7-updated --real-ref-alleles --recode vcf --chr 9 --out /scratch/richards/tomoko.nakanishi/09.COVID19/05.BQC/01.genotypeQC/v5.0/bqc19-v5.0-qc7-updated-chr9
plink --bfile /scratch/richards/tomoko.nakanishi/09.COVID19/05.BQC/01.genotypeQC/v5.0/bqc19-v5.0-qc7-updated --real-ref-alleles --make-bed --chr 10 --out /scratch/richards/tomoko.nakanishi/09.COVID19/05.BQC/01.genotypeQC/v5.0/bqc19-v5.0-qc7-updated-chr10
plink --bfile /scratch/richards/tomoko.nakanishi/09.COVID19/05.BQC/01.genotypeQC/v5.0/bqc19-v5.0-qc7-updated --real-ref-alleles --recode vcf --chr 10 --out /scratch/richards/tomoko.nakanishi/09.COVID19/05.BQC/01.genotypeQC/v5.0/bqc19-v5.0-qc7-updated-chr10
plink --bfile /scratch/richards/tomoko.nakanishi/09.COVID19/05.BQC/01.genotypeQC/v5.0/bqc19-v5.0-qc7-updated --real-ref-alleles --make-bed --chr 11 --out /scratch/richards/tomoko.nakanishi/09.COVID19/05.BQC/01.genotypeQC/v5.0/bqc19-v5.0-qc7-updated-chr11
plink --bfile /scratch/richards/tomoko.nakanishi/09.COVID19/05.BQC/01.genotypeQC/v5.0/bqc19-v5.0-qc7-updated --real-ref-alleles --recode vcf --chr 11 --out /scratch/richards/tomoko.nakanishi/09.COVID19/05.BQC/01.genotypeQC/v5.0/bqc19-v5.0-qc7-updated-chr11
plink --bfile /scratch/richards/tomoko.nakanishi/09.COVID19/05.BQC/01.genotypeQC/v5.0/bqc19-v5.0-qc7-updated --real-ref-alleles --make-bed --chr 12 --out /scratch/richards/tomoko.nakanishi/09.COVID19/05.BQC/01.genotypeQC/v5.0/bqc19-v5.0-qc7-updated-chr12
plink --bfile /scratch/richards/tomoko.nakanishi/09.COVID19/05.BQC/01.genotypeQC/v5.0/bqc19-v5.0-qc7-updated --real-ref-alleles --recode vcf --chr 12 --out /scratch/richards/tomoko.nakanishi/09.COVID19/05.BQC/01.genotypeQC/v5.0/bqc19-v5.0-qc7-updated-chr12
plink --bfile /scratch/richards/tomoko.nakanishi/09.COVID19/05.BQC/01.genotypeQC/v5.0/bqc19-v5.0-qc7-updated --real-ref-alleles --make-bed --chr 13 --out /scratch/richards/tomoko.nakanishi/09.COVID19/05.BQC/01.genotypeQC/v5.0/bqc19-v5.0-qc7-updated-chr13
plink --bfile /scratch/richards/tomoko.nakanishi/09.COVID19/05.BQC/01.genotypeQC/v5.0/bqc19-v5.0-qc7-updated --real-ref-alleles --recode vcf --chr 13 --out /scratch/richards/tomoko.nakanishi/09.COVID19/05.BQC/01.genotypeQC/v5.0/bqc19-v5.0-qc7-updated-chr13
plink --bfile /scratch/richards/tomoko.nakanishi/09.COVID19/05.BQC/01.genotypeQC/v5.0/bqc19-v5.0-qc7-updated --real-ref-alleles --make-bed --chr 14 --out /scratch/richards/tomoko.nakanishi/09.COVID19/05.BQC/01.genotypeQC/v5.0/bqc19-v5.0-qc7-updated-chr14
plink --bfile /scratch/richards/tomoko.nakanishi/09.COVID19/05.BQC/01.genotypeQC/v5.0/bqc19-v5.0-qc7-updated --real-ref-alleles --recode vcf --chr 14 --out /scratch/richards/tomoko.nakanishi/09.COVID19/05.BQC/01.genotypeQC/v5.0/bqc19-v5.0-qc7-updated-chr14
plink --bfile /scratch/richards/tomoko.nakanishi/09.COVID19/05.BQC/01.genotypeQC/v5.0/bqc19-v5.0-qc7-updated --real-ref-alleles --make-bed --chr 15 --out /scratch/richards/tomoko.nakanishi/09.COVID19/05.BQC/01.genotypeQC/v5.0/bqc19-v5.0-qc7-updated-chr15
plink --bfile /scratch/richards/tomoko.nakanishi/09.COVID19/05.BQC/01.genotypeQC/v5.0/bqc19-v5.0-qc7-updated --real-ref-alleles --recode vcf --chr 15 --out /scratch/richards/tomoko.nakanishi/09.COVID19/05.BQC/01.genotypeQC/v5.0/bqc19-v5.0-qc7-updated-chr15
plink --bfile /scratch/richards/tomoko.nakanishi/09.COVID19/05.BQC/01.genotypeQC/v5.0/bqc19-v5.0-qc7-updated --real-ref-alleles --make-bed --chr 16 --out /scratch/richards/tomoko.nakanishi/09.COVID19/05.BQC/01.genotypeQC/v5.0/bqc19-v5.0-qc7-updated-chr16
plink --bfile /scratch/richards/tomoko.nakanishi/09.COVID19/05.BQC/01.genotypeQC/v5.0/bqc19-v5.0-qc7-updated --real-ref-alleles --recode vcf --chr 16 --out /scratch/richards/tomoko.nakanishi/09.COVID19/05.BQC/01.genotypeQC/v5.0/bqc19-v5.0-qc7-updated-chr16
plink --bfile /scratch/richards/tomoko.nakanishi/09.COVID19/05.BQC/01.genotypeQC/v5.0/bqc19-v5.0-qc7-updated --real-ref-alleles --make-bed --chr 17 --out /scratch/richards/tomoko.nakanishi/09.COVID19/05.BQC/01.genotypeQC/v5.0/bqc19-v5.0-qc7-updated-chr17
plink --bfile /scratch/richards/tomoko.nakanishi/09.COVID19/05.BQC/01.genotypeQC/v5.0/bqc19-v5.0-qc7-updated --real-ref-alleles --recode vcf --chr 17 --out /scratch/richards/tomoko.nakanishi/09.COVID19/05.BQC/01.genotypeQC/v5.0/bqc19-v5.0-qc7-updated-chr17
plink --bfile /scratch/richards/tomoko.nakanishi/09.COVID19/05.BQC/01.genotypeQC/v5.0/bqc19-v5.0-qc7-updated --real-ref-alleles --make-bed --chr 18 --out /scratch/richards/tomoko.nakanishi/09.COVID19/05.BQC/01.genotypeQC/v5.0/bqc19-v5.0-qc7-updated-chr18
plink --bfile /scratch/richards/tomoko.nakanishi/09.COVID19/05.BQC/01.genotypeQC/v5.0/bqc19-v5.0-qc7-updated --real-ref-alleles --recode vcf --chr 18 --out /scratch/richards/tomoko.nakanishi/09.COVID19/05.BQC/01.genotypeQC/v5.0/bqc19-v5.0-qc7-updated-chr18
plink --bfile /scratch/richards/tomoko.nakanishi/09.COVID19/05.BQC/01.genotypeQC/v5.0/bqc19-v5.0-qc7-updated --real-ref-alleles --make-bed --chr 19 --out /scratch/richards/tomoko.nakanishi/09.COVID19/05.BQC/01.genotypeQC/v5.0/bqc19-v5.0-qc7-updated-chr19
plink --bfile /scratch/richards/tomoko.nakanishi/09.COVID19/05.BQC/01.genotypeQC/v5.0/bqc19-v5.0-qc7-updated --real-ref-alleles --recode vcf --chr 19 --out /scratch/richards/tomoko.nakanishi/09.COVID19/05.BQC/01.genotypeQC/v5.0/bqc19-v5.0-qc7-updated-chr19
plink --bfile /scratch/richards/tomoko.nakanishi/09.COVID19/05.BQC/01.genotypeQC/v5.0/bqc19-v5.0-qc7-updated --real-ref-alleles --make-bed --chr 20 --out /scratch/richards/tomoko.nakanishi/09.COVID19/05.BQC/01.genotypeQC/v5.0/bqc19-v5.0-qc7-updated-chr20
plink --bfile /scratch/richards/tomoko.nakanishi/09.COVID19/05.BQC/01.genotypeQC/v5.0/bqc19-v5.0-qc7-updated --real-ref-alleles --recode vcf --chr 20 --out /scratch/richards/tomoko.nakanishi/09.COVID19/05.BQC/01.genotypeQC/v5.0/bqc19-v5.0-qc7-updated-chr20
plink --bfile /scratch/richards/tomoko.nakanishi/09.COVID19/05.BQC/01.genotypeQC/v5.0/bqc19-v5.0-qc7-updated --real-ref-alleles --make-bed --chr 21 --out /scratch/richards/tomoko.nakanishi/09.COVID19/05.BQC/01.genotypeQC/v5.0/bqc19-v5.0-qc7-updated-chr21
plink --bfile /scratch/richards/tomoko.nakanishi/09.COVID19/05.BQC/01.genotypeQC/v5.0/bqc19-v5.0-qc7-updated --real-ref-alleles --recode vcf --chr 21 --out /scratch/richards/tomoko.nakanishi/09.COVID19/05.BQC/01.genotypeQC/v5.0/bqc19-v5.0-qc7-updated-chr21
plink --bfile /scratch/richards/tomoko.nakanishi/09.COVID19/05.BQC/01.genotypeQC/v5.0/bqc19-v5.0-qc7-updated --real-ref-alleles --make-bed --chr 22 --out /scratch/richards/tomoko.nakanishi/09.COVID19/05.BQC/01.genotypeQC/v5.0/bqc19-v5.0-qc7-updated-chr22
plink --bfile /scratch/richards/tomoko.nakanishi/09.COVID19/05.BQC/01.genotypeQC/v5.0/bqc19-v5.0-qc7-updated --real-ref-alleles --recode vcf --chr 22 --out /scratch/richards/tomoko.nakanishi/09.COVID19/05.BQC/01.genotypeQC/v5.0/bqc19-v5.0-qc7-updated-chr22
plink --bfile /scratch/richards/tomoko.nakanishi/09.COVID19/05.BQC/01.genotypeQC/v5.0/bqc19-v5.0-qc7-updated --real-ref-alleles --make-bed --chr 23 --out /scratch/richards/tomoko.nakanishi/09.COVID19/05.BQC/01.genotypeQC/v5.0/bqc19-v5.0-qc7-updated-chr23
plink --bfile /scratch/richards/tomoko.nakanishi/09.COVID19/05.BQC/01.genotypeQC/v5.0/bqc19-v5.0-qc7-updated --real-ref-alleles --recode vcf --chr 23 --out /scratch/richards/tomoko.nakanishi/09.COVID19/05.BQC/01.genotypeQC/v5.0/bqc19-v5.0-qc7-updated-chr23
rm TEMP*