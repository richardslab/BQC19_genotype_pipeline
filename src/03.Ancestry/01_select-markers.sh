path=/home/richards/tomoko.nakanishi/scratch/09.COVID19/05.BQC/01.genotypeQC/v4.0
data=/home/richards/tomoko.nakanishi/09.COVID19/src/05.BQC/01.genotypeQC/data/v4.0/

#select hapmap3 SNPs from hgdp and 1kg dataset
for chrom in {1..22} X Y
do
plink --vcf /project/richards/public/hgdp_tgp/gnomad.genomes.v3.1.hgdp_1kg_subset.chr${chrom}.vcf.bgz \
--extract ../../data/hapmap3.list \
--maf 0.01 \
--mac 5 \
--geno 0.05 \
--hwe 1e-6 \
--exclude range ../../data/LdRegion-AbecasisHg38.txt \
--snps-only just-acgt \
--make-bed --out $path/hgdp_1kg_${chrom} --threads 40 
awk '($5~/[AT]/ && $6 ~/[AT]/) || ($5~/[GC]/ && $6~/[GC]/) { print $2 }' $path/hgdp_1kg_${chrom}.bim > $path/hgdp_1kg_${chrom}.palin
plink --bfile $path/hgdp_1kg_${chrom} --exclude $path/hgdp_1kg_${chrom}.palin --make-bed --out $path/hgdp_1kg_${chrom}.rev
done

for chrom in {1..22}
do
echo $path/hgdp_1kg_${chrom}.rev.bed $path/hgdp_1kg_${chrom}.rev.bim $path/hgdp_1kg_${chrom}.rev.fam >> $path/hgdp_1kg.mergelst
done

plink --bfile $path/hgdp_1kg_1.rev --merge-list $path/hgdp_1kg.mergelst --make-bed --out ../../data/hgdp_1kg

for chrom in {1..22} X Y
do
rm $path/hgdp_1kg_${chrom}.*
done

cut -f 2 $path/hgdp_1kg.bim > $path/hgdp_1kg.snps

#select intersect SNPs in BQC19
plink --bfile $path/bqc19-v4.0-qc2 \
--make-bed \
--out $path/bqc19-v4.0-qc3

plink --bfile $path/bqc19-v4.0-qc3 \
--extract $path/hgdp_1kg.snps \
--indep-pairwise 1000 kb 5 0.1 \
--out $path/bqc19-v4.0-qc3.hgdp_1kg

plink --bfile $path/bqc19-v4.0-qc3 \
--extract $path/bqc19-v4.0-qc3.hgdp_1kg.prune.in \
--make-bed \
--out $path/bqc19-v4.0-qc4.hgdp_1kg

awk '{print $2}' $path/bqc19-v4.0-qc4.hgdp_1kg.bim > $path/bqc19-v4.0-qc4.hgdp_1kg.snp

plink --bfile ../../data/hgdp_1kg \
--a1-allele $path/bqc19-v4.0-qc4.hgdp_1kg.bim 5 2 \
--extract $path/bqc19-v4.0-qc4.hgdp_1kg.snp \
--exclude $path/bqc19_hgdp_1kg-merge.missnp \
--make-bed \
--out $path/hgdp_1kg.bqc19-v4.0-qc4

echo "$path/bqc19-v4.0-qc4.hgdp_1kg.bed $path/bqc19-v4.0-qc4.hgdp_1kg.bim $path/bqc19-v4.0-qc4.hgdp_1kg.fam" > $path/bqc19-v4.0-qc4.list

plink --bfile $path/hgdp_1kg.bqc19-v4.0-qc4 --merge-list $path/bqc19-v4.0-qc4.list \
--make-bed --out $path/bqc19_hgdp_1kg

king -b $path/hgdp_1kg.bqc19-v4.0-qc4.bed --unrelated --degree 2

plink --bfile ../../data/hgdp_1kg \
--a1-allele $path/bqc19-v4.0-qc4.hgdp_1kg.bim 5 2 \
--extract $path/bqc19-v4.0-qc4.hgdp_1kg.snp \
--exclude $path/bqc19_hgdp_1kg-merge.missnp \
--remove kingunrelated_toberemoved.txt \
--make-bed \
--out $path/hgdp_1kg.bqc19-v4.0-qc4

plink --bfile $path/hgdp_1kg.bqc19-v4.0-qc4 --merge-list $path/bqc19-v4.0-qc4.list \
--remove kingunrelated_toberemoved.txt \
--make-bed --out $path/bqc19_hgdp_1kg

flashpca --bfile $path/hgdp_1kg.bqc19-v4.0-qc4 --outload $path/hgdp_1kg.bqc19-v4.0-loadings.pop --outmeansd $path/hgdp_1kg.bqc19-v4.0-meansd.pop --suffix .pop --ndim 10
flashpca --bfile $path/bqc19-v4.0-qc4.hgdp_1kg --project --inmeansd $path/hgdp_1kg.bqc19-v4.0-meansd.pop --outproj $path/bqc19-v4.0-qc4.projections --inload $path/hgdp_1kg.bqc19-v4.0-loadings.pop -v

mv *.pop $path
