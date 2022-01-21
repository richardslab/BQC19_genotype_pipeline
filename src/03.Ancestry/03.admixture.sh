path=/scratch/richards/tomoko.nakanishi/09.COVID19/05.BQC/01.genotypeQC/v5.0
data=/project/richards/tomoko.nakanishi/09.COVID19/data/05.BQC/01.genotype/v5.0/


for K in {5..7}
do
admixture --cv $path/bqc19_hgdp_1kg.bed $K -j40 | tee $path/log${K}.out
done
