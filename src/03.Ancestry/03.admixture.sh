path=/home/richards/tomoko.nakanishi/scratch/09.COVID19/05.BQC/01.genotypeQC/v4.0
data=/home/richards/tomoko.nakanishi/09.COVID19/src/05.BQC/01.genotypeQC/data/v4.0/


for K in {5..8}
do
admixture --cv $path/bqc19_hgdp_1kg.bed $K -j40 | tee $path/log${K}.out
done
