path=/scratch/richards/tomoko.nakanishi/09.COVID19/05.BQC/01.genotypeQC/v5.0/
data=/project/richards/tomoko.nakanishi/09.COVID19/data/05.BQC/01.genotype/v5.0

awk '$2 ~ /CRCHUM/ {print $1"_"$1,$2,$3}' $data/01.batch/sampleid.studyid.sex.map > $data/09.CHUM/CHUM_sample_list

for i in {1..23}
do
bcftools view -S <(awk '{print $1}' $data/09.CHUM/CHUM_sample_list) --force-samples $data/07.TopMed/chr${i}.dose.vcf.gz -Oz -o $path/09.CHUM/${i}.vcf.gz &
done

