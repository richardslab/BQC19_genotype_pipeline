PRSDIR=/home/richards/tomoko.nakanishi/scratch/09.COVID19/03.ICDA/04.PRS
GENODIR=/home/richards/tomoko.nakanishi/09.COVID19/data/05.BQC/01.genotype/v4.0/
OUTDIR=/home/richards/tomoko.nakanishi/09.COVID19/data/05.BQC/01.genotype/v4.0/09.CHUM

for chrom in {1..22}
do
plink2 --bgen ${GENODIR}/07.TopMed/${chrom}.bgen 'ref-first' \
--sample ${GENODIR}/07.TopMed/1-23.sample \
--memory 4000 \
--threads 40 \
--score <(awk 'NR != 1{print "chr"$3":"toupper($4)":"toupper($5), toupper($4),$8}NR != 1{print "chr"$3":"toupper($5)":"toupper($4), toupper($4),$8}' ${PRSDIR}/5e-4_0.7_B2_meta_filtered.rsid.txt) 'cols=+scoresums' \
--out ${OUTDIR}/all_${chrom}
done
