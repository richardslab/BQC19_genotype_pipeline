path=/project/richards/tomoko.nakanishi/09.COVID19/data/05.BQC/01.genotype/v5.0/07.TopMed

qctool2 -filetype "vcf" -g $path/chr$PBS_ARRAYID.dose.vcf.gz -og $path/$PBS_ARRAYID.bgen -vcf-genotype-field GP -bgen-bits "8" -os $path/$PBS_ARRAYID.sample
