#!/bin/bash
#PBS -N plink
#PBS -o logs/
#PBS -e logs/
#PBS -l mem=70G
#PBS -l vmem=70G
#PBS -l walltime=5:00:00
#PBS -l nodes=1:ppn=1
##PBS -t 1-23
# Below is where to put your job code:

cd /home/richards/tomoko.nakanishi/09.COVID19/data/05.BQC/01.genotype/v4.0
plink --bfile /home/richards/tomoko.nakanishi/09.COVID19/src/05.BQC/01.genotypeQC/data/v4.0/06.SampleQC/bqc19-v4.0-qc6 \
--keep 01.batch/JGH.geno.sample \
--memory 50000 \
--make-bed --out 10.JGH/bqc19_jgh_v4.0
