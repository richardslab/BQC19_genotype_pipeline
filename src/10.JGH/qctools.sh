#!/bin/bash
#PBS -N qctools
#PBS -o logs/
#PBS -e logs/
#PBS -l mem=5G
#PBS -l vmem=5G
#PBS -l walltime=10:00:00
#PBS -l nodes=1:ppn=1
#PBS -t 1-23
# Below is where to put your job code:

cd /home/richards/tomoko.nakanishi/09.COVID19/data/05.BQC/01.genotype/v4.0
qctool2 -g 07.TopMed/$PBS_ARRAYID.bgen \
-incl-samples 01.batch/JGH.sample \
-og 10.JGH/${PBS_ARRAYID}.bgen -bgen-bits "8" \
-os 10.JGH/${PBS_ARRAYID}.sample
