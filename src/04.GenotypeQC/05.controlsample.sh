plink --bed $path/bqc19-v4.0-raw2.bed --bim $path/bqc19-v4.0-qc1.bim --fam $path/bqc19-v4.0-raw2.fam \
--keep <(grep CEPH1463 $path/bqc19-v4.0-raw2.fam | awk -F" " '{print $1,$2}') --freqx --missing --out $path/CEPH1463

plink --bed $path/bqc19-v4.0-raw2.bed --bim $path/bqc19-v4.0-qc1.bim --fam $path/bqc19-v4.0-raw2.fam \
--keep <(grep NA24385 $path/bqc19-v4.0-raw2.fam | awk -F" " '{print $1,$2}') --freqx --missing --out $path/NA24385

