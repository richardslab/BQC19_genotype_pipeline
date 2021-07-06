# Download scripts
wget https://imputation.sanger.ac.uk/www/plink2ensembl.txt
wget https://www.well.ox.ac.uk//~wrayner/tools/HRC-1000G-check-bim-v4.3.0.zip
unzip HRC-1000G-check-bim-v4.3.0.zip
# Download dbSNP
curl 'https://bravo.sph.umich.edu/freeze5/hg38/download/all' -H 'Accept-Encoding: gzip, deflate, br' -H 'Cookie: _ga=GA1.2.2125178096.1598536256; _gat_gtag_UA_73910830_2=1; _gid=GA1.2.918640196.1598536256; remember_token="tomoco_n@kuhp.kyoto-u.ac.jp|3fe2cac6c8ac02e2554720148800da830285c0692e93aa6f3fd415f06a997577f5d6ca0b69608e9bbf47f6c2d62453613322601c0a8185276a0e72679a97ebb3"' --compressed > bravo-dbsnp-all.vcf.gz
# VCF can be converted to an HRC formatted reference legend
wget https://www.well.ox.ac.uk//~wrayner/tools/CreateTOPMed.zip
unzip CreateTOPMed.zip
perl CreateTOPMed.pl -i bravo-dbsnp-all.vcf.gz -o ALL.TOPMed_freeze5_hg38_dbSNP.PASS.sites.tab.gz
bgzip -d ALL.TOPMed_freeze5_hg38_dbSNP.PASS.sites.tab.gz
