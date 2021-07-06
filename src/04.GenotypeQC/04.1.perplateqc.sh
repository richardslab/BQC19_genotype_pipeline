path=/home/richards/tomoko.nakanishi/scratch/09.COVID19/05.BQC/01.genotypeQC/v4.0

paste $path/bqc-b{A,B,C,D,E,F,G,H,I,J,K,L,M,N,O,P,Q,R,S,T,U,V,W,X,Y,Z,a,b,c,d}-qc3-EUR.frq \
| tail -n +2 \
| awk 'BEGIN { OFS=" "; print "CHR SNP A1 A2 A.MAF B.MAF C.MAF D.MAF E.MAF F.MAF G.MAF H.MAF I.MAF J.MAF K.MAF L.MAF M.MAF N.MAF O.MAF P.MAF Q.MAF R.MAF S.MAF T.MAF U.MAF V.MAF W.MAF X.MAF Y.MAF Z.MAF a.MAF b.MAF c.MAF d.MAF"} \
{ print $1,$2,$3,$4,$5,$11,$17,$23,$29,$35,$41,$47,$53,$59,$65,$71,$77,$83,$89,$95,$101,$107,$113,$119,$125,$131,$137,$143,$149,$155,$161,$167,$173,$179,$185,$191 }' \
> $path/bqc-bAtod-qc3-EUR.frq

cut -f23 $path/bqc-bAtod-qc3-EUR.frq -d " " | head

for b in {A..Z} {a..d}; {
PREFIX=$path/bqc-b${b}-qc3-EUR
ls ${PREFIX}.bed ${PREFIX}.bim ${PREFIX}.fam
# All
plink \
--bfile ${PREFIX} \
--exclude <(cat $path/het.fisher.remove.b${b}.snp_id $path/hwe.remove.b${b}.snp_id) \
--make-bed \
--out $path/bqc-b${b}-qc4-EUR
}

cut -f1 $path/het.fisher.remove.b*.snp_id | sort | uniq | wc -l #96
cut -f1 $path/hwe.remove.b*.snp_id | sort | uniq | wc -l #11

rm $path/mergelist 
for b in {A..Z} {a..d};{ 
PREFIX=$path/bqc-b${b}-qc4-EUR; echo $PREFIX >> $path/mergelist 
}
plink --merge-list $path/mergelist --make-bed --out $path/bqc-v4.0-qc4-EUR
