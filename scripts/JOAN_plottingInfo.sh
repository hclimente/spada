gene=IGF2BP3
nTx=uc003swf.2
tTx=uc003swg.2

~/SmartAS/Pipeline/scripts/JOAN_plottingInfo.R $gene $nTx $tTx

cut -f7 ~/Desktop/plottingInfo_$gene.tsv | cut -d"|" -f2 >tmp.tsv
paste ~/Desktop/plottingInfo_$gene.tsv tmp.tsv >tmp2.tsv
awk -F'\t' 'BEGIN {OFS="\t"} {print $1,$2,$12,$3,$4,$5,$6,$8,$9,$10,$11}' tmp2.tsv >~/Desktop/plottingInfo_$gene.tsv 
rm tmp.tsv tmp2.tsv