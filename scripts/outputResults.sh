#!/bin/sh

while read line
do
    parameter=`echo $line | cut -d"=" -f1`
    value=`echo $line | cut -d"=" -f2`

    if [ "$parameter" == "minCandidateExpression" ]
    	then
    	minExpresion=$value
    elif [ "$parameter" == "minPSI" ]
    	then
    	minPSI=$value
    fi
done < "Parameters.cfg"

resultsDir=~/testResults/GENECODE19
outFolder=e"$minExpresion"p"$minPSI"

mkdir -p $resultsDir/$outFolder

cp candidates_normal.gff $resultsDir/$outFolder/candidates_normal.e"$minExpresion"p"$minPSI".gff
cp candidates_tumor.gff $resultsDir/$outFolder/candidates_tumor.e"$minExpresion"p"$minPSI".gff
cp candidateInteractions.tsv $resultsDir/$outFolder/candidateInteractions.e"$minExpresion"p"$minPSI".tsv
cp candidates_normal.top.gff $resultsDir/$outFolder/candidates_normal.top.e"$minExpresion"p"$minPSI".gff
cp candidates_tumor.top.gff $resultsDir/$outFolder/candidates_tumor.top.e"$minExpresion"p"$minPSI".gff

cd $resultsDir/$outFolder/

tar -czvf $outFolder.tar.gz *