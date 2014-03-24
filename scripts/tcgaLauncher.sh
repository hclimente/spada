#!/bin/bash

fileList=$1
origin=$2

while read line
do
	echo === $line ===
	if [[ $origin -eq "ext" ]]; then
		if [[ ! -f Data/TCGA/External/"$line"_expressedGenes.lst ]]; then
			echo Data/TCGA/External/"$line"_expressedGenes.lst
			continue
		fi
		cat michal.cfg >$line.cfg
   		echo tag1=$line >>$line.cfg
		echo external=Data/TCGA/External/$line >>$line.cfg
   		Pipeline/SmartAS.py -f $line.cfg 2>&1

   		rm $line.cfg
    else
		Pipeline/scripts/standarizeInput.py TCGA $line
		echo Input standarized
		Pipeline/SmartAS.py -f Data/Input/TCGA/$line/config.cfg 2>&1
	fi

done < "$fileList"

#Data/Input/TCGA_tags.txt