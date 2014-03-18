#!/bin/bash

while read tag
do
	echo === $tag ===
	Pipeline/scripts/standarizeInput.py TCGA $tag
	echo Input standarized
	Pipeline/SmartAS.py -f Data/Input/TCGA/$tag/config.cfg 2>&1

done < $1

#Data/Input/TCGA_tags.txt