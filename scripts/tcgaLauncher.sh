#!/bin/bash

while read tag
do
	echo === $tag ===
	Pipeline/scripts/standarizeInput.py TCGA $tag
	echo Input standarized
	Pipeline/SmartAS.py -f Data/Input/TCGA/$tag/config.cfg

done < Data/Input/TCGA_tags.txt
