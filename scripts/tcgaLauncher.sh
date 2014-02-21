#!/bin/bash

while read tag
do
	echo === $tag ===
	Pipeline/standarizeInput.py TCGA $tag
	echo Input standarized
	Pipeline/SmartAS.py -f Data/Input/TCGA/$tag/config.cfg &>$tag.log &
	wait $!

done < Data/Input/TCGA_tags.txt