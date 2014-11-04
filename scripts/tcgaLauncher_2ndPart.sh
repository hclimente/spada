#!/bin/bash

#Data/Input/TCGA_tags.txt
fileList=$1
initialStep=$2

while read fullTag
do

	echo $fullTag 
	cancerTag=`echo $fullTag | cut -d'-' -f1`

	echo paired

	echo initial-step=$initialStep >$fullTag.cfg
	echo minimum-expression=-1.0 >>$fullTag.cfg
	echo tag=$fullTag >>$fullTag.cfg
	echo specific-drivers=Data/"$fullTag"Drivers.txt >>$fullTag.cfg

	SmartAS.py -f $fullTag.cfg
	rm $fullTag.cfg
	continue

	echo unpaired

	echo initial-step=$initialStep >u_$fullTag.cfg
	echo minimum-expression=-1.0 >>u_$fullTag.cfg
	echo tag=$fullTag >>u_$fullTag.cfg
	echo specific-drivers=Data/TCGA/specificDrivers/"$cancerTag"Drivers.txt >>u_$fullTag.cfg
	echo unpaired-replicates=Yes >>u_$fullTag.cfg
	echo specific-drivers=Data/"$fullTag"Drivers.txt >>u_$fullTag.cfg	
	
	SmartAS.py -f u_$fullTag.cfg
	rm u_$fullTag.cfg

done < "$fileList"
