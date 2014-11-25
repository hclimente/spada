#!/bin/bash

#Data/Input/TCGA_tags.txt
fileList=$1

function launchQ {

	source ~/.bashrc

	tag=$1

	echo '#!/bin/sh' >$tag.sh
	echo '# SmartAS launch' >>$tag.sh
	echo '#$ -q sbi' >>$tag.sh
	echo '#$ -cwd' >>$tag.sh
	echo "#$ -e /sbi/users/hectorc/SmartAS_createGeneNetwork/esmartas_$tag.txt" >>$tag.sh
	echo "#$ -o /sbi/users/hectorc/SmartAS_createGeneNetwork/osmartas_$tag.txt" >>$tag.sh
	echo "#$ -V" >>$tag.sh
	echo "#$ -N $tag" >>$tag.sh

	echo "/sbi/users/hectorc/SmartAS_createGeneNetwork/Pipeline/SmartAS.py -f $tag.cfg" >>$tag.sh

	qsub $tag.sh

}

while read fullTag
do

	echo $fullTag 
	cancerTag=`echo $fullTag | cut -d'-' -f1`

	#echo paired

	#echo initial-step=3  >$fullTag.cfg
	#echo minimum-expression=-1.0 >>$fullTag.cfg
	#echo tag=$fullTag >>$fullTag.cfg
	#echo specific-drivers=Data/"$fullTag"Drivers.txt >>$fullTag.cfg
	#echo working-directory=/sbi/users/hectorc/SmartAS_createGeneNetwork >>$fullTag.cfg

	#launchQ $fullTag &

	echo unpaired

	echo initial-step=3  >u_$fullTag.cfg
	echo minimum-expression=-1.0 >>u_$fullTag.cfg
	echo tag=$fullTag >>u_$fullTag.cfg
	echo specific-drivers=Data/TCGA/specificDrivers/"$cancerTag"Drivers.txt >>u_$fullTag.cfg
	echo unpaired-replicates=Yes >>u_$fullTag.cfg
	echo specific-drivers=Data/"$fullTag"Drivers.txt >>u_$fullTag.cfg
	echo working-directory=/sbi/users/hectorc/SmartAS_createGeneNetwork >>u_$fullTag.cfg
	
	launchQ u_$fullTag &

done < "$fileList"
