#!/bin/bash

#Data/Input/TCGA_tags.txt
fileList=$1
initialStep=$2

function launchQ {

	source ~/.bashrc

	tag=$1
	echo '#!/bin/sh' >$tag.sh
	echo '# SmartAS import data' >>$tag.sh
	echo '#$ -q normal' >>$tag.sh
	echo '#$ -cwd' >>$tag.sh
	echo "#$ -e esmartas_$tag.txt" >>$tag.sh
	echo "#$ -o osmartas_$tag.txt" >>$tag.sh
	echo "#$ -V" >>$tag.sh
	echo "#$ -N $tag" >>$tag.sh

	echo "Pipeline/SmartAS.py -f $tag.cfg" >>$tag.sh

	qsub $tag.sh

}

while read fullTag
do

	echo $fullTag 
	cancerTag=`echo $fullTag | cut -d'-' -f1`

	echo paired

	echo initial-step=$initialStep  >$fullTag.cfg
	echo minimum-expression=-1.0 >>$fullTag.cfg
	echo tag=$fullTag >>$fullTag.cfg
	echo specific-drivers=Data/"$fullTag"Drivers.txt >>$fullTag.cfg
	echo working-directory=/data/users/hector >>$fullTag.cfg

	launchQ $fullTag &

	echo unpaired
	mkdir -p Results/TCGA/u_"$fullTag"_mE-1.0

	echo initial-step=$initialStep  >u_$fullTag.cfg
	echo minimum-expression=-1.0 >>u_$fullTag.cfg
	echo tag=$fullTag >>u_$fullTag.cfg
	echo specific-drivers=Data/TCGA/specificDrivers/"$cancerTag"Drivers.txt >>u_$fullTag.cfg
	echo unpaired-replicates=Yes >>u_$fullTag.cfg
	echo specific-drivers=Data/"$fullTag"Drivers.txt >>u_$fullTag.cfg
	echo working-directory=/data/users/hector >>u_$fullTag.cfg
	
	launchQ u_$fullTag &

done < "$fileList"
