#!/bin/bash

#Data/Input/TCGA_tags.txt
fileList=$1
action=$2

function launchQ {

	source ~/.bashrc

	tag=$2
	echo '#!/bin/sh' >$tag.sh
	echo "# SmartAS $tag" >>$tag.sh
	echo '#$ -q normal' >>$tag.sh
	echo '#$ -cwd' >>$tag.sh
	echo "#$ -e esmartas_$tag.txt" >>$tag.sh
	echo "#$ -o osmartas_$tag.txt" >>$tag.sh
	echo "#$ -V" >>$tag.sh
	echo "#$ -N $tag" >>$tag.sh

	echo "Pipeline/gSmartAS.py -f $tag.cfg" >>$tag.sh

	qsub $tag.sh

}

while read fullTag
do

	echo $fullTag 
	cfgFile="$fullTag"_$action.cfg

	cancerTag=`echo $fullTag | cut -d'-' -f1`

	echo initial-step=$action  >$cfgFile
	echo tag=$fullTag >>$cfgFile
	echo unpaired-replicates=Yes >>$cfgFile
	echo working-directory=/data/users/hector >>$cfgFile

	launchQ $fullTag "$fullTag"_$action

done < "$fileList"
