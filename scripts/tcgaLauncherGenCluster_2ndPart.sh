#!/bin/bash

#Data/Input/TCGA_tags.txt
fileList=$1
action=$2
models=$3

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

	echo "Pipeline/SmartAS.py -f $tag.cfg" >>$tag.sh

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
	if [[ "$models" == "all-switches" ]]; then
		echo $models >>$cfgFile
	fi

	rm esmartas_"$fullTag"_$action.txt osmartas_"$fullTag"_$action.txt

	launchQ $fullTag "$fullTag"_$action
	#Pipeline/SmartAS.py -f "$fullTag"_$action.cfg

done < "$fileList"
