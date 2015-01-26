#!/bin/bash

#Data/Input/TCGA_tags.txt
fileList=$1

function launchQ {

	source ~/.bashrc

	thisTag="$1"_$2
	prevTag="$1"_$3

	echo '#!/bin/sh' >$thisTag.sh
	echo '# SmartAS launch' >>$thisTag.sh
	echo '#$ -q sbi' >>$thisTag.sh
	echo '#$ -cwd' >>$thisTag.sh
	echo "#$ -e /sbi/users/hectorc/SmartAS_experimental/esmartas_$thisTag.txt" >>$thisTag.sh
	echo "#$ -o /sbi/users/hectorc/SmartAS_experimental/osmartas_$thisTag.txt" >>$thisTag.sh
	echo "#$ -V" >>$thisTag.sh
	echo "#$ -N $thisTag" >>$thisTag.sh

	echo "/sbi/users/hectorc/SmartAS_experimental/Pipeline/gSmartAS.py -f $thisTag.cfg" >>$thisTag.sh
	echo "tar -zcvf $fullTag.tar.gz testResults/TCGA/$fullTag" >>$thisTag.sh

	qsub -hold_jid $prevTag $thisTag.sh

}

function printFile {

	fullTag=$1
	action=$2
	cfgFile="$fullTag"_$action.cfg
	cancerTag=`echo $fullTag | cut -d'-' -f1`

	echo initial-step=$action  >$cfgFile
	echo tag=$fullTag >>$cfgFile
	echo specific-drivers=Data/TCGA/specificDrivers/"$cancerTag"Drivers.txt >>$cfgFile
	echo unpaired-replicates=Yes >>$cfgFile
	echo specific-drivers=Data/"$fullTag"Drivers.txt >>$cfgFile
	echo working-directory=/sbi/users/hectorc/SmartAS_experimental >>$cfgFile

}

while read fullTag
do

	echo $fullTag 
	printFile $fullTag get-switches
	launchQ $fullTag get-switches

	# printFile $fullTag get-relevant-switches
	# launchQ $fullTag get-relevant-switches get-switches

 	# printFile $fullTag neighborhood-analysis
 	# launchQ $fullTag neighborhood-analysis get-relevant-switches

 	# printFile $fullTag experimental-network-analysis
 	# launchQ $fullTag experimental-network-analysis neighborhood-analysis

done < "$fileList"
