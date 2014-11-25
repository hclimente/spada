#!/bin/bash

#Data/Input/TCGA_tags.txt
fileList=$1

function launchQ {

	source ~/.bashrc

	thisTag=$1$2
	prevTag=$1$3

	echo '#!/bin/sh' >$thisTag.sh
	echo '# SmartAS launch' >>$thisTag.sh
	echo '#$ -q sbi' >>$thisTag.sh
	echo '#$ -cwd' >>$thisTag.sh
	echo "#$ -e /sbi/users/hectorc/SmartAS_experimental/esmartas_$thisTag.txt" >>$thisTag.sh
	echo "#$ -o /sbi/users/hectorc/SmartAS_experimental/osmartas_$thisTag.txt" >>$thisTag.sh
	echo "#$ -V" >>$thisTag.sh
	echo "#$ -N $thisTag" >>$thisTag.sh

	echo "/sbi/users/hectorc/SmartAS_experimental/Pipeline/SmartAS.py -f $thisTag.cfg" >>$thisTag.sh

	qsub -hold_jid $prevTag $thisTag.sh

}

function printFile {

	fullTag = $1
	action = $2
	cancerTag=`echo $fullTag | cut -d'-' -f1`

	echo initial-step=$action  >$fullTag$action.cfg
	echo tag=$fullTag >>$fullTag$action.cfg
	echo specific-drivers=Data/TCGA/specificDrivers/"$cancerTag"Drivers.txt >>$fullTag$action.cfg
	echo unpaired-replicates=Yes >>$fullTag$action.cfg
	echo specific-drivers=Data/"$fullTag"Drivers.txt >>$fullTag$action.cfg
	echo working-directory=/sbi/users/hectorc/SmartAS_experimental >>$fullTag$action.cfg

}

while read fullTag
do

	echo $fullTag 
	printFile $fullTag get-switches
	launchQ $fullTag get-switches

	printFile $fullTag get-relevant-switches
	launchQ $fullTag get-relevant-switches get-switches

    printFile $fullTag predicted-network-analysis
    launchQ $fullTag predicted-network-analysis get-relevant-switches

    printFile $fullTag neighborhood-analysis
    launchQ $fullTag neighborhood-analysis predicted-network-analysis

    printFile $fullTag experimental-network-analysis
    launchQ $fullTag experimental-network-analysis neighborhood-analysis

done < "$fileList"
