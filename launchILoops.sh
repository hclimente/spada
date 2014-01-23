#!/bin/bash

cd $1

counter=0
declare -a pidQueue

for transcript in `ls input/ | egrep ENST`
do
	for configFile in `ls input/$transcript | egrep *net`
	do

		echo "/soft/devel/python-2.7/bin/python /sbi/programs/iLoops_devel/iLoops.py \
		-f input/ExpressedTranscripts.fasta \
		-q input/$transcript/$configFile \
		-j output/$configFile \
		-x $configFile.xml \
		-g all \
		-n 25 \
		-Q sbi \
		-c 1,5,6,7,8,9,10,11,12,13,14,15,20,30,40,50 \
		-v" &

	 	pidQueue+=($!)

	done
	counter="1"
	if [[ counter -ge 1 ]]; then
		break
	fi
done

for job in ${pidQueue[*]}
do
	while :
	do
		finish=`ps --pid $job | grep -v "$job" | wc -l`
		if [[ finish -eq 0 ]]; then
			break
		else
			echo "Awaiting for completion of iLoops jobs."
			sleep 900
		fi
	done
done