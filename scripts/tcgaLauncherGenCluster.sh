#!/bin/bash

#Data/Input/TCGA_tags.txt
fileList=$1

function getSpecificDrivers {
	cancerType=$1
	tissue=`python -c "dTissue = {\"urinary_tract\": [\"blca\",\"blca-epithelial\",\"blca-mesenchymal\"],\"breast\": [\"brca\",\"brca-basal\",\"brca-her2\",\"brca-luminala\",\"brca-luminalb\"],\"large_intestine\": [\"coad\",\"coad-hypermutated\",\"coad-hypomutated\"],\"upper_aerodigestive_tract\": [\"hnsc\"],\"kidney\": [\"kich\",\"kirc\",\"kirp\"],\"liver\": [\"lihc\"],\"lung\": [\"luad\",\"lusc\",\"lusc-basal\",\"lusc-classical\",\"lusc-primitive\",\"lusc-secretory\"],\"prostate\": [\"prad\"],\"\": [\"tcga\"],\"thyroid\": [\"thca\"]}; print [ x for x in dTissue if '$cancerType' in dTissue[x] ][0]"`
	subtissue=`python -c "dSubtissue = {\"bladder\": [\"blca\",\"blca-epithelial\",\"blca-mesenchymal\"],\"colon\": [\"coad\",\"coad-hypermutated\",\"coad-hypomutated\"],\"head_neck\": [\"hnsc\"],\"\": [\"brca\",\"brca-basal\",\"brca-her2\",\"brca-luminala\",\"brca-luminalb\",\"kich\",\"kirc\",\"kirp\",\"lihc\",\"luad\",\"lusc\",\"lusc-basal\",\"lusc-classical\",\"lusc-primitive\",\"lusc-secretory\",\"prad\",\"thca\",\"tcga\"]};print [ x for x in dSubtissue if '$cancerType' in dSubtissue[x] ][0]"`
	histology=`python -c "dHistology = {\"carcinoma\": [\"blca\",\"blca-epithelial\",\"blca-mesenchymal\",\"brca\",\"brca-basal\",\"brca-her2\",\"brca-luminala\",\"brca-luminalb\",\"hnsc\",\"kich\",\"kirc\",\"kirp\",\"lihc\",\"thca\"],\"\": [\"coad\",\"coad-hypermutated\",\"coad-hypomutated\",\"luad\",\"lusc\",\"lusc-basal\",\"lusc-classical\",\"lusc-primitive\",\"lusc-secretory\",\"prad\",\"tcga\"]};print [ x for x in dHistology if '$cancerType' in dHistology[x] ][0]"`
	subhistology=`python -c "dSubhistology = {\"\": [\"blca\",\"brca\",\"thca\",\"blca-epithelial\",\"blca-mesenchymal\",\"tcga\"],\"HER-positive_carcinoma\": [\"brca-her2\"],\"luminal_A_carcinoma\": [\"brca-luminala\"],\"luminal_B_carcinoma\": [\"brca-luminalb\"],\"basal_(triple-negative)_carcinoma\": [\"brca-basal\"],\"adenocarcinoma\": [\"coad\",\"coad-hypermutated\",\"coad-hypomutated\",\"luad\",\"prad\"],\"chromophobe_renal_cell_carcinoma\": [\"kich\"],\"clear_cell_renal_cell_carcinoma\": [\"kirc\"],\"papillary_renal_cell_carcinoma\": [\"kirp\"],\"hepatocellular_carcinoma\": [\"lihc\"],\"squamous_cell_carcinoma\": [\"hnsc\",\"lusc\",\"lusc-basal\",\"lusc-classical\",\"lusc-primitive\",\"lusc-secretory\"]};print [ x for x in dSubhistology if '$cancerType' in dSubhistology[x] ][0]"`
	awk '$27 == "primary"' /genomics/users/hector/SmartAS/Data/CosmicCompleteExport_v70_100814.tsv >kk.tmp
	if [[ $tissue -ne "" ]]
		then
		awk -v k=$tissue '$8==k' kk.tmp >kk1.tmp
	else
		mv kk.tmp kk1.tmp
	fi
	if [[ $subTissue -ne "" ]]
		then
		awk -v k=$subTissue '$8==k' kk1.tmp >kk2.tmp
	else
		mv kk1.tmp kk2.tmp
	fi
	if [[ $histology -ne "" ]]
		then
		awk -v k=$histology '$8==k' kk2.tmp >kk3.tmp
	else
		mv kk2.tmp kk3.tmp
	fi
	if [[ $subhistology -ne "" ]]
		then
		awk -v k=$subhistology '$8==k' kk3.tmp >kk4.tmp
	else
		mv kk3.tmp kk4.tmp
	fi

	cut -f1  kk4.tmp | sed 's/_ENST.\+//' | sort | uniq -c | sort -nr | head -n30 | sed 's/^[^A-Z]\+//' >Data/TCGA/specificDrivers/"$cancerType"Drivers.txt

	rm kk4.tmp 
}

function launchQ {

	source ~/.bashrc

	tag=$1
	echo '#!/bin/sh' >$tag.sh
	echo '# SmartAS import data' >>$tag.sh
	echo '#$ -q normal' >>$tag.sh
	echo '#$ -cwd' >>$tag.sh
	echo "#$ -e /data/users/hector/esmartas_$tag.1.txt" >>$tag.sh
	echo "#$ -o /data/users/hector/osmartas_$tag.1.txt" >>$tag.sh

	echo "/data/users/hector/Pipeline/SmartAS.py -f $tag.cfg" >>$tag.sh

	qsub -V -N $tag $tag.sh

	echo '#!/bin/sh' >$tag.2.sh
	echo '# SmartAS getCandidates' >>$tag.2.sh
	echo '#$ -q normal' >>$tag.2.sh
	echo '#$ -cwd' >>$tag.2.sh
	echo "#$ -e /data/users/hector/esmartas_$tag.2.txt" >>$tag.2.sh
	echo "#$ -o /data/users/hector/osmartas_$tag.2.txt" >>$tag.2.sh

	echo "/data/users/hector/Pipeline/SmartAS.py -f $tag.cfg" >>$tag.2.sh

	qsub -V -N $tag.2 -hold_jid $tag  $tag.2.sh

}

while read fullTag
do

	echo $fullTag 
	cancerTag=`echo $fullTag | cut -d'-' -f1`
	getSpecificDrivers $fullTag

	echo paired
	mkdir -p Results/TCGA/"$fullTag"_mE-1.0

	echo initial-step=0  >$fullTag.cfg
	echo minimum-expression=-1.0 >>$fullTag.cfg
	echo tag=$fullTag >>$fullTag.cfg
	echo specific-drivers=Data/"$fullTag"Drivers.txt >>$fullTag.cfg
	echo working-directory=/data/users/hector >>$fullTag.cfg

	launchQ $fullTag &

	echo unpaired
	mkdir -p Results/TCGA/u_"$fullTag"_mE-1.0

	echo initial-step=0  >u_$fullTag.cfg
	echo minimum-expression=-1.0 >>u_$fullTag.cfg
	echo tag=$fullTag >>u_$fullTag.cfg
	echo specific-drivers=Data/TCGA/specificDrivers/"$cancerTag"Drivers.txt >>u_$fullTag.cfg
	echo unpaired-replicates=Yes >>u_$fullTag.cfg
	echo specific-drivers=Data/"$fullTag"Drivers.txt >>u_$fullTag.cfg
	echo working-directory=/data/users/hector >>u_$fullTag.cfg
	
	launchQ u_$fullTag &

done < "$fileList"
