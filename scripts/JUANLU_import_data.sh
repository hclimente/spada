# put files in correct format
for file in `ls *_SmartAS.txt`
do
    newFile=`echo $file | sed 's/\.txt/_new\.txt/'`
    echo -e "Gene_name\ttIsoform\tnIsoform" >$newFile
    sed 1d $file | while read line
    do
        geneName=`echo -e $line | sed 's/ .*//' `
        echo $line | sed "s/ $geneName[#;]/\t/g" >>$newFile
    done
done

# after scp the folder, import data into smartas
for file in `ls juanluAnalyses5/*_pairs_isoforms.txt`
do

    knsur=`echo $file | sed 's/.\+\///' | sed 's/_.*//' | tr '[:upper:]' '[:lower:]'`
    part=`echo $file  | sed 's/.\+\///' | sed 's/[^_]\+_//' | sed 's/_.\+//' ` 

    cfgFile=`echo juanluAnalyses5/$knsur"_"$part".cfg"`

    mkdir -p testResults/TCGA/"$knsur"_$part/logs

    echo "initial-step=get-switches" >$cfgFile
    echo "unpaired-replicates=Yes" >>$cfgFile
    echo "tag="$knsur"_"$part >>$cfgFile
    echo "specific-drivers=Data/TCGA/specificDrivers/"$knsur"Drivers.txt" >>$cfgFile
    echo "working-directory=/data/users/hector" >>$cfgFile
    echo "external-switches="$file >>$cfgFile
    if [ "$knsur" = "ov" -o "$knsur" = "skcm" -o "$knsur" = "er+" -o "$knsur" = "er-" -o "$knsur" = "mitf-" -o "$knsur" = "mitf+" ]
    then
        echo "parent-tag=brca" >>$cfgFile
    else
        echo "parent-tag="$knsur >>$cfgFile
    fi

    Pipeline/SmartAS.py -f $cfgFile

done

# run structural analysis (x2)
for file in `ls juanluAnalyses5/*_pairs_isoforms.txt`
do

    knsur=`echo $file | sed 's/.\+\///' | sed 's/_.*//' | tr '[:upper:]' '[:lower:]'`
    part=`echo $file  | sed 's/.\+\///' | sed 's/[^_]\+_//' | sed 's/_.\+//' ` 

    cfgFile=`echo juanluAnalyses5/$knsur"_"$part"_structural_analysis.cfg"`

    echo "initial-step=get-relevant-switches" >$cfgFile
    echo "unpaired-replicates=Yes" >>$cfgFile
    echo "tag="$knsur"_"$part >>$cfgFile
    echo "specific-drivers=Data/TCGA/specificDrivers/"$knsur"Drivers.txt" >>$cfgFile
    echo "working-directory=/data/users/hector" >>$cfgFile

    Pipeline/SmartAS.py -f $cfgFile

done

# check analyses
mkdir -p /projects_rg/TCGA/users/hector/juanluAnalysis/v5/
for a in `echo 'prosite interpro iupred anchor'`
do
    grep Gene testResults/TCGA/*_*/structural_analysis/"$a"_analysis.tsv | head -n1 | sed 's/testResults.\+tsv://' >/projects_rg/TCGA/users/hector/juanluAnalysis/v5/"$a"_analysis.tsv
    grep -v Nothing testResults/TCGA/*_*/structural_analysis/"$a"_analysis.tsv | grep -v tsv:Gene | sed 's/testResults\/TCGA\///' | sed 's/\/structural_analysis\/'$a'_analysis\.tsv//' >>/projects_rg/TCGA/users/hector/juanluAnalysis/v5/"$a"_analysis.tsv
done

# get-i3d-broken-interactions
ls testResults/TCGA | grep _ >juanlu_tags.txt
Pipeline/scripts/tcgaLauncherGenCluster_2ndPart.sh juanlu_tags.txt get-i3d-broken-interactions

# copy candidate list
mkdir -p /projects_rg/TCGA/users/hector/juanluAnalysis/candidateList/v1/
while read a
do
    cp testResults/TCGA/$a/candidateList_v6.tsv /projects_rg/TCGA/users/hector/juanluAnalysis/candidateList/v1/$a.candidateList.tsv
done <juanlu_tags.txt