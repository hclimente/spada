# # structural
# grep -v ^Cancer ~/testResults/TCGA/*/result_summary/structural_summary_onlyModels.tsv | cut -d':' -f2- >~/temp/structural_summary$tag.onlyModels.tsv
# grep -v ^Cancer ~/testResults/TCGA/*/result_summary/structural_loops_onlyModels.tsv | cut -d':' -f2- >~/temp/structural_loops$tag.onlyModels.tsv

# echo -e "Patient\tCancer" >~/temp/hallmark_info/patient_annotation.txt
# grep -v Patient ~/temp/hallmark_info/HALLMARK_DNA_REPAIR.tsv | cut -f1,4 | sort -k2 | uniq >>~/temp/hallmark_info/patient_annotation.txt

function getVersion {

    analysis_type=$1
    copyFlag=""

    if [ "`ls /projects_rg/TCGA/users/hector/SmartAS/$analysis_type`" == "" ]
    then
        last=0
        version=1
    else
        last=`ls /projects_rg/TCGA/users/hector/SmartAS/$analysis_type | sed 's/v//' | sort -n | uniq | tail -n1`
        version=$(($last + 1))
    fi

}

function copyFile {

    analysis_type=$1
    f=$2

    dest=/projects_rg/TCGA/users/hector/SmartAS/$analysis_type/v$version

    if [ ! -d "$dest" ]; then
        mkdir $dest
        mkdir $dest/figures
        mkdir $dest/tables
        rm ~/TCGA_analysis/$analysis_type
        ln -s $dest ~/TCGA_analysis/$analysis_type
    fi

    # if the filename includes a folder substructure, copy it
    if [[ "`echo $f | sed 's/[^\/]\+//g'`" != "" ]]
        then
        finalDest=$dest/`echo $f | sed 's/[^\/]\+$//'`
        if [ ! -d "$finalDest" ]; then
            mkdir -p $finalDest
        fi
    else
        finalDest=$dest
    fi

    cp ~/temp/$f $finalDest

}

function checkFile {

    analysis_type=$1
    f=$2

    if [ -e /projects_rg/TCGA/users/hector/SmartAS/$analysis_type/v$last/$f ]
        then
        diffOut=`diff /projects_rg/TCGA/users/hector/SmartAS/$analysis_type/v$last/$f ~/temp/$f`
    else
        diffOut="UNEXISTANT"
    fi

    if [[ "$last" == "0" || "$diffOut" != "" ]]
        then
        copyFlag="TRUE"
    fi

}

cancerTypes='brca coad hnsc kich kirc kirp lihc luad lusc prad thca'
mkdir ~/temp

########### SWITCHES GENERAL INFO ###########

getVersion switches 

# candidateList
for knsur in $cancerTypes
do
    cp ~/testResults/TCGA/$knsur/candidateList_smartas.tsv ~/temp/$knsur.candidateList.tsv
    checkFile switches $knsur.candidateList.tsv
done

# switches.tsv
analyses='CDS_change CDS_study UTR_change'
for a in $analyses
do
    cat ~/testResults/TCGA/????/result_summary/switches_onlyModels.tsv | head -n1 >~/temp/$a.tsv
    grep $a ~/testResults/TCGA/????/result_summary/switches_onlyModels.tsv | cut -d':' -f2 >>~/temp/$a.tsv
    checkFile switches $a.tsv
done

#analyses='Driver_D0_enrichment Driver_D0_patients Driver_D0_patients_relevant Driver_D1_enrichment Driver_D1_patients Driver_D1_patients_relevant Driver_relevance_enrichment Driver_d1_relevance_enrichment'
analyses="d0_enrichment d0_relevant_enrichment d1_enrichment d1_relevant_enrichment"
for a in $analyses
do
    grep $a ~/testResults/TCGA/????/result_summary/switches_onlyModels.tsv | cut -d':' -f2 >~/temp/$a.tsv
    checkFile switches $a.tsv
done

# exons.tsv
echo -e "Cancer\tRandom\tSwitch\tOrigin\tType\tLength\tCDSLength\tCDSRelativeSize\tPosition\tKeepOrf" >~/temp/exons.tsv
ls ~/testResults/TCGA/????/result_summary/exons_onlyModels.tsv | xargs -n 1 tail -n +2 >>~/temp/exons.tsv
checkFile switches exons.tsv

cut -f1,2 ~/temp/exons.tsv | sort | uniq -c | sed 's/^ \+//' | sed 's/ /\t/g' >~/temp/exonsPerSwitch.tsv
checkFile switches exonsPerSwitch.tsv

echo -e "Cancer\tRandom\tGene\tSymbol\tnTranscript\ttTranscript\tTag\tOrfChange\tnormalSegment\ttumorSegment" >~/temp/exons_new.tsv
grep -v ^Cancer ~/testResults/TCGA/????/result_summary/exons_new_onlyModels.tsv | cut -d':' -f2 >>~/temp/exons_new.tsv
checkFile switches exons_new.tsv

# isoform length
echo -e "Cancer\tRandom\tnIsoLength\ttIsoLength\tnIsoSpecificLength\ttIsoSpecificLength" >~/temp/isoform_length.tsv
ls ~/testResults/TCGA/????/result_summary/isoform_length_onlyModels.tsv | xargs -n 1 tail -n +2 >>~/temp/isoform_length.tsv
checkFile switches isoform_length.tsv

if [[ "$copyFlag" != "" ]]
    then
    for knsur in $cancerTypes
    do
        copyFile switches $knsur.candidateList.tsv
    done

    # switches.tsv
    analyses='CDS_change CDS_study UTR_change'
    for a in $analyses
    do
        copyFile switches $a.tsv
    done

    analyses="d0_enrichment d0_relevant_enrichment d1_enrichment d1_relevant_enrichment"
    for a in $analyses
    do
        copyFile switches $a.tsv
    done

    # exons.tsv
    copyFile switches exons.tsv
    copyFile switches exonsPerSwitch.tsv
    copyFile switches exons_new.tsv
    # isoform length
    copyFile switches isoform_length.tsv
fi

########### STRUCTURAL ANALYSIS ###########

analyses='anchor_analysis iupred_analysis interpro_analysis prosite_analysis structural_summary'
getVersion structural_analysis 

for knsur in $cancerTypes
do
    for a in $analyses
    do
        cp ~/testResults/TCGA/$knsur/$analysis_type/"$a"_random.tsv ~/temp/$knsur."$a"_random.tsv
        cp ~/testResults/TCGA/$knsur/$analysis_type/$a.tsv ~/temp/$knsur.$a.tsv
        checkFile structural_analysis $knsur."$a"_random.tsv
        checkFile structural_analysis $knsur.$a.tsv
    done
done

echo -e "Cancer\tSymbol\tAnnotation\tPartner\tPartnerAnnotation\tWhatsHappening\tInteractionAffection\tSequenceCover\tPartnerCover\tGene\tnormalTranscript\ttumorTranscript\tUniprot\tPartnerUniprot\tSequenceInformation\tIsoformSpecific\tPDB\tpymolCommands"  >~/temp/i3d_broken.tsv
grep -v ^Cancer ~/testResults/TCGA/????/i3d/i3d_broken.tsv | cut -d':' -f2 >>~/temp/i3d_broken.tsv
checkFile structural_analysis i3d_broken.tsv

echo -e "Cancer\tGene\tSymbol\tnTx\ttTx\tRandom\tAnalysis\tWhatsHappenning\tFeature\tDriver\tASDriver\tDriverType" >~/temp/structural_features.onlyModels.tsv
grep -v ^Cancer ~/testResults/TCGA/????/result_summary/structural_features_onlyModels.tsv | cut -d':' -f2- >>~/temp/structural_features.onlyModels.tsv
checkFile structural_analysis structural_features.onlyModels.tsv

echo -e "Gene\tSymbol\tNormalTranscript\tTumorTranscript\tWhat\tFeature\tnormalReps\ttumorReps\tnMacroScore\tnMicroScore\tnJaccard\ttMacroScore\ttMicroScore\ttJaccard" >~/temp/interpro_analysis.tsv
grep -v ^Gene ~/testResults/TCGA/????/structural_analysis/interpro_analysis.tsv | cut -d':' -f2- | sort | uniq >>~/temp/interpro_analysis.tsv
checkFile structural_analysis interpro_analysis.tsv

echo -e "Gene\tSymbol\tNormalTranscript\tTumorTranscript\tWhat\tFeature\tnormalReps\ttumorReps\tnMacroScore\tnMicroScore\tnJaccard\ttMacroScore\ttMicroScore\ttJaccard" >~/temp/prosite_analysis.tsv
grep -v ^Gene ~/testResults/TCGA/????/structural_analysis/prosite_analysis.tsv | cut -d':' -f2- | sort | uniq >>~/temp/prosite_analysis.tsv
checkFile structural_analysis prosite_analysis.tsv

if [[ "$copyFlag" != "" ]]
    then

    dest=/projects_rg/TCGA/users/hector/SmartAS/structural_analysis/v$version

    mkdir $dest
    mkdir $dest/figures
    mkdir $dest/tables
    rm ~/TCGA_analysis/structural_analysis
    ln -s $dest ~/TCGA_analysis/structural_analysis
    ln -s /projects_rg/TCGA/users/hector/SmartAS/EduardPorta-Interactions/ $dest
    
    for knsur in $cancerTypes
    do
        for a in $analyses
        do
            copyFile structural_analysis $knsur."$a"_random.tsv
            copyFile structural_analysis $knsur.$a.tsv
        done
    done
    copyFile structural_analysis i3d_broken.tsv
    copyFile structural_analysis structural_features.onlyModels.tsv
    copyFile structural_analysis interpro_analysis.tsv
    copyFile structural_analysis prosite_analysis.tsv
fi

########### NEIGHBORHOOD ANALYSIS ###########

analyses='canonical_pathways hallmarks go_biological_process oncogenic_signatures'
genesubgroups="all relevant"
getVersion neighborhood_analysis 

for a in $analyses
do
    for g in $genesubgroups
    do
        echo -e "GeneSet\tCancer\tpval\tqval\tNormalizedAverageAffection\tSwitchingGenes\tOR\tsg\tsng\tnsg\tnsng" >~/temp/"$a"_$g.txt
        grep -v ^GeneSet testResults/TCGA/????/neighborhood_analysis/"$a"_"$g"_onlyModels.txt  | cut -d':' -f2- >>~/temp/"$a"_$g.txt
        checkFile neighborhood_analysis "$a"_$g.txt
    done
done

hallmark_files=`ls testResults/TCGA/????/mutations/hallmark_info/*tsv | cut -d'/' -f6 | sort | uniq`
mkdir ~/temp/hallmark_info

for hll_file in $hallmark_files
do
    echo -e "Patient\tGene\tState\tCancer" >~/temp/hallmark_info/$hll_file
    for knsur in $cancerTypes
    do
        grep -v ^patient testResults/TCGA/$knsur/mutations/hallmark_info/$hll_file | cut -d':' -f2- | awk -v knsur=$knsur 'BEGIN{OFS="\t";} {print $0,knsur }' >>~/temp/hallmark_info/$hll_file
    done
    checkFile neighborhood_analysis hallmark_info/$hll_file
done

if [[ "$copyFlag" != "" ]]
    then
    for a in $analyses
    do
        for g in $genesubgroups
        do
            copyFile neighborhood_analysis "$a"_$g.txt
        done
    done

    for hll_file in $hallmark_files
    do
        copyFile neighborhood_analysis hallmark_info/$hll_file
    done
fi

########### MUTATIONS ###########
subsetTypes='all functional'

getVersion mutations

for s in $subsetTypes
do
    for t in $subsetTypes
    do
        echo -e "Cancer\tGene\tSymbol\tnTx\ttTx\tms\tm\ts\tn\tp_me\tpadj_me\tp_o\tpadj_o\tJaccard" >~/temp/gene_"$s"_mutations_"$t"_switches.txt
        cat ~/testResults/TCGA/????/mutations/gene_"$s"_mutations_"$t"_switches.txt >>~/temp/gene_"$s"_mutations_"$t"_switches.txt
        checkFile mutations gene_"$s"_mutations_"$t"_switches.txt

        echo -e "Cancer\tGene\tSymbol\tnTx\ttTx\tHallmark\tms\tm\ts\tn\tp\tpadj\tGeneset_genes" >~/temp/geneset_"$s"_mutations_"$t"_switches.txt
        cat ~/testResults/TCGA/????/mutations/geneset_"$s"_mutations_"$t"_switches.txt >>~/temp/geneset_"$s"_mutations_"$t"_switches.txt
        checkFile mutations geneset_"$s"_mutations_"$t"_switches.txt

        echo -e "Cancer\tGene\tSymbol\tnTx\ttTx\tHallmark\tms\tm\ts\tn\tp\tpadj\tGeneset_genes" >~/temp/geneset_"$s"_mutations_"$t"_switches_onlyDrivers.txt
        cat ~/testResults/TCGA/????/mutations/geneset_"$s"_mutations_"$t"_switches_onlyDrivers.txt >>~/temp/geneset_"$s"_mutations_"$t"_switches_onlyDrivers.txt
        checkFile mutations geneset_"$s"_mutations_"$t"_switches_onlyDrivers.txt

        echo -e "Cancer\tGene\tSymbol\tnTx\ttTx\tms\tm\ts\tn\tp\tpadj" >~/temp/pannegative_"$s"_mutations_"$t"_switches.txt
        cat ~/testResults/TCGA/????/mutations/pannegative_"$s"_mutations_"$t"_switches.txt >>~/temp/pannegative_"$s"_mutations_"$t"_switches.txt
        checkFile mutations pannegative_"$s"_mutations_"$t"_switches.txt

    done
done

# proteome features
echo -e "Cancer\tGene\tSymbol\tTranscript\tAnalysis\tFeature\tn\tFeatureLength\tStart\tEnd" >~/temp/proteome_features.txt
grep -v ^Cancer testResults/TCGA/????/mutations/proteome_features.txt | cut -d':' -f2- >>~/temp/proteome_features.txt
checkFile mutations proteome_features.txt

# proteome mutations
echo -e "Cancer\tGene\tSymbol\tTranscript\tAnalysis\tFeature\tn\tType\tPatient" >~/temp/proteome_mutations.txt
grep -v ^Cancer testResults/TCGA/????/mutations/proteome_mutations.txt | cut -d':' -f2- >>~/temp/proteome_mutations.txt
checkFile mutations proteome_mutations.txt

# proteome information
echo -e "Cancer\tGene\tSymbol\tTranscript\tTPM\tProteinLength" >~/temp/proteome_information.txt
grep -v ^Cancer testResults/TCGA/????/mutations/proteome_information.txt | cut -d':' -f2- >>~/temp/proteome_information.txt
checkFile mutations proteome_information.txt

# switches features
echo -e "Cancer\tGene\tSymbol\tTranscript\tAnalysis\tFeature\tn\tFeatureLength\tStart\tEnd" >~/temp/switch_features.txt
grep -v ^Cancer testResults/TCGA/????/mutations/switch_features.txt | cut -d':' -f2- >>~/temp/switch_features.txt
checkFile mutations switch_features.txt

# switches mutations
echo -e "Cancer\tGene\tSymbol\tTranscript\tAnalysis\tFeature\tn\tType\tPatient" >~/temp/switch_mutations.txt
grep -v ^Cancer testResults/TCGA/????/mutations/switch_mutations.txt | cut -d':' -f2- >>~/temp/switch_mutations.txt
checkFile mutations switch_mutations.txt

# switches information
echo -e "Cancer\tGene\tSymbol\tTranscript\tTPM\tProteinLength" >~/temp/switch_information.txt
grep -v ^Cancer testResults/TCGA/????/mutations/switch_information.txt | cut -d':' -f2- >>~/temp/switch_information.txt
checkFile mutations switch_information.txt

if [[ "$copyFlag" != "" ]]

    then
    for s in $subsetTypes
    do
        for t in $subsetTypes
        do
            copyFile mutations gene_"$s"_mutations_"$t"_switches.txt
            copyFile mutations geneset_"$s"_mutations_"$t"_switches.txt
            copyFile mutations geneset_"$s"_mutations_"$t"_switches_onlyDrivers.txt
            copyFile mutations pannegative_"$s"_mutations_"$t"_switches.txt
        done
    done
    
    copyFile mutations proteome_features.txt
    copyFile mutations proteome_mutations.txt
    copyFile mutations proteome_information.txt
    copyFile mutations switch_features.txt
    copyFile mutations switch_mutations.txt
    copyFile mutations switch_information.txt

    dest=/projects_rg/TCGA/users/hector/SmartAS/mutations/v$version
    ln -s /projects_rg/TCGA/users/hector/SmartAS/comet/ $dest

fi

rm -r ~/temp

#Pipeline/scripts/PLOT_analyzeResults.R