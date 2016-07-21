# # structural
# grep -v ^Cancer ~/smartas/analyses/*/result_summary/structural_summary_onlyModels.tsv | cut -d':' -f2- >/structural_summary$tag.onlyModels.tsv
# grep -v ^Cancer ~/smartas/analyses/*/result_summary/structural_loops_onlyModels.tsv | cut -d':' -f2- >/structural_loops$tag.onlyModels.tsv

# echo -e "Patient\tCancer" >/hallmark_info/patient_annotation.txt
# grep -v Patient /hallmark_info/HALLMARK_DNA_REPAIR.tsv | cut -f1,4 | sort -k2 | uniq >>/hallmark_info/patient_annotation.txt

function getVersion {

    analysis_type=$1
    copyFlag=""

    if [ "`ls smartas/projects_rg/$analysis_type`" == "" ]
    then
        last=0
        version=1
    else
        last=`ls smartas/projects_rg/$analysis_type | sed 's/v//' | sort -n | uniq | tail -n1`
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
        rm ~/smartas/results/$analysis_type
        ln -s $dest ~/smartas/results/$analysis_type
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

    cp ~/smartas/notebook/data/$analysis_type/$f $finalDest

}

function checkFile {

    analysis_type=$1
    f=$2

    if [ -e smartas/projects_rg/$analysis_type/v$last/$f ]
        then
        diffOut=`diff smartas/projects_rg/$analysis_type/v$last/$f ~/smartas/notebook/data/$analysis_type/$f`
    else
        diffOut="UNEXISTANT"
    fi

    if [[ "$last" == "0" || "$diffOut" != "" ]]
        then
        copyFlag="TRUE"
    fi

}

cancerTypes='brca coad hnsc kich kirc kirp lihc luad lusc prad thca'

########### SWITCHES GENERAL INFO ###########
getVersion switches
rm -r ~/smartas/notebook/data/switches
mkdir ~/smartas/notebook/data/switches

# candidateList
for knsur in $cancerTypes
do
    cp ~/smartas/analyses/$knsur/candidateList_info.tsv ~/smartas/notebook/data/switches/$knsur.candidateList.tsv
    checkFile switches $knsur.candidateList.tsv
done

# switches.tsv
echo -e "Cancer\tAnalysis\tBoth\tOnly_nIso\tOnly_tIso\tNone\tRandom_Both\tRandom_Only_nIso\tRandom_Only_tIso\tRandom_None" >~/smartas/notebook/data/switches/$a.tsv
grep CDS_study ~/smartas/analyses/????/result_summary/switches_onlyModels.tsv | cut -d':' -f2 >>~/smartas/notebook/data/switches/CDS_study.tsv
checkFile switches CDS_study.tsv

analyses='CDS_change UTR_change'
for a in $analyses
do
    echo -e "Cancer\tAnalysis\tYes\tNo\tRandom_Yes\tRandom_No\tp\tOR" >~/smartas/notebook/data/switches/$a.tsv
    grep $a ~/smartas/analyses/????/result_summary/switches_onlyModels.tsv | cut -d':' -f2 >>~/smartas/notebook/data/switches/$a.tsv
    checkFile switches $a.tsv
done

analyses="d0_enrichment d1_enrichment d0_functional_enrichment d1_functional_enrichment"
for a in $analyses
do
    echo -e "Cancer\tAnalysis\tFS\tFNS\tNFS\tNFNS\tp.me\tOR" >~/smartas/notebook/data/switches/$a.tsv
    grep $a ~/smartas/analyses/????/result_summary/switches_onlyModels.tsv | cut -d':' -f2 >~/smartas/notebook/data/switches/$a.tsv
    checkFile switches $a.tsv
done

# exons.tsv
echo -e "Cancer\tRandom\tSwitch\tOrigin\tType\tLength\tCDSLength\tCDSRelativeSize\tPosition\tKeepOrf" >~/smartas/notebook/data/switches/exons.tsv
ls ~/smartas/analyses/????/result_summary/exons_onlyModels.tsv | xargs -n 1 tail -n +2 >>~/smartas/notebook/data/switches/exons.tsv
checkFile switches exons.tsv

cut -f1,2 ~/smartas/notebook/data/switches/exons.tsv | sort | uniq -c | sed 's/^ \+//' | sed 's/ /\t/g' >~/smartas/notebook/data/switches/exonsPerSwitch.tsv
checkFile switches exonsPerSwitch.tsv

echo -e "Cancer\tRandom\tGene\tSymbol\tnTranscript\ttTranscript\tTag\tOrfChange\tnormalSegment\ttumorSegment" >~/smartas/notebook/data/switches/exons_new.tsv
grep -v ^Cancer ~/smartas/analyses/????/result_summary/exons_new_onlyModels.tsv | cut -d':' -f2 >>~/smartas/notebook/data/switches/exons_new.tsv
checkFile switches exons_new.tsv

# isoform length
echo -e "Cancer\tRandom\tnIsoLength\ttIsoLength\tnIsoSpecificLength\ttIsoSpecificLength" >~/smartas/notebook/data/switches/isoform_length.tsv
ls ~/smartas/analyses/????/result_summary/isoform_length_onlyModels.tsv | xargs -n 1 tail -n +2 >>~/smartas/notebook/data/switches/isoform_length.tsv
checkFile switches isoform_length.tsv

# evidence of driverness
echo -e "Tumor\tGeneId\tSymbol\tNormal_transcript\tTumor_transcript\tRecurrence\tAffects_mutated_feature\tMutual_exclusion\tCoocurrence\tPPI" >~/smartas/notebook/data/switches/driverEvidence.tsv
grep -v ^Tumor ~/smartas/analyses/????/candidateList_driverEvidence.tsv | cut -d':' -f2 >>~/smartas/notebook/data/switches/driverEvidence.tsv
checkFile switches driverEvidence.tsv

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

    analyses="d0_enrichment d1_enrichment d0_functional_enrichment d1_functional_enrichment"
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

    # driverness
    copyFile switches driverEvidence.tsv

    #smartas/pipeline/figures/switches_plots.r

fi

########### STRUCTURAL ANALYSIS ###########
analyses='anchor_analysis iupred_analysis interpro_analysis prosite_analysis structural_summary'
getVersion structural_analysis 
rm -r ~/smartas/notebook/data/structural_analysis
mkdir ~/smartas/notebook/data/structural_analysis
cp -r /projects_rg/TCGA/users/hector/SmartAS/structural_analysis/pfam_go_term ~/smartas/notebook/data/structural_analysis/
cp -r /projects_rg/TCGA/users/hector/SmartAS/eporta/comet_output/ ~/smartas/notebook/data/structural_analysis/
cp /projects_rg/TCGA/users/hector/SmartAS/eporta/Switched_interactions_consensus.txt /projects_rg/TCGA/users/hector/SmartAS/eporta/num_patients_with_each_switch.txt ~/smartas/notebook/data/structural_analysis/

for knsur in $cancerTypes
do
    for a in $analyses
    do
        cp ~/smartas/analyses/$knsur/$analysis_type/"$a"_random.tsv ~/smartas/notebook/data/structural_analysis/$knsur."$a"_random.tsv
        cp ~/smartas/analyses/$knsur/$analysis_type/$a.tsv ~/smartas/notebook/data/structural_analysis/$knsur.$a.tsv
        checkFile structural_analysis $knsur."$a"_random.tsv
        checkFile structural_analysis $knsur.$a.tsv
    done
done

echo -e "Cancer\tGene\tSymbol\tnTx\ttTx\tRandom\tAnalysis\tWhatsHappenning\tFeature\tDriver\tASDriver\tDriverType" >~/smartas/notebook/data/structural_analysis/structural_features.onlyModels.tsv
grep -v ^Cancer ~/smartas/analyses/????/result_summary/structural_features_onlyModels.tsv | cut -d':' -f2- >>~/smartas/notebook/data/structural_analysis/structural_features.onlyModels.tsv
checkFile structural_analysis structural_features.onlyModels.tsv

echo -e "Gene\tSymbol\tNormalTranscript\tTumorTranscript\tWhat\tFeature\tnormalReps\ttumorReps\tnMacroScore\tnMicroScore\tnJaccard\ttMacroScore\ttMicroScore\ttJaccard" >~/smartas/notebook/data/structural_analysis/interpro_analysis.tsv
grep -v ^Gene ~/smartas/analyses/????/structural_analysis/interpro_analysis.tsv | cut -d':' -f2- | sort | uniq >>~/smartas/notebook/data/structural_analysis/interpro_analysis.tsv
checkFile structural_analysis interpro_analysis.tsv

echo -e "Gene\tSymbol\tNormalTranscript\tTumorTranscript\tWhat\tFeature\tnormalReps\ttumorReps\tnMacroScore\tnMicroScore\tnJaccard\ttMacroScore\ttMicroScore\ttJaccard" >~/smartas/notebook/data/structural_analysis/prosite_analysis.tsv
grep -v ^Gene ~/smartas/analyses/????/structural_analysis/prosite_analysis.tsv | cut -d':' -f2- | sort | uniq >>~/smartas/notebook/data/structural_analysis/prosite_analysis.tsv
checkFile structural_analysis prosite_analysis.tsv

if [[ "$copyFlag" != "" ]]
    then

    dest=/projects_rg/TCGA/users/hector/SmartAS/structural_analysis/v$version

    mkdir $dest
    mkdir $dest/figures
    mkdir $dest/tables
    rm ~/smartas/results/structural_analysis
    ln -s $dest ~/smartas/results/structural_analysis
    ln -s ~/smartas/projects_rg/EduardPorta-Interactions/ $dest
    
    for knsur in $cancerTypes
    do
        for a in $analyses
        do
            copyFile structural_analysis $knsur."$a"_random.tsv
            copyFile structural_analysis $knsur.$a.tsv
        done
    done

    copyFile structural_analysis structural_features.onlyModels.tsv
    copyFile structural_analysis interpro_analysis.tsv
    copyFile structural_analysis prosite_analysis.tsv

    smartas/pipeline/figures/structural_analysis.r
fi

########### NEIGHBORHOOD ANALYSIS ###########

analyses='canonical_pathways hallmarks go_biological_process oncogenic_signatures'
genesubgroups="all functional"
getVersion switches 

for a in $analyses
do
    for g in $genesubgroups
    do
        echo -e "GeneSet\tCancer\tpval\tqval\tNormalizedAverageAffection\tSwitchingGenes\tOR\tsg\tsng\tnsg\tnsng" >~/smartas/notebook/data/switches/"$a"_$g.txt
        grep -v ^GeneSet smartas/analyses/????/switches/"$a"_"$g"_onlyModels.txt  | cut -d':' -f2- >>~/smartas/notebook/data/switches/"$a"_$g.txt
        checkFile switches "$a"_$g.txt
    done
done

hallmark_files=`ls smartas/analyses/????/mutations/hallmark_info/*tsv | cut -d'/' -f6 | sort | uniq`
#mkdir ~/smartas/notebook/data/switches/hallmark_info

for hll_file in $hallmark_files
do
    echo -e "Patient\tGene\tState\tCancer" >~/smartas/notebook/data/switches/hallmark_info/$hll_file
    for knsur in $cancerTypes
    do
        grep -v ^patient smartas/analyses/$knsur/mutations/hallmark_info/$hll_file | cut -d':' -f2- | awk -v knsur=$knsur 'BEGIN{OFS="\t";} {print $0,knsur }' >>~/smartas/notebook/data/switches/hallmark_info/$hll_file
    done
    checkFile switches hallmark_info/$hll_file
done

if [[ "$copyFlag" != "" ]]
    then

    for a in $analyses
    do
        for g in $genesubgroups
        do
            copyFile switches "$a"_$g.txt
        done
    done

    for hll_file in $hallmark_files
    do
        copyFile switches hallmark_info/$hll_file
    done

    smartas/pipeline/figures/neighborhoods_plots.r
fi

########### MUTATIONS ###########
getVersion mutations
rm -r ~/smartas/notebook/data/mutations
mkdir ~/smartas/notebook/data/mutations

# WES ME analysis
echo -e "Tumor\tGeneId\tSymbol\tNormal_transcript\tTumor_transcript\tMS\tM\tS\tN\tp\tOR\tSwitched\tMutated" >~/smartas/notebook/data/mutations/gene_functional_mutations_all_switches.txt
grep -v ^Tumor smartas/analyses/????/mutations/gene_functional_mutations_all_switches.txt | cut -d':' -f2- >>~/smartas/notebook/data/mutations/gene_functional_mutations_all_switches.txt
checkFile mutations gene_functional_mutations_all_switches.txt

# WGS coocurrence analysis
echo -e "GeneId\tSymbol\tNormal_transcript\tTumor_transcript\tMS\tM\tS\tN\tp.o" >~/smartas/notebook/data/mutations/gene_wgs_mutations_all_switches.txt
grep -v ^Tumor smartas/analyses/????/mutations/gene_wgs_mutations_all_switches.txt | cut -d':' -f2- >>~/smartas/notebook/data/mutations/gene_wgs_mutations_all_switches.txt
checkFile mutations gene_wgs_mutations_all_switches.txt

# proteome features
echo -e "Cancer\tGene\tSymbol\tTranscript\tAnalysis\tFeature\tn\tFeatureLength\tStart\tEnd" >~/smartas/notebook/data/mutations/proteome_features.txt
grep -v ^Cancer smartas/analyses/????/mutations/proteome_features.txt | cut -d':' -f2- >>~/smartas/notebook/data/mutations/proteome_features.txt
checkFile mutations proteome_features.txt

# proteome mutations
echo -e "Cancer\tGene\tSymbol\tTranscript\tAnalysis\tFeature\tn\tType\tPatient" >~/smartas/notebook/data/mutations/proteome_mutations.txt
grep -v ^Cancer smartas/analyses/????/mutations/proteome_mutations.txt | cut -d':' -f2- >>~/smartas/notebook/data/mutations/proteome_mutations.txt
checkFile mutations proteome_mutations.txt

# proteome information
echo -e "Cancer\tGene\tSymbol\tTranscript\tTPM\tProteinLength\tasEvidence" >~/smartas/notebook/data/mutations/proteome_information.txt
grep -v ^Cancer smartas/analyses/????/mutations/proteome_information.txt | cut -d':' -f2- >>~/smartas/notebook/data/mutations/proteome_information.txt
checkFile mutations proteome_information.txt

# switches features
echo -e "Cancer\tGene\tSymbol\tTranscript\tAnalysis\tFeature\tn\tFeatureLength\tStart\tEnd" >~/smartas/notebook/data/mutations/switch_features.txt
grep -v ^Cancer smartas/analyses/????/mutations/switch_features.txt | cut -d':' -f2- >>~/smartas/notebook/data/mutations/switch_features.txt
checkFile mutations switch_features.txt

# switches mutations
echo -e "Cancer\tGene\tSymbol\tTranscript\tAnalysis\tFeature\tn\tType\tPatient" >~/smartas/notebook/data/mutations/switch_mutations.txt
grep -v ^Cancer smartas/analyses/????/mutations/switch_mutations.txt | cut -d':' -f2- >>~/smartas/notebook/data/mutations/switch_mutations.txt
checkFile mutations switch_mutations.txt

# switches information
echo -e "Cancer\tGene\tSymbol\tTranscript\tTPM\tProteinLength\tasEvidence" >~/smartas/notebook/data/mutations/switch_information.txt
grep -v ^Cancer smartas/analyses/????/mutations/switch_information.txt | cut -d':' -f2- >>~/smartas/notebook/data/mutations/switch_information.txt
checkFile mutations switch_information.txt

# wgs mutations
echo -e "Tumor\tGene\tSymbol\tPatient\tPosition\tReference\tVariant" >~/smartas/notebook/data/mutations/wgs_mutations.txt
grep -v ^Tumor smartas/analyses/????/mutations/wgs_mutations.txt | cut -d':' -f2- >>~/smartas/notebook/data/mutations/wgs_mutations.txt
checkFile mutations wgs_mutations.txt

# wes mutations
echo -e "Tumor\tGene\tSymbol\tTranscript\tPatient\tStart\tEnd\tType\tMedianExpression" >~/smartas/notebook/data/mutations/wes_mutations.txt
grep -v ^Tumor smartas/analyses/????/mutations/wes_mutations.txt | cut -d':' -f2- >>~/smartas/notebook/data/mutations/wes_mutations.txt
checkFile mutations wes_mutations.txt

# me with top drivers
for i in $(seq 1 10)
do
    echo -e "Tumor\tGeneId\tSymbol\tNormal_transcript\tTumor_transcript\tMS\tM\tS\tN\tp.me\tadjp.me" >~/smartas/notebook/data/mutations/pannegative_mutual_exclusion.top_"$i"_drivers.txt
    grep -v ^Tumor smartas/analyses/????/mutations/pannegative_mutual_exclusion.top_"$i"_drivers.txt | cut -d':' -f2- >>~/smartas/notebook/data/mutations/pannegative_mutual_exclusion.top_"$i"_drivers.txt
    checkFile mutations pannegative_mutual_exclusion.top_"$i"_drivers.txt
done

echo -e "Tumor\tGeneId\tSymbol\tNormal_transcript\tTumor_transcript\tMS\tM\tS\tN\tp.me\tadjp.me" >~/smartas/notebook/data/mutations/pannegative_mutual_exclusion.top_all_drivers.txt
grep -v ^Tumor smartas/analyses/????/mutations/pannegative_mutual_exclusion.top_all_drivers.txt | cut -d':' -f2- >>~/smartas/notebook/data/mutations/pannegative_mutual_exclusion.top_all_drivers.txt
checkFile mutations pannegative_mutual_exclusion.top_all_drivers.txt

# top drivers
echo -e "Tumor\tGeneId\tSymbol\tSamples" >~/smartas/notebook/data/mutations/driver_mutation_number.txt
grep -v ^Tumor smartas/analyses/????/mutations/driver_mutation_number.txt | cut -d':' -f2- >>~/smartas/notebook/data/mutations/driver_mutation_number.txt
checkFile mutations driver_mutation_number.txt

echo -e "Tumor\tGeneId\tSymbol\tNormal_transcript\tTumor_transcript\tDriver\tDriverSymbol\tPathway\tDistance\tMS\tM\tS\tN\tp.me" >~/smartas/notebook/data/mutations/mutual_exclusion_top_drivers.txt
grep -v ^Tumor smartas/analyses/????/mutations/mutual_exclusion_top_drivers.txt | cut -d':' -f2- >>~/smartas/notebook/data/mutations/mutual_exclusion_top_drivers.txt
checkFile mutations mutual_exclusion_top_drivers.txt

if [[ "$copyFlag" != "" ]]
    then
    
    copyFile mutations gene_functional_mutations_all_switches.txt
    copyFile mutations gene_wgs_mutations_all_switches.txt
    copyFile mutations mutual_exclusion_top_drivers.txt
    
    copyFile mutations proteome_features.txt
    copyFile mutations proteome_mutations.txt
    copyFile mutations proteome_information.txt
    copyFile mutations switch_features.txt
    copyFile mutations switch_mutations.txt
    copyFile mutations switch_information.txt
    copyFile mutations wgs_mutations.txt
    copyFile mutations wes_mutations.txt

    for i in $(seq 1 10)
    do
        copyFile mutations pannegative_mutual_exclusion.top_"$i"_drivers.txt
    done

    copyFile mutations pannegative_mutual_exclusion.top_all_drivers.txt
    copyFile mutations driver_mutation_number.txt

    dest=smartas/projects_rg/mutations/v$version
    ln -s smartas/projects_rg/comet/ $dest

    smartas/pipeline/figures/mutations_plots.r

fi

# COPY PANCANCER DATA
~/smartas/pipeline/scripts/create_macrotable.R
rm -r ~/smartas/notebook/data/pancancer
cp -r smartas/analyses/pancancer ~/smartas/notebook/data

smartas/pipeline/figures/summary_plots.r 