cancerTypes='brca coad hnsc kich kirc kirp lihc luad lusc prad thca'

echo "######################"
echo "##     SWITCHES     ##"
echo "######################"

echo Delete previous
rm -r ~/smartas/notebook/data/switches
mkdir ~/smartas/notebook/data/switches

echo candidateList
for knsur in $cancerTypes
do
    cp ~/smartas/analyses/$knsur/candidateList_info.tsv ~/smartas/notebook/data/switches/$knsur.candidateList.tsv
done

echo switches
echo -e "Cancer\tAnalysis\tBoth\tOnly_nIso\tOnly_tIso\tNone\tRandom_Both\tRandom_Only_nIso\tRandom_Only_tIso\tRandom_None" >~/smartas/notebook/data/switches/$a.tsv
grep CDS_study ~/smartas/analyses/????/result_summary/switches_onlyModels.tsv | cut -d':' -f2 >>~/smartas/notebook/data/switches/CDS_study.tsv

analyses='CDS_change UTR_change'
for a in $analyses
do
    echo -e "Cancer\tAnalysis\tYes\tNo\tRandom_Yes\tRandom_No\tp\tOR" >~/smartas/notebook/data/switches/$a.tsv
    grep $a ~/smartas/analyses/????/result_summary/switches_onlyModels.tsv | cut -d':' -f2 >>~/smartas/notebook/data/switches/$a.tsv
done

analyses="d0_enrichment d1_enrichment d0_functional_enrichment d1_functional_enrichment"
for a in $analyses
do
    echo -e "Cancer\tAnalysis\tFS\tFNS\tNFS\tNFNS\tp.me\tOR" >~/smartas/notebook/data/switches/$a.tsv
    grep $a ~/smartas/analyses/????/result_summary/switches_onlyModels.tsv | cut -d':' -f2 >~/smartas/notebook/data/switches/$a.tsv
done

echo exons
echo -e "Cancer\tRandom\tSwitch\tOrigin\tType\tLength\tCDSLength\tCDSRelativeSize\tPosition\tKeepOrf" >~/smartas/notebook/data/switches/exons.tsv
ls ~/smartas/analyses/????/result_summary/exons_onlyModels.tsv | xargs -n 1 tail -n +2 >>~/smartas/notebook/data/switches/exons.tsv

cut -f1,2 ~/smartas/notebook/data/switches/exons.tsv | sort | uniq -c | sed 's/^ \+//' | sed 's/ /\t/g' >~/smartas/notebook/data/switches/exonsPerSwitch.tsv

echo -e "Cancer\tRandom\tGene\tSymbol\tnTranscript\ttTranscript\tTag\tOrfChange\tnormalSegment\ttumorSegment" >~/smartas/notebook/data/switches/exons_new.tsv
grep -v ^Cancer ~/smartas/analyses/????/result_summary/exons_new_onlyModels.tsv | cut -d':' -f2 >>~/smartas/notebook/data/switches/exons_new.tsv

echo isoform length
echo -e "Cancer\tRandom\tnIsoLength\ttIsoLength\tnIsoSpecificLength\ttIsoSpecificLength" >~/smartas/notebook/data/switches/isoform_length.tsv
ls ~/smartas/analyses/????/result_summary/isoform_length_onlyModels.tsv | xargs -n 1 tail -n +2 >>~/smartas/notebook/data/switches/isoform_length.tsv

echo evidence of driverness
echo -e "Tumor\tGeneId\tSymbol\tNormal_transcript\tTumor_transcript\tRecurrence\tAffects_mutated_feature\tPPI\tPannegative\tDriverME" >~/smartas/notebook/data/switches/driverEvidence.tsv
grep -v ^Tumor ~/smartas/analyses/????/candidateList_driverEvidence.tsv | cut -d':' -f2 >>~/smartas/notebook/data/switches/driverEvidence.tsv

echo "######################"
echo "##   NEIGHBORHOOD   ##"
echo "######################"

analyses='canonical_pathways hallmarks go_biological_process oncogenic_signatures'
genesubgroups="all functional"

for a in $analyses
do
    for g in $genesubgroups
    do
        echo -e "GeneSet\tCancer\tpval\tqval\tNormalizedAverageAffection\tSwitchingGenes\tOR\tsg\tsng\tnsg\tnsng" >~/smartas/notebook/data/switches/"$a"_$g.txt
        grep -v ^GeneSet smartas/analyses/????/neighborhood_analysis/"$a"_"$g"_onlyModels.txt  | cut -d':' -f2- >>~/smartas/notebook/data/switches/"$a"_$g.txt
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
done

echo "######################"
echo "##    STRUCTURAL    ##"
echo "######################"

echo Delete previous
rm -r ~/smartas/notebook/data/structural_analysis
mkdir ~/smartas/notebook/data/structural_analysis

echo Copy Eduard\'s files
cp -r /projects_rg/TCGA/users/hector/SmartAS/structural_analysis/pfam_go_term ~/smartas/notebook/data/structural_analysis/
cp -r /projects_rg/TCGA/users/hector/SmartAS/eporta/comet_output/ ~/smartas/notebook/data/structural_analysis/
cp -r /projects_rg/TCGA/users/hector/SmartAS/eporta/estimate/ ~/smartas/notebook/data/switches/
cp /projects_rg/TCGA/users/hector/SmartAS/eporta/Switched_interactions_consensus.txt /projects_rg/TCGA/users/hector/SmartAS/eporta/num_patients_with_each_switch.txt ~/smartas/notebook/data/structural_analysis/

analyses='anchor_analysis iupred_analysis interpro_analysis prosite_analysis structural_summary'

echo Structural analyses
for knsur in $cancerTypes
do
    for a in $analyses
    do
        cp ~/smartas/analyses/$knsur/structural_analysis/"$a"_random.tsv ~/smartas/notebook/data/structural_analysis/$knsur."$a"_random.tsv
        cp ~/smartas/analyses/$knsur/structural_analysis/$a.tsv ~/smartas/notebook/data/structural_analysis/$knsur.$a.tsv
    done
done

echo structural_features
echo -e "Cancer\tGene\tSymbol\tnTx\ttTx\tRandom\tAnalysis\tWhatsHappenning\tFeature\tDriver\tASDriver\tDriverType" >~/smartas/notebook/data/structural_analysis/structural_features.onlyModels.tsv
grep -v ^Cancer ~/smartas/analyses/????/result_summary/structural_features_onlyModels.tsv | cut -d':' -f2- >>~/smartas/notebook/data/structural_analysis/structural_features.onlyModels.tsv

echo interpro_analysis
echo -e "Gene\tSymbol\tNormalTranscript\tTumorTranscript\tWhat\tFeature\tnormalReps\ttumorReps\tnMacroScore\tnMicroScore\tnJaccard\ttMacroScore\ttMicroScore\ttJaccard" >~/smartas/notebook/data/structural_analysis/interpro_analysis.tsv
grep -v ^Gene ~/smartas/analyses/????/structural_analysis/interpro_analysis.tsv | cut -d':' -f2- | sort | uniq >>~/smartas/notebook/data/structural_analysis/interpro_analysis.tsv

echo prosite_analysis
echo -e "Gene\tSymbol\tNormalTranscript\tTumorTranscript\tWhat\tFeature\tnormalReps\ttumorReps\tnMacroScore\tnMicroScore\tnJaccard\ttMacroScore\ttMicroScore\ttJaccard" >~/smartas/notebook/data/structural_analysis/prosite_analysis.tsv
grep -v ^Gene ~/smartas/analyses/????/structural_analysis/prosite_analysis.tsv | cut -d':' -f2- | sort | uniq >>~/smartas/notebook/data/structural_analysis/prosite_analysis.tsv

echo "######################"
echo "##     MUTATIONS    ##"
echo "######################"

echo Delete previous
rm -r ~/smartas/notebook/data/mutations
mkdir ~/smartas/notebook/data/mutations

echo WES ME analysis
echo -e "Tumor\tGeneId\tSymbol\tNormal_transcript\tTumor_transcript\tMS\tM\tS\tN\tp\tOR\tSwitched\tMutated" >~/smartas/notebook/data/mutations/gene_functional_mutations_all_switches.txt
grep -v ^Tumor smartas/analyses/????/mutations/gene_functional_mutations_all_switches.txt | cut -d':' -f2- >>~/smartas/notebook/data/mutations/gene_functional_mutations_all_switches.txt

echo WGS coocurrence analysis
echo -e "Tumor\tGeneId\tSymbol\tNormal_transcript\tTumor_transcript\tMS\tM\tS\tN\tp.o" >~/smartas/notebook/data/mutations/gene_wgs_mutations_all_switches.txt
grep -v ^Tumor smartas/analyses/????/mutations/gene_wgs_mutations_all_switches.txt | cut -d':' -f2- >>~/smartas/notebook/data/mutations/gene_wgs_mutations_all_switches.txt

echo proteome features
echo -e "Cancer\tGene\tSymbol\tTranscript\tAnalysis\tFeature\tn\tFeatureLength\tStart\tEnd" >~/smartas/notebook/data/mutations/proteome_features.txt
grep -v ^Cancer smartas/analyses/????/mutations/proteome_features.txt | cut -d':' -f2- >>~/smartas/notebook/data/mutations/proteome_features.txt

echo proteome mutations
echo -e "Cancer\tGene\tSymbol\tTranscript\tAnalysis\tFeature\tn\tType\tPatient" >~/smartas/notebook/data/mutations/proteome_mutations.txt
grep -v ^Cancer smartas/analyses/????/mutations/proteome_mutations.txt | cut -d':' -f2- >>~/smartas/notebook/data/mutations/proteome_mutations.txt

echo proteome information
echo -e "Cancer\tGene\tSymbol\tTranscript\tTPM\tProteinLength\tasEvidence" >~/smartas/notebook/data/mutations/proteome_information.txt
grep -v ^Cancer smartas/analyses/????/mutations/proteome_information.txt | cut -d':' -f2- >>~/smartas/notebook/data/mutations/proteome_information.txt

echo switches features
echo -e "Cancer\tGene\tSymbol\tTranscript\tAnalysis\tFeature\tn\tFeatureLength\tStart\tEnd" >~/smartas/notebook/data/mutations/switch_features.txt
grep -v ^Cancer smartas/analyses/????/mutations/switch_features.txt | cut -d':' -f2- >>~/smartas/notebook/data/mutations/switch_features.txt

echo switches mutations
echo -e "Cancer\tGene\tSymbol\tTranscript\tAnalysis\tFeature\tn\tType\tPatient" >~/smartas/notebook/data/mutations/switch_mutations.txt
grep -v ^Cancer smartas/analyses/????/mutations/switch_mutations.txt | cut -d':' -f2- >>~/smartas/notebook/data/mutations/switch_mutations.txt

echo switches information
echo -e "Cancer\tGene\tSymbol\tTranscript\tTPM\tProteinLength\tasEvidence" >~/smartas/notebook/data/mutations/switch_information.txt
grep -v ^Cancer smartas/analyses/????/mutations/switch_information.txt | cut -d':' -f2- >>~/smartas/notebook/data/mutations/switch_information.txt

echo wgs mutations
echo -e "Tumor\tGene\tSymbol\tPatient\tPosition\tReference\tVariant" >~/smartas/notebook/data/mutations/wgs_mutations.txt
grep -v ^Tumor smartas/analyses/????/mutations/wgs_mutations.txt | cut -d':' -f2- >>~/smartas/notebook/data/mutations/wgs_mutations.txt

echo wes mutations
echo -e "Tumor\tGene\tSymbol\tTranscript\tPatient\tStart\tEnd\tType\tMedianExpression" >~/smartas/notebook/data/mutations/wes_mutations.txt
grep -v ^Tumor smartas/analyses/????/mutations/wes_mutations.txt | cut -d':' -f2- >>~/smartas/notebook/data/mutations/wes_mutations.txt

echo me with top drivers
for i in $(seq 1 10)
do
    echo -e "Tumor\tGeneId\tSymbol\tNormal_transcript\tTumor_transcript\tMS\tM\tS\tN\tp.me\tadjp.me" >~/smartas/notebook/data/mutations/pannegative_mutual_exclusion.top_"$i"_drivers.txt
    grep -v ^Tumor smartas/analyses/????/mutations/pannegative_mutual_exclusion.top_"$i"_drivers.txt | cut -d':' -f2- >>~/smartas/notebook/data/mutations/pannegative_mutual_exclusion.top_"$i"_drivers.txt
done

echo -e "Tumor\tGeneId\tSymbol\tNormal_transcript\tTumor_transcript\tMS\tM\tS\tN\tp.me\tadjp.me" >~/smartas/notebook/data/mutations/pannegative_mutual_exclusion.top_all_drivers.txt
grep -v ^Tumor smartas/analyses/????/mutations/pannegative_mutual_exclusion.top_all_drivers.txt | cut -d':' -f2- >>~/smartas/notebook/data/mutations/pannegative_mutual_exclusion.top_all_drivers.txt

echo top drivers
echo -e "Tumor\tGeneId\tSymbol\tSamples" >~/smartas/notebook/data/mutations/driver_mutation_number.txt
grep -v ^Tumor smartas/analyses/????/mutations/driver_mutation_number.txt | cut -d':' -f2- >>~/smartas/notebook/data/mutations/driver_mutation_number.txt

echo -e "Tumor\tGeneId\tSymbol\tNormal_transcript\tTumor_transcript\tDriver\tDriverSymbol\tPathway\tDistance\tMS\tM\tS\tN\tp.me" >~/smartas/notebook/data/mutations/mutual_exclusion_top_drivers.txt
grep -v ^Tumor smartas/analyses/????/mutations/mutual_exclusion_top_drivers.txt | cut -d':' -f2- >>~/smartas/notebook/data/mutations/mutual_exclusion_top_drivers.txt

echo "######################"
echo "##     PANCANCER    ##"
echo "######################"

echo Creating macrotable
~/smartas/pipeline/scripts/create_macrotable.R

echo Delete previous and copying new one
rm -r ~/smartas/notebook/data/pancancer
cp -r smartas/analyses/pancancer ~/smartas/notebook/data

echo "######################"
echo "##      SUMMARY     ##"
echo "######################"

echo Creating additional summary plots
smartas/pipeline/figures/summary_plots.r 