switches=$1

if [ "$switches" == "all" ]
    then
    tag=""
elif [ "$switches" == "relevant" ]
    then
    tag="_relevant"
elif [ "$switches" == "nonrelevant" ]
    then
    tag="_nonrelevant"
else
    echo Error
    exit
fi

# switches.tsv
grep CDS_change ~/testResults/TCGA/*/result_summary/switches_onlyModels.tsv | cut -d':' -f2 >~/TCGA_analysis/CDS_change$tag.tsv
grep CDS_study ~/testResults/TCGA/*/result_summary/switches_onlyModels.tsv | cut -d':' -f2 >~/TCGA_analysis/CDS_study$tag.tsv
grep UTR_change ~/testResults/TCGA/*/result_summary/switches_onlyModels.tsv | cut -d':' -f2 >~/TCGA_analysis/UTR_change$tag.tsv

grep Driver_D0_enrichment ~/testResults/TCGA/*/result_summary/switches_onlyModels.tsv | cut -d':' -f2 >~/TCGA_analysis/Driver_D0_enrichment$tag.tsv
grep Driver_D0_patients[^_] ~/testResults/TCGA/*/result_summary/switches_onlyModels.tsv | cut -d':' -f2 >~/TCGA_analysis/Driver_D0_patients$tag.tsv
grep Driver_D0_patients_relevant ~/testResults/TCGA/*/result_summary/switches_onlyModels.tsv | cut -d':' -f2 >~/TCGA_analysis/Driver_D0_patientsRelevant$tag.tsv
grep Driver_D1_enrichment ~/testResults/TCGA/*/result_summary/switches_onlyModels.tsv | cut -d':' -f2 >~/TCGA_analysis/Driver_D1_enrichment$tag.tsv
grep Driver_D1_patients[^_] ~/testResults/TCGA/*/result_summary/switches_onlyModels.tsv | cut -d':' -f2 >~/TCGA_analysis/Driver_D1_patients$tag.tsv
grep Driver_D1_patients_relevant ~/testResults/TCGA/*/result_summary/switches_onlyModels.tsv | cut -d':' -f2 >~/TCGA_analysis/Driver_D1_patientsRelevant$tag.tsv
grep Driver_relevance_enrichment ~/testResults/TCGA/*/result_summary/switches_onlyModels.tsv | cut -d':' -f2 >~/TCGA_analysis/Driver_D0_relevance$tag.tsv
grep Driver_d1_relevance_enrichment ~/testResults/TCGA/*/result_summary/switches_onlyModels.tsv | cut -d':' -f2 >~/TCGA_analysis/Driver_D1_relevance$tag.tsv

# exons.tsv
ls ~/testResults/TCGA/*/result_summary/exons_onlyModels.tsv | xargs -n 1 tail -n +2 >~/TCGA_analysis/exons$tag.tsv
cut -f1,2 ~/TCGA_analysis/exons$tag.tsv | sort | uniq -c | sed 's/^ \+//' | sed 's/ /\t/g' >~/TCGA_analysis/exonsPerSwitch$tag.tsv

echo -e "Cancer\tRandom\tGene\tSymbol\tnTranscript\ttTranscript\tTag\tOrfChange\tnormalSegment\ttumorSegment" >~/TCGA_analysis/exons_new$tag.tsv
grep -v ^Cancer ~/testResults/TCGA/*/result_summary/exons_new_onlyModels.tsv | cut -d':' -f2 >>~/TCGA_analysis/exons_new$tag.tsv
echo -e "Cancer\tSymbol\tAnnotation\tPartner\tPartnerAnnotation\tWhatsHappening\tInteractionAffection\tSequenceCover\tPartnerCover\tGene\tnormalTranscript\ttumorTranscript\tUniprot\tPartnerUniprot\tSequenceInformation\tIsoformSpecific\tPDB\tpymolCommands"  >~/TCGA_analysis/i3d_broken$tag.tsv
grep -v ^Cancer ~/testResults/TCGA/*/i3d/i3d_broken.tsv | cut -d':' -f2 >>~/TCGA_analysis/i3d_broken$tag.tsv
cp ~/TCGA_analysis/i3d_broken$tag.tsv /projects_rg/TCGA/users/hector/

# isoform length
echo -e "Cancer\tRandom\tnIsoLength\ttIsoLength\tnIsoSpecificLength\ttIsoSpecificLength" >~/TCGA_analysis/isoform_length$tag.tsv
ls ~/testResults/TCGA/*/result_summary/isoform_length_onlyModels.tsv | xargs -n 1 tail -n +2 >>~/TCGA_analysis/isoform_length$tag.tsv

# structural
grep -v ^Cancer ~/testResults/TCGA/*/result_summary/structural_summary_onlyModels.tsv | cut -d':' -f2- >~/TCGA_analysis/structural_summary$tag.onlyModels.tsv
grep -v ^Cancer ~/testResults/TCGA/*/result_summary/structural_loops_onlyModels.tsv | cut -d':' -f2- >~/TCGA_analysis/structural_loops$tag.onlyModels.tsv
grep -v ^Cancer ~/testResults/TCGA/*/result_summary/structural_features_onlyModels.tsv | cut -d':' -f2- >~/TCGA_analysis/structural_features$tag.onlyModels.tsv

# mutations
echo -e "Gene\tSymbol\tCancer\tp_me\tpadj_me\tp_o\tpadj_o\tms\tm\ts\tn" >~/TCGA_analysis/gene_functional_mutations_all_switches.txt
cat ~/testResults/TCGA/*/mutations/gene_functional_mutations_all_switches.txt >>~/TCGA_analysis/gene_functional_mutations_all_switches.txt
echo -e "Gene\tSymbol\tCancer\tp_me\tpadj_me\tp_o\tpadj_o\tms\tm\ts\tn" >~/TCGA_analysis/gene_all_mutations_all_switches.txt
cat ~/testResults/TCGA/*/mutations/gene_all_mutations_all_switches.txt >>~/TCGA_analysis/gene_all_mutations_all_switches.txt
echo -e "Gene\tSymbol\tCancer\tp_me\tpadj_me\tp_o\tpadj_o\tms\tm\ts\tn" >~/TCGA_analysis/gene_functional_mutations_functional_switches.txt
cat ~/testResults/TCGA/*/mutations/gene_functional_mutations_functional_switches.txt >>~/TCGA_analysis/gene_functional_mutations_functional_switches.txt
echo -e "Gene\tSymbol\tCancer\tp_me\tpadj_me\tp_o\tpadj_o\tms\tm\ts\tn" >~/TCGA_analysis/gene_all_mutations_functional_switches.txt
cat ~/testResults/TCGA/*/mutations/gene_all_mutations_functional_switches.txt >>~/TCGA_analysis/gene_all_mutations_functional_switches.txt

echo -e "Gene\tSymbol\tCancer\tHallmark\tp\tpadj\tms\tm\ts\tn\tGeneset_genes" >~/TCGA_analysis/geneset_functional_mutations_all_switches.txt
cat ~/testResults/TCGA/*/mutations/geneset_functional_mutations_all_switches.txt >>~/TCGA_analysis/geneset_functional_mutations_all_switches.txt
echo -e "Gene\tSymbol\tCancer\tHallmark\tp\tpadj\tms\tm\ts\tn\tGeneset_genes" >~/TCGA_analysis/geneset_all_mutations_all_switches.txt
cat ~/testResults/TCGA/*/mutations/geneset_all_mutations_all_switches.txt >>~/TCGA_analysis/geneset_all_mutations_all_switches.txt
echo -e "Gene\tSymbol\tCancer\tHallmark\tp\tpadj\tms\tm\ts\tn\tGeneset_genes" >~/TCGA_analysis/geneset_functional_mutations_functional_switches.txt
cat ~/testResults/TCGA/*/mutations/geneset_functional_mutations_functional_switches.txt >>~/TCGA_analysis/geneset_functional_mutations_functional_switches.txt
echo -e "Gene\tSymbol\tCancer\tHallmark\tp\tpadj\tms\tm\ts\tn\tGeneset_genes" >~/TCGA_analysis/geneset_all_mutations_functional_switches.txt
cat ~/testResults/TCGA/*/mutations/geneset_all_mutations_functional_switches.txt >>~/TCGA_analysis/geneset_all_mutations_functional_switches.txt

echo -e "Gene\tSymbol\tCancer\tHallmark\tp\tpadj\tms\tm\ts\tn\tGeneset_genes" >~/TCGA_analysis/geneset_functional_mutations_all_switches_onlyDrivers.txt
cat ~/testResults/TCGA/*/mutations/geneset_functional_mutations_all_switches_onlyDrivers.txt >>~/TCGA_analysis/geneset_functional_mutations_all_switches_onlyDrivers.txt
echo -e "Gene\tSymbol\tCancer\tHallmark\tp\tpadj\tms\tm\ts\tn\tGeneset_genes" >~/TCGA_analysis/geneset_all_mutations_all_switches_onlyDrivers.txt
cat ~/testResults/TCGA/*/mutations/geneset_all_mutations_all_switches_onlyDrivers.txt >>~/TCGA_analysis/geneset_all_mutations_all_switches_onlyDrivers.txt
echo -e "Gene\tSymbol\tCancer\tHallmark\tp\tpadj\tms\tm\ts\tn\tGeneset_genes" >~/TCGA_analysis/geneset_functional_mutations_functional_switches_onlyDrivers.txt
cat ~/testResults/TCGA/*/mutations/geneset_functional_mutations_functional_switches_onlyDrivers.txt >>~/TCGA_analysis/geneset_functional_mutations_functional_switches_onlyDrivers.txt
echo -e "Gene\tSymbol\tCancer\tHallmark\tp\tpadj\tms\tm\ts\tn\tGeneset_genes" >~/TCGA_analysis/geneset_all_mutations_functional_switches_onlyDrivers.txt
cat ~/testResults/TCGA/*/mutations/geneset_all_mutations_functional_switches_onlyDrivers.txt >>~/TCGA_analysis/geneset_all_mutations_functional_switches_onlyDrivers.txt

echo -e "Gene\tSymbol\tCancer\tp\tpadj\tms\tm\ts\tn" >~/TCGA_analysis/pannegative_functional_mutations_all_switches.txt
cat ~/testResults/TCGA/*/mutations/pannegative_functional_mutations_all_switches.txt >>~/TCGA_analysis/pannegative_functional_mutations_all_switches.txt
echo -e "Gene\tSymbol\tCancer\tp\tpadj\tms\tm\ts\tn" >~/TCGA_analysis/pannegative_all_mutations_all_switches.txt
cat ~/testResults/TCGA/*/mutations/pannegative_all_mutations_all_switches.txt >>~/TCGA_analysis/pannegative_all_mutations_all_switches.txt
echo -e "Gene\tSymbol\tCancer\tp\tpadj\tms\tm\ts\tn" >~/TCGA_analysis/pannegative_functional_mutations_functional_switches.txt
cat ~/testResults/TCGA/*/mutations/pannegative_functional_mutations_functional_switches.txt >>~/TCGA_analysis/pannegative_functional_mutations_functional_switches.txt
echo -e "Gene\tSymbol\tCancer\tp\tpadj\tms\tm\ts\tn" >~/TCGA_analysis/pannegative_all_mutations_functional_switches.txt
cat ~/testResults/TCGA/*/mutations/pannegative_all_mutations_functional_switches.txt >>~/TCGA_analysis/pannegative_all_mutations_functional_switches.txt

cp ~/TCGA_analysis/gene_functional_mutations_all_switches.txt ~/TCGA_analysis/geneset_functional_mutations_all_switches.txt ~/TCGA_analysis/gene_all_mutations_all_switches.txt ~/TCGA_analysis/geneset_all_mutations_all_switches.txt ~/TCGA_analysis/gene_functional_mutations_functional_switches.txt ~/TCGA_analysis/geneset_functional_mutations_functional_switches.txt ~/TCGA_analysis/gene_all_mutations_functional_switches.txt ~/TCGA_analysis/geneset_all_mutations_functional_switches.txt ~/TCGA_analysis/geneset_functional_mutations_all_switches_onlyDrivers.txt ~/TCGA_analysis/geneset_all_mutations_all_switches_onlyDrivers.txt ~/TCGA_analysis/geneset_functional_mutations_functional_switches_onlyDrivers.txt ~/TCGA_analysis/geneset_all_mutations_functional_switches_onlyDrivers.txt ~/TCGA_analysis/pannegative_all_mutations_all_switches.txt ~/TCGA_analysis/pannegative_all_mutations_functional_switches.txt ~/TCGA_analysis/pannegative_functional_mutations_all_switches.txt ~/TCGA_analysis/pannegative_functional_mutations_functional_switches.txt /projects_rg/TCGA/users/hector/mutation_comparison

# neighborhoods
echo -e "GeneSet\tCancer\tpval\tqval\tNormalizedAverageAffection\tSwitchingGenes\tOR\tsg\tsng\tnsg\tnsng" >~/TCGA_analysis/canonical_pathways.txt
grep -v ^GeneSet testResults/TCGA/*/neighborhood_analysis/canonical_pathways_onlyModels.txt  | cut -d':' -f2- >>~/TCGA_analysis/canonical_pathways.txt
echo -e "GeneSet\tCancer\tpval\tqval\tNormalizedAverageAffection\tSwitchingGenes\tOR\tsg\tsng\tnsg\tnsng" >~/TCGA_analysis/hallmarks.txt
grep -v ^GeneSet testResults/TCGA/*/neighborhood_analysis/hallmarks_onlyModels.txt  | cut -d':' -f2- >>~/TCGA_analysis/hallmarks.txt
echo -e "GeneSet\tCancer\tpval\tqval\tNormalizedAverageAffection\tSwitchingGenes\tOR\tsg\tsng\tnsg\tnsng" >~/TCGA_analysis/go_biological_process.txt
grep -v ^GeneSet testResults/TCGA/*/neighborhood_analysis/go_biological_process_onlyModels.txt  | cut -d':' -f2- >>~/TCGA_analysis/go_biological_process.txt
echo -e "GeneSet\tCancer\tpval\tqval\tNormalizedAverageAffection\tSwitchingGenes\tOR\tsg\tsng\tnsg\tnsng" >~/TCGA_analysis/oncogenic_signatures.txt
grep -v ^GeneSet testResults/TCGA/*/neighborhood_analysis/oncogenic_signatures_onlyModels.txt  | cut -d':' -f2- >>~/TCGA_analysis/oncogenic_signatures.txt

cancerTypes='brca coad hnsc kich kirc kirp lihc luad lusc prad thca'

for knsur in $cancerTypes
do
    cp ~/testResults/TCGA/$knsur/candidateList_v6.tsv ~/TCGA_analysis/$knsur.candidateList$tag.tsv
    cp ~/testResults/TCGA/$knsur/candidateList_v6.tsv /projects_rg/TCGA/users/hector/switches/$knsur.candidateList$tag.tsv

    # protein_centrality.tsv nIso_length.tsv tIso_length.tsv
    # cp ~/testResults/TCGA/$knsur/result_summary/protein_centrality.tsv ~/TCGA_analysis/$knsur.protein_centrality$tag.tsv

    # interpro anchor iupred prosite
    cp ~/testResults/TCGA/$knsur/structural_analysis/anchor_analysis_random.tsv ~/TCGA_analysis/$knsur.anchor_analysis_random$tag.tsv
    cp ~/testResults/TCGA/$knsur/structural_analysis/iupred_analysis_random.tsv ~/TCGA_analysis/$knsur.iupred_analysis_random$tag.tsv
    cp ~/testResults/TCGA/$knsur/structural_analysis/interpro_analysis_random.tsv ~/TCGA_analysis/$knsur.interpro_analysis_random$tag.tsv
    cp ~/testResults/TCGA/$knsur/structural_analysis/prosite_analysis_random.tsv ~/TCGA_analysis/$knsur.prosite_analysis_random$tag.tsv
    cp ~/testResults/TCGA/$knsur/structural_analysis/anchor_analysis.tsv ~/TCGA_analysis/$knsur.anchor_analysis$tag.tsv
    cp ~/testResults/TCGA/$knsur/structural_analysis/iupred_analysis.tsv ~/TCGA_analysis/$knsur.iupred_analysis$tag.tsv
    cp ~/testResults/TCGA/$knsur/structural_analysis/interpro_analysis.tsv ~/TCGA_analysis/$knsur.interpro_analysis$tag.tsv
    cp ~/testResults/TCGA/$knsur/structural_analysis/prosite_analysis.tsv ~/TCGA_analysis/$knsur.prosite_analysis$tag.tsv
done


hallmark_files=`ls testResults/TCGA/*/mutations/hallmark_info/*tsv | cut -d'/' -f6 | sort | uniq`
mkdir ~/TCGA_analysis/hallmark_info

for hll_file in $hallmark_files
do
    rm ~/TCGA_analysis/hallmark_info/$hll_file
    echo -e "Patient\tGene\tState\tCancer" >~/TCGA_analysis/hallmark_info/$hll_file
    for knsur in $cancerTypes
    do
        grep -v ^patient testResults/TCGA/$knsur/mutations/hallmark_info/$hll_file | cut -d':' -f2- | awk -v knsur=$knsur 'BEGIN{OFS="\t";} {print $0,knsur }' >>~/TCGA_analysis/hallmark_info/$hll_file
    done
done

echo -e "Patient\tCancer" >~/TCGA_analysis/hallmark_info/patient_annotation.txt
grep -v Patient ~/TCGA_analysis/hallmark_info/HALLMARK_DNA_REPAIR.tsv | cut -f1,4 | sort -k2 | uniq >>~/TCGA_analysis/hallmark_info/patient_annotation.txt
