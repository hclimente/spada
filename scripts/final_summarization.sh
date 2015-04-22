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

# isoform length
ls ~/testResults/TCGA/*/result_summary/isoform_length_onlyModels.tsv | xargs -n 1 tail -n +2 >~/TCGA_analysis/isoform_length$tag.tsv

# structural
grep -v ^Cancer ~/testResults/TCGA/*/result_summary/structural_summary_onlyModels.tsv | cut -d':' -f2- >~/TCGA_analysis/structural_summary$tag.onlyModels.tsv
grep -v ^Cancer ~/testResults/TCGA/*/result_summary/structural_loops_onlyModels.tsv | cut -d':' -f2- >~/TCGA_analysis/structural_loops$tag.onlyModels.tsv
grep -v ^Cancer ~/testResults/TCGA/*/result_summary/structural_features_onlyModels.tsv | cut -d':' -f2- >~/TCGA_analysis/structural_features$tag.onlyModels.tsv

cancerTypes='brca coad hnsc kich kirc kirp lihc luad lusc prad thca'

for knsur in $cancerTypes
do
	cp ~/testResults/TCGA/$knsur/candidateList_v5.tsv ~/TCGA_analysis/$knsur.candidateList$tag.tsv
	cp ~/testResults/TCGA/$knsur/candidateList_v5.tsv /projects_rg/TCGA/users/hector/switches/$knsur.candidateList$tag.tsv

	# protein_centrality.tsv nIso_length.tsv tIso_length.tsv
	cp ~/testResults/TCGA/$knsur/result_summary/protein_centrality.tsv ~/TCGA_analysis/$knsur.protein_centrality$tag.tsv

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