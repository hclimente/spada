mkdir -p testResults/TCGA/analysis/

#num paired patients
while read line; do number=$((`ls Data/Input/TCGA/$line | wc -l`/2)) ; echo -e $line $number; done <Data/Input/TCGA_bigTags.lst >testResults/TCGA/analysis/numPatients.txt

#accuracy info
tail -n2 testResults/TCGA/*/result_summary/accuracy.tsv | grep -v '==>' | grep -v '^$' | grep kmeans | cut -f2-6 >testResults/TCGA/analysis/kmeans.txt
tail -n2 testResults/TCGA/*/result_summary/accuracy.tsv | grep -v '==>' | grep -v '^$' | grep hclust | cut -f2-6 >testResults/TCGA/analysis/hclust.txt

#Switch info
wc -l testResults/TCGA/*/candidateList_v3.tsv | sed 's/\ \+/\t/g' | cut -f2,3 | cut -d'/' -f1,3 | cut -d'_' -f1 | sed 's/Results\///' | grep -v total | awk '{print $2 " " $1-1}' >testResults/TCGA/analysis/switchNum.txt


grep [^_]CDS_study testResults/TCGA/*/result_summary/switches.tsv | cut -d':' -f2 | cut -f1,3-7 >testResults/TCGA/analysis/CDS_study.txt
grep [^_]CDS_change testResults/TCGA/*/result_summary/switches.tsv | cut -d':' -f2 | cut -f1,3-5 >testResults/TCGA/analysis/CDS_change.txt
grep [^_]UTR_change testResults/TCGA/*/result_summary/switches.tsv | cut -d':' -f2 | cut -f1,3-5 >testResults/TCGA/analysis/UTR_change.txt
grep [^_]Driver_affection testResults/TCGA/*/result_summary/switches.tsv | cut -d':' -f2 | cut -f1,3-5 >testResults/TCGA/analysis/Driver_affection.txt

#structural info
tail -n1 testResults/TCGA/*/result_summary/structural_loops.tsv | grep -v '==>' | grep -v '^$' >testResults/TCGA/analysis/loops.txt

grep functional_change testResults/TCGA/*/result_summary/structural_function.tsv | cut -d':' -f2 | cut -f1,3-5  >testResults/TCGA/analysis/functional_change.txt
grep change_type testResults/TCGA/*/result_summary/structural_function.tsv | cut -d':' -f2 | cut -f1,3-7  >testResults/TCGA/analysis/change_type.txt
grep motifs testResults/TCGA/*/result_summary/structural_function.tsv | cut -d':' -f2 | cut -f1,3-4  >testResults/TCGA/analysis/specificFeatures.txt

grep disordered_change testResults/TCGA/*/result_summary/structural_disorder.tsv | cut -d':' -f2 | cut -f1,3-5 >testResults/TCGA/analysis/disordered_change.txt
grep mean_length testResults/TCGA/*/result_summary/structural_disorder.tsv | cut -d':' -f2 | cut -f1,3 >testResults/TCGA/analysis/disordered_meanLength.txt